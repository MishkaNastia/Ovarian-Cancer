import pandas as pd
from cobra.flux_analysis import single_reaction_deletion
import pandas as pd
import re


def rxn_essential(model):
    df_res_rxn = single_reaction_deletion(model, processes=10)
    df_essential=df_res_rxn[df_res_rxn["status"]!="optimal"] # remove non essential
    df_essential["rxn_ids"]=df_essential.ids.apply(lambda x: list(x)[0] if isinstance(x, (set, frozenset)) else x) # remove {} around ids

    return df_essential



def clean_gene_name(name: str) -> str:
    """
    Clean up a gene name so it matches the model's format.
    Removes extra info in parentheses, spaces, etc.
    Example: 'AOC3 (8639)' → 'AOC3'
    """
    cleaned = re.sub(r"\s*\(.*\)$", "", str(name))  # remove text in parentheses
    return cleaned.strip().upper()  # remove spaces, make uppercase for consistency

#IS REACTION DEFINED IN MEDIA
def is_defined_inmedia(rxn_id):
    """
    Check if a reaction is already defined in the media conditions.
    
    Parameters
    ----------
    rxn_id : str
        Reaction ID (e.g., 'EX_glc__D_e')

    Returns
    -------
    bool
        True if the reaction is marked as an exchange reaction
    """
    return rxn_id.startswith("EX")

def is_essential_reaction(df, rxn):
    """
    Check if a reaction ID is present in the dataframe of essential reactions.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing at least a 'rxn_ids' column.
    rxn : str
        Reaction ID to check.
    
    Returns
    -------
    bool
        True if rxn is in df['rxn_ids'], else False.
    """
    return rxn in df['rxn_ids'].values

#DO WE HAVE TRANSCIPTOMICS FOR THE REACTION GENES
def has_transcriptomics(reaction, CCLE_genes):
    """
    Checks whether transcriptomics data are available for the gene(s)
    associated with a given reaction.
    
    Parameters
    ----------
    reaction : cobra.Reaction
        Reaction object from the COBRA model.
    CCLE_expression : pandas.DataFrame
        Transcriptomics dataset where columns are gene IDs (e.g., ENSG IDs).
    
    Returns
    -------
    bool
        True if at least one gene in the reaction is found in CCLE_expression, else False.
    """
    
    # Get all gene IDs linked to this reaction
    reaction_genes = [gene.name for gene in reaction.genes]
    
    
    # Check if any gene in the reaction exists in the expression dataframe
    for g in reaction_genes:
        if g in CCLE_genes:
            return True  # we have transcriptomics for at least one gene
    
    # If none of the genes matched
    return False

def open_bounds(rxn):
    """
    Set reaction bounds fully open depending on reversibility.
    """
    if rxn.reversibility:
        rxn.bounds = (-1000, 1000)
    else:
        rxn.bounds = (0, 1000)

def clean_bounds(rxn, tol=1e-6):
    lb, ub = rxn.lower_bound, rxn.upper_bound

    # Round tiny bounds to zero
    if abs(lb) < tol:
        lb = 0.0
    if abs(ub) < tol:
        ub = 0.0

    # If still inverted, collapse to zero flux
    if lb > ub:
        lb, ub = 0.0, 0.0

    rxn.bounds = (lb, ub)

def classify_rule(rxn):
    rule = rxn.gene_reaction_rule.lower()
    if "and" in rule and "or" not in rule:
        return "and_rule"
    elif "or" in rule and "and" not in rule:
        return "or_rule"
    elif "and" not in rule and "or" not in rule and rule != "":
        return "one_gene"
    else:
        return None
    
def calculate_new_bounds(rxn, rule_type, cell_line, CCLE_name_map, CCLE_expression, tol=1e-6):
    matched_cols = [CCLE_name_map[g.name] for g in rxn.genes if g.name in CCLE_name_map]
    if not matched_cols:
        return

    expr_values = CCLE_expression.loc[cell_line, matched_cols].astype(float).tolist()

    if rule_type == "one_gene":
        E = expr_values[0]
    elif rule_type == "or_rule":
        E = sum(expr_values)
    elif rule_type == "and_rule":
        E = min(expr_values)
    else:
        return

    # Ensure E is valid and non-negative
    if not pd.notnull(E) or E < 0:
        return

    # Determine new bounds
    if rxn.reversibility:
        lb = -E
        ub = E
    else:
        lb = 0.0
        ub = E

    # --- Critical part: clean tiny noise BEFORE assigning ---
    if abs(lb) < tol: lb = 0.0
    if abs(ub) < tol: ub = 0.0
    if lb > ub:
        lb, ub = (0.0, 0.0)

    # --- Assign both bounds at once (avoids ValueError) ---
    rxn.bounds = (lb, ub)