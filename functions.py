import pandas as pd
from cobra.io.json import load_json_model
import numpy as np
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import single_reaction_deletion
import math
from tqdm import tqdm

from cobra.io import read_sbml_model

def rxn_essential_df(model):
    df_res_rxn = single_reaction_deletion(model, processes=10)
    df_essential=df_res_rxn[df_res_rxn["status"]!="optimal"] # remove non essential
    df_essential["rxn_ids"]=df_essential.ids.apply(lambda x: list(x)[0] if isinstance(x, (set, frozenset)) else x) # remove {} around ids

    return df_essential


def rxn_essential(model):
    
    df_essential = rxn_essential(model)

    def change_bounds(rxn_id):
        if rxn_id in model.reactions:
            rxn= model.reactions.get_by_id(rxn_id)
            lb,ub= rxn.lower_bound, rxn.upper_bound
            M=1000

            if ub <= 0 and lb < 0:
                rxn_bounds=(-M, 0.0)
            elif lb >= 0 and ub > 0:
                rxn_bounds=(0.0, +M)
            else:
                rxn_bounds=(-M, +M)
        
            
            rxn.bounds=rxn_bounds
        else:
            print("error reaction doesn't exist")
    
    df_essential["rxn_ids"].apply(change_bounds)
