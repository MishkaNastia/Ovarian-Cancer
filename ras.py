#load the libraries
import pandas as pd
import re
import numpy as np
from progressbar import ProgressBar,Bar,Percentage
from scanpy import AnnData
from cobra.flux_analysis.variability import find_essential_reactions,find_essential_genes
import scanpy.external as sce

"""
Class to compute the RAS values

"""


class RAS_computation:

    def __init__(self,adata,model):
                                                       
        self._logic_operators = ['and', 'or', '(', ')']
        self.val_nan = np.nan

        # Build the dictionary for the GPRs
        df_reactions = pd.DataFrame(index=[reaction.id for reaction in model.reactions])
        gene_rules=[reaction.gene_reaction_rule for reaction in model.reactions]   
        
        gene_rules=[el.replace("OR","or").replace("AND","and").replace("(","( ").replace(")"," )") for el in gene_rules]        
        df_reactions['rule'] = gene_rules
        df_reactions = df_reactions.reset_index()
        df_reactions = df_reactions.groupby('rule').agg(lambda x: sorted(list(x)))
        
        self.dict_rule_reactions = df_reactions.to_dict()['index']

        # build useful structures for RAS computation
        self.model = model
        self.count_adata = adata.copy()
        self.genes = self.count_adata.var.index.intersection([gene.id for gene in model.genes])
        
        #check if there is one gene at least 
        if len(self.genes)==0:
            print("ERROR: no gene of the count matrix is in the metabolic model!")
            print(" are you sure that the gene annotation is the same for the model and the count matrix?")
            return -1
        
        self.cell_ids = list(self.count_adata.obs.index.values)
        self.count_df_filtered = self.count_adata.to_df().T.loc[self.genes]
 
    def compute(self,
                or_expression=np.nansum,    # type of operation to do in case of an or expression (max, sum, mean)
                and_expression=np.nanmin,   # type of operation to do in case of an and expression(min, sum)
                drop_na_rows=True,          # if True remove the nan rows of the ras  matrix
                drop_duplicates=False,      # if true, remove duplicates rows
                regexp=re.compile(r"\([a-zA-Z0-9-.:\s]+\)"),  # regular expression inside a parenthesis
                print_progressbar=True,     # if True, print the progress bar
                add_count_metadata=True,    # if True add metadata of cells in the ras adata
                add_met_metadata=True,      # if True add metadata from the metabolic model (gpr and compartments of reactions)
                add_essential_reactions=False,
                add_essential_genes=False
                ):

        self.or_function = or_expression
        self.and_function = and_expression
        
        ras_df = pd.DataFrame(index=range(len(self.dict_rule_reactions)), columns=self.cell_ids)
        ras_df.loc[:, :] = self.val_nan
        
        if print_progressbar:
            pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=len(self.dict_rule_reactions)).start()
            i = 0
        
        # for loop on reactions
        ind = 0       
        for rule, reaction_ids in self.dict_rule_reactions.items():
            if len(rule) != 0:
                # there is one gene at least in the formula
                rule_split = rule.split()
                rule_split_elements = list(filter(lambda x: x not in self._logic_operators, rule_split))  # remove of all logical operators
                rule_split_elements = list(np.unique(rule_split_elements))                                # genes in formula
                
                # which genes are in the count matrix?                
                genes_in_count_matrix = list(set([el for el in rule_split_elements if el in self.genes]))
                genes_notin_count_matrix = list(set([el for el in rule_split_elements if el not in self.genes]))


                if len(genes_in_count_matrix) > 0: #there is at least one gene in the count matrix
                     if len(rule_split) == 1:
                         #one gene --> one reaction
                         ras_df.iloc[ind] = self.count_df_filtered.loc[genes_in_count_matrix]
                     else:                        
                        # more genes in the formula
                        lista = re.findall(regexp, rule)
                        if len(lista) == 0:
                             #or/and sequence
                             matrix = self.count_df_filtered.loc[genes_in_count_matrix].values
                             if len(genes_notin_count_matrix) > 0:
                                matrix = np.vstack([matrix, [self.val_nan for el in self.cell_ids]])

                             if 'or' in rule_split: 
                                ras_df.iloc[ind] = self.or_function(matrix, axis=0)
                             else:
                                ras_df.iloc[ind] = self.and_function(matrix, axis=0)
                        else:
                            # ho almeno una tonda
                            data = self.count_df_filtered.loc[genes_in_count_matrix]  # dataframe of genes in the GPRs
                            genes = data.index
                            j = 0
                             
                            for cellid in self.cell_ids:    #for loop on the cells
                                lista_cell = lista.copy()
                                rule_cell = rule
                                 
                                while len(lista_cell) > 0:
                                    #
                                    for el in lista_cell:
                                        #print(el[1:-1])
                                        value = self._evaluate_expression(el[1:-1].split(), data[cellid], genes)
                                        rule_cell = rule_cell.replace(el, str(value))   
                                    lista_cell = re.findall(regexp, rule_cell)      
         
                                ras_df.iloc[ind, j] = self._evaluate_expression(rule_cell.split(), data[cellid], genes)
                                j=j+1
      
            ind = ind+1
            #update percentage
            if print_progressbar:
                pbar.update(i+1)
                i = i+1
        
        if print_progressbar:
            pbar.finish()
        
        ras_df=ras_df.astype("float")    
        ras_df['REACTIONS'] = [reaction_ids for rule,reaction_ids in self.dict_rule_reactions.items()]
        
        reactions_common = pd.DataFrame()
        reactions_common["REACTIONS"] = ras_df['REACTIONS']
        reactions_common["proof2"] = ras_df['REACTIONS']
        reactions_common = reactions_common.explode('REACTIONS')
        reactions_common = reactions_common.set_index("REACTIONS")

        ras_df = ras_df.explode("REACTIONS")
        ras_df = ras_df.set_index("REACTIONS")

        if drop_na_rows:
            ras_df = ras_df.dropna(how="all")
            
        if drop_duplicates:
            ras_df = ras_df.drop_duplicates()
        
        #create AnnData structure for RAS
        ras_adata = AnnData(ras_df.T)

        #add metadata
        if add_count_metadata:
            ras_adata.var["common_gprs"] = reactions_common.loc[ras_df.index]
            ras_adata.var["common_gprs"] = ras_adata.var["common_gprs"].apply(lambda x: ",".join(x))
            for el in self.count_adata.obs.columns:
                ras_adata.obs["countmatrix_"+el]=self.count_adata.obs[el]

        if add_met_metadata:
            if len(self.model.compartments)>0:
                  ras_adata.var['compartments']=[list(self.model.reactions.get_by_id(reaction).compartments) for reaction in ras_adata.var.index]  
                  ras_adata.var['compartments']=ras_adata.var["compartments"].apply(lambda x: ",".join(x))
            
            ras_adata.var['GPR rule'] = [self.model.reactions.get_by_id(reaction).gene_reaction_rule for reaction in ras_adata.var.index]

        if add_essential_reactions:            
            essential_reactions=find_essential_reactions(self.model)
            essential_reactions=[el.id for el in essential_reactions]            
            ras_adata.var['essential reactions']=["yes" if el in essential_reactions else "no" for el in ras_adata.var.index]
        
        if add_essential_genes:
            essential_genes=find_essential_genes(self.model)
            essential_genes=[el.id for el in essential_genes]
            ras_adata.var['essential genes']=[" ".join([gene for gene in genes.split()  if gene in essential_genes]) for genes in ras_adata.var["GPR rule"]]

        
        return ras_adata


    def _check_number(self,value):
      try:
        float(value)
        return True
      except ValueError:
        return False

    def _evaluate_expression(self, rule_split, values_cellid, genes):
        
        #ci sono per forza solo or
        rule_split2 = list(filter(lambda x: x != "or" and x!="and", rule_split))   

        values = list()
        i=0
        for el in rule_split2:
             if self._check_number(el):
                 values.append(float(el))
             elif el in genes:
                 values.append(values_cellid[el])
             else:
                 values.append(self.val_nan)
                 i=i+1
                 
        if i==len(rule_split2):
            return self.val_nan
        if "or" in rule_split:
            #or sequence
            return self.or_function(values)
        else:
            #and sequence
            return self.and_function(values)
