#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dic 13 15:03:20 2023

@author: lizy
"""

import pandas as pd
import numpy as np
import math
import os.path
import shutil

os.makedirs('genes_folder', exist_ok=True)
path = "./genes_folder/"

def sum_pl(my_list, condition):  #select a list with the file name and condition
    '''
    Calculate the proteomic load sum of any list of E. coli genes, in 1 of the 11 Carbon source conditions (Glucose, Fructose, etc)    '''
    eco_db4691_final = "~/Notebooks/LIZ/Ntbooks_finales/infiles/allgen4691_fitfg_newname.csv"
    proteomic = "~/Notebooks/LIZ/Ntbooks_finales/infiles/proteomicf_2346.csv"
    
    my_list = path + my_list
    all_db =  pd.read_csv(eco_db4691_final)
    total_proteo = all_db.sum(numeric_only=True)
    df = pd.read_csv(my_list)
    df2 = pd.read_csv(proteomic)
    print('Total proteome =', total_proteo[condition], 'fg')
    print(f'in {condition} condition', '\n')
    
    lista2 = list(df["genes"])
    matching_genes = all_db[all_db["gene"].isin(lista2)]
    PL = matching_genes[condition].sum(numeric_only=True)
    tgenes = len(lista2)
    mw_prote = df2[df2["gene"].isin(lista2)]
    genes_withpl = len(mw_prote)

    print(tgenes, 'genes were found in the list')
    print(genes_withpl, 'with PL')
 ##   print(mw_prote)
    
    print('Total proteomic load in the list is: {:.6f} fg'.format(PL))
    proteome_fraction = (PL / total_proteo[condition]) * 100
    print('Corresponding to {:.2f}% of total proteome'.format(proteome_fraction))
    shutil.rmtree(path)
    
    
synonim = pd.read_csv('/home/notebooks/Notebooks/LIZ/Regulon12/GeneProductAllIdentifiersSet.tsv',sep="\t", comment='#')
tbl_synon = synonim[['2)geneName','6)geneSynonyms']]
tbl_synon.columns = ["gene","synonim_genes"]
syn_sel = tbl_synon.drop_duplicates()
 

def search_synonyms(genes):
    '''
    This function take any list of E. coli genes or file and find their synomym genes   '''
    
    result = []
    flag = 0

    for indice, row in syn_sel.iterrows():
        gene_name_reg12 = row['gene']
        synonym = row['synonim_genes'].split(',')

        for item in synonym:
            if item.startswith('b') and item[1:].isdigit():
                bnum = item

        if flag == 1:
            break

        for gene in genes:
            if gene == gene_name_reg12 or gene in synonym:
                result.append((gene_name_reg12, gene))
                break

    return result


def find_TFs(my_path):
    #Datasets
    essential_genes = "~/Notebooks/LIZ/Ntbooks_finales/infiles/Essential_Genes_Glucose.csv"
    proteomic_data = "~/Notebooks/LIZ/Ntbooks_finales/infiles/proteomicf_2346.csv"
    tabla_TFs =  pd.read_csv("~/Notebooks/LIZ/Ntbooks_finales/infiles/RealPB_FTs_regulon12_212ne.csv")
    network = pd.read_csv('~/Notebooks/LIZ/Ntbooks_finales/infiles/NetworkRegulatorGene.tsv', sep='\t')

    network.columns = ["id_TF","TF","tf_gname","id_Target","Target","Regulation","confidence"]
    network = network[['TF','Target','Regulation','confidence']]
    proteomic_data =  pd.read_csv(proteomic_data)
    proteomic_data.rename(columns={"gene": "Target"}, inplace=True)
    ess_gen = pd.read_csv(essential_genes)
    ess_gen.columns  = ["Target"]
    ess_gen["Essential"] = ["Yes"]*len(ess_gen["Target"])
    net_fgs = network.merge(ess_gen, on="Target", how="left").drop_duplicates().fillna("No")
    
    df = pd.read_csv(my_path)
    lista2 = list(df["genes"])
    print(len(lista2), 'genes were found in the list')
    matching_genes = net_fgs[net_fgs["Target"].isin(lista2)]
    ngenes = len(matching_genes)
    print(ngenes, 'interactions created')

    ##cuantos TF regulan mi lista
    tf = np.unique(matching_genes["TF"])
    longitud = len(tf)
    print(longitud, 'TFs regulate the total network')
    
    #cuantos son esenciales y cuantos no esenciales?
    lista3= [] 
    for g in lista2:
        if isinstance(g, str):
            res = g[0].upper() + g[1:]
        else:
            res = str(g).capitalize()
        lista3.append(res)

    TF_nes = tabla_TFs[tabla_TFs["TF"].isin(lista3)]
    nTF = len(TF_nes)
    print(nTF, 'are TFs no essential', '\n')
    
    #
    TF_f = list(TF_nes.TF)
    if TF_f:  
        print(TF_f, '\n')
        
        
def gene_category(gen, condition=None):
    all_genes = "~/Notebooks/LIZ/Ntbooks_finales/infiles/all_genes_categorized.csv"
    df = pd.read_csv(all_genes)
    if condition:
        df = df[df['gene'] == gen]
        df = df[df['condition'] == condition]  
        df2 = df[['category']]
        df2.set_index('category')
        
        if len(df2) == 0:
            print(gen, "is not in the fitness base")
            return
        else:
            print(gen, "categorized in", condition, 'is a gene', )
            
        return df2
        
    else:
        # Buscar el gen en todas las condiciones
        df = df[df['gene'] == gen]
        df2 = df[['condition','category']]
        if len(df2) == 0:
            print(gen, "is not in the fitness base")
        else:
            print(gen, "categorized in all conditions")
        return df2