import cobra
import pdb
import sys
import json
import os.path
import glob, os
sys.path.append('../')

from util import *
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import pandas as pd

# calculate molecular signature
def init_with_atom_feature(metList):
    molecular_signature = dict()
    all_features = dict()
    all_features[0] = set()
    for met in metList:
        if met in missing: continue
        if met == 'C00080': continue # this is a proton, ignore it
        if met == 'C00282': continue # this is h2, ignore it
        try:
            # molecular_signature = dict()
            mol = Chem.MolFromMolFile('../../novoStoic_data/KEGG/' + met + '.mol')
            G = mol_to_nx(mol)

            # assign feature string                             
            met_features, G = init_atomfeature(G)
            all_features[0] = all_features[0].union(met_features)
            # pdb.set_trace()
            nx.write_yaml(G, './radius_0/' + met + '.yaml')
        except Exception as e:
            print met
    df = pd.DataFrame(sorted(all_features[0]),columns=None)
    df.to_csv('moieties_0.csv',header=None,index=False)

def get_prime(metList,radius):
    for met in metList:
      if met in missing: continue
      if met == 'C00080': continue # this is just a proton, ignore it
      if met == 'C00282': continue # this is just h2, ignore it

      try:      
        # at radius = 0, moieties
        # radius = 0
        G = nx.read_yaml('./radius_' +str(radius) + '/' + met + '.yaml') 
        product_list = pd.read_csv('./moieties_' + str(radius) + '.csv',\
                                                    header=None)[0].tolist()
        # pdb.set_trace()
        if radius >= 1:
            product_list = map(int, product_list)
        G = processPrm(G,radius,product_list)
        nx.write_yaml(G,'./radius_' +str(radius) + '/' + met + '.yaml')
      except Exception as e:
          print met
          raise e

def get_product(metList,radius):
    molecular_signature = dict()
    all_features = dict()
    all_features[0] = set()
    for met in metList:
        if met in missing: continue
        if met == 'C00080': continue # this is just a proton, ignore it
        if met == 'C00282': continue # this is just h2, ignore it
        try:
            # pdb.set_trace()
            G = nx.read_yaml('./radius_' +str(radius-1) + '/' + met + '.yaml') 

            # assign feature string                             
            met_features, G = processPrd(G,radius)
            all_features[0] = all_features[0].union(met_features)
            # pdb.set_trace()
            nx.write_yaml(G,'./radius_' +str(radius) + '/' + met + '.yaml')
        except Exception as e:
            print met
            raise e
    df = pd.DataFrame(sorted(all_features[0]),columns=None)
    df.to_csv('moieties_'+str(radius) +'.csv',header=None,index=False)

# get molecular signature: count the number of moieties
def calmolecularsignature(metList,radius):
    molecular_signature = dict()
    for met in metList:
        if met in missing: continue
        if met == 'C00080': continue # this is just a proton, ignore it
        if met == 'C00282': continue # this is just h2, ignore it
        G = nx.read_yaml('./radius_' +str(radius) + '/' + met + '.yaml')
        # reactantsGraph.append(primesGraph)
        cardinality_dict = cal_cardinality(G,radius)
        # reactantsCardinality.append(cardinality_dict)
        
        molecular_signature[met] = cardinality_dict
    molsigna_df = pd.DataFrame.from_dict(molecular_signature).fillna(0)
    molsigna_df.to_csv('molecular_signature_'+str(radius) + '.csv',index=True)
    return molecular_signature

# calculate reaction rule
def get_rxn_rule(radius):
    reaction_dict = json.load(open('./optstoic_v3_reduced_Sji.json','rb'))
    molsigna_df = pd.read_csv('molecular_signature_'+ str(radius)+'.csv',index_col=0)
    
    rule_df = pd.DataFrame(index=molsigna_df.index)
    for rid,value in reaction_dict.iteritems():
        # skip the reactions with missing metabolites
        b = value.keys()
        if not set(missing).isdisjoint(b): continue

        rule_df[rid] = 0
        for met,stoic in value.iteritems():
            if met == 'C00080' or met == 'C00282': continue # hydogen is zero
            rule_df[rid] += molsigna_df[met]*stoic

    rule_df.to_csv('rule_'+ str(radius)+'.csv',index=True)

def remove_duplicate(radius):
    rule1 = pd.read_csv('rule_'+str(radius) + '.csv',index_col=0).T.drop_duplicates().T
    rule1.to_csv('rule_'+str(radius) + '_noduplic.csv',index=True)

def get_molsig_exmetab(met):
    mol = Chem.MolFromMolFile('../../novoStoic_data/KEGG/' + met + '.mol')
    filename = './exchange/' + met + '.yaml'

    G = mol_to_nx(mol)
    # assign feature string                             
    met_features,G = init_atomfeature(G)
    molecular_signature = dict()

    for radius in range(0,2):
        product_list = pd.read_csv('./moieties_' + str(radius) + '.csv',\
                                                    header=None)[0].tolist()
        # pdb.set_trace()
        processPrm(G,radius,product_list)
        met_features, G = processPrd(G,radius+1)
    nx.write_yaml(G,filename)

    for radius in range(0,2):
        cardinality_dict = cal_cardinality(G,radius)
        molecular_signature[met] = cardinality_dict
        molsigna_df = pd.DataFrame.from_dict(molecular_signature).fillna(0)
        molsigna_df.to_csv('./exchange/molsig_'+ met + '_' + str(radius) +\
                             '.csv',index=True)

if __name__ == '__main__':
    os.chdir("../../novoStoic_data/KEGG/") # list of mol file for metabolites
                                           # from database (e.g. KEGG)
    metList = [file.replace('.mol','') for file in glob.glob("*.mol")]

    ### reprime

    #Step 1: Identification of moieties for all metabolites in the MetRxn database
    init_with_atom_feature(metList)
    get_prime(metList,0)
    for radius in range(0,3):
        get_product(metList,radius)
        get_prime(metList,radius)
    
    for radius in range(0,3):
    # Step 2: Determination of the molecular signature of each metabolite.
        calmolecularsignature(metList,radius)

    # Step 3: Inference of the associated reaction rule for each reaction in MetRxn.
        get_rxn_rule(radius)
        remove_duplicate(radius)





