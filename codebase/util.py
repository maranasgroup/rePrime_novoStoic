##########################################
##########################################
#########                 ################
######### rePrime util ###################
#########                 ################
##########################################
##########################################
import networkx as nx
import random
from sortedcontainers import SortedSet
from collections import defaultdict, Counter
import cPickle as pickle
import os 
import matplotlib.pyplot as plt
import csv
from rdkit.Chem import AllChem, rdmolops
#SanitizeMol, RemoveHs, Kekulize
import pdb

# the level of reaction rules
maxRadius = 3

# main function for rePrime
def cal_cardinality(G,radius):
    maxRadius = 3
    cardinality_dict = dict()
    
    prmId = "prm"+str(radius)
    prmDict = nx.get_node_attributes(G, prmId)

    for i,prime in prmDict.iteritems():
        if prime not in cardinality_dict.keys():
            cardinality_dict[prime] = 1
            # cardinality_dict[prmId][prime][prmId] = 1
        else:
            cardinality_dict[prime] = cardinality_dict[prime] + 1
    # print cardinality_dict
    return cardinality_dict

def cal_rxn_rule(reactantsCardinality,productsCardinality,maxRadius):
    rxn_rules = []
    for radius in range(maxRadius):
        for cardi in reactantsCardinality:
            rxn_rules
    return

def loadPrimes():
    f = open("primes.txt")
    ret = []    
    for line in f:
        a = line.split()
        for i in a:
            ret.append(int(i))
    return ret        

def init(G):
    nodeCnt = len(G.nodes())
    
    # test using atomic number as atom features instead of the feature string
    # C N O P S  to 2, 3, 5, 7, 11
    all_features = set()       
    for i in range(nodeCnt):
        value = G.node[i]['atomic_num']
        G.node[i]['prd0'] = value
        all_features.add(value)
    # pdb.set_trace()
    # save to pickle file
    # p_Lkup[0] = sorted(all_features)
    return all_features,G

def init_atomfeature(G):
    nodeCnt = len(G.nodes())
    
    # test using atomic number as atom features instead of the feature string
    # degree=atom.GetDegree(),
    # explicit_valence=atom.GetExplicitValence(),  
    # atomic_num=atom.GetAtomicNum(),
    # num_implicit_hs=atom.GetNumExplicitHs()
    all_features = set()       
    for i in range(nodeCnt):
        degree = G.node[i]['degree']
        explicit_valence = G.node[i]['explicit_valence']
        atomic_num = G.node[i]['atomic_num']
        num_implicit_hs = G.node[i]['num_implicit_hs']
        if degree < 10:
            value = str(degree) + '-' + str(explicit_valence) + '-0' \
                    + str(atomic_num) + '-' + str(num_implicit_hs)
        else:
             value = str(degree) + '-' + 'explicit_valence' + '-' \
                    + str(atomic_num) + '-' + str(num_implicit_hs)           
        G.node[i]['prd0'] = value
        all_features.add(value)
    # pdb.set_trace()
    # save to pickle file
    # p_Lkup[0] = sorted(all_features)
    return all_features,G

def processPrm(G,radius,product_list): 
    primes = loadPrimes()   

    prdId = "prd"+str(radius)
    prmId = "prm"+str(radius)

    # assign prime number based on lexco grphic oder
    prdDict = nx.get_node_attributes(G, prdId)

    prdPrm = {v : primes[product_list.index(prdDict.get(v))] \
                for v in prdDict}       
    nx.set_node_attributes(G, prmId, prdPrm)    
    return G
        
def processPrd(G,radius):    
    prm_cur = "prm"+str(radius-1)
    prdId_nxt = "prd"+str(radius)
    prevPrm = nx.get_node_attributes(G, prm_cur)

    # claculate prime product
    prdList = {v : reduce(lambda x,y:x*y, map(lambda x: prevPrm[x], \
        nx.all_neighbors(G, v)),1)*prevPrm[v]**2 for v in G.nodes()}
    nx.set_node_attributes(G, prdId_nxt, prdList)

    all_features = set(prdList.values())
    return all_features,G

# def prdLkup(rw = "ab",prdLkup = defaultdict(dict)):
#     if rw == "w":
#         pickle.dump(prdLkup, open("prd.lkup","wb"))        
#     else:
#         prdLkup = pickle.load(open("prd.lkup","rb"))
        
#     return prdLkup


# create a graph based on the rxns after added the pseudo nodes
# this function is from https://github.com/dakoner/keras-molecules
def mol_to_nx(mol):
    G = nx.Graph()

    # pdb.set_trace()
    # print len(mol.GetAtoms())
    for atom in mol.GetAtoms():
        # print G.nodes()
        # if atom.GetIdx() == 78:
        #     pdb.set_trace()
        G.add_node(atom.GetIdx(),
                   degree=atom.GetDegree(),
                   explicit_valence=atom.GetExplicitValence(),  
                   atomic_num=atom.GetAtomicNum(),
                   num_implicit_hs=atom.GetNumImplicitHs(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   is_aromatic=atom.GetIsAromatic())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
    return G