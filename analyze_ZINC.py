from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import MolDrawing
from rdkit.Chem.Draw import IPythonConsole
import numpy as np
import pickle

def read_file(file_name):
  smiles_list = []
  with open(file_name,'r') as file:
    for smiles in file:
      smiles_list.append(smiles)

  return smiles_list

def get_probs(smarts_list,smiles_list,ring=False):
  import collections
  bonds = []
  probs = collections.OrderedDict()

  for smarts in smarts_list:
    probs[smarts] = 0

  number_of_molecules = 0
  tot = 0
  for smiles in smiles_list:
    #print smiles
    number_of_molecules += 1
    mol = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol)
    for smarts in smarts_list:
      matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts),uniquify=ring)
      num_bonds = len(matches) 
      probs[smarts] += num_bonds
      tot += num_bonds
   
  tot = 0
  probs2 = collections.OrderedDict()
  for key in probs:
    if probs[key] > 0:
      #print key, probs[key]
      tot += probs[key]
      probs2[key] = probs[key]
            
  return tot, probs2

def get_p(probs):
  p = []
  for key in probs:
    p.append(float(probs[key])/tot)
    
  return p

def get_rxn_smarts_make_rings(probs):
  X = {'[#6R': 'X4', '[#7R': 'X3'}
  rxn_smarts = []
  for key in probs:
    tokens = key.split(']')

    smarts = ''
    if '=' in key:
      smarts += tokens[0][:-1] + X[tokens[0]] + ';!R:1]'
    else:
      smarts += tokens[0][:-1] + ';!R:1]=,'
    smarts += tokens[2][:-1] + ';!R:2]>>'
    smarts += '[*:1]1' + tokens[1] + '][*:2]1'

    rxn_smarts.append(smarts)
    
  return rxn_smarts

def get_rxn_smarts_rings(probs):
  X = {'[#6R': 'X4', '[#7R': 'X3'}
  rxn_smarts = []
  for key in probs:
    tokens = key.split(']')

    smarts = ''
    if '=' in key:
      smarts += tokens[0] + X[tokens[0]] + ';!r6;!r7;!R2:1]'
    else:
      smarts += tokens[0] + ';!r6;!r7;!R2:1]'

    smarts += tokens[2] + ';!r6;!r7:2]>>'
    smarts += '[*:1]' + tokens[1] + '][*:2]'

    rxn_smarts.append(smarts)
    
  return rxn_smarts

def get_rxn_smarts(probs):
  rxn_smarts = []
  for key in probs:
    smarts = ''
    tokens = key.split(']')
    smarts = tokens[0]
    if '-' in key and '#16' not in key:
      smarts += ';!H0:1]>>[*:1]'
    if '=' in key and '#16' not in key:
      smarts += ';!H1;!H0:1]>>[*:1]'
    if ']#[' in key:
      smarts += ';H3:1]>>[*:1]'
    if '#16' in key:
      smarts += ':1]>>[*:1]'
      
    smarts += tokens[-2] + ']'
    rxn_smarts.append(smarts)
    
  return rxn_smarts

file_name = '1000.smi'

elements = ['#5','#6','#7','#8','#9','#14','#15','#16','#17','#35','#53']
bonds = ['-','=','#']

smiles_list = read_file(file_name)

smarts = []
for element in elements:
  smarts.append('['+element+']')

print get_probs(smarts,smiles_list)

smarts = []
for element in elements:
  smarts.append('['+element+'R]')

tot_Ratoms,probs_Ratoms = get_probs(smarts,smiles_list)
print tot_Ratoms,probs_Ratoms

R_elements = []
for key in probs_Ratoms:
  R_elements.append(key)
  
print R_elements

smarts = []

for i,e1 in enumerate(R_elements):
  for e2 in R_elements:
    for j,e3 in enumerate(R_elements):
      if j >= i:
        sm_s = e1 + '-' + e2 + '-' + e3
        if sm_s not in smarts:
          smarts.append(sm_s)
      sm_d = e1 + '=' + e2 + '-' + e3
      if sm_d not in smarts:
        smarts.append(sm_d)
        
print len(smarts),smarts

tot,probs = get_probs(smarts,smiles_list,ring=True)

print tot,probs

import operator
sorted_x = sorted(probs.items(), key=operator.itemgetter(1), reverse=True)
count = 0
for i in range(len(sorted_x)):
  print sorted_x[i][0],sorted_x[i][1]
  if '=' in sorted_x[i][0]: count += sorted_x[i][1]

print count

rxn_smarts_rings = get_rxn_smarts_rings(probs)
print rxn_smarts_rings

rxn_smarts_make_rings = get_rxn_smarts_make_rings(probs)
print rxn_smarts_make_rings

p_rings = get_p(probs)
print p_rings

pickle.dump(p_rings,open('p_ring.p','wb')) 
pickle.dump(rxn_smarts_rings,open('rs_ring.p','wb'))
pickle.dump(rxn_smarts_make_rings,open('rs_make_ring.p','wb'))

smarts = []

for bond in bonds:
  for element1 in elements:
    for element2 in elements:
      smarts.append('['+element1+']'+bond+'['+element2+';!R]')

print len(smarts)
tot,probs = get_probs(smarts,smiles_list)
print tot, probs
p = get_p(probs)
print p
pickle.dump(p,open('p1.p','wb')) 
rxn_smarts = get_rxn_smarts(probs)
print rxn_smarts
pickle.dump(rxn_smarts,open('r_s1.p','wb')) 

