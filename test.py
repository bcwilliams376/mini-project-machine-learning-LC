import tmap as tm
import numpy as np
from map4 import MAP4Calculator
from rdkit import Chem
import csv
from sklearn.model_selection import train_test_split, GridSearchCV

def bring():
    smiles = []
    file = open("SMILES.txt") 
    for line in file.readlines(): #reads each SMILES string from txt file and adds to the array
        smiles.append(line)
    file.close()
    return smiles

def send(list, file):
    wtr = csv.writer(open (file, 'w'), delimiter=';', lineterminator='\n')
    for x in list : 
        wtr.writerow ([x])
    return "done"


smiles_list = bring() 

dim = 1024

MAP4 = MAP4Calculator(dimensions=dim)
ENC = tm.Minhash(dim)

fingerprint_list = []

for item in smiles_list:
    smiles_a = item
    mol_a = Chem.MolFromSmiles(smiles_a)
    map4_a = MAP4.calculate(mol_a)
    fingerprint_list.append(map4_a)

distances = []
fingerprint_list_length = len(fingerprint_list)

for x in range(fingerprint_list_length-1): 
    map1 = fingerprint_list[x]
    map2 = fingerprint_list[x+1]
    distances.append(ENC.get_distance(map1, map2))

print(distances)



exp_kernel = []
for distance in distances:
    k = np.exp(-(distance)**2/1**2)
    exp_kernel.append(k)

print(send(exp_kernel, "Exponential_kernel.csv"))

linear_kernel = []
for distance in distances:
    k = 1 - (distance)**2
    linear_kernel.append(k)

print(send(linear_kernel, "Linear_kernel.csv"))
    


   


