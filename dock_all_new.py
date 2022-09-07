import os
import math
import subprocess
import numpy as np
import pandas as pd
import pubchempy as pcp
from Bio import PDB
from Bio.PDB import *
from multiprocessing import Pool

global compound_list
global pdb_list
compound_list=np.loadtxt('target_ligand.csv',dtype=str,encoding="utf-8", delimiter='\n')
#pdb_list=np.loadtxt('pdb.csv',dtype=str,encoding="utf-8", delimiter='\n')
pdb_list=['1s78']

def docking(receptorligand):
    receptor = receptorligand[0]
    ligand = receptorligand[1].replace(' ','_')
    energy = 0.0
    if os.path.isfile(str(ligand+receptor)+'.log')==False:
        parser = PDB.PDBParser()
        struct = parser.get_structure(receptor,'pdb'+str(receptor)+'.pdbqt')
        atom_cord = []
        for atom in struct.get_atoms():
            atom_cord.append(atom.get_coord())
        atom_cord = np.asarray(atom_cord)
        size_x, size_y, size_z = atom_cord.max(axis=0)-atom_cord.mean(axis=0)
        size_x, size_y, size_z = 2*math.ceil(size_x), 2*math.ceil(size_y), 2*math.ceil(size_z)
        if size_x > 126:
            size_x = 126
        if size_y > 126:
            size_y = 126
        if size_z > 126:
            size_z = 126
        center_x, center_y, center_z = atom_cord.mean(axis=0)
        cmd = 'vina --receptor pdb'+str(receptor)+'.pdbqt --center_x '+str(center_x)+' --center_y '+str(center_y)+' --center_z '+str(center_z)+' --size_x '+str(size_x)+' --size_y '+str(size_y)+' --size_z '+str(size_z)+' --out '+str(ligand+receptor)+'_out.pdbqt --log '+str(ligand+receptor)+'.log --ligand '+str(ligand)+'.pdbqt'
        subprocess.run(cmd)
        try:
            energy = open(str(ligand+receptor)+'.log',mode='r').read().split('   1         ')[1].split('      0.000      0.000')[0]
        except:
            pass

        try:
            energy = open(str(ligand+receptor)+'.log',mode='r').read().split('   1        ')[1].split('      0.000      0.000')[0]
        except:
            pass
    else:
        try:
            energy = open(str(ligand+receptor)+'.log',mode='r').read().split('   1         ')[1].split('      0.000      0.000')[0]
        except:
            pass

        try:
            energy = open(str(ligand+receptor)+'.log',mode='r').read().split('   1        ')[1].split('      0.000      0.000')[0]
        except:
            pass

    return (ligand, receptor), float(energy)



if __name__ == '__main__':

    for b in compound_list:
        if os.path.isfile(b+'.pdb')==False:
                print(b+'not exist')
        if os.path.isfile(b+'.pdb')==True:
            print(b+'.pdb already exists...')
    for a in pdb_list:
        if os.path.isfile('pdb'+a+'.pdb')==False:
            print(a+' does not exists')
        elif os.path.isfile('pdb'+a+'.pdb')==True:
            print(a+' already exists')

    for b in compound_list:
        cmd=r'"C:\Program Files (x86)\MGLTools-1.5.6\python.exe" "C:\Program Files (x86)\MGLTools-1.5.6\Lib\site-packages\AutoDockTools\Utilities24\prepare_ligand4.py" -l "'+str(b.replace(' ','_'))+'.pdb" -v'
        subprocess.run(cmd)
    for p in pdb_list:
        cmd=r'"C:\Program Files (x86)\MGLTools-1.5.6\python.exe" "C:\Program Files (x86)\MGLTools-1.5.6\Lib\site-packages\AutoDockTools\Utilities24\prepare_receptor4.py" -r "pdb'+str(p)+'.pdb"'
        subprocess.run(cmd)
    
    df = pd.DataFrame(index=compound_list,columns=pdb_list)
    

    with Pool(1) as p:
        a = p.map(docking, [(x, y) for x in pdb_list for y in compound_list])
    for b in a:
        df.loc[[b[0][0]],[b[0][1]]] = b[1]
    df.to_excel('result.xlsx')
