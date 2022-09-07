import os
import math
import time
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


pdb_list=np.loadtxt('pdb.csv',dtype=str,encoding="utf-8", delimiter='\n')


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

    os.system('vina_split --input '+str(ligand+receptor)+'_out.pdbqt')
    os.system('obabel '+str(ligand+receptor)+'_out_ligand_1.pdbqt -opdb -m')
    time.sleep(2)
    os.system('pdb_merge pdb'+str(receptor)+'.pdb '+str(ligand+receptor)+r'_out_ligand_1.pdb>merged_'+str(ligand+receptor)+'.pdb')
    print('pdb_merge pdb'+str(receptor)+'.pdb '+str(ligand+receptor)+r'_out_ligand_1.pdb>merged_'+str(ligand+receptor)+'.pdb')
    time.sleep(2)
    os.system('pdb_tidy merged_'+str(ligand+receptor)+r'.pdb>tidy_merged_'+str(ligand+receptor)+'.pdb')
    time.sleep(2)
    os.system('plip -f tidy_merged_'+str(ligand+receptor)+'.pdb -ypt -o tidy_merged_'+str(ligand+receptor)+' --name result_'+str(ligand+receptor))
    os.system('plip -f tidy_merged_'+str(ligand+receptor)+'.pdb -ypt --dnareceptor -o tidy_merged_'+str(ligand+receptor)+' --name result_'+str(ligand+receptor))
    
    return (ligand.replace('_',' '), receptor), float(energy)



if __name__ == '__main__':

    for b in compound_list:
        if os.path.isfile(b.replace(' ','_')+'.pdb')==False:
                print(b.replace(' ','_')+'.sdf downloaded...')
        elif os.path.isfile(b.replace(' ','_')+'.pdb')==True:
            print(b.replace(' ','_')+'.sdf already exists...')
    for a in pdb_list:
        if os.path.isfile('pdb'+a+'.pdb')==False:
            print(a+' does not exists')
        elif os.path.isfile('pdb'+a+'.pdb')==True:
            print(a+' already exists')
    specific_list=[('4ey7','compound1 1'),
                ('4zvm','compound1 5'),
                ('3qx3','compound1 8'),
                ('4ey7','compound2 2'),
                ('4ey7','compound2 7'),
                ('3qx3','compound2 3'),
                ('3qx3','compound2 9'),
                ('4zvm','compound2 17'),
                ('4zvm','Doxorubicin'),
                ('3qx3','Quizartinib'),
                ('4ey7','Tacrine')
                   ]
    with Pool(1) as p:
        #a = p.map(docking, [(x, y) for x in pdb_list for y in compound_list])
        a = p.map(docking, specific_list)
    
