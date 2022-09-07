import pubchempy as pcp
import numpy as np
import os

compound_list=np.loadtxt('target_ligand.csv',dtype=str,encoding="utf-8", delimiter='\n')
print(compound_list)
print(len(compound_list))
for b in compound_list:
    
    if os.path.isfile(b.replace(' ','_')+'.sdf')==False:
        print(b)
        print(pcp.get_compounds(b, 'name', record_type='2d'))
        '''
        for compound in pcp.get_compounds(b, 'name', record_type='3d'):
            print(compound)
            pcp.download('SDF', b.replace(' ','_')+'.sdf', compound.cid, 'cid', overwrite=True, record_type='3d')
            print(b.replace(' ','_')+'.sdf downloaded...')
    elif os.path.isfile(b.replace(' ','_')+'.sdf')==True:
        print(b.replace(' ','_')+'.sdf already exists')'''