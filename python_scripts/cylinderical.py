import numpy as np
import os
import shutil
test_name = 'cylinder_belgium'
data_path = '/home/hasitha/Desktop/data/fabric/'+test_name
old_path  = '/home/hasitha/Desktop/data/fabric/full_belgium'
index     = np.loadtxt('/home/hasitha/Desktop/data/fabric/'+test_name+'/qualified_id.txt')

print(len(index))

for i in np.arange(len(index)):
    if i % 100 == 0: print("copy Morphologies: ", i)
    idx = int(index[i])
    shutil.copy2(old_path+'/Morphologies/morph_'+str(idx+1)+'.dat', \
                 data_path+'/Morphologies/morph_'+str(i+1)+'.dat')
"""
for i in np.arange(len(index)):
    if i % 100 == 0: print("copy Polyhedrons: ", i)
    idx = int(index[i])
    shutil.copy2(old_path+'/Polyhedrons/poly_'+str(idx+1)+'.dat', \
                 data_path+'/Polyhedrons/poly_'+str(i+1)+'.dat')
"""
with open(old_path+'/InitState/positions_full_belgium.dat', 'r') as file:
    positions = file.readlines()
with open(old_path+'/InitState/rotations_full_belgium.dat', 'r') as file:
    rotations = file.readlines()

positions_file  = data_path+'/InitState/positions_'+test_name+'.dat'
rotations_file  = data_path+'/InitState/rotations_'+test_name+'.dat'

with open(positions_file, 'w') as file:
    for i in np.arange(len(index)):
        file.write(positions[int(index[i])])

with open(rotations_file, 'w') as file:
    for i in np.arange(len(index)):
        file.write(rotations[int(index[i])])
