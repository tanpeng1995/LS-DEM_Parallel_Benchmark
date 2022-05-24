import numpy as np
import os
from shutil import copyfile

testName = 'smallSample_100'
filePath = '/home/hasitha/Desktop/data/fabric/' + testName
posPath  = filePath + '/InitState/positions_' + testName + '.dat'
rotPath  = filePath + '/InitState/rotations_' + testName + '.dat'
morPath  = filePath + '/Morphologies/'

saveName = 'smallSample_400'
N = 4

# do positions
positions = []
with open('/home/hasitha/Desktop/data/fabric/' + testName + '/InitState/positions_' + testName + '.dat', 'r') as file:
    temps = file.readlines()
    for temp in temps:
        temp = temp.split()
        positions.append([float(temp[0]), float(temp[1]), float(temp[2])])


with open('/home/hasitha/Desktop/data/fabric/' + saveName + '/InitState/positions_' + saveName + '.dat', 'w') as file:
    for i in np.arange(N):
        for j in np.arange(N):
            for k in np.arange(N):
                for pos in positions:
                    file.write('{:.3f} {:.3f} {:.3f}\n'.format(pos[0]+i*100,pos[1]+j*100,pos[2]+k*100))

# do rotations
rotations = []
with open('/home/hasitha/Desktop/data/fabric/' + testName + '/InitState/rotations_' + testName + '.dat', 'r') as file:
    rotations = file.readlines()

with open('/home/hasitha/Desktop/data/fabric/' + saveName + '/InitState/rotations_' + saveName + '.dat', 'w') as file:
    for i in np.arange(N):
        for j in np.arange(N):
            for k in np.arange(N):
                for rot in rotations:
                    file.write(rot)

# do morph files and poly files
for i in np.arange(N):
    for j in np.arange(N):
        for k in np.arange(N):
            for ii in np.arange(1,75):
                src = '/home/hasitha/Desktop/data/fabric/' + testName + '/Morphologies/morph_' + str(ii) + '.dat'
                dst = '/home/hasitha/Desktop/data/fabric/' + saveName + '/Morphologies/morph_' + str(ii+(N**2*i+N*j+k)*74) + '.dat'
                copyfile(src, dst)

# do morph files and poly files
for i in np.arange(N):
    for j in np.arange(N):
        for k in np.arange(N):
            for ii in np.arange(1,75):
                src = '/home/hasitha/Desktop/data/fabric/' + testName + '/Polyhedrons/poly_' + str(ii) + '.dat'
                dst = '/home/hasitha/Desktop/data/fabric/' + saveName + '/Polyhedrons/poly_' + str(ii+(N**2*i+N*j+k)*74) + '.dat'
                copyfile(src, dst)
