import numpy as np
import os

print("hello world")

ng = 127900
testName = '7-6B-2'
deleteList = []
filePath = '/home/hasitha/Desktop/data/fabric/'+testName
#filePath = '/global/scratch/users/tanpeng/LSDEM_Modelling/Input/'+testName
maxMass = 0.0
maxId   = -1
"""
for i in np.arange(1,ng+1):
    if i % 100 == 0: print('have read {} files'.format(i))
    with open(filePath+'/Morphologies/'+'morph_'+str(i)+'.dat', 'r') as file:
        mass = float(file.readline())
        if mass > maxMass:
            maxMass = mass
            maxId   = i
        file.readline()
        file.readline()
        file.readline()
        mark = False
        nodes = file.readline().split()
        for j in np.arange(len(nodes)):
            if len(nodes[j]) > 10:
                mark = True
                break
        #if 'nan' in nodes or '-nan' in nodes or mark or mass < 400:
        if mass < 1000:
            deleteList.append(i)
"""

with open(filePath+'/InitState/rotations_'+testName+'.dat', 'r') as file:
    for i in np.arange(1,ng+1):
        if i % 1000 == 0: print('processed {} grains'.format(i))
        rot = file.readline().split()
        flag = True
        for j in rot:
            if float(j) != 0:
                flag = False
                break
        if 'nan' in rot or '-nan' in rot or 'inf' in rot or '-inf' in rot or flag:
            deleteList.append(i)

#pos = np.loadtxt(filePath+'/InitState/positions_'+testName+'.dat')
deleteList = sorted(deleteList)[::-1]
print(deleteList)
print(len(deleteList))
#print('max mass is {}, max ID is {}'.format(maxMass, maxId))

"""
posPath  = filePath+'/InitState/positions_'+testName+'.dat'
rotPath  = filePath+'/InitState/rotations_'+testName+'.dat'
bk_pos   = filePath+'/InitState/positions_'+testName+'.backup'
bk_rot   = filePath+'/InitState/rotations_'+testName+'.backup'

#backup
for i in np.arange(len(deleteList)):
    os.rename(filePath+'/Morphologies/morph_'+str(ng-i)+'.dat', \
              filePath+'/Morphologies/morph_'+str(deleteList[i])+'.dat')


for i in np.arange(len(deleteList)):
    os.rename(filePath+'/Polyhedrons/poly_'+str(ng-i)+'.dat', \
              filePath+'/Polyhedrons/poly_'+str(deleteList[i])+'.dat')


os.rename(posPath, bk_pos)
os.rename(rotPath, bk_rot)
#readfiles
with open(bk_pos, 'r') as file:
    positions = file.readlines()
    for idx in deleteList:
        positions[idx-1] = positions.pop()

with open(bk_rot, 'r') as file:
    rotations = file.readlines()
    for idx in deleteList:
        rotations[idx-1] = rotations.pop()
#savefile
with open(posPath, "w") as file:
    for position in positions:
        file.write(position)

with open(rotPath, "w") as file:
    for rotation in rotations:
        file.write(rotation)
"""
