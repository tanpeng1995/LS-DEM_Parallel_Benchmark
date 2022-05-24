import numpy as np
import os

print("hello world")

ng = 127800
testName = '7-6B_full'
filePath = '/home/hasitha/Desktop/data/fabric/'+testName
maxMass = 0.0
maxId   = -1
deleteList = []

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
        if mass < 400:
            deleteList.append(i)

deleteList = deleteList[::-1]
print(deleteList)
print(len(deleteList))
print('max mass is {}, max ID is {}'.format(maxMass, maxId))

"""
pos = np.loadtxt(filePath+'/results/positions_'+testName+'.dat')
last = pos[-ng:]
#dim 0
sublist1 = list(np.where(np.abs(last[:,0]) > 450.)[0])
print(len(sublist1))
#dim 1
sublist2 = list(np.where(np.abs(last[:,1]) > 450.)[0])
print(len(sublist2))
#dim 2
sublist3 = list(np.where(np.abs(last[:,2]) > 1330.)[0])
print(len(sublist3))

deleteList = sublist1 + sublist2 + sublist3
deleteList = np.array(deleteList) + 1
deleteList = list(set(deleteList))
deleteList = sorted(deleteList)
deleteList = deleteList[::-1]
#print(deleteList)
print(len(deleteList))
"""


posPath  = filePath+'/InitState/positions_'+testName+'.dat'
rotPath  = filePath+'/InitState/rotations_'+testName+'.dat'
bk_pos   = filePath+'/InitState/positions_'+testName+'.backup'
bk_rot   = filePath+'/InitState/rotations_'+testName+'.backup'

#backup

for i in np.arange(len(deleteList)):
    os.rename(filePath+'/Morphologies/morph_'+str(ng-i)+'.dat',
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
