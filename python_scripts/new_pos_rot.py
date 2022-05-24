import numpy as np
import os

testName = 'topography_long'
ng = 853
posName = '/home/hasitha/Desktop/data/fabric/'+testName+'/results/positions_'+testName+'.dat'
rotName = '/home/hasitha/Desktop/data/fabric/'+testName+'/results/rotations_'+testName+'.dat'

posSave = '/home/hasitha/Desktop/data/fabric/'+testName+'/InitState/positions_'+testName+'.dat'
rotSave = '/home/hasitha/Desktop/data/fabric/'+testName+'/InitState/rotations_'+testName+'.dat'

bkPos = '/home/hasitha/Desktop/data/fabric/'+testName+'/InitState/positions_'+testName+'.bk'
bkRot = '/home/hasitha/Desktop/data/fabric/'+testName+'/InitState/rotations_'+testName+'.bk'

os.rename(posSave, bkPos)
os.rename(rotSave, bkRot)

with open(posName, 'r') as file:
    allPositions = file.readlines()

with open(rotName, 'r') as file:
    allRotations = file.readlines()

ns = int(len(allPositions)/ng)

N = 9
positions = allPositions[(N-1)*ng:N*ng]
rotations = allRotations[(N-1)*ng:N*ng]

with open(posSave, 'w') as file:
    for pos in positions:
        file.write(pos)

with open(rotSave, 'w') as file:
    for rot in rotations:
        file.write(rot)
