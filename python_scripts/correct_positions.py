import numpy as np
import os

testName = 'parallel_benchmark_3'
deleteList = []
filePath = '/home/hasitha/Desktop/data/fabric/'+testName+'/InitState/positions_'+testName+'.dat'
bkPath   = '/home/hasitha/Desktop/data/fabric/'+testName+'/InitState/positions_'+testName+'.bk'

os.rename(filePath, bkPath)

radius = 75.
new_pos= []
with open(bkPath, 'r') as file:
    positions = file.readlines()
    for position in positions:
        pos = position.split()
        pos[0] = float(pos[0])+radius
        pos[1] = float(pos[1])+radius
        pos[2] = float(pos[2])-50.
        pos = '{:.3f} {:.3f} {:.3f}\n'.format(pos[0],pos[1],pos[2])
        new_pos.append(pos)

with open(filePath, 'w') as file:
    for pos in new_pos:
        file.write(pos)
