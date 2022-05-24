import numpy as np

testName = 'lowres_cylinder'
ng = 2594 #19330
fileName = '/home/hasitha/Desktop/data/fabric/'+testName+'/results/positions_'+testName+'.dat'

print(fileName)

height = []
with open(fileName, 'r') as file:
    allPositions = file.readlines()

ns = int(len(allPositions)/ng)
for i in np.arange(ns):
    positions = allPositions[i*ng:(i+1)*ng]
    maxHeight = 0.0
    for position in positions:
        pos = float(position.split()[2])
        if pos > maxHeight: maxHeight = pos
    height.append(maxHeight)

print(height)
print(len(height))
