from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

'''
Plots a contact map
'''

filename = "best/0best.map"
if (len(sys.argv) > 1):
    filename = sys.argv[1]

#Build the matrix
f = open(filename, 'r')
lines = f.readlines()
matrix = []
maxvalue = -1
for line in lines:
    line = line[0:len(line) - 1]
    split = line.split(',')
    newRow = []
    for num in split:
        if int(num) > maxvalue:
            maxvalue = int(num)
        newRow.append(int(num))
    matrix.append(newRow)
print("Loaded matrix: ")
print(matrix)


for i in range(0, len(matrix)):
    for j in range(0, len(matrix)):
        if j in range(i - 1, i + 2):
            matrix[i][j] = -10
fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(matrix, cmap='coolwarm', vmin=-1, vmax=maxvalue)
fig.colorbar(cax)
plt.show()