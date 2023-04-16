from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os


def load_map(filename):
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
    f.close()
    return matrix
    


def save_matrix(matrix, filename):
    f = open(filename, 'w')
    for i in range(0, 20):
        line = ''
        for j in range(0, 20):
            if j == 0:
                line += str(matrix[i][j])
            else:
                line += ',' + str(matrix[i][j])
        #print(line)
        f.write(line + '\n')

    f.close()

Type = "best"
Lim = 10000
Output_name = "res.sum.map"


for i in range(1, len(sys.argv)):
    arg = sys.argv[i]
    split = arg.split(':')
    if split[0] == '-t':
        Type = split[1]
    elif split[0] == '-o':
        Output_name = split[1]
    elif split[0] == '-l':
        Lim = int(split[1])

print("Type: ", Type, "Lim: ", Lim, "output filename: ", Output_name)
dir_path = os.path.dirname(os.path.realpath(__file__))
directory = dir_path + '/' + Type + '/'
print("Looking in directory", directory)
dirlist = os.listdir(directory)
#print(dirlist)

filelist = []
for name in dirlist:
    if 'map' in name:
        num = int(name.split('_')[0])
        if num < Lim:
            filelist.append(name)

print("Files found: ", filelist)
#Build the result matrix
resultMatrix = []
for i in range(0, 20):
    row = []
    for j in range(0, 20):
        row.append(0)
    resultMatrix.append(row)

for name in filelist:
    matrix = load_map(directory + name)
    for i in range(0, 20):
        for j in range(0, 20):
            resultMatrix[i][j] += matrix[i][j]
save_matrix(resultMatrix, Output_name)

