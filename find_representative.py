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

sumMap = "res.sum.map"
threshold_modifier = 2
Type = "best"
outfile = "r_coords.out"

for i in range(1, len(sys.argv)):
    arg = sys.argv[i]
    split = arg.split(':')
    if split[0] == '-t':
        Type = split[1]
    elif split[0] == '-s':
        sumMap = split[1]
    elif split[0] == '-th':
        threshold_modifier = int(split[1])
    elif split[0] == '-o':
        outfile = split[1]

print("Source: ", sumMap, "threshold: ", threshold_modifier, "Type: ", Type, "outfile: ", outfile)

summed = load_map(sumMap)
print(summed)


max_value = -1
for i in range(0, 20):
    for j in range(0, 20):
        if summed[i][j] > max_value:
            max_value = summed[i][j]
print("Max found: ", max_value)
threshold = max_value - threshold_modifier
print("So will consider all entries with value", threshold, "or higher" )

contact_list = []
for i in range(0, 20):
    for j in range(i, 20):
        if threshold <= summed[i][j]:
            contact_list.append((i, j))
print(contact_list)


dir_path = os.path.dirname(os.path.realpath(__file__))
directory = dir_path + '/' + Type + '/'
print("Looking in directory", directory)
dirlist = os.listdir(directory)
#print(dirlist)

filelist = []
for name in dirlist:
    if 'map' in name:
        filelist.append(name)

representatives = []
for name in filelist:
    map = load_map(directory + name)
    hasAll = True
    for pair in contact_list:
        if map[pair[0]][pair[1]] != 1:
            hasAll = False
            break
    if hasAll == True:
        representatives.append(name)

print("Representatives found: ", representatives)

f = open(outfile, "w")
for i in range(0, len(contact_list)):
    f.write(str(contact_list[i][0]) + ',' + str(contact_list[i][1]) + '\n')
f.close()

        

