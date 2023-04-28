from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt 
from  matplotlib.colors import *
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import sys
import os

LENGTH = 20
POLAR = [0, 4, 7, 8, 9, 10, 12, 13, 14, 15, 19]

def plot_grid(filename, title=None):
    if title == None:
        title = filename
    data = pd.read_csv(filename, ',')
    fig = plt.figure()
    print(type(fig))
    print(fig)
    ax = plt.axes(projection='3d')
    xData = []
    yData = []
    zData = []
    colors = []
    size = []
    for index in range(0, len(data)):
        row = data.loc[index]
        xData.append(row['x'])
        yData.append(row['y'])
        zData.append(row['z'])
        if row['hydrophilic'] == True:
            colors.append("green")
        else:
            colors.append("red")
        size.append(20)
    ax.set_title(title)
    ax.scatter3D(xData, yData, zData, '--o', s = size, c = colors)
    for index in range(0, len(data)):
        row = data.loc[index]
        ax.text(row['x'], row['y'], row['z'], row['acid'] + ' ' + str(index))
    ax.plot3D(xData, yData, zData, '--')
    colors = ["green", "red"]
    labels = ["Hydrophilic", "Hydrophobic"]
    patches = [ mpatches.Patch(color=colors[i], label=labels[i] ) for i in range(0, len(colors)) ]
    fig.legend(handles=patches)



    


def get_files(best_or_worst, map_or_grid):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    directory = dir_path + '/' + best_or_worst + '/'
    filelist = os.listdir(directory)
    ret = []
    for name in filelist:
        if map_or_grid in name:
            ret.append(dir_path + '/' + best_or_worst + '/' + name)
    return ret


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


def plot_map(filename, matrix=None, title=None):
    if title == None:
        title = filename
    if matrix == None:
        matrix = load_map(filename)
    max_value = -1
    for i in range(0, LENGTH):
        for j in range(0, LENGTH):
            if matrix[i][j] > max_value:
                max_value = matrix[i][j]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(matrix, cmap='coolwarm', vmin=0, vmax=max_value)
    ax.set_title(title)
    fig.colorbar(cax)

def plot_map_contact_type(filename, title=None):
    if title == None:
        title = filename
    matrix = load_map(filename)
    colors = ["white", "blue", "purple", "red", "black"]
    for i in range(0, LENGTH):
        for j in range(0, LENGTH):
            if j in range(i - 1, i + 2):
                matrix[i][j] = 4
            if matrix[i][j] == 1:
                if i in POLAR and j in POLAR:
                    matrix[i][j] = 1
                elif i not in POLAR and j not in POLAR:
                    matrix[i][j] = 3
                else:
                    matrix[i][j] = 2
    labels = ["No Contact", "Polar-Polar", "Hydrophobic-Polar", "Hydrophobic-Hydrophobic", "Contact Not Possible"]
    patches = [ mpatches.Patch(color=colors[i], label=labels[i] ) for i in range(0, len(colors)) ]
    cmap = LinearSegmentedColormap.from_list("mycmap", colors)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(matrix, cmap=cmap, vmin=0, vmax=len(colors) - 1)
    ax.set_title(title)
    fig.legend(handles=patches)
    #fig.colorbar(cax)


def sum_matrices(filelist, num):
    newFilelist = []
    for filename in filelist:
        split = filename.split('/')
        name = split[len(split) - 1]
        split = name.split('_')
        if int(split[0]) <= num:
            newFilelist.append(filename)
    print("sum matrices: files found: ", newFilelist)
    master_matrix = []
    for i in range(0, LENGTH):
        newRow = []
        for j in range(0, LENGTH):
            newRow.append(0)
        master_matrix.append(newRow)
    
    for filename in newFilelist:
        matrix = load_map(filename)
        for i in range(0, LENGTH):
            for j in range(0, LENGTH):
                master_matrix[i][j] += matrix[i][j]

    return master_matrix

def find_max(matrix):
    max_value = -1
    for i in range(0, LENGTH):
        for j in range(0, LENGTH):
            if matrix[i][j] > max_value:
                max_value = matrix[i][j]
    return max_value

def cell_color(index):
    if index in POLAR:
        return '\\cellcolor{red}'
    else:
        return '\\cellcolor{blue}'

ACIDS = ["ASN", "LEU", "TYR", "ILE", "GLN", "TRP", "LEU", "LYS", "ASP", "GLY", "GLY", "PRO", "SER", "SER", "GLY", "ARG", "PRO", "PRO", "PRO", "SER"]
def find_representative(filelist, num, sum_matrix, threshold=1):
    newFilelist = []
    for filename in filelist:
        split = filename.split('/')
        name = split[len(split) - 1]
        split = name.split('_')
        if int(split[0]) <= num:
            newFilelist.append(filename)
    
    contact_pairs = []
    max_value = find_max(sum_matrix)
    print("max: ", max_value)
    for i in range(0, LENGTH):
        for j in range(i, LENGTH):
            if max_value - threshold <= sum_matrix[i][j]:
                contact_pairs.append((i, j))
    print("Contact pairs found: ", contact_pairs)
    for pair in contact_pairs:
        line = cell_color(pair[0]) + ' ' + ACIDS[pair[0]] + ' (' + str(pair[0]) + ')' + ' & ' +  cell_color(pair[1]) + ' ' + ACIDS[pair[1]] + ' (' + str(pair[1]) + ') \\\\'
        print(line)
        print("\\hline")
    representatives = []
    for filename in newFilelist:
        matrix = load_map(filename)
        shouldAdd = True
        for pair in contact_pairs:
            if matrix[pair[0]][pair[1]] != 1:
                shouldAdd = False
                break
        if shouldAdd == True:
            representatives.append(filename)
    print("representatives: ", representatives)
    return (representatives[0], contact_pairs)
            
    
def plot_representation(filename, pairs, title):
    if title == None:
        title = filename
    data = pd.read_csv(filename, ',')
    fig = plt.figure()
    print(type(fig))
    print(fig)
    ax = plt.axes(projection='3d')
    xData = []
    yData = []
    zData = []
    colors = []
    size = []
    for index in range(0, len(data)):
        row = data.loc[index]
        xData.append(row['x'])
        yData.append(row['y'])
        zData.append(row['z'])
        if row['hydrophilic'] == True:
            colors.append("green")
        else:
            colors.append("red")
        size.append(20)
    ax.set_title(title)
    ax.scatter3D(xData, yData, zData, '--o', s = size, c = colors)
    for index in range(0, len(data)):
        row = data.loc[index]
        #ax.text(row['x'], row['y'], row['z'], row['acid'] + ' ' + str(index))
    ax.plot3D(xData, yData, zData, '--', alpha=.3)
    colors = ["green", "red"]
    labels = ["Hydrophilic", "Hydrophobic"]
    patches = [ mpatches.Patch(color=colors[i], label=labels[i] ) for i in range(0, len(colors)) ]

    
    for pair in pairs:
        contact_x = []
        contact_y = []
        contact_z = []
        for index in pair:
            row = data.loc[index]
            contact_x.append(row['x'])
            contact_y.append(row['y'])
            contact_z.append(row['z'])
            ax.text(row['x'], row['y'], row['z'], row['acid'] + ' ' + str(index))
        ax.plot3D(contact_x, contact_y, contact_z, '-', color='orange')


    fig.legend(handles=patches)

    



NUM = 10

best_grid_name = "best/0_best.grid"
worst_grid_name = "worst/0_worst.grid"
best_map_name = "best/0_best.map"
worst_map_name = "worst/0_worst.map"
print(get_files("worst", "grid"))
plot_grid(best_grid_name,title="Best Structure")
plot_grid(worst_grid_name, title="Worst Structure")


plot_map_contact_type(best_map_name, "Best Map")
plot_map_contact_type(worst_map_name, "Worst Map")

BEST_GRIDS = get_files("best", "grid")
BEST_MAPS = get_files("best", "map")

sum_matrix = sum_matrices(BEST_MAPS, NUM)
map_name = "Frequency of Contacts in the Top " + str(NUM) + " Structures"
plot_map(None, matrix=sum_matrix, title=map_name)

(rep, pairs) = find_representative(BEST_MAPS, NUM, sum_matrix, threshold=2)
print("Rep: ", rep)
split = rep.split('.')
rep = split[0] + '.grid'
plot_representation(rep, pairs, "Best Contacts")

plt.show()