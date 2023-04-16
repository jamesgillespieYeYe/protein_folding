from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

filename = "grid.out.txt"
r_mode = False
r_coords_file = None
if (len(sys.argv) > 1):
    filename = sys.argv[1]
if (len(sys.argv) > 2):
    if sys.argv[2] == "r_mode":
        r_mode = 1
        r_coords_file = sys.argv[3]





print("Reading grid from", filename)
print("rmode: ", r_mode, "r_coords_file", r_coords_file)
data = pd.read_csv(filename, ',')
print(data)


fig = plt.figure()
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
ax.set_title(filename)
ax.scatter3D(xData, yData, zData, '--o', s = size, c = colors)
for index in range(0, len(data)):
    row = data.loc[index]
    ax.text(row['x'], row['y'], row['z'], row['acid'] + ' ' + str(index))
ax.plot3D(xData, yData, zData, '--')




plt.show()