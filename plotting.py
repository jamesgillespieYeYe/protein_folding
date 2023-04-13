from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


data = pd.read_csv("grid.out.txt", ',')
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

ax.scatter3D(xData, yData, zData, '--o', s = size, c = colors)
for index in range(0, len(data)):
    row = data.loc[index]
    ax.text(row['x'], row['y'], row['z'], row['acid'] + ' ' + str(index))
ax.plot3D(xData, yData, zData, '--')




plt.show()