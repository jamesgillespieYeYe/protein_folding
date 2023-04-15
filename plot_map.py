from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

filename = "best/0best.map"
if (len(sys.argv) > 1):
    filename = sys.argv[1]


print("Reading grid from", filename)
data = pd.read_csv(filename, ',')
print(data)

matrix = []
for i in range(0, len(data)):
    newRow = []
    row = data.loc[i]
    for j in range(1, len(row)):
        newRow.append(row[j])
    matrix.append(newRow)
print(matrix)

for i in range(0, len(matrix)):
    for j in range(0, len(matrix)):
        if j in range(i - 1, i + 2):
            matrix[i][j] = -1
fig = plt.figure()
ax = fig.add_subplot(111)
ax.matshow(matrix, cmap='coolwarm', vmin=-1, vmax=1)
ticks = np.arange(0,len(data.columns) - 1,1)
ax.set_xticks(ticks)
ax.set_xticklabels(data.columns[1:len(data.columns)])
ax.set_yticks(ticks)
ax.set_yticklabels(data.columns[1:len(data.columns)])
plt.show()