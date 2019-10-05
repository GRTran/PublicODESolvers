import numpy as np
import matplotlib.pyplot as plt
import re
from mpl_toolkits.mplot3d import Axes3D

def ConvertData(datafile):
    t = open(datafile,'r')
    r = t.readlines()
    out = []
    for rx in r:
        data = re.split(' ',rx)

        temp =[float(data[0]), float(data[1]), float(data[2]), float(data[3])]
        out = out + [temp]
    return out


def PlotData(input_matrix):
    fig = plt.figure()
    print(len(input_matrix))
    ax = plt.axes(projection='3d')
    for count in range(0,len(input_matrix)):
        ax.scatter(input_matrix[count][1],input_matrix[count][2], input_matrix[count][3])
    plt.show()

dataout = ConvertData('lorenz_output.txt')
print(dataout)
PlotData(dataout)
