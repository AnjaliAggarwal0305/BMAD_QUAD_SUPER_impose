from ROOT import *
from math import *
from sys import stdin
import numpy as np 
import os
import random
import matplotlib.pyplot as plt

def readData(path):
    
    data = [[],[],[]]
    with open(path) as f:
        for _ in xrange(9):
            next(f)
        for line in f:
            line = line.split()
            if line[1] != "WF_B" and line[1] != "WF_E":
                
                data[0].append(line[5])
                data[1].append(line[7])
                data[2].append(line[3])
            
        
        x = np.asarray(map(lambda l: float(l),data[0]))
        y = np.asarray(map(lambda l: float(l),data[1]))
        position = np.asarray(map(lambda l: float(l),data[2]))
    return x,y,position

path ='./filename_plus_seedID_number.txt'

x,y,position = readData(path)

plt.plot(position,x, 'r')
plt.plot(position,y, 'b')
plt.title('COSY orbit plot')
#plt.title('particle tracking')
plt.xlabel('position in meter')
plt.ylabel('closed orbit')
plt.legend(('x','y'),loc='upper right')

plt.show()
