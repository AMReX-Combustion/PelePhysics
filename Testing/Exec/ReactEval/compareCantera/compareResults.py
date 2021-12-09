import sys
sys.path.append('../utils')
import fileio as io
import numpy as np
import matplotlib.pyplot as plt
from plotsUtil import *

filePP = '../PPreaction.txt'
fileCT = 'canteraSim/output/CanteraReaction_dodecane.txt'

A1 = io.readMultiColFile(filePP)
A3 = io.readMultiColFile(fileCT)

fig=plt.figure()
plt.plot(A3[:,0],A3[:,1],'x',linewidth=3, color='b',label='Cantera')
plt.plot(A1[:,0],A1[:,1],linewidth=3,color='k',label='PelePhysics')
prettyLabels('time[s]','T[K]',14)
plotLegend()
plt.show()
