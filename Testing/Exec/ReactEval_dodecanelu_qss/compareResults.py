import sys
sys.path.append('utils')
import fileio as io
import numpy as np
import matplotlib.pyplot as plt
from plotsUtil import *

filePP_FDJ = 'PPreaction_FDJ.txt'
filePP_AJ = 'PPreaction_AJ.txt'
filePP_GMRES = 'PPreaction_GMRES.txt'

A_FDJ = io.readMultiColFile(filePP_FDJ)
A_AJ = io.readMultiColFile(filePP_AJ)
A_GMRES = io.readMultiColFile(filePP_GMRES)

fig=plt.figure()
plt.plot(A_FDJ[:,0],A_FDJ[:,1],linewidth=6, color='r',label='Finite Diff Jac')
plt.plot(A_AJ[:,0],A_AJ[:,1],linewidth=2, color='k',label='Analytical Jac')
plt.plot(A_GMRES[:,0],A_GMRES[:,1],'x', color='k',label='GMRES')
prettyLabels('time[s]','T[K]',14)
plotLegend()
plt.show()
