import sys
sys.path.append('../utils')
import fileio as io
import numpy as np
import matplotlib.pyplot as plt
from plotsUtil import *

filePP = '../PPreaction.txt'
A1 = io.readMultiColFile(filePP)
ignitionDelayPP = (2/20000)*np.amin(np.argwhere(A1[:,1]>(A1[0,1]+400)))
print("Ignition delay PP = ", ignitionDelayPP)
