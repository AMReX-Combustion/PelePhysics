import numpy as np
import sys

f = open(sys.argv[1],'r')
lines = f.readlines()
f.close()

time=[]
for line in lines:
    time.append(float(line))

time = np.array(time)

print("time =  %.3f +/- %.3f" % (np.mean(time), np.std(time)))


