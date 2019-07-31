import pandas as pd
import matplotlib.pyplot as plt
import os
import sys


## Usage: python plot_wfn.py (wfn number)

file = 'wfn'+str(sys.argv[-1])

df = pd.read_csv(file,comment='#',sep='\s+',header=None)
e = df[0]
exact = df[1]
pseudo = df[2]
proj = df[3]

plt.figure()

plt.plot(e,exact,'r-',label='Exact')
plt.plot(e,pseudo,'b-',label='Pseudo')
plt.plot(e,proj,'g-',label='Projector')

## Some plotting features
plt.legend(loc='upper left')
plt.xlabel('Energy',fontsize=16)
plt.ylabel('Log Deriv',fontsize=16)
plt.minorticks_on()

## Show you the plot
plt.show()
