import pandas as pd
import matplotlib.pyplot as plt
import os

## Gather log deriv files in directory
## We only care about derivs up to l_max - 1
files = os.listdir('./')
log_derivs = []
for file in files:
    if file[:4] == 'logd':
        log_derivs.append(file)


## Begin figure
plt.figure()

## Possible orbitals and colors
orbital_list = ['s','p','d']
color_list = ['r','b','g']
color_index = 0

## For each l, plot exact and pseudo derivs
for file in log_derivs[:-1]:
    df = pd.read_table(file,sep='\s+',header=None)
    e = df[0]
    log_pseudo = df[2]
    log_exact = df[1]
    plt.plot(e,log_exact,color_list[color_index],label='Exact '+orbital_list[color_index])
    plt.plot(e,log_pseudo,color_list[color_index],linestyle='dashed',label='Pseudo '+orbital_list[color_index])
    color_index += 1

## Some plotting features
plt.legend(prop={'size':12},loc='lower left')
plt.xlabel('Energy (Ry)',fontsize=16)
plt.ylabel('Log Deriv',fontsize=16)
plt.minorticks_on()

plt.xlim(-1,4)
plt.ylim(-150,100)

## Show you the plot
plt.savefig('Logderivs.png')
