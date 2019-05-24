import pandas as pd
import sys
import numpy as np


## Determine accuracy of each final solution from a completed Dakota run
## Best solution: minimize logderiv difference and lattice constant differences 
## To ensure accuracy among all structures, also minimize standard deviation
## Usage: python check_solns.py (number of variables to check)

num_obj = int(sys.argv[-1])
df = pd.read_table('finaldata1.dat',sep='\s+',header=None)
print ('\n',df,'\n')
i = 0
obj_list = []
index = len(df.keys()) - 1
for i in range(num_obj):
    obj_list.append(df[index])
    index -= 1
obj_list = np.array(obj_list)
obj_list = np.transpose(obj_list)
solutions = []
index = 1
indices = []
for obj in obj_list:
    net_diff = 0
    std_dev = np.std(obj[1:])
    for value in obj:
        net_diff += value**2
    net_diff = (net_diff + std_dev**2)**0.5
    print (index-1, ':', net_diff)
    solutions.append(net_diff)
    indices.append(index)
    index += 1
min_diff = min(solutions)
for (a,b) in zip(indices,solutions):
    if b == min_diff:
        print ('\nBest Solution: ', a-1, b,'\n')
    else:
        pass
