import math
import os
import pandas as pd
import sys
import numpy as np
from shutil import copyfile
import sys


### Post analysis script for finding best solution
### Usage: python get_best_soln.py (number of variables in dakota run)
### Be sure to run in file containing dakota_tabular.dat

num_var = sys.argv[-1]
begin_obj_fn = int(num_var)+2
df = pd.read_table('dakota_tabular.dat',sep='\s+',header=None)
var_matrix = np.array(df.transpose()[:begin_obj_fn])
obj_fn_matrix = np.array(df.transpose()[begin_obj_fn:])
new_matrix = []
for row in obj_fn_matrix.transpose()[1:]:
    new_matrix.append([float(value) for value in row])
obj_fn_matrix = new_matrix

###for (var_list,obj_fn_list) in zip(var_matrix,obj_fn_matrix):
index = 0
for obj_fn_list in obj_fn_matrix:
    if float(obj_fn_list[0]) == 100.0:
        index += 1
        continue
    if 'Best_Solution' not in os.listdir('./'):
        os.mkdir('./Best_Solution')
    if 'results.out' in os.listdir('./Best_Solution/'):
        last_results_df = pd.read_table('./Best_Solution/results.out',sep='\s+',header=None)
        last_obj_fn_list = [float(value) for value in list(last_results_df[0])]
        index = 1
        for obj_fn in obj_fn_list:
            last_max = float(np.loadtxt('./Best_Solution/Max_Error_'+str(index)))
            if obj_fn > last_max:
                os.remove('./Best_Solution/Max_Error_'+str(index))
                f = open('./Best_Solution/Max_Error_'+str(index),'w+')
                f.write(str(obj_fn))
                f.close()
            index += 1
    else:
        index = 1
        for obj_fn in obj_fn_list:
            f = open('./Best_Solution/Max_Error_'+str(index),'w+')
            f.write(str(obj_fn))
            f.close()
            index += 1
    index = 1
    norm_obj_fn_list = []
    for obj_fn in obj_fn_list:
        max_value = float(np.loadtxt('./Best_Solution/Max_Error_'+str(index)))
        norm_obj_fn_list.append(obj_fn/max_value)
        index += 1
    rms_error = 0
    for obj_fn in norm_obj_fn_list:
        rms_error += obj_fn**2
    rms_error = math.sqrt(rms_error/len(norm_obj_fn_list))
    if 'results.out' in os.listdir('./Best_Solution/'):
        index = 1
        last_norm_obj_fn_list = []
        for obj_fn in last_obj_fn_list:
            max_value = float(np.loadtxt('./Best_Solution/Max_Error_'+str(index)))
            last_norm_obj_fn_list.append(obj_fn/max_value)
            index += 1
        last_rms_error = 0
        for obj_fn in last_norm_obj_fn_list:
            last_rms_error += obj_fn**2
        last_rms_error = math.sqrt(last_rms_error/len(last_norm_obj_fn_list))
    else:
        last_rms_error = 999999999.0
    if rms_error < last_rms_error:
        files_in_dir = os.listdir('./Best_Solution/')
        files_to_del = []
        for file in files_in_dir:
            if 'Max_Error' not in file:
                files_to_del.append(file)
        for filename in files_to_del:
            os.remove('./Best_Solution/'+filename)
        f = open('./Best_Solution/results.out','w+')
        for value in obj_fn_list:
            f.write(str(value)+'\n')
        f.close()
        f = open('./Best_Solution/rms_error','w+')
        f.write(str(rms_error))
        f.close()
        f = open('./Best_Solution/atompaw_var','w+')
        for var in var_matrix.transpose()[index]:
            f.write(str(var)+'\n')
        f.close()
    index += 1

