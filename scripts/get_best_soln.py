import math
import os
import pandas as pd
import sys
import numpy as np
from shutil import copyfile
import sys


df = pd.read_table('dakota_tabular.dat',sep='\s+',header=None)
result_matrix = np.array(df)[1:].transpose()
num_objs = int(sys.argv[-1])
slice_index = -1
obj_fn_matrix = []
for i in range(num_objs):
    obj_fn_matrix.append(result_matrix[slice_index])
    slice_index -= 1
obj_fn_matrix.reverse()
new_matrix = []
for obj_fn_list in obj_fn_matrix:
    new_list = []
    for value in obj_fn_list:
        if round(float(value),1) != 100.0:
            new_list.append(float(value))
    new_matrix.append(new_list)
obj_fn_matrix = np.array(new_matrix)
min_values = []
max_values = []
for obj_fn_list in obj_fn_matrix:
    min_values.append(min(obj_fn_list))
    max_values.append(max(obj_fn_list))
obj_fn_matrix = obj_fn_matrix.transpose()
rms_errors = []
for result_set in obj_fn_matrix:
    result_index = 0
    norm_errors = []
    for value in result_set:
        utopia = min_values[result_index]
        nadir = max_values[result_index]
        norm_errors.append( abs( (value - utopia)/(nadir - utopia) ) )
        result_index += 1
    rms = 0
    for value in norm_errors:
        rms += value
    rms = rms/len(norm_errors)
    rms_errors.append(rms)
min_rms = min(rms_errors)
min_index = np.argmin(rms_errors)
final_result = obj_fn_matrix[min_index]
for value in final_result:
    print(round(value,5))
