import numpy as np
import sys


template_dir = '/scr/szymansk/gaopaw/Elem_Templates'
elem_list = sys.argv[1:]
vars = []
var_labels = []
for elem in elem_list:
    with open(template_dir+'/BOUNDS') as var_bounds:
        var_list = var_bounds.readlines()
        index = 0
        for line in var_list:
            if elem+':' in line:
                elem_vars = var_list[index+1]
                elem_var_labels = var_list[index+2]
                break
            index += 1
    vars.append(elem_vars)
    elem_var_labels = ['DAKOTA_'+elem+'_'+label for label in elem_var_labels.split()]
    elem_var_labels = ' '.join(elem_var_labels)
    var_labels.append(elem_var_labels)

init_pts = np.loadtxt(vars).flatten()
num_vars = len(init_pts)
var_labels = ' '.join(var_labels)
lower_bounds = [round(0.85*value,3) for value in init_pts]
upper_bounds = [round(1.15*value,3) for value in init_pts]
init_pts = ' '.join([str(value) for value in init_pts])
lower_bounds = ' '.join([str(value) for value in lower_bounds])
upper_bounds = ' '.join([str(value) for value in upper_bounds])

with open('dakota.in') as dakota_input:
    orig_dakota = dakota_input.readlines()

new_dakota = []

for line in orig_dakota:
    new_line = line
    if 'continuous_design' in line:
        new_line = '    continuous_design =  '+str(num_vars)+'\n'
    if 'initial_point' in line:
        new_line = '    initial_point =  '+init_pts+'\n'
    if 'lower_bounds' in line:
        new_line = '    lower_bounds =  '+lower_bounds+'\n'
    if 'upper_bounds' in line:
        new_line = '    upper_bounds =  '+upper_bounds+'\n'
    if 'descriptors' in line:
        new_line = '    descriptors =  '+var_labels+'\n'
    new_dakota.append(new_line)

with open('test.in','w+') as dakota_input:
    for line in new_dakota:
        dakota_input.write(line)
