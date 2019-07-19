import numpy as np
import sys
import json
from types import SimpleNamespace


template_dir = '/scr/szymansk/gaopaw/Elem_Templates'

def main():
    with open('input.json') as input:
        input_settings = json.load(input,object_hook=lambda d: SimpleNamespace(**d))
    cmpd_list = input_settings.compounds
    element_list = []
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        element_list.extend(parse_elems(formula))
    element_list = unique(element_list)
    vars = []
    var_labels = []
    for elem in element_list:
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
    temp_pts = []
    for elem_set in vars:
        temp_pts.append([float(value) for value in elem_set.split()])
    init_pts = np.concatenate(temp_pts)
    num_vars = len(init_pts)
    var_labels = ' '.join(var_labels)
    lower_bounds = [round(0.9*value,3) for value in init_pts]
    check_lower = []
    for value in lower_bounds:
        if value < 0:
            check_lower.append(0.0)
        else:
            check_lower.append(value)
    lower_bounds = check_lower
    upper_bounds = [round(1.1*value,3) for value in init_pts]
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
    with open('dakota.in','w+') as dakota_input:
        for line in new_dakota:
            dakota_input.write(line)


def unique(value_list):
    """
    Get list of unique elements to be tested
    """
    try: ## if list of numbers
        value_list = [round(float(value),3) for value in value_list]
    except ValueError: ## if list of strings
        list = [str(value) for value in value_list]
    unique_list = []
    for value in value_list:
        if value not in unique_list:
            unique_list.append(value)
    return unique_list

def parse_elems(formula):
    """
    Parse compound formula to obtain constituent elements
    """
    letters_only = ''.join([letter for letter in formula if not letter.isdigit()])
    index = -1
    elems = []
    for letter in letters_only:
        if letter.isupper():
            elems.append(letter)
            index += 1
        else:
            elems[index] += letter
    return elems

if __name__=='__main__':
    main()

