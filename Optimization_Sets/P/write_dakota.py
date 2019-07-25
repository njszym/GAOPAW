import numpy as np
import json
from types import SimpleNamespace


template_dir = '/home/szymansk/Elem_Templates'

def main():
    """
    Parse information from input.json and update dakota.in accordingly
    """
    with open('input.json') as input:
        input_settings = json.load(input,object_hook=lambda d: SimpleNamespace(**d))
    cmpd_list = input_settings.compounds
    update_num_obj(cmpd_list)
    update_vars(cmpd_list)


def update_num_obj(cmpd_list):
    """
    Update total number of objective functions in dakota.in
    """
    element_list = []
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        element_list.extend(parse_elems(formula))
    element_list = unique(element_list)
    num_elems = len(element_list)
    if num_elems == 1 and len(cmpd.keys()) == 1:
        if element_list[0] in ['N','P']:
            print('\n'+'Number of objective functions: 2\n')
        else:
            print('\n'+'Number of objective functions: 3\n')
        return
    cmpd_diff_dict = form_cmpd_dict(cmpd_list)
    num_properties = 0
    for cmpd in cmpd_diff_dict.keys():
        for lat_type in cmpd_diff_dict[cmpd].keys():
            for property in cmpd_diff_dict[cmpd][lat_type].keys():
                num_properties += 1
    num_obj_fns = 0
    for elem in element_list:
        if elem in ['N','P']:
            num_obj_fns += 2
        else:
            num_obj_fns += 3
    num_obj_fns += num_properties
    with open('dakota.in') as dakota_input:
        orig_dakota = dakota_input.readlines()
    new_dakota = []
    for line in orig_dakota:
        new_line = line
        if 'num_objective_functions' in line:
            new_line = '    num_objective_functions =  '+str(num_obj_fns)+'\n'
        new_dakota.append(new_line)
    with open('dakota.in','w+') as dakota_input:
        for line in new_dakota:
            dakota_input.write(line)

def update_vars(cmpd_list):
    """
    Update variable bounds in dakota.in
    """
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
        elem_var_labels = ['"'+'DAKOTA_'+elem+'_'+label+'"' for label in elem_var_labels.split()]
        elem_var_labels = ' '.join(elem_var_labels)
        var_labels.append(elem_var_labels)
    temp_pts = []
    for elem_set in vars:
        temp_pts.append([float(value) for value in elem_set.split()])
    init_pts = np.concatenate(temp_pts)
    num_vars = len(init_pts)
    var_labels = ' '.join(var_labels)
    lower_bounds = [round(0.95*value,3) for value in init_pts] ## May be tuned
    check_lower = []
    for value in lower_bounds:
        if value < 0:
            check_lower.append(0.0)
        else:
            check_lower.append(value)
    lower_bounds = check_lower
    upper_bounds = [round(1.05*value,3) for value in init_pts] ## May be tuned
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

def form_cmpd_dict(cmpd_list):
    """
    Constructs empty dictionary of correct length for testing
    """
    cmpd_diff_dict = {}
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        lat_type = cmpd['lattice_type']
        property_list = [property for property in cmpd if property not in ['formula','lattice_type']]
        cmpd_diff_dict[formula] = {}
        cmpd_diff_dict[formula][lat_type] = {}
        for property in property_list:
            cmpd_diff_dict[formula][lat_type][property] = {}
    return cmpd_diff_dict

def test_cmpd_list(cmpd_list,cmpd_diff_dict,cmpd_template_dir):
    """
    Perform property tests for all compounds
    """
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        lat_type = cmpd['lattice_type']
        property_list = [property for property in cmpd if property not in ['formula','lattice_type']]
        for property in property_list:
            ae_value = cmpd[property]
            cmpd_diff_dict[formula][lat_type][property], error_check = test_property(formula,lat_type,property,ae_value,cmpd_template_dir)
            if error_check:
                bad_run(cmpd_diff_dict)
                return cmpd_diff_dict, True
    return cmpd_diff_dict, False

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

