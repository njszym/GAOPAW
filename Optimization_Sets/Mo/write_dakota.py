from gaopaw import *


def main():
    """
    Parse information from input.json and update dakota.in accordingly
    """
    with open('input.json') as input:
        input_settings = json.load(input,object_hook=lambda d: SimpleNamespace(**d))
    template_dir = input_settings.directories.elem_template_dir
    cmpd_list = input_settings.compounds
    update_num_obj(cmpd_list, template_dir)
    update_vars(cmpd_list, template_dir)
    update_labels(cmpd_list, template_dir)


def update_num_obj(cmpd_list, template_dir):
    """
    Update total number of objective functions in dakota.in
    """
    element_list = []
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        element_list.extend(parse_elems(formula))
    element_list = list(set(element_list))
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

def update_vars(cmpd_list, template_dir):
    """
    Update variable bounds in dakota.in
    """
    element_list = []
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        element_list.extend(parse_elems(formula))
    element_list = list(set(element_list))
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

def update_labels(cmpd_list, template_dir):
    """
    Write labels for objective functions.
    """
    element_list = []
    for cmpd in cmpd_list:
        element_list.extend(parse_elems(cmpd.formula))
    element_list = list(set(element_list))
    num_objs = get_num_objs(cmpd_list, element_list)
    cmpd_dict = form_cmpd_dict(cmpd_list)
    elem_dict = form_element_dict(element_list, template_dir)
    obj_fn_dict = merge_dicts(elem_dict, cmpd_dict)
    label_list = dict_to_list(obj_fn_dict)[1]
    label_string = ' '.join("""'%s'""" % label for label in label_list)
    with open('dakota.in') as dakota_input:
        orig_dakota = dakota_input.readlines()
    new_dakota = []
    index = 0
    for line in orig_dakota:
        if 'responses' in line:
            response_index = index
        index += 1
    index = 0
    for line in orig_dakota:
        new_line = line
        if 'descriptors' in line and index > response_index:
            new_line = '    descriptors =  '+label_string+'\n'
        new_dakota.append(new_line)
        index += 1
    with open('dakota.in','w+') as dakota_input:
        for line in new_dakota:
            dakota_input.write(line)



if __name__=='__main__':
    main()

