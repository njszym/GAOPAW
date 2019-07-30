from gaopaw import *


def main():
    """
    Get best solution from dakota_tabular.dat
    """
    with open('input.json') as input:
        input_settings = json.load(input, object_hook=lambda d: SimpleNamespace(**d))
    template_dir = input_settings.directories.elem_template_dir
    cmpd_list = input_settings.compounds
    element_list = []
    for cmpd in cmpd_list:
        element_list.extend(parse_elems(cmpd.formula))
    element_list = list(set(element_list))
    num_objs = get_num_objs(cmpd_list, element_list)
    cmpd_dict = form_cmpd_dict(cmpd_list)
    elem_dict = form_element_dict(element_list, template_dir)
    obj_fn_dict = merge_dicts(elem_dict, cmpd_dict)
    label_list = dict_to_list(obj_fn_dict)[1]
    df = pd.read_table('dakota_tabular.dat',sep='\s+',header=None)
    col_headers = np.array(df)[0][2:]
    var_headers = col_headers[:-num_objs]
    obj_headers = col_headers[-num_objs:]
    value_table = np.array(df)[1:][:, 2:]
    var_table = value_table[:, :-num_objs]
    obj_table = value_table[:, -num_objs:]
    good_vars, good_objs = [], []
    for (var_list, obj_list) in zip(var_table, obj_table):
        if '100' not in obj_list:
            good_vars.append([float(value) for value in var_list])
            good_objs.append([float(value) for value in obj_list])
    min_values, max_values = [], []
    for all_obj_values in np.array(good_objs).transpose():
        min_values.append(min(all_obj_values))
        max_values.append(max(all_obj_values))
    mae_list = []
    for result_set in good_objs:
        obj_index = 0
        normalized_err = []
        for value in result_set:
            utopia = min_values[obj_index]
            nadir = max_values[obj_index]
            normalized_err.append( abs( (value - utopia)/(nadir - utopia) ) )
            obj_index += 1
        mae_list.append(np.average(normalized_err))
    best_result_set = good_objs[np.argmin(mae_list)]
    best_var_set = good_vars[np.argmin(mae_list)]
    for (label, var) in zip(var_headers, best_var_set):
        print(label, var)
    for (label, obj) in zip(obj_headers, best_result_set):
        print(label, obj)



if __name__=='__main__':
    main()


