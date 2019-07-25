from gaopaw import *


def main():
    """
    First test all constituent elements for compounds given in input.json, 
    ensuring pseudopotential generation and associated QE runs proceed without error, 
    then test specified properties and optimize using genetical algorithm.
    """
    working_dir = os.getcwd()
    num_obj_fns = parse_num_objs(working_dir)
    with open(os.path.join(working_dir, os.pardir, 'input.json')) as input:
        input_settings = json.load(input, object_hook=lambda d: SimpleNamespace(**d))
    elem_template_dir = input_settings.directories.elem_template_dir
    if hasattr(input_settings.directories, 'cmpd_template_dir'):
        cmpd_template_dir = input_settings.directories.cmpd_template_dir
    cmpd_list = input_settings.compounds
    element_list = []
    for cmpd in cmpd_list:
        element_list.extend(parse_elems(cmpd.formula))
    element_list = list(set(element_list))
    assert num_obj_fns == get_num_objs(cmpd_list, element_list), \
        'Wrong number of objective functions specified, should be %s' \
             % get_num_objs(cmpd_list, element_list)
    elem_diff_dict, error_check = test_element_list(element_list, elem_template_dir)
    if error_check:
        bad_run(num_obj_fns)
        return
    if len(cmpd_list) == 1 and len(element_list) == 1 and len(vars(cmpd_list[0])) == 1:
        update_obj_file(elem_diff_dict)
        obj_fn_list = dict_to_list(elem_diff_dict)[0]
        update_best_result(obj_fn_list)
        update_dakota(elem_diff_dict)
        return
    cmpd_diff_dict = form_cmpd_dict(cmpd_list)
    cmpd_diff_dict = merge_dicts(elem_diff_dict, cmpd_diff_dict)
    cmpd_diff_dict, error_check = test_cmpd_list(cmpd_list, 
        cmpd_diff_dict, cmpd_template_dir, elem_template_dir)
    if error_check:
        bad_run(num_obj_fns)
        return
    update_obj_file(cmpd_diff_dict)
    obj_fn_list = dict_to_list(cmpd_diff_dict)[0]
    update_best_result(obj_fn_list)
    update_dakota(cmpd_diff_dict)

if __name__=='__main__':
    main()
