from gaopaw import *


def main():
    """
    First test all constituent elements for compounds given in input.json, 
    ensuring pseudopotential generation and associated QE runs proceed without error, 
    then test specified properties and optimize using genetical algorithm.
    """
<<<<<<< HEAD
    gp_run = Runner()
    gp_run.get_paws()
    elem_objs, err_check = gp_run.test_element_list()
    if err_check:
        gp_run.bad_run()
        return
    if gp_run.is_elem:
        gp_run.update_obj_file(elem_objs)
        gp_run.update_dakota(elem_objs)
    cmpd_objs, err_check = gp_run.test_cmpd_list()
    if err_check:
        gp_run.bad_run()
        return
    all_objs = gp_run.merge_dicts(elem_objs, cmpd_objs)
    gp_run.update_obj_file(all_objs)
    gp_run.update_dakota(all_objs)
=======
    working_dir = os.getcwd()
    num_obj_fns = parse_num_objs(working_dir)
    with open(os.path.join(working_dir, os.pardir, 'input.json')) as input:
        input_settings = json.load(input, object_hook=lambda d: SimpleNamespace(**d))
    elem_template_dir = input_settings.directories.elem_template_dir
    if hasattr(input_settings.directories, 'cmpd_template_dir'):
        cmpd_template_dir = input_settings.directories.cmpd_template_dir
    if hasattr(input_settings.directories, 'include_paw'):
        include_paw = input_settings.directories.include_paw
        for paw_elem in include_paw:
            copyfile(os.path.join(cmpd_template_dir,'%s.GGA-PBE-paw.UPF' % paw_elem), 
                './%s.GGA-PBE-paw.UPF' % paw_elem)
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
        bad_run()
        return
    if len(cmpd_list) == 1 and len(element_list) == 1 and len(vars(cmpd_list[0])) == 1:
        update_obj_file(elem_diff_dict)
        update_dakota(elem_diff_dict)
        return
    cmpd_diff_dict = form_cmpd_dict(cmpd_list)
    cmpd_diff_dict = merge_dicts(elem_diff_dict, cmpd_diff_dict)
    cmpd_diff_dict, error_check = test_cmpd_list(cmpd_list, 
        cmpd_diff_dict, cmpd_template_dir, elem_template_dir)
    if error_check:
        bad_run()
        return
    update_obj_file(cmpd_diff_dict)
    update_dakota(cmpd_diff_dict)
>>>>>>> 018c78d0f5e5af81763f348bb96ab4f5366c7b45

if __name__=='__main__':
    main()
