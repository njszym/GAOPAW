from gaopaw import *


def main():
    """
    Main function which Dakota executes to carry out evaluations.
    Parse input from gaopaw.yaml contained in a subfolder.
    For each element: write AtomPAW input, run AtomPAW to generate .UPF,
    write and run QE input, parse specified properties to be tested/optimized,
    compare with known AE properties and update objective functions accordingly.
    Logarithmic derivatives (arctan) of exact and pseudized partial waves are compared.
    PAWs may be optimized based on the following properties:
    lattice constants, atomic positions, net magnetization, individual magnetic moments,
    band gaps, bulk moduli, phonon frequencies, and delta-factors.
    """
    working_dir = os.getcwd()
    num_obj_fns = parse_num_objs(working_dir)
    with open(working_dir+'/../input.json') as input:
        input_settings = json.load(input,object_hook=lambda d: SimpleNamespace(**d))
    template_settings = input_settings.directories[0].__dict__
    elem_template_dir = template_settings['elem_template_dir']
    if 'cmpd_template_dir' in template_settings.keys():
        cmpd_template_dir = template_settings['cmpd_template_dir']
    cmpd_list = input_settings.compounds
    element_list = []
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        element_list.extend(parse_elems(formula))
    element_list = unique(element_list)
    assert num_obj_fns == get_num_objs(cmpd_list,element_list), \
        'Wrong number of objective functions specified, should be '+str(get_num_objs(cmpd_list,element_list))
    elem_diff_dict, error_check = test_element_list(element_list,elem_template_dir)
    if error_check:
        bad_run(num_obj_fns)
        return
    if len(element_list) == 1 and len(cmpd.keys()) == 1:
        update_dakota(elem_diff_dict)
        update_best_result(elem_diff_dict)
        return
    cmpd_diff_dict = form_cmpd_dict(cmpd_list)
    cmpd_diff_dict = merge_dicts(elem_diff_dict,cmpd_diff_dict)
    cmpd_diff_dict, error_check = test_cmpd_list(cmpd_list,cmpd_diff_dict,cmpd_template_dir)
    if error_check:
        bad_run(num_obj_fns)
        return
    update_dakota(cmpd_diff_dict)
    update_best_result(cmpd_diff_dict)

if __name__=='__main__':
    main()
