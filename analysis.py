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
    working_dir = sys.argv[-3]
    elem_template_dir = '/scr/szymansk/gaopaw/Elem_Templates'
    cmpd_template_dir = working_dir+'/../Input/'
    ## Will probably define these as global vars and remove from function args
    with open(working_dir+'/../input.json') as input:
        input_settings = json.load(input,object_hook=lambda d: SimpleNamespace(**d))
    cmpd_list = input_settings.compounds
    element_list = []
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        element_list.extend(parse_elems(formula))
    element_list = unique(element_list)
    elem_diff_dict, error_check = test_element_list(element_list,elem_template_dir)
    if len(element_list) == 1 and len(cmpd.keys()) == 1:
        if error_check:
            bad_run(elem_diff_dict)
            return
        else:
            update_dakota(elem_diff_dict)
            return
    cmpd_diff_dict = form_cmpd_dict(cmpd_list)
    cmpd_diff_dict.update(elem_diff_dict)
    if error_check:
        bad_run(cmpd_diff_dict)
        return
    cmpd_diff_dict, error_check = test_cmpd_list(cmpd_list,cmpd_diff_dict,cmpd_template_dir)
    if error_check:
        bad_run(cmpd_diff_dict)
        return
    update_dakota(cmpd_diff_dict)

if __name__=='__main__':
    main()
