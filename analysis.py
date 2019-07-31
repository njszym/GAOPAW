from gaopaw import *


def main():
    """
    First test all constituent elements for compounds given in input.json, 
    ensuring pseudopotential generation and associated QE runs proceed without error, 
    then test specified properties and optimize using genetical algorithm.
    """
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

if __name__=='__main__':
    main()
