from gaopaw import *


def main():
    """
    First test all constituent elements for compounds given in input.json, 
    ensuring pseudopotential generation and associated QE runs proceed without error,
    then test specified properties and optimize using genetical algorithm.
    """
    gp_run = Runner(input_dir='current', test_paw=True)
    elem_objs, err_check = gp_run.testElementList()
    err_mssg = \
        'The pseudopotential caused issues, check log for more details'
    if err_check:
        raise ValueError(err_mssg)
    cmpd_objs, err_check = gp_run.testCmpdList()
    if err_check:
        raise ValueError(err_mssg)
    all_objs = gp_run.mergeDicts(elem_objs, cmpd_objs)
    gp_run.updateObjFile(all_objs)

if __name__=='__main__':
    main()

