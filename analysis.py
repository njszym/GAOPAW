from gaopaw import *


def main():
    """
    First test all constituent elements for compounds given in input.json, 
    ensuring pseudopotential generation and associated QE runs proceed without error, 
    then test specified properties and optimize using genetical algorithm.
    """
    gp_run = Runner()
    gp_run.getPaws()
    elem_objs, err_check = gp_run.testElementList()
    if err_check:
        gp_run.badRun()
        return
    if gp_run.is_elem:
        gp_run.updateObjFile(elem_objs)
        gp_run.updateDakota(elem_objs)
        return
    cmpd_objs, err_check = gp_run.testCmpdList()
    if err_check:
        gp_run.badRun()
        return
    all_objs = gp_run.mergeDicts(elem_objs, cmpd_objs)
    gp_run.updateObjFile(all_objs)
    gp_run.updateDakota(all_objs)

if __name__=='__main__':
    main()
