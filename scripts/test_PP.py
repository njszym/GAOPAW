from gaopaw import *


UPF_dir = '/home/szymansk/GBRV_Tests/PAWs'

def main():
    """
    First test all constituent elements for compounds given in input.json, 
    ensuring pseudopotential generation and associated QE runs proceed without error,
    then test specified properties and optimize using genetical algorithm.
    """
    working_dir = os.getcwd()
    with open(working_dir+'/input.json') as input:
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
    element_list = list(set(element_list))
    elem_diff_dict, error_check = test_elem_structs(element_list,elem_template_dir,UPF_dir)
    if error_check:
        raise ValueError('Your pseudopotential caused some issues')
        return
    if len(element_list) == 1 and len(cmpd.keys()) == 1:
        update_obj_file(elem_diff_dict)
        return
    cmpd_diff_dict = form_cmpd_dict(cmpd_list)
    cmpd_diff_dict = merge_dicts(elem_diff_dict,cmpd_diff_dict)
    cmpd_diff_dict, error_check = test_cmpd_list(cmpd_list,cmpd_diff_dict,cmpd_template_dir, elem_template_dir)
    if error_check:
        raise ValueError('Your pseudopotential caused some issues')
        return
    update_obj_file(cmpd_diff_dict)

def test_elem_structs(elem_list,template_dir,UPF_dir):
    """
    Only structures are tested (no log derivs) usin pre-made .UPF files in UPF_dir
    """
    elem_diff_dict = {}
    elemental_data = get_element_info(template_dir)
    for elem in elem_list:
        assert elem in elemental_data.keys(), \
            'No AE data available for your element: '+elem
        elem_diff_dict[elem] = {}
        elem_diff_dict[elem]['elemental'] = {}
        if elem in ['N','P']:
            if elem == 'N':
                elem_diff_dict[elem]['SC'] = {}
                elem_diff_dict[elem]['SC']['atomic_positions'] = {}
            if elem == 'P':
                elem_diff_dict[elem]['ortho'] = {}
                elem_diff_dict[elem]['ortho']['lattice_constant'] = {}
        else:
            for lat_type in ['FCC','BCC']:
                elem_diff_dict[elem][lat_type] = {}
                elem_diff_dict[elem][lat_type]['lattice_constant'] = {}
    for elem in elem_list:
        os.mkdir(elem)
        copyfile(UPF_dir+'/'+elem+'.GGA-PBE-paw.UPF','./'+elem+'.GGA-PBE-paw.UPF')
        with fileutils.chdir(elem):
            copyfile('../'+elem+'.GGA-PBE-paw.UPF','./'+elem+'.GGA-PBE-paw.UPF')
            if elem in ['N','P']:
                if elem == 'N':
                    run_qe(elem,'SC','relax',template_dir)
                    if not check_convergence(elem,'SC','relax'):
                        return elem_diff_dict, True
                    elem_diff_dict[elem]['SC']['atomic_positions'] = \
                         compare_atoms(elem,'SC',template_dir)
                if elem == 'P':
                    run_qe(elem,'ortho','relax',template_dir)
                    if not check_convergence(elem,'ortho','relax'):
                        return elem_diff_dict, True
                    ae_lat = elemental_data[elem]['ortho']
                    elem_diff_dict[elem]['ortho']['lattice_constant'] = \
                        compare_lat(ae_lat,elem,'ortho')
            else:
                for lat_type in ['FCC','BCC']:
                    run_qe(elem,lat_type,'relax',template_dir)
                    if not check_convergence(elem,lat_type,'relax'):
                        return elem_diff_dict, True
                    ae_lat = elemental_data[elem][lat_type]
                    elem_diff_dict[elem][lat_type]['lattice_constant'] = \
                        compare_lat(ae_lat,elem,lat_type)
    return elem_diff_dict, False


if __name__=='__main__':
    main()

