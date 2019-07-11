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
    with open(working_dir+'/../gaopaw.yaml') as f:
        input_settings = yaml.load(f)
    element_list = input_settings['elements']
    template_dir = input_settings['template_dir']
    lat_type_list = input_settings['lattice_type']
    lat_const_list = input_settings['lattice_constant']
    num_tests = 0
    if 'cmpd_formula' in input_settings.keys():
        cmpd_formula_list = input_settings['cmpd_formula']
        test_cmpds = True
        cmpd_lat_type_list = input_settings['cmpd_lattice_type']
        param_list = ['test_atomic_positions','test_magnetization','test_magnetic_moment','test_gap','test_bulk','test_delta','test_phonon','test_lattice']
        test_param_list = [cmpd_formula_list,cmpd_lat_type_list]
        for param in param_list:
            if param in input_settings.keys():
                test_param_list.append(input_settings[param])
                for boolean_val in input_settings[param]:
                    if boolean_val == True:
                        num_tests += 1
            else:
                test_param_list.append([False]*len(cmpd_formula_list))
    diff_list = []
    lanthanides = ['Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']
    check_error = False
    for (elem,lat_type,lat_const) in zip(element_list,lat_type_list,lat_const_list):
        if check_error == False:
            if elem not in os.listdir('.'):
                os.mkdir(elem)
            write_atompaw_input(elem, template_dir)
            copyfile('./'+elem+'.atompaw.in',elem+'/'+elem+'.atompaw.in')
            os.chdir(elem)
            run_atompaw(elem)
            if check_UPF() == True:
                if elem not in lanthanides:
                    run_QE(elem,lat_type,'relax',template_dir)
                    if check_convergence(elem,lat_type,'relax') == True and check_error == False:
                        lat_diff = compare_lat(lat_const,elem,lat_type)
                        diff_list.append(lat_diff)
                        copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
                        if elem == 'N':
                            copyfile('N.SC.relax.out','../N.SC.relax.out')
                        if elem == 'P':
                            copyfile('P.ortho.relax.in','../P.ortho.relax.in')
                    else:
                        check_error = True
                else:
                    copyfile(template_dir+'/N.GGA-PBE-paw.UPF','./N.GGA-PBE-paw.UPF')
                    run_QE(elem+'N','RS','relax',template_dir)
                    if check_convergence(elem+'N','RS','relax') == True and check_error == False:
                        diff_list.append(compare_lat(lat_const,elem+'N','RS'))
                        copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
                        copyfile(elem+'N.RS.relax.out','../'+elem+'N.RS.relax.out')
                    else:
                        check_error = True
            else:
                check_error = True
            os.chdir('../')
    if len(diff_list) == len(lat_type_list):
        if 'PAWs' in input_settings.keys():
            PAW_list = input_settings['PAWs']
            UPF_list = [elem_name+'.GGA-PBE-paw.UPF' for elem_name in PAW_list]
            for UPF in UPF_list:
                copyfile(template_dir+'/'+UPF,'./'+UPF)
        if test_cmpds == True:
            cmpd_index = 0
            error_check = []
            for (cmpd,cmpd_lat_type,test_atoms,test_mag,test_mag_mom,test_gap,test_bulk,test_delta,test_phonon,test_lat) in zip(*test_param_list):
                if test_atoms == True:
                    run_QE(cmpd,cmpd_lat_type,'relax',template_dir)
                    if check_convergence(cmpd,cmpd_lat_type,'relax') == True and True not in error_check:
                        atom_diff = compare_atoms(cmpd,cmpd_lat_type,template_dir)
                        diff_list.append(atom_diff)
                    else:
                        error_check.append(True)
                if test_mag == True:
                    run_QE(cmpd,cmpd_lat_type,'relax',template_dir)
                    if check_convergence(cmpd,cmpd_lat_type,'relax') == True and True not in error_check:
                        run_QE(cmpd,cmpd_lat_type,'scf',template_dir)
                        if check_convergence(cmpd,cmpd_lat_type,'scf') == True and True not in error_check:
                            AE_mag = float(input_settings['magnetization'][cmpd_index])
                            QE_mag = float(get_mag(cmpd,cmpd_lat_type))
                            diff_list.append(abs(QE_mag-AE_mag)/AE_mag)
                        else:
                            error_check.append(True)
                    else:
                        error_check.append(True)
                if test_mag_mom == True:
                    run_QE(cmpd,cmpd_lat_type,'relax',template_dir)
                    if check_convergence(cmpd,cmpd_lat_type,'relax') == True and True not in error_check:
                        run_QE(cmpd,cmpd_lat_type,'scf',template_dir)
                        if check_convergence(cmpd,cmpd_lat_type,'scf') == True and True not in error_check:
                            mag_mom_diff = compare_mag_mom(cmpd,cmpd_lat_type,template_dir)
                            diff_list.append(mag_mom_diff)
                        else:
                            error_check.append(True)
                    else:
                        error_check.append(True)
                if test_gap == True:
       	            run_QE(cmpd,cmpd_lat_type,'relax',template_dir)
                    if check_convergence(cmpd,cmpd_lat_type,'relax') == True and True not in error_check:
                        run_QE(cmpd,cmpd_lat_type,'scf',template_dir)
                        if check_convergence(cmpd,cmpd_lat_type,'scf') == True and True not in error_check:
                            AE_gap = input_settings['band_gap'][cmpd_index]
                            QE_gap = get_gap(cmpd,cmpd_lat_type)
                            diff_list.append(abs(QE_gap-AE_gap)/AE_gap)
                        else:
                            error_check.append(True)
                    else:
                        error_check.append(True)
                if test_delta == True:
                    run_QE(cmpd,cmpd_lat_type,'relax',template_dir)
                    if check_convergence(cmpd,cmpd_lat_type,'relax') == True and True not in error_check:
                        run_scale_lat(cmpd,cmpd_lat_type,template_dir)
                        V0, QE_bulk, B_prime = get_bulk(cmpd,cmpd_lat_type)
                        QE_EOS_data, AE_EOS_data = read_eos(cmpd,cmpd_lat_type,template_dir)
                        delta_factor = calcDelta(QE_EOS_data,AE_EOS_data,[cmpd],False)
                        diff_list.append(delta_factor[0])
                    else:
                        error_check.append(True)
                if test_phonon == True:
                    run_QE(cmpd,cmpd_lat_type,'relax',template_dir)
                    if check_convergence(cmpd,cmpd_lat_type,'relax') == True and True not in error_check:
                        run_QE(cmpd,cmpd_lat_type,'scf',template_dir)
                        if check_convergence(cmpd,cmpd_lat_type,'scf') == True and True not in error_check:
                            run_phonon(cmpd,cmpd_lat_type,template_dir)
                            phonon_diff = compare_phonon(cmpd,cmpd_lat_type,template_dir)
                            if phonon_diff == 'bad_run':
                                error_check.append(True)
                            else:
                                diff_list.append(phonon_diff)
                        else:
                           diff_list.append(phonon_diff)
                    else:
                        error_check.append(True)
                if test_bulk == True:
                    run_QE(cmpd,cmpd_lat_type,'relax',template_dir)
                    if check_convergence(cmpd,cmpd_lat_type,'relax') == True and True not in error_check:
                        run_scale_lat(cmpd,cmpd_lat_type,template_dir)
                        V0, QE_bulk, B_prime = get_bulk(cmpd,cmpd_lat_type)
                        AE_bulk = input_settings['bulk_modulus'][cmpd_index]
                        bulk_diff = abs(AE_bulk-QE_bulk)/AE_bulk
                        diff_list.append(bulk_diff)
                    else:
                        error_check.append(True)
                if test_lat == True:
                    run_QE(cmpd,cmpd_lat_type,'relax',template_dir)
                    if check_convergence(cmpd,cmpd_lat_type,'relax') == True and True not in error_check:
                        AE_lat = input_settings['cmpd_lattice_constant'][cmpd_index]
                        lat_diff = compare_lat(AE_lat,cmpd,cmpd_lat_type)
                        diff_list.append(lat_diff)
                    else:
                        error_check.append(True)
                cmpd_index += 1
            if True not in error_check:
                update_best_result(diff_list)
                update_dakota(element_list,diff_list)
            else:
                for i in range(num_tests):
                    lat_type_list.append('placeholder')
                bad_run(element_list,lat_type_list)
        else:
            update_best_result(diff_list)
            update_dakota(element_list,diff_list)
    else:
        for i in range(num_tests):
            lat_type_list.append('placeholder')
        bad_run(element_list,lat_type_list)


if __name__=='__main__':
    main()
