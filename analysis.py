import pandas as pd
from scipy.signal import argrelextrema
from scipy.optimize import curve_fit as cf
import shutil
import os
import sys
sys.path.insert(0, '/scr/szymansk/dakota/share/dakota/Python')
import dakota.interfacing as di
from schrodinger.utils import subprocess
from schrodinger.utils import imputils
from schrodinger.utils import fileutils
from schrodinger.application.matsci import property_names as pnames
from schrodinger.application.matsci.nano import xtal
import math
import yaml
import numpy as np
from shutil import copyfile

def main():
    """
    Main function which Dakota executes to carry out evaluations.
    Parse input from gaopaw.yaml contained in a subfolder.
    For each element: write AtomPAW input, run AtomPAW and generate .UPF,
    write QE input, run QE relaxation, parse equilibrium lattice constant,
    compare with AE lattice constant and update objective function.
    Logarithmic derivatives (arctan) of exact and pseudized partial waves are compared.
    May also study binaries and ternaries, comparing the following:
    lattice constants, atomic positions, magnetic moments, ...
    """
    working_dir = sys.argv[-3]
    with open(working_dir+'/../gaopaw.yaml') as f:
        input_settings = yaml.load(f)
    element_list = input_settings['elements']
    template_dir = input_settings['template_dir']
    lat_type_list = input_settings['lattice_type']
    lat_const_list = input_settings['lattice_constant']
    num_tests = []
    try:
        cmpd_formula_list = input_settings['cmpd_formula']
        try:
            test_atoms_list = input_settings['test_atomic_positions']
            for TF_val in test_atoms_list:
                if TF_val == True:
                    num_tests.append('placeholder')
        except:
            test_atoms_list = [False]*len(cmpd_formula_list)
        try:
            test_mag_list = input_settings['test_magnetization']
            for TF_val in test_mag_list:
                if TF_val == True:
                    num_tests.append('placeholder')
        except:
            test_mag_list = [False]*len(cmpd_formula_list)
        try:
            test_mag_mom_list = input_settings['test_magnetic_moment']
            for TF_val in test_mag_mom_list:
                if TF_val == True:
                    num_tests.append('placeholder')
        except:
            test_mag_mom_list = [False]*len(cmpd_formula_list)
        try:
            test_gap_list = input_settings['test_gap']
            for TF_val in test_gap_list:
                if TF_val == True:
                    num_tests.append('placeholder')
        except:
            test_gap_list = [False]*len(cmpd_formula_list)
        try:
            test_bulk_list = input_settings['test_bulk']
            for TF_val in test_bulk_list:
                if TF_val == True:
                    num_tests.append('placeholder')
        except:
            test_bulk_list = [False]*len(cmpd_formula_list)
        try:
            test_delta_list = input_settings['test_delta']
            for TF_val in test_delta_list:
                if TF_val == True:
                    num_tests.append('placeholder')
        except:
            test_delta_list = [False]*len(cmpd_formula_list)
        try:
            test_phonon_list = input_settings['test_phonon']
            for TF_val in test_phonon_list:
                if TF_val == True:
                    num_tests.append('placeholder')
        except:
            test_phonon_list = [False]*len(cmpd_formula_list)
        try:
            test_lat_list = input_settings['test_lattice']
            for TF_val in test_lat_list:
                if TF_val == True:
                    num_tests.append('placeholder')
        except:
            test_lat_list = [False]*len(cmpd_formula_list)
        try:
            test_coh_list = input_settings['test_coh_energy']
            for TF_val in test_coh_list:
                if TF_val == True:
                    num_tests.append('placeholder')
        except:
            test_coh_list = [False]*len(cmpd_formula_list)
    except:
        pass
    lat_diff_list = []
    lanthanides = ['Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']
    check_error = False
    for (elem,lat_type,lat_const) in zip(element_list,lat_type_list,lat_const_list):
        if check_error == False:
            if elem not in lanthanides:
                if elem not in os.listdir('.'):
                    os.mkdir(elem)
                write_atompaw_input(elem, template_dir)
                copyfile('./'+elem+'.atompaw.in',elem+'/'+elem+'.atompaw.in')
                os.chdir(elem)
                run_atompaw(elem)
                if check_UPF() == True:
                    write_QE_input(elem,lat_type,'relax',template_dir)
                    run_QE(elem,lat_type,'relax')
                    try:
                        if lat_type in ['SC','RS','BCC','FCC','ZB','diamond','CsCl']:
                            QE_lat = get_lattice_constant(elem,lat_type)
                            AE_lat = lat_const
                            if check_convergence(elem,lat_type,'relax') == True:
                                lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                                copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
                                if elem == 'N':
                                    copyfile('N.SC.relax.out','../N.SC.relax.out')
                            else:
                                check_error = True
                        if lat_type == 'tetrag':
                            QE_a, QE_c = get_lattice_constant(elem,lat_type)
                            AE_lat = lat_const
                            lat_diff = compare_lat(AE_lat[0],QE_a)
                            lat_diff += compare_lat(AE_lat[1],QE_c)
                            lat_diff = lat_diff/2.
                            if check_convergence(elem,lat_type,'relax') == True:
                                lat_diff_list.append(lat_diff)
                                copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
                            else:
                                check_error = True
       	                if lat_type == 'ortho':
                            QE_a, QE_b, QE_c = get_lattice_constant(elem,lat_type)
                            AE_lat = lat_const
                            lat_diff = compare_lat(AE_lat[0],QE_a)
       	                    lat_diff +=	compare_lat(AE_lat[1],QE_b)
                            lat_diff += compare_lat(AE_lat[2],QE_c)
       	               	    lat_diff = lat_diff/3.
                            if check_convergence(elem,lat_type,'relax') == True:
                                lat_diff_list.append(lat_diff)
                                copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
                                if elem == 'P':
                                    copyfile('P.ortho.relax.out','../P.ortho.relax.out')
                                    copyfile('P.ortho.relax.in','../P.ortho.relax.in')
                            else:
       	       	       	        check_error = True
                    except:
                        check_error = True
                else:
                    check_error = True
                os.chdir('../')
            else:
                if elem not in os.listdir('.'):
                    os.mkdir(elem)
                write_atompaw_input(elem, template_dir)
                copyfile('./'+elem+'.atompaw.in',elem+'/'+elem+'.atompaw.in')
                os.chdir(elem)
                run_atompaw(elem)
                if check_UPF() == True:
                    copyfile(template_dir+'/N.GGA-PBE-paw.UPF','./N.GGA-PBE-paw.UPF')
                    write_QE_input(elem+'N','RS','relax',template_dir)
                    run_QE(elem+'N','RS','relax')
                    try:
                        QE_lat = get_lattice_constant(elem+'N','RS')
                        AE_lat = lat_const
                        if check_convergence(elem+'N','RS','relax') == True:
                            lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                            copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
                            copyfile(elem+'N.RS.relax.out','../'+elem+'N.RS.relax.out')
                        else:
                            check_error = True
                    except:
                        check_error = True
                else:
                    check_error = True
                os.chdir('../')
    if len(lat_diff_list) == len(lat_type_list):
        try:
            cmpd_index = 0
            cmpd_formula_list = input_settings['cmpd_formula']
            cmpd_lat_type_list = input_settings['cmpd_lattice_type']
            for (cmpd,cmpd_lat_type,test_atoms,test_mag,test_mag_mom,test_gap,test_bulk,test_delta,test_phonon,test_lat,test_coh) in zip(cmpd_formula_list,cmpd_lat_type_list,test_atoms_list,test_mag_list,test_mag_mom_list,test_gap_list,test_bulk_list,test_delta_list,test_phonon_list,test_lat_list,test_coh_list):
                try:
                    if test_atoms == True:
                        if str(cmpd+'.'+cmpd_lat_type+'.relax.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                            run_QE(cmpd,cmpd_lat_type,'relax')
                        if check_convergence(cmpd,cmpd_lat_type,'relax') == True:
                            atom_diff = compare_atoms(cmpd,cmpd_lat_type,template_dir)
                            lat_diff_list.append(atom_diff)
                       	else:
                            lat_type_list.append('bad_run')
                    if test_mag == True:
                        if str(cmpd+'.'+cmpd_lat_type+'.relax.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                            run_QE(cmpd,cmpd_lat_type,'relax')
                       	if str(cmpd+'.'+cmpd_lat_type+'.scf.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'scf',template_dir)
                            update_structure(cmpd,cmpd_lat_type,'scf')
                            run_QE(cmpd,cmpd_lat_type,'scf')
                       	AE_mag = float(input_settings['magnetization'][cmpd_index])
                       	if check_convergence(cmpd,cmpd_lat_type,'scf') == True:
                            QE_mag = float(get_mag(cmpd,cmpd_lat_type))
                            lat_diff_list.append(abs(QE_mag-AE_mag)/AE_mag)
                       	else:
                            lat_type_list.append('bad_run')
                    if test_mag_mom == True:
                       	if str(cmpd+'.'+cmpd_lat_type+'.relax.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                            run_QE(cmpd,cmpd_lat_type,'relax')
                        if str(cmpd+'.'+cmpd_lat_type+'.scf.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'scf',template_dir)
                            update_structure(cmpd,cmpd_lat_type,'scf')
                            run_QE(cmpd,cmpd_lat_type,'scf')
                       	if check_convergence(cmpd,cmpd_lat_type,'scf') == True:
                            mag_mom_diff = compare_mag_mom(cmpd,cmpd_lat_type,template_dir)
                            lat_diff_list.append(mag_mom_diff)
                       	else:
                            lat_type_list.append('bad_run')
                    if test_gap == True:
       	               	if str(cmpd+'.'+cmpd_lat_type+'.relax.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                            run_QE(cmpd,cmpd_lat_type,'relax')
       	               	if str(cmpd+'.'+cmpd_lat_type+'.scf.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'scf',template_dir)
                            update_structure(cmpd,cmpd_lat_type,'scf')
                            run_QE(cmpd,cmpd_lat_type,'scf')
                       	AE_gap = input_settings['band_gap'][cmpd_index]
                        if check_convergence(cmpd,cmpd_lat_type,'scf') == True:
                            QE_gap = get_gap(cmpd,cmpd_lat_type)
                            lat_diff_list.append(abs(QE_gap-AE_gap)/AE_gap)
                       	else:
                            lat_type_list.append('bad_run')
                    if test_coh == True:
                       	if str(cmpd+'.'+cmpd_lat_type+'.relax.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                            run_QE(cmpd,cmpd_lat_type,'relax')
                       	if check_convergence(cmpd,cmpd_lat_type,'relax') == True:
                            AE_coh = input_settings['cohesive_energy'][cmpd_index]
                            QE_coh = get_coh_energy(unique_elem_list,cmpd_lat_type,template_dir) ## deprecated
                            lat_diff_list.append(abs(QE_coh-AE_coh)/AE_coh)
                       	else:
                            lat_type_list.append('bad_run')
                    if test_delta == True:
                       	if str(cmpd+'.'+cmpd_lat_type+'.relax.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                            run_QE(cmpd,cmpd_lat_type,'relax')
                       	if check_convergence(cmpd,cmpd_lat_type,'relax') == True:
                            run_scale_lat(cmpd,cmpd_lat_type,template_dir)
                            V0, QE_bulk, B_prime = get_bulk(cmpd,cmpd_lat_type)
                            QE_EOS_data, AE_EOS_data = read_eos(cmpd,cmpd_lat_type,template_dir)
                            delta_factor = calcDelta(QE_EOS_data,AE_EOS_data,[cmpd],False)
                            lat_diff_list.append(delta_factor)
                       	else:
                            lat_type_list.append('bad_run')
                    if test_phonon == True: ## Testing required
       	               	if str(cmpd+'.'+cmpd_lat_type+'.relax.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                            run_QE(cmpd,cmpd_lat_type,'relax')
                       	if str(cmpd+'.'+cmpd_lat_type+'.scf.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'scf',template_dir)
                            update_structure(cmpd,cmpd_lat_type,'scf')
                            run_QE(cmpd,cmpd_lat_type,'scf')
                       	if check_convergence(cmpd,cmpd_lat_type,'scf') == True:
                            copyfile(template_dir+'/phonon.in','./phonon.in')
                            if 'phonon.save' in os.listdir('.'):
                                os.remove('phonon.out')
                                shutil.rmtree('phonon.save')
                                os.remove('phonon.save.qegz')
                            os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py ph.x phonon.in -input_save '+cmpd+'.'+cmpd_lat_type+'.scf.save.qegz -MPICORES 4')
                            copyfile(template_dir+'/dynmat.in','./dynmat.in')
                            if 'dynmat.save' in os.listdir('.'):
                                os.remove('dynmat.out')
                                shutil.rmtree('dynmat.save')
                                os.remove('dynmat.save.qegz')
                            os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py dynmat.x dynmat.in -input_save phonon.save.qegz -MPICORES 4')
                            phonon_diff = compare_phonon(cmpd,cmpd_lat_type,template_dir)
                            if phonon_diff == 'bad_run':
                               	lat_type_list.append('bad_run')
                            else:
                               	lat_diff_list.append(phonon_diff)
                       	else:
                            lat_type_list.append('bad_run')
                    if test_bulk == True:
                       	if str(cmpd+'.'+cmpd_lat_type+'.relax.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                            run_QE(cmpd,cmpd_lat_type,'relax')
                       	if check_convergence(cmpd,cmpd_lat_type,'relax') == True:
                            run_scale_lat(cmpd,cmpd_lat_type,template_dir)
                            V0, QE_bulk, B_prime = get_bulk(cmpd,cmpd_lat_type)
                            AE_bulk = input_settings['bulk_modulus'][cmpd_index]
                            bulk_diff = abs(AE_bulk-QE_bulk)/AE_bulk
                            lat_diff_list.append(bulk_diff)
                       	else:
                            lat_type_list.append('bad_run')
                    if test_lat == True:
                       	if str(cmpd+'.'+cmpd_lat_type+'.relax.out') not in os.listdir('.'):
                            write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                            run_QE(cmpd,cmpd_lat_type,'relax')
                       	if check_convergence(cmpd,cmpd_lat_type,'relax') == True:
                            if cmpd_lat_type in ['SC','FCC','BCC','ZB','per','RS','diamond','CsCl','HH']:
                               	QE_lat = get_lattice_constant(cmpd,cmpd_lat_type)
                               	AE_lat = input_settings['cmpd_lattice_constant'][cmpd_index]
                               	lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                            if cmpd_lat_type in ['tetrag','hex']:
                               	QE_a, QE_c = get_lattice_constant(cmpd,cmpd_lat_type)
                               	AE_lat = input_settings['cmpd_lattice_constant'][cmpd_index]
                               	lat_diff = compare_lat(AE_lat[0],QE_a)
                               	lat_diff += compare_lat(AE_lat[1],QE_c)
                               	lat_diff = lat_diff/2.
                               	lat_diff_list.append(lat_diff)
                            if cmpd_lat_type == 'ortho':
                               	QE_a, QE_b, QE_c = get_lattice_constant(cmpd,cmpd_lat_type)
                               	AE_lat = input_settings['cmpd_lattice_constant'][cmpd_index]
       	       	       	       	lat_diff = compare_lat(AE_lat[0],QE_a)
                               	lat_diff += compare_lat(AE_lat[1],QE_b)
                               	lat_diff += compare_lat(AE_lat[2],QE_c)
                                lat_diff = lat_diff/3.
                               	lat_diff_list.append(lat_diff)
                            if cmpd_lat_type == 'rhomb':
                               	QE_a, QE_angle = get_lattice_constant(cmpd,cmpd_lat_type)
                               	AE_lat = input_settings['cmpd_lattice_constant'][cmpd_index]
       	       	       	       	lat_diff = compare_lat(AE_lat[0],QE_a)
                               	lat_diff += compare_lat(AE_lat[1],QE_angle)
       	       	       	       	lat_diff = lat_diff/2.
                               	lat_diff_list.append(lat_diff)
                            if cmpd_lat_type == 'monoclin':
                               	QE_a, QE_b, QE_c, QE_angle = get_lattice_constant(cmpd,cmpd_lat_type)
                               	AE_lat = input_settings['cmpd_lattice_constant'][cmpd_index]
                               	lat_diff = compare_lat(AE_lat[0],QE_a)
                               	lat_diff += compare_lat(AE_lat[1],QE_b)
                               	lat_diff += compare_lat(AE_lat[2],QE_c)
                               	lat_diff += compare_lat(AE_lat[3],QE_angle)
                               	lat_diff = lat_diff/4.
                               	lat_diff_list.append(lat_diff)
                            if cmpd_lat_type == 'triclin':
                               	QE_a, QE_b, QE_c, QE_angle_1, QE_angle_2, QE_angle_3 = get_lattice_constant(cmpd,cmpd_lat_type)
                               	AE_lat = input_settings['cmpd_lattice_constant'][cmpd_index]
                               	lat_diff = compare_lat(AE_lat[0],QE_a)
                               	lat_diff += compare_lat(AE_lat[1],QE_b)
                               	lat_diff += compare_lat(AE_lat[2],QE_c)
                               	lat_diff += compare_lat(AE_lat[3],QE_angle_1)
                               	lat_diff += compare_lat(AE_lat[4],QE_angle_2)
                               	lat_diff += compare_lat(AE_lat[5],QE_angle_3)
                               	lat_diff = lat_diff/6.
                                lat_diff_list.append(lat_diff)
                        else:
                       	    lat_type_list.append('bad_run')
                except:
                    lat_type_list.append('bad_run')
                cmpd_index += 1
            if 'bad_run' not in lat_type_list:
                update_dakota(element_list,lat_diff_list)
                update_best_result()
            else:
                lat_type_list = [value for value in lat_type_list if value != 'bad_run']
                for i in range(len(num_tests)):
                    lat_type_list.append('placeholder')
                bad_run(element_list,lat_type_list)
        except:
            update_dakota(element_list,lat_diff_list)
            update_best_result()
    else:
        for i in range(len(num_tests)):
            lat_type_list.append('placeholder')
        bad_run(element_list,lat_type_list)

def check_UPF():
    """
    Check if a .UPF file was succesfully created by AtomPAW
    """
    files = os.listdir('./')
    check = False
    for file in files:
        if file[-3:] == 'UPF':
            check = True
    return check

def bad_run(element_list,lat_type_list):
    """
    If something went wrong with the run, e.g., no .UPF file created or
    if running QE raised an error, set the objective function to 100
    """
    params, results = di.read_parameters_file('params.in','results.out')
    unique_elem_list = unique(element_list)
    for index in range(len(unique_elem_list)):
        label = 'obj_fn_'+str(index+1)
        results[label].function = 100.0
    add_index = len(unique_elem_list)+1
    for index in range(len(lat_type_list)):
        label = 'obj_fn_'+str(index+add_index)
        results[label].function = 100
    results.write()

def compare_lat(AE_lat,QE_lat):
    """
    Compute difference between AE and PAW lattice constants
    """
    return abs(AE_lat - QE_lat)/AE_lat

def update_dakota(element_list,lat_diff_list):
    """
    Set the parameters and results files to be used by Dakota
    The objective function is equal to the difference between the lattice
    constants of AE calculations and PAW calculations performed here
    """
    params, results = di.read_parameters_file('params.in','results.out')
    unique_elem_list = unique(element_list)
    for (index,elem) in zip(range(len(unique_elem_list)),unique_elem_list):
        os.chdir(elem)
        label = 'obj_fn_'+str(index+1)
        results[label].function = compare_log()
        os.chdir('../')
    add_index = len(unique_elem_list)+1
    for index in range(len(lat_diff_list)):
        label = 'obj_fn_'+str(index+add_index)
        results[label].function = lat_diff_list[index]
    results.write()

def write_atompaw_input(elem,template_path):
    """
    Write AtomPAW input file based on some template specified in template_path
    """
    env = os.environ.copy()
    env['PATH'] = '/scr/fonari/dakota/bin:' + env['PATH']
    env['LD_LIBRARY_PATH'] = '/scr/fonari/dakota/bin:' + env['LD_LIBRARY_PATH']
    template_file = os.path.join(template_path, elem+'.atompaw.template')
    new_input_file = elem+'.atompaw.in'
    subprocess.check_call(['run', 'dprepro.py', 'params.in', template_file, new_input_file], env=env)

def run_atompaw(elem):
    """
    Run AtomPAW using (elem).atompaw.in
    """
    env = os.environ.copy()
    env['PATH'] = '/scr/szymansk/atompaw-4.1.0.5/src:' + env['PATH']
    with open(elem+'.atompaw.in','r') as input_fin, open('log_atompaw', 'w') as log_fout: 
        subprocess.call(['atompaw'], stdin=input_fin, stdout=log_fout, env=env)

def write_QE_input(elem,lat_type,calc_type,template_path):
    """
    Write QE input file based on some template specified in template_path
    """
    template_file = os.path.join(template_path, elem+'.'+lat_type+'.'+calc_type+'.template')
    new_input_file = elem+'.'+lat_type+'.'+calc_type+'.in'
    shutil.copy(template_file,new_input_file)

def run_QE(elem,lat_type,calc_type):
    """
    Run QE relaxation using (elem).(calc_type).in
    """
    os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py pw.x '+elem+'.'+lat_type+'.'+calc_type+'.in -MPICORES 4')

def get_lattice_constant(elem,lat_type):
    """
    Get relaxed lattice constant from QE run
    """
    qe_reader_path = os.path.join(fileutils.get_mmshare_scripts_dir(),'periodic_dft_gui_dir', 'qe2mae.py')
    qe_reader_mod = imputils.import_module_from_file(qe_reader_path)
    qe_reader = qe_reader_mod.QEOutputReader(elem+'.'+lat_type+'.relax.out')
    struct = qe_reader.structs[qe_reader.final_struct_id]
    cparams = xtal.get_chorus_properties(struct)
    params = xtal.get_params_from_chorus(cparams)
    if lat_type == 'FCC' or lat_type == 'RS' or lat_type == 'ZB' or lat_type == 'HH' or lat_tye == 'diamond':
        return math.sqrt(2)*params[0]
    if lat_type == 'BCC':
        return (2./3.)*math.sqrt(3)*params[0]
    if lat_type == 'per' or lat_type == 'SC' or lat_type == 'CsCl':
        return params[0]
    if lat_type == 'tetrag':
        unique_lat_list = sorted(unique(params[:3]))
        return unique_lat_list[0], unique_lat_list[1]
    if lat_type == 'ortho':
        unique_lat_list = sorted(params[:3])
        return unique_lat_list[0], unique_lat_list[1], unique_lat_list[2]
    if lat_type == 'hex':
        unique_lat_list = sorted(unique(params[:3]))
        return unique_lat_list[0], unique_lat_list[1]
    if lat_type == 'rhomb':
        lat = params[0]
        angle = params[4]
        return lat, angle
    if lat_type == 'monoclin':
        unique_lat_list = sorted(unique(params[:3]))
        for value in params[3:]:
            if value > 90.01 or value < 89.99:
                angle = value
        return unique_lat_list[0], unique_lat_list[1], unique_lat_list[2], angle
    if lat_type == 'triclin':
        unique_lat_list = sorted(unique(params[:3]))
        unique_angle_list = sorted(unique(params[3:]))
        all_params = unique_lat_list + unique_angle_list
        return all_params

def check_convergence(elem,lat_type,calc_type):
    """
    Check if the QE run converged
    """
    check = True
    with open(elem+'.'+lat_type+'.'+calc_type+'.out') as f:
        output = f.readlines()
    for line in output:
        if 'convergence NOT' in line:
            check = False
    return check

def compare_log():
    """"
    Compare arctan of the logarithmic derivatives of the pseudized and 
    exact wavefunctions produced by AtomPAW...want to minimize...
    Comparing arctan works better than using log derivs explicitly.
    Note that columns in logderiv.l correspond to:
    energy, exact logderiv, pseudized logderv, exact arctan, pseudized arctan
    """
    files = os.listdir('./')
    log_derivs = []
    for file in files:
        if file[:4] == 'logd':
            log_derivs.append(file)
    sum_log = 0
    total_diff = 0
    for file in log_derivs[:-1]:
        df = pd.read_table(file,sep='\s+',header=None)
        e = df[0]
        log_exact = df[3]
        log_pseudo = df[4]
        sum_log += sum([abs(value) for value in log_exact])
        diff = []
        for (ps, ex) in zip(log_pseudo,log_exact):
            diff.append(abs(ps-ex))
        net_diff = sum(diff)
        total_diff += net_diff
    return total_diff/sum_log

def unique(list): 
    """
    Get list of unique elements to be tested
    """
    unique_list = []   
    for x in list:
        try:
            value = round(float(x),3)
            if value not in unique_list:
                unique_list.append(value)
        except:
            if x not in unique_list: 
                unique_list.append(x) 
    return unique_list

def compare_atoms(elem,lat_type,template_path):
    """
    Compare atomic positions of QE-relaxed structure
    and those of the AE-relaxed structure...
    """
    df_AE = pd.read_table(template_path+'/AE_Struct.'+elem+'.'+lat_type,sep='\s+',header=None)
    df_AE = df_AE.drop(0,1)
    df_AE = df_AE.transpose()
    distance_list = []
    qe_reader_path = os.path.join(fileutils.get_mmshare_scripts_dir(),'periodic_dft_gui_dir', 'qe2mae.py')
    qe_reader_mod = imputils.import_module_from_file(qe_reader_path)
    qe_reader = qe_reader_mod.QEOutputReader(elem+'.'+lat_type+'.relax.out')
    struct = qe_reader.structs[qe_reader.final_struct_id]
    df_QE = struct.getXYZ()
    for index in range(len(df_AE.keys())):
        AE_position = np.array(df_AE[index])
        QE_position = np.array(df_QE[index])
        distance = np.linalg.norm(AE_position-QE_position)
        distance_list.append(distance)
    return sum(distance_list)

def get_mag(elem,lat_type):
    """
    Parse QE output (scf run) to obtain total magnetization
    """
    with open(elem+'.'+lat_type+'.scf.out') as f:
        lines = f.readlines()
    mag = []
    for line in lines:
        if 'absolute' in line.split():
            mag.append(line.split()[3])
    return float(mag[-1])

def compare_mag_mom(elem,lat_type,template_path):
    """
    Parse QE output (scf run) to obtain individual
    magnetic moments of atoms in given structure.
    Compare these with AE magnetic moments in AE_mag.
    """
    with open(elem+'.'+lat_type+'.scf.out') as f:
        lines = f.readlines()
    QE_mag_mom = []
    for line in lines:
        if 'magn:' in line.split():
            QE_mag_mom.append(line.split()[5])
    QE_mag_mom = [float(value) for value in mag_mom]
    with open(template_path+'/AE_mag.'+elem+'.'+lat_type) as f:
        lines = f.readlines()
    AE_mag_mom = []
    for line in lines:
        try:
            AE_mag_mom.append(float(line[:-1]))
        except:
            pass
    rel_diff = []
    for (QE,AE) in zip(QE_mag_mom,AE_mag_mom):
        rel_diff.append(abs((QE-AE)/AE))
    net_diff = sum(rel_diff)/len(rel_diff)
    return float(net_diff)

def get_gap(elem,lat_type):
    """
    Parse QE output (scf run) to obtain band gap.
    Note that unoccupied bands need to be included
    and occupations need to be fixed in the scf run.
    """
    with open(elem+'.'+lat_type+'.scf.out') as f:
        lines = f.readlines()
    for line in lines:
        if 'highest' and 'lowest' in line.split():
            band_gap = (float(line.split()[7]) - float(line.split()[6]))
    return float(band_gap)

def birch_murnaghan(V, V0, B0, B0_prime, E0):
    """
    3rd order Birch-Murnaghan equation of state, in the energy-volume form
    """
    V = np.array(V)
    return E0 + 9 * V0 * B0 / 16. * (
        ((V0 / V) ** (2 / 3.) - 1) ** 3 * B0_prime +
        ((V0 / V) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (V0 / V) ** (2 / 3.)))

def get_bulk(elem,lat_type):
    """
    Reads in energy-volume data from E_V.txt and calculates:
    equilibrium volume, bulk modulus, and dB/dP...
    i.e., fit to Birch Murnaghan equation of state
    """
    df = pd.read_table('E_V.txt',sep='\s+',header=None)
    E = list(df[0])
    V = list(df[1])
    V = np.array([0.14818453429566825*value for value in V]) ## bohr^3 to A^3
    Y = np.array([13.6056980659*value for value in E]) ## Ry to eV
    initial_parameters = [V.mean(), 2.5, 4, Y.mean()]
    fit_eqn = eval('birch_murnaghan')
    popt, pcov = cf(fit_eqn, V, Y, initial_parameters)

    with open(elem+'.'+lat_type+'.relax.in') as f:
        lines = f.readlines()
    for line in lines:
        if 'nat=' in line:
            num_atoms = float(line.split('=')[1][:-1])
    volume = popt[0]/num_atoms ## A^3/atom
    bulk = popt[1]*160.2 ## GPa
    B_prime = popt[2] ## Dimensionless
    f = open('QE_EOS.txt','w+')
    f.write(str(volume)+' '+str(bulk)+' '+str(B_prime))
    f.close()
    return float(volume), float(bulk), float(B_prime)

def run_scale_lat(elem,lat_type,template_path):
    """
    Read in relaxed cell parameters from (elem).(lat_type).relax.out,
    scale this lattice constant from 94% to 106% (7 values created),
    write input and run QE at each value, write corresponding volumes
    and energies into E_V.txt (units of Bohr^3 and Ry^3)
    """
    scale_num = [0.94,0.96,0.98,1.0,1.02,1.04,1.06]
    with open(template_path+'/'+elem+'.'+lat_type+'.relax.template') as f:
        lines = f.readlines()
    relax_file = elem+'.'+lat_type+'.relax.in'
    UPF_files = []
    files_in_folder = os.listdir('.')
    for file in files_in_folder:
        if file[-3:] == 'UPF':
            UPF_files.append(file)
    energies = []
    volumes = []
    folder_index = 1
    for value in scale_num:
        folder = elem+'_'+str(folder_index)
        new_cell_params = scale_cell(elem,lat_type,value)
        new_cell_matrix = np.matrix(new_cell_params)
        volumes.append(np.linalg.det(new_cell_matrix))
        os.mkdir(folder)
        for file in UPF_files:
            copyfile(file,folder+'/'+file)
        copyfile(relax_file,folder+'/'+relax_file)
        copyfile(relax_file[:-2]+'out',folder+'/'+relax_file[:-2]+'out')
        os.chdir(folder)
        update_structure(elem,lat_type,'relax')
        write_cell(elem,lat_type,new_cell_params)
        os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py pw.x '+relax_file+' -MPICORES 4')
        with open(relax_file[:-2]+'out') as f:
            out_lines = f.readlines()
        temp_energies = []
        for line in out_lines:
            if '!    total energy              =' in line:
                temp_energies.append(line.split()[4])
        energies.append(temp_energies[-1])
        os.chdir('../')
        folder_index += 1
    f = open('E_V.txt','w+')
    for (e,v) in zip(energies,volumes):
        f.write(str(e)+' '+str(v)+'\n')
    f.close()

def read_eos(elem,lat_type,template_path):
    """
    Read in QE and AE equilibrium volume, bulk modulus, and dB/dP
    from QE_EOS.txt and AE_EOS.txt
    """
    with open('QE_EOS.txt') as f:
        lines = f.readlines()
    QE_data = [float(value) for value in lines[0].split()]
    with open(template_path+'/AE_EOS.'+elem+'.'+lat_type) as f:
        lines = f.readlines()
    AE_data = [float(value) for value in lines[0].split()]
    data_f = {'element': [elem], 'V0': [QE_data[0]], 'B0': [QE_data[1]], 'BP': [QE_data[2]]}
    data_w = {'element': [elem], 'V0': [AE_data[0]], 'B0': [AE_data[1]], 'BP': [AE_data[2]]}
    return data_f, data_w

def calcDelta(data_f, data_w, eloverlap, useasymm):
    """
    Calculate the Delta using the data in data_f, data_w on
    element in eloverlap...taken from Delta Package...
    data_f: QE data (dict) calculated with PAWs
    data_w: AE data (dict) calculated with WIEN2k
    eloverlap: names (list) of elements/compounds
    useasymm: set to False for now...
    """
    v0w = np.zeros(len(eloverlap))
    b0w = np.zeros(len(eloverlap))
    b1w = np.zeros(len(eloverlap))
    v0f = np.zeros(len(eloverlap))
    b0f = np.zeros(len(eloverlap))
    b1f = np.zeros(len(eloverlap))
    elw = list(data_w['element'])
    elf = list(data_f['element'])
    for i in range(len(eloverlap)):
        searchnr = elw.index(eloverlap[i])
        v0w[i] = data_w['V0'][searchnr]
        b0w[i] = data_w['B0'][searchnr] * 10.**9. / 1.602176565e-19 / 10.**30.
        b1w[i] = data_w['BP'][searchnr]
        searchnr = elf.index(eloverlap[i])
        v0f[i] = data_f['V0'][searchnr]
        b0f[i] = data_f['B0'][searchnr] * 10.**9. / 1.602176565e-19 / 10.**30.
        b1f[i] = data_f['BP'][searchnr]
    vref = 30.
    bref = 100. * 10.**9. / 1.602176565e-19 / 10.**30.
    if useasymm:
        Vi = 0.94 * v0w
        Vf = 1.06 * v0w
    else:
        Vi = 0.94 * (v0w + v0f) / 2.
        Vf = 1.06 * (v0w + v0f) / 2.
    a3f = 9. * v0f**3. * b0f / 16. * (b1f - 4.)
    a2f = 9. * v0f**(7./3.) * b0f / 16. * (14. - 3. * b1f)
    a1f = 9. * v0f**(5./3.) * b0f / 16. * (3. * b1f - 16.)
    a0f = 9. * v0f * b0f / 16. * (6. - b1f)
    a3w = 9. * v0w**3. * b0w / 16. * (b1w - 4.)
    a2w = 9. * v0w**(7./3.) * b0w / 16. * (14. - 3. * b1w)
    a1w = 9. * v0w**(5./3.) * b0w / 16. * (3. * b1w - 16.)
    a0w = 9. * v0w * b0w / 16. * (6. - b1w)
    x = [0, 0, 0, 0, 0, 0, 0]
    x[0] = (a0f - a0w)**2
    x[1] = 6. * (a1f - a1w) * (a0f - a0w)
    x[2] = -3. * (2. * (a2f - a2w) * (a0f - a0w) + (a1f - a1w)**2.)
    x[3] = -2. * (a3f - a3w) * (a0f - a0w) - 2. * (a2f - a2w) * (a1f - a1w)
    x[4] = -3./5. * (2. * (a3f - a3w) * (a1f - a1w) + (a2f - a2w)**2.)
    x[5] = -6./7. * (a3f - a3w) * (a2f - a2w)
    x[6] = -1./3. * (a3f - a3w)**2.
    y = [0, 0, 0, 0, 0, 0, 0]
    y[0] = (a0f + a0w)**2 / 4.
    y[1] = 3. * (a1f + a1w) * (a0f + a0w) / 2.
    y[2] = -3. * (2. * (a2f + a2w) * (a0f + a0w) + (a1f + a1w)**2.) / 4.
    y[3] = -(a3f + a3w) * (a0f + a0w) / 2. - (a2f + a2w) * (a1f + a1w) / 2.
    y[4] = -3./20. * (2. * (a3f + a3w) * (a1f + a1w) + (a2f + a2w)**2.)
    y[5] = -3./14. * (a3f + a3w) * (a2f + a2w)
    y[6] = -1./12. * (a3f + a3w)**2.
    Fi = np.zeros_like(Vi)
    Ff = np.zeros_like(Vf)
    Gi = np.zeros_like(Vi)
    Gf = np.zeros_like(Vf)
    for n in range(7):
        Fi = Fi + x[n] * Vi**(-(2.*n-3.)/3.)
        Ff = Ff + x[n] * Vf**(-(2.*n-3.)/3.)
        Gi = Gi + y[n] * Vi**(-(2.*n-3.)/3.)
        Gf = Gf + y[n] * Vf**(-(2.*n-3.)/3.)
    Delta = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi))
    Deltarel = 100. * np.sqrt((Ff - Fi) / (Gf - Gi))
    if useasymm:
        Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
                 / v0w / b0w * vref * bref
    else:
        Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
                 / (v0w + v0f) / (b0w + b0f) * 4. * vref * bref
    return Delta ## Just return Delta-factor for now

def compare_phonon(elem,lat_type,template_path):
    """
    Parse optical phonon frequencies from QE run
    and compare with AE frequencies
    """
    with open('dynmat.out') as f:
        lines = f.readlines()
    index = 0
    for line in lines:
        if 'mode' in line:
            mode_line = index
        index += 1
    try:
        freq_index = mode_line + 1
    except:
        return 'bad_run'
    freq = []
    check = True
    while check == True:
        try:
            freq.append(lines[freq_index].split()[2])
            freq_index += 1
        except:
            check = False
    QE_freq = sorted([float(value) for value in freq[3:]])
    with open(template_path+'/AE_freq.'+elem+'.'+lat_type) as f:
        lines = f.readlines()
    AE_freq = []
    for line in lines:
        try:
            AE_freq.append(float(line[:-1]))
        except:
            pass
    AE_freq = sorted(AE_freq)
    rel_diff = []
    for (QE,AE) in zip(QE_freq,AE_freq):
        rel_diff.append(abs((QE-AE)/AE))
    net_diff = sum(rel_diff)/len(rel_diff)
    return net_diff

def update_structure(elem,lat_type,calc_type):
    """
    Parse equilibrium structure from completed relaxation
    and update the corresponding calculation input file.
    """
    with open(elem+'.'+lat_type+'.relax.out') as f:
        lines = f.readlines()
    index = 0
    for line in lines:
        if 'ATOMIC_POSITIONS' in line:
            start = index+1
            if 'alat' in line:
                coord_type = 'alat'
       	    if 'angstrom' in line:
                coord_type = 'angstrom'
       	    if 'crystal' in line:
                coord_type = 'crystal'
       	    if 'bohr' in line:
                coord_type = 'bohr'
        index += 1
    coords = []
    for line in lines[start:]:
        if 'End' in line:
            break
        else:
            coords.append(line)
    coords_header = 'ATOMIC_POSITIONS '+coord_type+'\n'
    index = 0
    for line in lines:
        if 'CELL_PARAMETERS' in line:
            cell_index = [index+1,index+2,index+3]
            split_line = line.split('=')
            try:
                alat = float(split_line[1][1:-2])
            except:
                alat = 1.88973 ## A to bohr
        index += 1
    vectors = []
    for i in cell_index:
        v = lines[i].split()
        v = [float(value) for value in v]
        v = [alat*value for value in v]
        vectors.append(v)
    cell_header = 'CELL_PARAMETERS bohr\n'
    v1 = str(vectors[0][0])+' '+str(vectors[0][1])+' '+str(vectors[0][2])+'\n'
    v2 = str(vectors[1][0])+' '+str(vectors[1][1])+' '+str(vectors[1][2])+'\n'
    v3 = str(vectors[2][0])+' '+str(vectors[2][1])+' '+str(vectors[2][2])+'\n'
    with open(elem+'.'+lat_type+'.'+calc_type+'.in') as f:
        lines = f.readlines()
    orig_struct = []
    line_index = 0
    for line in lines:
        if 'nat=' in line:
            natoms = int(line.split('=')[1][:-1])
    while line_index < len(lines):
        line = lines[line_index]
        if ('ATOMIC_POSITIONS' not in line) and ('CELL_PARAMETERS' not in line):
            if ('ibrav' in line) or ('celldm' in line):
                if 'ibrav' in line:
                    orig_struct.append('  ibrav=0\n')
            else:
                orig_struct.append(line)
            line_index += 1
        else:
            if 'CELL_PARAMETERS' in line:
                line_index += 4
            if 'ATOMIC_POSITIONS' in line:
                jump = natoms+1
                line_index += jump
    f = open(elem+'.'+lat_type+'.'+calc_type+'.in','w+')
    for line in orig_struct:
        f.write(line)
    f.write(coords_header)
    for line in coords:
        f.write(line)
    f.write(cell_header)
    f.write(v1)
    f.write(v2)
    f.write(v3)
    f.close()

def scale_cell(elem,lat_type,scale_factor):
    """
    Scale cell volume according to scale_factor
    """
    with open(elem+'.'+lat_type+'.relax.out') as f:
        lines = f.readlines()
    index = 0
    for line in lines:
        if 'CELL_PARAMETERS' in line:
            cell_index = [index+1,index+2,index+3]
            split_line = line.split('=')
            if 'alat' in line:
                alat = float(split_line[1][1:-2])
            if 'angstrom' in line:
                alat = 1.88973 ## A to Bohr
            if 'bohr' in line:
                alat = 1.00
        index += 1
    vectors = []
    for i in cell_index:
        v = lines[i].split()
        v = [float(value) for value in v]
        v = np.array([alat*value for value in v])
        vectors.append(v)
    M = np.matrix(vectors)
    i, j = [0,0]
    scaled_M = [[0,0,0],[0,0,0],[0,0,0]]
    for vector in np.array(M):
        j = 0
        for value in vector:
            scaled_M[i][j] = value*(scale_factor**(1./3.))
            j += 1
        i += 1
    return np.array(scaled_M)

def write_cell(elem,lat_type,cell):
    """
    Write given cell to QE relaxation input
    """
    vectors = np.array(cell)
    v1 = str(vectors[0][0])+' '+str(vectors[0][1])+' '+str(vectors[0][2])+'\n'
    v2 = str(vectors[1][0])+' '+str(vectors[1][1])+' '+str(vectors[1][2])+'\n'
    v3 = str(vectors[2][0])+' '+str(vectors[2][1])+' '+str(vectors[2][2])+'\n'
    cell_header = 'CELL_PARAMETERS bohr\n'
    with open(elem+'.'+lat_type+'.relax.in') as f:
        lines = f.readlines()
    orig_struct = []
    line_index = 0
    while line_index < len(lines):
        line = lines[line_index]
        if 'CELL_PARAMETERS' not in line:
            if ('ibrav' in line) or ('celldm' in line) or ('vc-relax' in line):
                if 'ibrav' in line:
                    orig_struct.append('  ibrav=0\n')
                if 'vc-relax' in line:
                    orig_struct.append("""  calculation='relax'\n""")
            else:
                orig_struct.append(line)
            line_index += 1
        else:
            line_index += 4
    f = open(elem+'.'+lat_type+'.relax.in','w+')
    for line in orig_struct:
        f.write(line)
    f.write(cell_header)
    f.write(v1)
    f.write(v2)
    f.write(v3)
    f.close()

def get_coh_energy(elem_list,lat_type,template_path):
    """
    Parse final energy from converged relaxation calculation
    for given cmpd, then for each constituent element in cmpd,
    perform scf run on the "atom in a box", obtain final energies
    per atom, and compute the difference, i.e., cohesive energy.
    Remains untested! Deprecated formatting (use of element_list)
    """
    cmpd = ''
    cmpd = cmpd.join(elem_list)
    with open(cmpd+'.'+lat_type+'.relax.out') as f:
        lines = f.readlines()
    index = 0
    for line in lines:
        if '!    total energy              =' in line:
            cmpd_energy = float(line.split()[4])
        if 'ATOMIC_POSITIONS' in line:
            start = index+1
        index += 1
    elems_in_cmpd = []
    for line in lines[start:]:
        if 'End' in line:
            break
        else:
            elems_in_cmpd.append(line.split()[0])
    formula = {}
    for elem in elems_in_cmpd:
        if elem in formula:
            formula[elem] += 1
        else:
            formula[elem] = 1
    elem_energies = {}
    for elem in elem_list:
        write_QE_input(elem,'atom','scf',template_path)
        run_QE(elem,'atom','scf')
        with open(elem+'.atom.scf.out') as f:
            if '!    total energy              =' in line:
                elem_energies[elem] = float(line.split()[4])
    sum_energies = 0
    for elem in elem_list:
        sum_energies += formula[elem]*elem_energies[elem]
    coh_energy = cmpd_energy - sum_energies
    return coh_energy

def update_best_result():
    """
    Parse dakota results and check overall fitness with
    respect to previous best solution. If current solution
    is better, replace old solution with current one.
    Note that fitness is normalized per the highest error
    for a given objective function.
    """
    UPF_files = []
    files_in_folder = os.listdir('.')
    for file in files_in_folder:
        if file[-3:] == 'UPF':
            UPF_files.append(file)
    atompaw_files = []
    for file in files_in_folder:
        if 'atompaw' in file:
            atompaw_files.append(file)
    if 'Best_Solution' not in os.listdir('../'):
        os.mkdir('../Best_Solution')
    while os.path.exists('../Best_Solution/WAIT'):
        time.sleep(1)
    f = open('../Best_Solution/WAIT','w+')
    f.write('wait to start until previous finishes')
    f.close()
    results_df = pd.read_table('results.out',sep='\s+',header=None)
    obj_fn_list = [float(value) for value in list(results_df[0])]
    if 'results.out' in os.listdir('../Best_Solution/'):
        last_results_df = pd.read_table('../Best_Solution/results.out',sep='\s+',header=None)
        last_obj_fn_list = [float(value) for value in list(last_results_df[0])]
        index = 1
        for obj_fn in obj_fn_list:
            last_max = float(np.loadtxt('../Best_Solution/Max_Error_'+str(index)))
            if obj_fn > last_max:
                os.remove('../Best_Solution/Max_Error_'+str(index))
                f = open('../Best_Solution/Max_Error_'+str(index),'w+')
                f.write(str(obj_fn))
                f.close()
            index += 1
    else:
        index = 1
        for obj_fn in obj_fn_list:
            f = open('../Best_Solution/Max_Error_'+str(index),'w+')
            f.write(str(obj_fn))
            f.close()
            index += 1
    index = 1
    norm_obj_fn_list = []
    for obj_fn in obj_fn_list:
        max_value = float(np.loadtxt('../Best_Solution/Max_Error_'+str(index)))
        norm_obj_fn_list.append(obj_fn/max_value)
        index += 1
    rms_error = 0
    for obj_fn in norm_obj_fn_list:
        rms_error += obj_fn**2
    rms_error = math.sqrt(rms_error/len(norm_obj_fn_list))
    if 'results.out' in os.listdir('../Best_Solution/'):
        index = 1
        last_norm_obj_fn_list = []
        for obj_fn in last_obj_fn_list:
            max_value = float(np.loadtxt('../Best_Solution/Max_Error_'+str(index)))
            last_norm_obj_fn_list.append(obj_fn/max_value)
            index += 1
        last_rms_error = 0
        for obj_fn in last_norm_obj_fn_list:
            last_rms_error += obj_fn**2
        last_rms_error = math.sqrt(last_rms_error/len(last_norm_obj_fn_list))
    else:
        last_rms_error = 999999999.0
    if rms_error < last_rms_error:
        files_in_dir = os.listdir('../Best_Solution/')
        files_to_del = []
        for file in files_in_dir:
            if ('Max_Error' not in file) and ('WAIT' not in file):
                files_to_del.append(file)
        for filename in files_to_del:
            os.remove('../Best_Solution/'+filename)
        copyfile('results.out','../Best_Solution/results.out')
        for (file_1,file_2) in zip(atompaw_files,UPF_files):
            copyfile(file_1,'../Best_Solution/'+file_1)
       	    copyfile(file_2,'../Best_Solution/'+file_2)
        f = open('../Best_Solution/rms_error','w+')
        f.write(str(rms_error))
        f.close()
    os.remove('../Best_Solution/WAIT')


if __name__=='__main__':
    main()


