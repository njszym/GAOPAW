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
    if len(element_list) == 4:
        test_binary = True
    else:
        test_binary = False
    if len(element_list) == 6:
       	test_ternary = True
    else:
       	test_ternary = False
    num_tests = []
    try:
        test_atoms = input_settings['test_atomic_positions']
        num_tests.append('placeholder')
    except:
        test_atoms = False
    try:
        test_mag = input_settings['test_magnetization']
        num_tests.append('placeholder')
    except:
        test_mag = False
    try:
        test_mag_mom = input_settings['test_magnetic_moment']
        num_tests.append('placeholder')
    except:
	test_mag_mom = False
    try:
        test_gap = input_settings['test_gap']
        num_tests.append('placeholder')
    except:
        test_gap = False
    try:
        test_bulk = input_settings['test_bulk']
        num_tests.append('placeholder')
    except:
        test_bulk = False
    try:
        test_delta = input_settings['test_delta']
        num_tests.append('placeholder')
    except:
        test_delta = False
    try:
        test_phonon = input_settings['test_phonon']
        num_tests.append('placeholder')
    except:
        test_phonon = False
    try:
        test_lat = input_settings['test_lattice']
        num_tests.append('placeholder')
    except:
        test_lat = False
    lat_diff_list = []
    lanthanides = ['Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']
    ## Not including La, for now
    if element_list[0] not in lanthanides:
        for (elem,lat_type,lat_const) in zip(element_list,lat_type_list,lat_const_list):
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
                    QE_lat = get_lattice_constant(elem,lat_type)
                    AE_lat = lat_const
                    if check_convergence(elem,lat_type,'relax') == True:
                        lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                        copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
                    else:
                        pass
                except:
                    pass
            else:
                pass
            os.chdir('../')
    else:
        elem = element_list[0]
        os.mkdir(elem)
        write_atompaw_input(elem, template_dir)
        copyfile('./'+elem+'.atompaw.in',elem+'/'+elem+'.atompaw.in')
        os.chdir(elem)
        run_atompaw(elem)
        if check_UPF() == True:
            lat_diff_list.append('placeholder')
            lat_diff_list.append('placeholder')
        else:
            pass
    if len(lat_diff_list) == len(lat_type_list):
        if element_list[0] in lanthanides:
            lat_diff_list = []
            copyfile(template_dir+'/'+element_list[1]+'.GGA-PBE-paw.UPF','./'+element_list[1]+'.GGA-PBE-paw.UPF')
        else:
            pass
        if test_binary == True:
            bin_lat_type = input_settings['binary_lattice_type']
            unique_elem_list = unique(element_list)
            cmpd = unique_elem_list[0]+unique_elem_list[1]
            if test_atoms == True:
                write_QE_input(cmpd,bin_lat_type,'relax',template_dir)
                run_QE(cmpd,bin_lat_type,'relax')
                if check_convergence(cmpd,bin_lat_type,'relax') == True:
                    atom_diff = compare_atoms(cmpd,bin_lat_type,template_dir)
                    lat_diff_list.append(atom_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_mag == True:
                write_QE_input(cmpd,bin_lat_type,'scf',template_dir)
                run_QE(cmpd,bin_lat_type,'scf')
                AE_mag = float(input_settings['magnetization'])
                if check_convergence(cmpd,bin_lat_type,'scf') == True:
                    QE_mag = float(get_mag(cmpd,bin_lat_type))
                    lat_diff_list.append(abs(QE_mag-AE_mag))
                else:
                    lat_type_list.append('bad_run')
            if test_mag_mom == True:
                write_QE_input(cmpd,bin_lat_type,'scf',template_dir)
                run_QE(cmpd,bin_lat_type,'scf')
                if check_convergence(cmpd,bin_lat_type,'scf') == True:
                    mag_mom_diff = compare_mag_mom(cmpd,bin_lat_type,template_dir)
                    lat_diff_list.append(mag_mom_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_gap == True:
                write_QE_input(cmpd,bin_lat_type,'scf',template_dir)
                run_QE(cmpd,bin_lat_type,'scf')
                AE_gap = input_settings['band_gap']
                if check_convergence(cmpd,bin_lat_type,'scf') == True:
                    QE_gap = get_gap(cmpd,bin_lat_type)
                    lat_diff_list.append(abs(QE_gap-AE_gap))
                else:
                    lat_type_list.append('bad_run')
            if test_delta == True:
                write_QE_input(cmpd,bin_lat_type,'relax',template_dir)
                run_QE(cmpd,bin_lat_type,'relax')
                if check_convergence(cmpd,bin_lat_type,'relax') == True:
                    num_atoms = input_settings['num_atoms']
                    run_scale_lat(cmpd,bin_lat_type,template_dir)
                    V0, QE_bulk, B_prime = get_bulk(num_atoms)
                    QE_EOS_data, AE_EOS_data = read_eos(cmpd,template_dir)
                    delta_factor = calcDelta(QE_EOS_data,AE_EOS_data,[cmpd],False)
                    lat_diff_list.append(delta_factor)
                else:
                    lat_type_list.append('bad_run')
            if test_phonon == True: ## Testing required
                write_QE_input(cmpd,bin_lat_type,'scf',template_dir)
                run_QE(cmpd,bin_lat_type,'scf')
                if check_convergence(cmpd,bin_lat_type,'scf') == True:
                    copyfile(template_dir+'/phonon.in','./phonon.in')
                    os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py ph.x phonon.in -input_save '+cmpd+'.'+bin_lat_type+'.scf.save.qegz -MPICORES 4')
                    copyfile(template_dir+'/dynmat.in','./dynmat.in')
                    os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py dynmat.x dynmat.in -input_save phonon.save.qegz -MPICORES 4')
                    phonon_diff = compare_phonon(template_dir)
                    lat_diff_list.append(phonon_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_bulk == True:
                write_QE_input(cmpd,bin_lat_type,'relax',template_dir)
                run_QE(cmpd,bin_lat_type,'relax')
                if check_convergence(cmpd,bin_lat_type,'relax') == True:
                    num_atoms = input_settings['num_atoms']
                    run_scale_lat(cmpd,bin_lat_type,template_dir)
                    V0, QE_bulk, B_prime = get_bulk(num_atoms)
                    AE_bulk = input_settings['bulk_modulus']
                    bulk_diff = abs(AE_bulk-QE_bulk)
                    lat_diff_list.append(bulk_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_lat == True:
                write_QE_input(cmpd,bin_lat_type,'relax',template_dir)
                run_QE(cmpd,bin_lat_type,'relax')
                if check_convergence(cmpd,bin_lat_type,'relax') == True:
                    QE_lat = get_lattice_constant(cmpd,bin_lat_type)
                    AE_lat = input_settings['binary_lattice_constant']
                    lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                else:
                    lat_type_list.append('bad_run')
        if test_ternary == True:
            tern_lat_type = input_settings['ternary_lattice_type']
            unique_elem_list = unique(element_list)
            cmpd = unique_elem_list[0]+unique_elem_list[1]+unique_elem_list[2]
            if test_atoms == True:
                write_QE_input(cmpd,tern_lat_type,'relax',template_dir)
                run_QE(cmpd,tern_lat_type,'relax')
                if check_convergence(cmpd,tern_lat_type,'relax') == True:
                    atom_diff = compare_atoms(cmpd,tern_lat_type,template_dir)
                    lat_diff_list.append(atom_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_mag == True:
                write_QE_input(cmpd,tern_lat_type,'scf',template_dir)
                run_QE(cmpd,tern_lat_type,'scf')
                AE_mag = float(input_settings['magnetization'])
                if check_convergence(cmpd,tern_lat_type,'scf') == True:
                    QE_mag = float(get_mag(cmpd,tern_lat_type))
                    lat_diff_list.append(abs(QE_mag-AE_mag))
                else:
                    lat_type_list.append('bad_run')
            if test_mag_mom == True:
                write_QE_input(cmpd,tern_lat_type,'scf',template_dir)
                run_QE(cmpd,tern_lat_type,'scf')
                if check_convergence(cmpd,tern_lat_type,'scf') == True:
                    mag_mom_diff = compare_mag_mom(cmpd,tern_lat_type,template_dir)
                    lat_diff_list.append(mag_mom_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_gap == True:
                write_QE_input(cmpd,tern_lat_type,'scf',template_dir)
                run_QE(cmpd,tern_lat_type,'scf')
                AE_gap = input_settings['band_gap']
                if check_convergence(cmpd,tern_lat_type,'scf') == True:
                    QE_gap = get_gap(cmpd,tern_lat_type)
                    lat_diff_list.append(abs(QE_gap-AE_gap))
                else:
                    lat_type_list.append('bad_run')
            if test_delta == True:
                write_QE_input(cmpd,tern_lat_type,'relax',template_dir)
                run_QE(cmpd,tern_lat_type,'relax')
                if check_convergence(cmpd,tern_lat_type,'relax') == True:
                    num_atoms = input_settings['num_atoms']
                    run_scale_lat(cmpd,tern_lat_type,template_dir)
                    V0, QE_bulk, B_prime = get_bulk(num_atoms)
                    QE_EOS_data, AE_EOS_data = read_eos(cmpd,template_dir)
                    delta_factor = calcDelta(QE_EOS_data,AE_EOS_data,[cmpd],False)
                    lat_diff_list.append(delta_factor)
                else:
                    lat_type_list.append('bad_run')
            if test_phonon == True: ## Testing required
                write_QE_input(cmpd,tern_lat_type,'scf',template_dir)
                run_QE(cmpd,tern_lat_type,'scf')
                if check_convergence(cmpd,tern_lat_type,'scf') == True:
                    copyfile(template_dir+'/phonon.in','./phonon.in')
                    os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py ph.x phonon.in -input_save '+cmpd+'.'+bin_lat_type+'.scf.save.qegz -MPICORES 4')
                    copyfile(template_dir+'/dynmat.in','./dynmat.in')
                    os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py dynmat.x dynmat.in -input_save phonon.save.qegz -MPICORES 4')
                    phonon_diff = compare_phonon(template_dir)
                    lat_diff_list.append(phonon_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_bulk == True:
                write_QE_input(cmpd,tern_lat_type,'relax',template_dir)
                run_QE(cmpd,tern_lat_type,'relax')
                if check_convergence(cmpd,tern_lat_type,'relax') == True:
                    num_atoms = input_settings['num_atoms']
                    run_scale_lat(cmpd,tern_lat_type,template_dir)
                    V0, QE_bulk, B_prime = get_bulk(num_atoms)
                    AE_bulk = input_settings['bulk_modulus']
                    bulk_diff = abs(AE_bulk-QE_bulk)
                    lat_diff_list.append(bulk_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_lat == True:
                write_QE_input(cmpd,tern_lat_type,'relax',template_dir)
                run_QE(cmpd,tern_lat_type,'relax')
                if check_convergence(cmpd,tern_lat_type,'relax') == True:
                    QE_lat = get_lattice_constant(cmpd,tern_lat_type)
                    AE_lat = input_settings['ternary_lattice_constant']
                    lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                else:
                    lat_type_list.append('bad_run')
        if test_binary == False and test_ternary == False:
            if test_atoms == True:
                cmpd = element_list[0]
                cmpd_lat_type = input_settings['test_lat_type']
                write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                run_QE(cmpd,cmpd_lat_type,'relax')
                if check_convergence(cmpd,cmpd_lat_type,'relax') == True:
                    atom_diff = compare_atoms(cmpd,'atoms',template_dir)
                    lat_diff_list.append(atom_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_mag == True:
                cmpd = element_list[0]
                cmpd_lat_type = input_settings['test_lat_type']
                write_QE_input(cmpd,cmpd_lat_type,'scf',template_dir)
                run_QE(cmpd,cmpd_lat_type,'scf')
                AE_mag = float(input_settings['magnetization'])
                if check_convergence(cmpd,cmpd_lat_type,'scf') == True:
                    QE_mag = float(get_mag(cmpd,cmpd_lat_type))
                    lat_diff_list.append(abs(QE_mag-AE_mag))
                else:
                    lat_type_list.append('bad_run')
            if test_mag_mom == True:
                cmpd = element_list[0]
                cmpd_lat_type = input_settings['test_lat_type']
                write_QE_input(cmpd,cmpd_lat_type,'scf',template_dir)
                run_QE(cmpd,cmpd_lat_type,'scf')
                if check_convergence(cmpd,cmpd_lat_type,'scf') == True:
                    mag_mom_diff = compare_mag_mom(cmpd,cmpd_lat_type,template_dir)
                    lat_diff_list.append(mag_mom_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_gap == True:
                cmpd = element_list[0]
                cmpd_lat_type = input_settings['test_lat_type']
                write_QE_input(cmpd,cmpd_lat_type,'scf',template_dir)
                run_QE(cmpd,cmpd_lat_type,'scf')
                AE_gap = input_settings['band_gap']
                if check_convergence(cmpd,cmpd_lat_type,'scf') == True:
                    QE_gap = get_gap(cmpd,cmpd_lat_type)
                    lat_diff_list.append(abs(QE_gap-AE_gap))
                else:
                    lat_type_list.append('bad_run')
            if test_phonon == True: ## Testing required
                cmpd = element_list[0]
                cmpd_lat_type = input_settings['test_lat_type']
                write_QE_input(cmpd,cmpd_lat_type,'scf',template_dir)
                run_QE(cmpd,cmpd_lat_type,'scf')
                if check_convergence(cmpd,cmpd_lat_type,'scf') == True:
                    copyfile(template_dir+'/phonon.in','./phonon.in')
                    os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py ph.x phonon.in -input_save '+cmpd+'.'+bin_lat_type+'.scf.save.qegz -MPICORES 4')
                    copyfile(template_dir+'/dynmat.in','./dynmat.in')
                    os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py dynmat.x dynmat.in -input_save phonon.save.qegz -MPICORES 4')
                    phonon_diff = compare_phonon(template_dir)
                    lat_diff_list.append(phonon_diff)
                else:
                    lat_type_list.append('bad_run')
            if test_bulk == True:
                num_atoms = input_settings['num_atoms']
                cmpd = element_list[0]
                cmpd_lat_type = input_settings['test_lat_type']
                write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                run_QE(cmpd,cmpd_lat_type,'relax')
                run_scale_lat(cmpd,cmpd_lat_type,template_dir)
                V0, QE_bulk, B_prime = get_bulk(num_atoms)
                AE_bulk = input_settings['bulk_modulus']
                bulk_diff = abs(AE_bulk-QE_bulk)
                lat_diff_list.append(bulk_diff)
            if test_delta == True:
                num_atoms = input_settings['num_atoms']
                cmpd = element_list[0]
                cmpd_lat_type = input_settings['test_lat_type']
                write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                run_QE(cmpd,cmpd_lat_type,'relax')
                if check_convergence(cmpd,cmpd_lat_type,'relax') == True:
                    run_scale_lat(cmpd,cmpd_lat_type,template_dir)
                    V0, QE_bulk, B_prime = get_bulk(num_atoms)
                    QE_EOS_data, AE_EOS_data = read_eos(cmpd,template_dir)
                    delta_factor = calcDelta(QE_EOS_data,AE_EOS_data,[cmpd],False)
                    lat_diff_list.append(delta_factor)
                else:
                    lat_type_list.append('bad_run')
            if test_lat == True:
                num_atoms = input_settings['num_atoms']
                cmpd = element_list[0]
                cmpd_lat_type = input_settings['test_lat_type']
                write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
                run_QE(cmpd,cmpd_lat_type,'relax')
                if check_convergence(cmpd,cmpd_lat_type,'relax') == True:
                    QE_lat = get_lattice_constant(cmpd,cmpd_lat_type)
                    AE_lat = input_settings['test_lattice_constant']
                    lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                else:
                    lat_type_list.append('bad_run')
        if 'bad_run' not in lat_type_list:
            update_dakota(element_list,lat_diff_list)
        else:
            bad_run(element_list,lat_type_list)
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
        else:
            pass
    return check

def bad_run(element_list,lat_type_list):
    """
    If something went wrong with the run, e.g., no .UPF file created or
    if running QE raised an error, set the objective function to 100
    """
    params, results = di.read_parameters_file('params.in','results.out')
    unique_elem_list = unique(element_list)
    for (index,elem) in zip(range(len(unique_elem_list)),unique_elem_list):
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
    return abs(AE_lat - QE_lat)

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
    if lat_type == 'FCC' or lat_type == 'RS' or lat_type == 'ZB':
        return math.sqrt(2)*params[0]
    if lat_type == 'BCC':
        return (2./3.)*math.sqrt(3)*params[0]
    if lat_type == 'per':
        return params[0]

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
        else:
            pass
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
        if x not in unique_list: 
            unique_list.append(x) 
    return unique_list

def compare_atoms(elem,lat_type,template_path):
    """
    Compare atomic positions of QE-relaxed structure
    and those of the AE-relaxed structure...
    """
    df_AE = pd.read_table(template_path+'/AE_Struct',sep='\s+',header=None)
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
        else:
            pass
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
        else:
            pass
    QE_mag_mom = [float(value) for value in mag_mom]
    with open(template_path+'/AE_mag') as f:
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
        else:
            pass
    return float(band_gap)

def birch_murnaghan(V, V0, B0, B0_prime, E0):
    """
    3rd order Birch-Murnaghan equation of state, in the energy-volume form
    """
    V = np.array(V)
    return E0 + 9 * V0 * B0 / 16. * (
        ((V0 / V) ** (2 / 3.) - 1) ** 3 * B0_prime +
        ((V0 / V) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (V0 / V) ** (2 / 3.)))

def get_bulk(num_atoms):
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
    volume = popt[0]/float(num_atoms) ## A^3/atom
    bulk = popt[1]*160.2 ## GPa
    B_prime = popt[2] ## Dimensionless
    f = open('QE_EOS.txt','w+')
    f.write(str(volume)+' '+str(bulk)+' '+str(B_prime))
    f.close()
    return float(volume), float(bulk), float(B_prime)

def run_scale_lat(elem,lat_type,template_path):
    """
    Read in relaxed lattice parameter from (elem).(lat_type).relax.out,
    scale this lattice constant from -1% to +1% (10 values created),
    write input and run QE at each value, write corresponding volumes
    and energies into E_V.txt (units of Bohr^3 and Ry^3)
    """
    with open(elem+'.'+lat_type+'.relax.out') as f:
        lines = f.readlines()
    volumes = []
    for line in lines:
        if 'volume' in line.split():
            if 'new' in line.split():
                volumes.append(line.split()[4])
            else:
                volumes.append(line.split()[3])
        else:
            pass
    final_vol = float(volumes[-1])
    if lat_type == 'FCC' or 'ZB' or 'diamond' or 'RS':
        lat_const = (final_vol*4.)**(1./3.)
    else: ## Need to implement other lattice types
        pass
    scale_num = [0.99,0.9922,0.9944,0.9966,0.9988,1.0,1.0022,1.0044,1.0066,1.0088,1.01]
    scaled_lat = [num*lat_const for num in scale_num]
    with open(template_path+'/'+elem+'.'+lat_type+'.scf.template') as f:
        lines = f.readlines()
    index = 0
    for line in lines:
        if 'celldm' in line:
            cell_index = index
        else:
            pass
        index += 1
    energies = []
    folder = 0
    scf_file = elem+'.'+lat_type+'.scf.in'
    UPF_files = []
    files_in_folder = os.listdir('.')
    for file in files_in_folder:
        if file[-3:] == 'UPF':
            UPF_files.append(file)
        else:
            pass
    for value in scaled_lat:
        os.mkdir(str(folder))
        for file in UPF_files:
            copyfile(file,str(folder)+'/'+file)
        os.chdir(str(folder))
        lines[cell_index] = '  celldm(1)='+str(value)+'\n'
        with open(scf_file,'w') as f:
            for line in lines:
                f.write(line)
        os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py pw.x '+scf_file+' -MPICORES 4')
        with open(scf_file[:-2]+'out') as f:
            out_lines = f.readlines()
        for line in out_lines:
            if '!    total energy              =' in line:
                energies.append(line.split()[4])
        os.chdir('../')
        folder += 1
    if lat_type == 'FCC' or 'ZB' or 'diamond' or 'RS':
        volumes = [(value**3.)/4. for value in scaled_lat]
    else:
        pass ## Need to implement other lattice types
    f = open('E_V.txt','w+')
    for (e,v) in zip(energies,volumes):
        f.write(str(e)+' '+str(v)+'\n')
    f.close()

def read_eos(elem,template_path):
    """
    Read in QE and AE equilibrium volume, bulk modulus, and dB/dP
    from QE_EOS.txt and AE_EOS.txt
    """
    with open('QE_EOS.txt') as f:
        lines = f.readlines()
    QE_data = [float(value) for value in lines[0].split()]
    with open(template_path+'/AE_EOS.txt') as f:
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

def compare_phonon(template_path):
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
    freq_index = mode_line + 1
    freq = []
    check = True
    while check == True:
        try:
            freq.append(lines[freq_index].split()[2])
            freq_index += 1
        except:
            check = False
    QE_freq = sorted([float(value) for value in freq[3:]])
    with open(template_path+'/AE_freq') as f:
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

if __name__=='__main__':
    main()


