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
    try:
        test_binary = input_settings['test_binary']
    except:
        test_binary = False
    try:
        test_ternary = input_settings['test_ternary']
    except:
        test_ternary = False
    try:
        test_atoms = input_settings['test_atomic_positions']
    except:
        test_atoms = False
    try:
        test_mag = input_settings['test_magnetic_moment']
    except:
        test_mag = False
    try:
        test_gap = input_settings['test_gap']
    except:
        test_gap = False
    lat_diff_list = []
    for (elem,lat_type,lat_const) in zip(element_list,lat_type_list,lat_const_list):
        if elem not in os.listdir('.'):
            os.mkdir(elem)
        write_atompaw_input(elem, template_dir)
        copyfile('./'+elem+'.atompaw.in',elem+'/'+elem+'.atompaw.in')
        os.chdir(elem)
        run_atompaw(elem)
        if check_UPF() == True:
            write_QE_input(elem,lat_type,template_dir)
            run_QE(elem,lat_type)
            try:
                QE_lat = get_lattice_constant(elem,lat_type)
                AE_lat = lat_const
                if check_convergence(elem,lat_type) == True:
                    lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                    copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
                else:
                    pass
            except:
                pass
        else:
            pass
        os.chdir('../')
    if len(lat_diff_list) == len(lat_type_list):
        if test_binary == True:
            bin_lat_type = input_settings['binary_lattice_type']
            unique_elem_list = unique(element_list)
            cmpd = unique_elem_list[0]+unique_elem_list[1]
            write_QE_input(cmpd,bin_lat_type,template_dir)
            run_QE(cmpd,bin_lat_type)
            if test_atoms == True:
                atom_diff = compare_atoms(cmpd,bin_lat_type,template_dir)
                lat_diff_list.append(atom_diff)
                if check_convergence(cmpd,bin_lat_type) == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append(bin_lat_type)
                    bad_run(element_list,lat_type_list)
            if test_mag == True:
                QE_mag = get_mag(cmpd,bin_lat_type)
                AE_mag = input_settings['magnetization']
                lat_diff_list.append(abs(QE_mag-AE_mag))
                if check_convergence(cmpd,bin_lat_type) == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append(bin_lat_type)
                    bad_run(element_list,lat_type_list)
            if test_gap == True:
                QE_gap = get_gap(cmpd,bin_lat_type)
                AE_gap = input_settings['band_gap']
                lat_diff_list.append(abs(QE_gap-AE_gap))
                if check_convergence(cmpd,bin_lat_type) == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append(bin_lat_type)
                    bad_run(element_list,lat_type_list)
            if test_atoms == False and test_mag == False and test_gap == False:
                QE_lat = get_lattice_constant(cmpd,bin_lat_type)
                AE_lat = input_settings['binary_lattice_constant']
                lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                if check_convergence(cmpd,bin_lat_type) == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append(bin_lat_type)
                    bad_run(element_list,lat_type_list)
        if test_ternary == True:
            tern_lat_type = input_settings['ternary_lattice_type']
            unique_elem_list = unique(element_list)
            cmpd = unique_elem_list[0]+unique_elem_list[1]+unique_elem_list[2]
            write_QE_input(cmpd,tern_lat_type,template_dir)
            run_QE(cmpd,tern_lat_type)
            if test_atoms == True:
                atom_diff = compare_atoms(cmpd,tern_lat_type,template_dir)
                lat_diff_list.append(atom_diff)
                if check_convergence(cmpd,tern_lat_type) == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append(tern_lat_type)
                    bad_run(element_list,lat_type_list)
            if test_mag == True:
                QE_mag = get_mag(cmpd,tern_lat_type)
                AE_mag = input_settings['magnetization']
                lat_diff_list.append(abs(QE_mag-AE_mag))
                if check_convergence(cmpd,tern_lat_type) == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append(tern_lat_type)
                    bad_run(element_list,lat_type_list)
            if test_gap == True:
                QE_gap = get_gap(cmpd,tern_lat_type)
                AE_gap = input_settings['band_gap']
                lat_diff_list.append(abs(QE_gap-AE_gap))
                if check_convergence(cmpd,tern_lat_type) == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append(tern_lat_type)
                    bad_run(element_list,lat_type_list)
            if test_atoms == False and test_mag == False and test_gap == False:
                QE_lat = get_lattice_constant(cmpd,tern_lat_type)
                AE_lat = input_settings['ternary_lattice_constant']
                lat_diff_list.append(compare_lat(AE_lat,QE_lat))
                if check_convergence(cmpd,tern_lat_type) == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append(tern_lat_type)
                    bad_run(element_list,lat_type_list)
        if test_binary == False and test_ternary == False:
            if test_atoms == True:
                cmpd = element_list[0]
                write_QE_input(cmpd,'atoms',template_dir)
                run_QE(cmpd,'atoms')
                atom_diff = compare_atoms(cmpd,'atoms',template_dir)
                lat_diff_list.append(atom_diff)
                if check_convergence(cmpd,'atoms') == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append('atoms')
                    bad_run(element_list,lat_type_list)
            if test_mag == True:
                cmpd = element_list[0]
                write_QE_input(cmpd,'mag',template_dir)
                run_QE(cmpd,'mag')
                QE_mag = get_mag(cmpd,'mag')
                AE_mag = input_settings['magnetization']
                lat_diff_list.append(abs(QE_mag-AE_mag))
                if check_convergence(cmpd,'mag') == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append('mag')
                    bad_run(element_list,lat_type_list)
            if test_gap == True:
                cmpd = element_list[0]
                write_QE_input(cmpd,'gap',template_dir)
                run_QE(cmpd,'gap')
                QE_gap = get_gap(cmpd,'gap')
                AE_gap = input_settings['band_gap']
                lat_diff_list.append(abs(QE_gap-AE_gap))
                if check_convergence(cmpd,'gap') == True:
                    update_dakota(element_list,lat_diff_list)
                else:
                    lat_type_list.append('gap')
                    bad_run(element_list,lat_type_list)
            if test_atoms == False and test_mag == False and test_gap == False:
                update_dakota(element_list,lat_diff_list)
    else:
        if test_binary == True:
            lat_type_list.append(input_settings['binary_lattice_type'])
            bad_run(element_list,lat_type_list)
        if test_ternary == True:
            lat_type_list.append(input_settings['ternary_lattice_type'])
            bad_run(element_list,lat_type_list)
        if test_binary == False and test_binary == False:
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

def write_QE_input(elem,lat_type,template_path):
    """
    Write QE input file based on some template specified in template_path
    """
    template_file = os.path.join(template_path, elem+'.'+lat_type+'.relax.template')
    new_input_file = elem+'.'+lat_type+'.relax.in'
    shutil.copy(template_file,new_input_file)

def run_QE(elem,lat_type):
    """
    Run QE relaxation using (elem).relax.in
    """
    os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py pw.x '+elem+'.'+lat_type+'.relax.in -MPICORES 4')

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

def check_convergence(elem,lat_type):
    """
    Check if the QE run converged
    """
    check = True
    with open(elem+'.'+lat_type+'.relax.out') as f:
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
    return mag[-1]

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
            gap = (float(line.split()[7]) - float(line.split()[6]))
        else:
            pass
    return gap

def birch_murnaghan(V, V0, B0, B0_prime, E0):
    """
    3rd order Birch-Murnaghan equation of state, in the energy-volume form
    """
    V = np.array(V)
    return E0 + 9 * V0 * B0 / 16. * (
        ((V0 / V) ** (2 / 3.) - 1) ** 3 * B0_prime +
        ((V0 / V) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (V0 / V) ** (2 / 3.)))

def get_bulk():
    """
    Reads in energy-volume data from E_V.txt and calculates bulk modulus
    """
    df = pd.read_table('E_V.txt',sep='\s+',header=None)
    E = list(df[0])
    V = list(df[1])
    V = np.array([0.14818453429566825*value for value in V]) ## bohr^3 to A^3
    Y = np.array([13.6056980659*value for value in E]) ## Ry to eV
    initial_parameters = [V.mean(), 2.5, 4, Y.mean()]
    fit_eqn = eval('birch_murnaghan')
    popt, pcov = cf(fit_eqn, V, Y, initial_parameters)
    bulk = popt[1]*160.2
    return bulk

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
            volumes.append(line.split()[3])
        else:
            pass
    final_vol = float(volumes[-1])
    if lat_type == 'FCC':
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
    for value in scaled_lat:
        os.mkdir(str(folder))
        os.chdir(str(folder))
        lines[cell_index] = '  celldm(1)='+str(value)+'\n'
        with open(scf_file,'w') as f:
            for line in lines:
                f.write(line)
        os.system('pw.x < '+scf_file+' > '+scf_file[:-2]+'out')
        with open(scf_file[:-2]+'out') as f:
            out_lines = f.readlines()
        for line in out_lines:
            if '!    total energy              =' in line:
                energies.append(line.split()[4])
        os.chdir('../')
        folder += 1
    if lat_type == 'FCC':
        volumes = [(value**3.)/4. for value in scaled_lat]
    else:
        pass ## Need to implement other lattice types
    f = open('E_V.txt','w+')
    for (e,v) in zip(energies,volumes):
        f.write(str(e)+' '+str(v)+'\n')
    f.close()


if __name__=='__main__':
    main()

