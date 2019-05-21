import pandas as pd
from scipy.signal import argrelextrema
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
    Logarithmic derivatives of exact and pseudized partial waves are also compared
    """
    working_dir = sys.argv[-3]
    with open(working_dir+'/../gaopaw.yaml') as f:
        input_settings = yaml.load(f)
    element_list = input_settings['elements']
    template_dir = input_settings['template_dir']
    lat_type_list = input_settings['lattice_type']
    lat_const_list = input_settings['lattice_constant']
    test_binary = input_settings['test_binary']
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
            unique_elem_list = unique(element_list)
            cmpd = unique_elem_list[0]+unique_elem_list[1]
            write_QE_input(cmpd,'RS',template_dir)
            run_QE(cmpd,'RS')
            QE_lat = get_lattice_constant(cmpd,'RS')
            AE_lat = input_settings['binary_lattice_constant']
            lat_diff_list.append(compare_lat(AE_lat,QE_lat))
            update_dakota(element_list,lat_diff_list)
        else:
            update_dakota(element_list,lat_diff_list)
    else:
        if test_binary == True:
            lat_type_list.append(input_settings['binary_lattice_type'])
            bad_run(element_list,lat_type_list)
        else:
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
    Compare the logarithmic derivatives of the pseudized and exact
    wavefunctions produced by AtomPAW...want to minimize
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
        log_pseudo = df[1]
        log_exact = df[2]
        sum_log += sum([abs(value) for value in log_exact])
        diff = []
        for (ps, ex) in zip(log_pseudo,log_exact):
            diff.append(abs(ps-ex))
        net_diff = sum(diff)
        total_diff += net_diff
    return total_diff/sum_log

def compare_log_peaks():
    """"
    Alternative method to compare log derivs...seems to converge a bit
    more slowly than compare_log()...
    Compare the energies at which local maxima occur in the log derivs
    of exact and pseudized partial waves...want to minimize difference
    If number of local maxima differs, a ghost state likely exists
    """
    files = os.listdir('./')
    log_derivs = []
    for file in files:
        if file[:4] == 'logd':
            log_derivs.append(file)
    net_diff = []
    total_diff = 0
    for file in log_derivs[:-1]:
        df = pd.read_table(file,sep='\s+',header=None)
        e = df[0]
        log_exact = np.array(df[1])
        log_pseudo = np.array(df[2])
        peaks_ind_ex = argrelextrema(log_exact,np.greater)[0]
        peaks_ind_ps = argrelextrema(log_pseudo,np.greater)[0]
        energy_diff = []
        if len(peaks_ind_ex) == len(peaks_ind_ps):
            for (p1, p2) in zip(peaks_ind_ex, peaks_ind_ps):
                energy_diff.append(abs(e[p1] - e[p2]))
        else:
            energy_diff.append(100.0)
        net_diff.append(sum(energy_diff))
    total_diff = sum(net_diff)
    return total_diff

def unique(list): 
    """
    Get list of unique elements to be tested
    """
    unique_list = []   
    for x in list:
        if x not in unique_list: 
            unique_list.append(x) 
    return unique_list

if __name__=='__main__':
    main()

