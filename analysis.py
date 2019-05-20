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
    ## Testing with multiples elements and/or lattice types is required
    ## Perhaps better to test all elem/lat, then sum up objectives and update dakota
    ## However, be careful, still need to detect errors, i.e., bad_run()
    for (elem,lat_type,lat_const) in zip(element_list,lat_type_list,lat_const_list):
        write_atompaw_input(elem, template_dir)
        run_atompaw(elem)
        if check_UPF() == True:
            write_QE_input(elem,lat_type,template_dir)
            run_QE(elem,lat_type)
            try:
                QE_lat = get_lattice_constant(elem,lat_type)
                AE_lat = lat_const
                if check_convergence(elem,lat_type) == True:
                    update_dakota(AE_lat,QE_lat)
                else:
                    bad_run()
            except:
                bad_run()
        else:
            bad_run()

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

def bad_run():
    """
    If something went wrong with the run, e.g., no .UPF file created or
    if running QE raised an error, set the objective function to 100
    """
    params, results = di.read_parameters_file('params.in','results.out')
    results['obj_fn_1'].function = 100
    results['obj_fn_2'].function = 100
    results.write()

def update_dakota(AE_lat,QE_lat):
    """
    Set the parameters and results files to be used by Dakota
    The objective function is equal to the difference between the lattice
    constants of AE calculations and PAW calculations performed here
    """
    params, results = di.read_parameters_file('params.in','results.out')
    results['obj_fn_1'].function = abs(AE_lat - QE_lat)
    results['obj_fn_2'].function = compare_log()
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
    if lat_type == 'FCC' or lat_type == 'RS':
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

#def compare_log():
#    """"
#    Compare the energies at which local maxima occur in the
#    exact and pseudized partial waves...want to minimize difference
#    If number of local maxima differs, a ghost state likely exists
#    Testing required, not sure if this works well
#    """
#    files = os.listdir('./')
#    log_derivs = []
#    for file in files:
#        if file[:4] == 'logd':
#            log_derivs.append(file)
#    net_diff = []
#    total_diff = 0
#    for file in log_derivs[:-1]:
#        df = pd.read_table(file,sep='\s+',header=None)
#        e = df[0]
#        log_exact = np.array(df[1])
#        log_pseudo = np.array(df[2])
#        peaks_ind_ex = argrelextrema(log_exact,np.greater)[0]
#        peaks_ind_ps = argrelextrema(log_pseudo,np.greater)[0]
#        energy_diff = []
#        if len(peaks_ind_ex) == len(peaks_ind_ps):
#            for (p1, p2) in zip(peaks_ind_ex, peaks_ind_ps):
#                energy_diff.append(abs(e[p1] - e[p2]))
#        else:
#            energy_diff.append(10.0)
#        net_diff.append(sum(energy_diff))
#    total_diff = sum(net_diff)
#    return total_diff

if __name__=='__main__':
    main()

