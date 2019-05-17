import pandas as pd
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


def main():

    write_atompaw_input(['Si'],'/scr/szymansk/my_opal/Si_Input')
    run_atompaw(['Si'])
    if check_UPF() == True:
        write_QE_input(['Si'],'/scr/szymansk/my_opal/Si_Input')
        run_QE(['Si'])
        try:
            QE_lat = get_lattice_constant(['Si'],'FCC')
            AE_lat = 3.857
            if check_convergence(['Si']) == True:
                parse_dakota(AE_lat,QE_lat)
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

def parse_dakota(AE_lat,QE_lat):
    """
    Set the parameters and results files to be used by Dakota
    The objective function is equal to the difference between the lattice
    constants of AE calculations and PAW calculations performed here
    """
    params, results = di.read_parameters_file('params.in','results.out')
    results['obj_fn_1'].function = abs(AE_lat - QE_lat)
    results['obj_fn_2'].function = compare_log()
    results.write()

def write_atompaw_input(element_list,template_path):
    """
    Write AtomPAW input file based on some template specified in template_path
    """
    env = os.environ.copy()
    env['PATH'] = '/scr/fonari/dakota/bin:' + env['PATH']
    env['LD_LIBRARY_PATH'] = '/scr/fonari/dakota/bin:' + env['LD_LIBRARY_PATH']
    for elem in element_list:
        template_file = os.path.join(template_path, elem+'.atompaw.template')
        new_input_file = elem+'.atompaw.in'
        subprocess.check_call(['run', 'dprepro.py', 'params.in', template_file, new_input_file], env=env)

def run_atompaw(element_list):
    """
    Run AtomPAW using (elem).atompaw.in
    """
    env = os.environ.copy()
    env['PATH'] = '/scr/szymansk/atompaw-4.1.0.5/src:' + env['PATH']
    for elem in element_list:
        with open(elem+'.atompaw.in','r') as input_fin, open('log_atompaw', 'w') as log_fout: 
            subprocess.call(['atompaw'], stdin=input_fin, stdout=log_fout, env=env)

def write_QE_input(element_list,template_path):
    """
    Write QE input file based on some template specified in template_path
    """
    for elem in element_list:
        template_file = os.path.join(template_path, elem+'.relax.template')
        new_input_file = elem+'.relax.in'
        shutil.copy(template_file,new_input_file)

def run_QE(element_list):
    """
    Run QE relaxation using (elem).relax.in
    """
    for elem in element_list:
        os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py pw.x '+elem+'.relax.in -MPICORES 4')

def get_lattice_constant(element_list,lat_type):
    """
    Get relaxed lattice constant from QE run
    """
    for elem in element_list:
        qe_reader_path = os.path.join(fileutils.get_mmshare_scripts_dir(),'periodic_dft_gui_dir', 'qe2mae.py')
        qe_reader_mod = imputils.import_module_from_file(qe_reader_path)
        qe_reader = qe_reader_mod.QEOutputReader(elem+'.relax.out')
        struct = qe_reader.structs[qe_reader.final_struct_id]
        cparams = xtal.get_chorus_properties(struct)
        params = xtal.get_params_from_chorus(cparams)
        if lat_type == 'FCC':
            return math.sqrt(2)*params[0]
        else: ## Need to add in BCC stuff
            pass

def check_convergence(element_list):
    """
    Check if the QE run converged
    """
    check = True
    for elem in element_list:
        with open(elem+'.relax.out') as f:
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

if __name__=='__main__':
    main()

