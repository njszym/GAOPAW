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


## Optimize the number of points in the logarithmic radial grid
## Decrease loggrid in intervals of 100
## At each value, run AtomPAW to generate .UPF and use this to run QE on FCC and BCC structures
## If differences between original and new values (energy and lattice constants) become too high, stop

def read_output(elem):
    with open(elem) as f:
        lines = f.readlines()
    for line in lines:
        if 'valence' in line.split():
            pseudized = float(line.split()[3])
        if 'Valence' in line.split():
            try:
                AE = float(line.split()[2])
            except:
                pass
    return [AE,pseudized]

def run_atompaw(elem):
    """
    Run AtomPAW using (elem).atompaw.in
    """
    env = os.environ.copy()
    env['PATH'] = '/scr/szymansk/atompaw-4.1.0.5/src:' + env['PATH']
    with open(elem+'.atompaw.in','r') as input_fin, open('log_atompaw', 'w') as log_fout: 
        subprocess.call(['atompaw'], stdin=input_fin, stdout=log_fout, env=env)

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

def get_walltime(elem,lat_type):
    pwscf_lines = []
    with open(elem+'.'+lat_type+'.relax.out') as f:
        lines = f.readlines()
    for line in lines:
        if 'PWSCF' in line.split():
            pwscf_lines.append(line.split())
    time = pwscf_lines[-1][2]
    return time

## Tolerances listed below may need to be tuned

elem = sys.argv[-1]
ediff = 0
etol = 1e-6
FCC_diff = 0
FCC_tol = 0.0001
BCC_diff = 0
BCC_tol = 0.0001
run_atompaw(elem)
benchmark_e = read_output(elem)[1]
run_QE(elem,'FCC')
benchmark_FCC = get_lattice_constant(elem,'FCC')
run_QE(elem,'BCC')
benchmark_BCC = get_lattice_constant(elem,'BCC')
differences = []

while (ediff < etol) and (FCC_diff < FCC_tol) and (BCC_diff < BCC_tol):
    filename = elem+'.atompaw.in'
    with open(filename) as f:
        lines = f.readlines()
    log_line = (lines[1]).split()
    index = 1
    for word in log_line:
        if word == 'loggrid':
            num_pts = float(log_line[index])
            break
        else:
            index += 1
    num_pts -= 100
    log_line[index] = str(num_pts)
    new_line = ''
    for word in log_line:
        var = word+' '
        new_line += var
    lines[1] = new_line+'\n'
    with open(filename,'w') as f:
        for line in lines:
            f.write(line)
    run_atompaw(elem)
    ediff = abs(benchmark_e - read_output(elem)[1])/benchmark_e
    run_QE(elem,'FCC')
    new_lat_FCC = get_lattice_constant(elem,'FCC')
    FCC_diff = abs(benchmark_FCC - new_lat_FCC)/benchmark_FCC
    FCC_walltime = get_walltime(elem,'FCC')
    run_QE(elem,'BCC')
    new_lat_BCC = get_lattice_constant(elem,'BCC')
    BCC_diff = abs(benchmark_BCC - new_lat_BCC)/benchmark_BCC
    BCC_walltime = get_walltime(elem,'BCC')
    differences.append([num_pts,ediff,FCC_diff,BCC_diff,FCC_walltime,BCC_walltime])

num_pts += 100
log_line[index] = str(num_pts)
new_line = ''
for word in log_line:
    var = word+' '
    new_line += var
lines[1] = new_line+'\n'
with open(filename,'w') as f:
    for line in lines:
        f.write(line)

print (differences)
