import time
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
import json
from types import SimpleNamespace


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

def bad_run(diff_dict):
    """
    If something went wrong with the run, e.g., no .UPF file created or
    if running QE raised an error, set the objective function to 100
    """
    params, results = di.read_parameters_file('params.in','results.out')
    label_index = 1
    for elem in diff_dict.keys():
        for lat_type in diff_dict[elem].keys():
            for property in diff_dict[elem][lat_type].keys():
                label = 'obj_fn_'+str(label_index)
                results[label].function = 100.0
                label_index += 1
    results.write()

def compare_lat(AE_lat,cmpd,cmpd_lat_type):
    """
    Compute difference between AE and PAW lattice constants.
    For cubic systems, a primitive cell is assumed.
    For all other lattice types, a conventional cell is assumed.
    """
    QE_lat = get_lattice_constant(cmpd,cmpd_lat_type)
    lat_diff = 0
    if cmpd_lat_type in ['SC','FCC','BCC','ZB','per','RS','diamond','CsCl','HH']:
        lat_diff = abs(QE_lat-AE_lat)/AE_lat
        return lat_diff
    else:
        for (QE,AE) in zip(QE_lat,AE_lat):
            lat_diff += abs(QE-AE)/AE
        avg_lat_diff = lat_diff/len(AE_lat)
        return avg_lat_diff

def update_dakota(diff_dict):
    """
    Set the parameters and results files to be used by Dakota
    The objective function is equal to the difference between the lattice
    constants of AE calculations and PAW calculations performed here
    """
    params, results = di.read_parameters_file('params.in','results.out')
    label_index = 1
    for elem in diff_dict.keys():
        for lat_type in diff_dict[elem].keys():
            for property in diff_dict[elem][lat_type].keys():
                label = 'obj_fn_'+str(label_index)
                results[label].function = diff_dict[elem][lat_type][property]
                label_index += 1
    results.write()

def write_atompaw_input(elem,template_dir):
    """
    Write AtomPAW input file based on some template specified in template_dir
    """
    env = os.environ.copy()
    env['PATH'] = '/scr/fonari/dakota/bin:' + env['PATH']
    env['LD_LIBRARY_PATH'] = '/scr/fonari/dakota/bin:' + env['LD_LIBRARY_PATH']
    template_file = os.path.join(template_dir, elem+'.atompaw.template')
    new_input_file = elem+'.atompaw.in'
    subprocess.check_call(['run', 'dprepro.py', 'params.in', template_file, new_input_file], env=env)

def run_atompaw(elem):
    """
    Run AtomPAW using (elem).atompaw.in
    """
    env = os.environ.copy()
    env['PATH'] = '/scr/szymansk/atompaw-4.1.0.6/src:' + env['PATH']
    with open(elem+'.atompaw.in','r') as input_fin, open('log_atompaw', 'w') as log_fout: 
        subprocess.call(['atompaw'], stdin=input_fin, stdout=log_fout, env=env)

def run_QE(elem,lat_type,calc_type,template_dir):
    """
    Write and run QE using elem.lat_type.calc_type in template_dir.
    """
    if str(elem+'.'+lat_type+'.'+calc_type+'.out') not in os.listdir('.'):
        template_file = os.path.join(template_dir, elem+'.'+lat_type+'.'+calc_type+'.template')
        new_input_file = elem+'.'+lat_type+'.'+calc_type+'.in'
        shutil.copy(template_file,new_input_file)
        if calc_type == 'scf':
            update_structure(elem,lat_type,'scf')
        os.system('$SCHRODINGER/run periodic_dft_gui_dir/runner.py pw.x '+elem+'.'+lat_type+'.'+calc_type+'.in -MPICORES 4')

def get_lattice_constant(elem,lat_type):
    """
    Get relaxed lattice constant from QE run.
    Note some tolerance is allowed.
    """
    qe_reader_path = os.path.join(fileutils.get_mmshare_scripts_dir(),'periodic_dft_gui_dir', 'qe2mae.py')
    qe_reader_mod = imputils.import_module_from_file(qe_reader_path)
    qe_reader = qe_reader_mod.QEOutputReader(elem+'.'+lat_type+'.relax.out')
    struct = qe_reader.structs[qe_reader.final_struct_id]
    cparams = xtal.get_chorus_properties(struct)
    params = xtal.get_params_from_chorus(cparams)
    if lat_type in ['FCC','ZB','RS','diamond','HH']:
        return math.sqrt(2)*params[0]
    if lat_type == 'BCC':
        return (2./3.)*math.sqrt(3)*params[0]
    if lat_type in ['per','SC','CsCl']:
        return params[0]
    if lat_type == 'tetrag':
        unique_lat_list = sorted(unique(params[:3]))
        return unique_lat_list[0], unique_lat_list[1]
    if lat_type == 'ortho':
        unique_lat_list = sorted(params[:3])
        return unique_lat_list[0], unique_lat_list[1], unique_lat_list[2]
    if lat_type == 'hex' or lat_type == 'WZ':
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
    with open(elem+'.'+lat_type+'.'+calc_type+'.out') as qe_output:
        for line in qe_output:
            if 'convergence NOT' in line:
                return False
            if 'S matrix not positive definite' in line:
                return False
            if 'stopping ...' in line:
                return False
    return True

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

def unique(value_list): 
    """
    Get list of unique elements to be tested
    """
    try: ## if list of numbers
        value_list = [round(float(value),3) for value in value_list]
    except ValueError: ## if list of strings
        list = [str(value) for value in value_list]
    unique_list = []
    for value in value_list:
        if value not in unique_list:
            unique_list.append(value)
    return unique_list

def compare_atoms(elem,lat_type,template_dir):
    """
    Compare atomic positions of QE-relaxed structure
    and those of the AE-relaxed structure...
    """
    df_AE = pd.read_table(template_dir+'/AE_Struct.'+elem+'.'+lat_type,sep='\s+',header=None)
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
    mag = []
    with open(elem+'.'+lat_type+'.scf.out') as qe_output:
        for line in qe_output:
            if 'absolute' in line.split():
                mag.append(line.split()[3])
    return float(mag[-1])

def compare_mag_mom(elem,lat_type,template_dir):
    """
    Parse QE output (scf run) to obtain individual
    magnetic moments of atoms in given structure.
    Compare these with AE magnetic moments in AE_mag.
    """
    qe_mag_mom = []
    with open(elem+'.'+lat_type+'.scf.out') as qe_output:
        for line in qe_output:
            if 'magn:' in line.split():
                qe_mag_mom.append(line.split()[5])
    qe_mag_mom = [float(value) for value in qe_mag_mom]
    ae_mag_mom = np.loadtxt(template_dir+'/AE_mag.'+elem+'.'+lat_type)
    rel_diff = []
    for (qe_val,ae_val) in zip(qe_mag_mom,ae_mag_mom):
        if float(ae_val) != 0.0:
            rel_diff.append(abs((qe_val-ae_val)/ae_val))
        else:
            pass
    net_diff = sum(rel_diff)/len(rel_diff)
    return float(net_diff)

def get_gap(elem,lat_type):
    """
    Parse QE output (scf run) to obtain band gap.
    Note that unoccupied bands need to be included
    and occupations need to be fixed in the scf run.
    """
    with open(elem+'.'+lat_type+'.scf.out') as qe_output:
        for line in qe_output:
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

    with open(elem+'.'+lat_type+'.relax.in') as qe_output:
        for line in qe_output:
            if 'nat=' in line:
                num_atoms = float(line.split('=')[1][:-1])
    volume = popt[0]/num_atoms ## A^3/atom
    bulk = popt[1]*160.2 ## GPa
    B_prime = popt[2] ## Dimensionless
    eos_file = open('QE_EOS.txt','w+')
    eos_file.write(str(volume)+' '+str(bulk)+' '+str(B_prime))
    eos_file.close()
    return float(volume), float(bulk), float(B_prime)

def run_scale_lat(elem,lat_type,template_dir):
    """
    Read in relaxed cell parameters from (elem).(lat_type).relax.out,
    scale this lattice constant from 94% to 106% (7 values created),
    write input and run QE at each value, write corresponding volumes
    and energies into E_V.txt (units of Bohr^3 and Ry^3)
    """
    scale_num = [0.94,0.96,0.98,1.0,1.02,1.04,1.06]
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
        all_energies = []
        with open(relax_file[:-2]+'out') as qe_output:
            for line in qe_output:
                if '!    total energy              =' in line:
                    all_energies.append(line.split()[4])
        energies.append(all_energies[-1])
        os.chdir('../')
        folder_index += 1
    ev_file = open('E_V.txt','w+')
    for (e,v) in zip(energies,volumes):
        ev_file.write(str(e)+' '+str(v)+'\n')
    ev_file.close()

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

def compare_phonon(elem,lat_type,template_dir):
    """
    Parse optical phonon frequencies from QE run
    and compare with AE frequencies
    """
    index = 0
    with open('dynmat.out') as dyn_output:
        for line in dyn_output:
            if 'mode' in line:
                mode_line = index
            index += 1
    with open(elem+'.'+lat_type+'.scf.in') as qe_input:
        for line in qe_input:
            if 'nat=' in line:
                num_atoms = int(line.split('=')[1][:-1])
    num_freq = 3*num_atoms
    freq_index = mode_line + 1
    freq = []
    for i in range(num_freq):
        freq.append(lines[freq_index].split()[2])
        freq_index += 1
    QE_freq = sorted([float(value) for value in freq])
    AE_freq = sorted(np.loadtxt(template_dir+'/AE_freq.'+elem+'.'+lat_type))
    rel_diff = []
    for (QE,AE) in zip(QE_freq,AE_freq):
        rel_diff.append(abs(QE-AE))
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
            try: ## if some scaling constant for cell
                alat = float(split_line[1][1:-2])
            except: ## else in angstroms
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
    with open(elem+'.'+lat_type+'.'+calc_type+'.in') as qe_input:
        lines = qe_input.readlines()
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
    index = 0
    with open(elem+'.'+lat_type+'.relax.out') as qe_output:
        for line in qe_output:
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
    with open(elem+'.'+lat_type+'.relax.in') as qe_input:
        lines = qe_input.readlines()
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

def run_phonon(cmpd,cmpd_lat_type,template_dir):
    """
    Run QE phonon calculations using ph.x and dynmat.x.
    """
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

def update_best_result(diff_dict):
    """
    Parse dakota results and check overall fitness with
    respect to previous best solution. If current solution
    is better, replace old solution with current one.
    Note that fitness is normalized per the highest error
    for a given objective function.
    """
    obj_fn_list = []
    obj_fn_labels = []
    for formula in diff_dict.keys():
        for lat_type in diff_dict[formula].keys():
            for property in diff_dict[formula][lat_type].keys():
                obj_fn_labels.append(formula+'_'+lat_type+'_'+property)
                obj_fn_list.append(diff_dict[formula][lat_type][property])
    f = open('OBJ_FN','w+')
    for (value,label) in zip(obj_fn_list,obj_fn_labels):
        value = round(float(value),6)
        if 'lattice_constant' in label:
            value = value*100
            f.write(label+':  '+str(value)+'%\n')
        if 'band_gap' in label:
            f.write(label+':  '+str(value)+' eV\n')
        if 'log' in label:
            f.write(label+':  '+str(value)+'\n')
    f.close()
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
            if 'Max_Error' not in file:
                files_to_del.append(file)
        for filename in files_to_del:
            os.remove('../Best_Solution/'+filename)
        copyfile('OBJ_FN','../Best_Solution/results.out')
        for (file_1,file_2) in zip(atompaw_files,UPF_files):
            copyfile(file_1,'../Best_Solution/'+file_1)
            copyfile(file_2,'../Best_Solution/'+file_2)
        f = open('../Best_Solution/rms_error','w+')
        f.write(str(rms_error))
        f.close()

def parse_elems(formula):
    """
    Parse compound formula to obtain constituent elements
    """
    letters_only = ''.join([letter for letter in formula if not letter.isdigit()])
    index = -1
    elems = []
    for letter in letters_only:
        if letter.isupper():
            elems.append(letter)
            index += 1
        else:
            elems[index] += letter

    return elems

def get_element_info():
    """
    Read in AE data for elemental lattice constants
    """
    elemental_data = {}
    df_FCC = pd.read_table(elem_template_dir+'/WIEN2k_FCC',sep='\s+',header=None)
    for (elem,lat_const) in zip(df_FCC[0],df_FCC[1]):
        elemental_data[elem] = {}
        elemental_data[elem]['FCC'] = lat_const
    df_BCC = pd.read_table(elem_template_dir+'/WIEN2k_BCC',sep='\s+',header=None)
    for (elem,lat_const) in zip(df_BCC[0],df_BCC[1]):
        elemental_data[elem]['BCC'] = lat_const
    elemental_data['N'] = {}
    elemental_data['N']['SC'] = 6.1902
    elemental_data['P'] = {}
    elemental_data['P']['ortho'] = [3.304659,4.573268,11.316935]
    return elemental_data

def test_element_list(elem_list,template_dir):
    """
    Test the atompaw generation, compare QE with AE lattice constants,
    and compare QE/AE log derivatives
    """
    elem_diff_dict = {}
    for elem in elem_list:
        elem_diff_dict[elem] = {}
        elem_diff_dict[elem]['elemental'] = {}
        elem_diff_dict[elem]['elemental']['log'] = {}
        for lat_type in ['FCC','BCC']:
            elem_diff_dict[elem][lat_type] = {}
            elem_diff_dict[elem][lat_type]['lattice_constant'] = {}
    for elem in elem_list:
        os.mkdir(elem)
        copyfile('params.in',elem+'/params.in')
        with fileutils.chdir(elem):
            write_atompaw_input(elem,template_dir)
            run_atompaw(elem)
            if not check_UPF():
                return elem_diff_dict, True
            copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
            elem_diff_dict[elem]['elemental']['log'] = compare_log()
            for lat_type in ['FCC','BCC']:
                run_QE(elem,lat_type,'relax',template_dir)
                if not check_convergence(elem,lat_type,'relax'):
                    return elem_diff_dict, True
                elemental_data = get_element_info()
                ae_lat = elemental_data[elem][lat_type]
                elem_diff_dict[elem][lat_type]['lattice_constant'] = compare_lat(ae_lat,elem,lat_type)
    return elem_diff_dict, False

def test_property(cmpd,lat_type,property,ae_data,template_dir):
    """
    General function to calculate a property with QE
    and compare with corresponding AE value
    """
    run_QE(cmpd,lat_type,'relax',template_dir)
    if not check_convergence(cmpd,lat_type,'relax'):
        return None, True
    if property == 'lattice_constant':
        return compare_lat(ae_data,cmpd,lat_type), False
    if property in ['eos','bulk_modulus']:
        run_scale_lat(cmpd,lat_type,template_dir)
        V0, QE_bulk, B_prime = get_bulk(cmpd,cmpd_lat_type)
        if property == 'eos':
            qe_data = np.loadtxt('QE_EOS.txt')
            qe_eos = {'element': [cmpd], 'V0': [qe_data[0]], 'B0': [qe_data[1]], 'BP': [qe_data[2]]}
            ae_eos = {'element': [cmpd], 'V0': [ae_data[0]], 'B0': [ae_data[1]], 'BP': [ae_data[2]]}
            delta_factor = calcDelta(qe_eos,ae_eos,[cmpd],False)
            return delta_factor, False
        if property == 'bulk_modulus':
            bulk_diff = abs(QE_bulk - ae_data)/ae_data
    if property == 'phonon_frequency': ## Probably just define AE freq in .json file
        run_QE(cmpd,lat_type,'scf',template_dir)
        if not check_convergence(cmpd,lat_type,'scf'):
            return None, True
        run_phonon(cmpd,lat_type,template_dir)
        phonon_diff = compare_phonon(cmpd,lat_type,template_dir)
        return phonon_diff, False
    if property == 'atomic_positions': ## Still use file?
        return compare_atoms(cmpd,lat_type,template_dir), False
    if property == 'band_gap':
        run_QE(cmpd,lat_type,'scf',template_dir)
        if not check_convergence(cmpd,lat_type,'scf'):
            return None, True
        qe_gap = get_gap(cmpd,lat_type)
        return abs(ae_data-qe_gap), False ## or maybe in eV?
    if property == 'magnetization':
        run_QE(cmpd,lat_type,'scf',template_dir)
        if not check_convergence(cmpd,lat_type,'scf'):
            return None, True
        qe_mag = get_mag(cmpd,lat_type)
        return abs(ae_mag-qe_mag)/ae_mag, False
    if property == 'magnetic_moment':
        run_QE(cmpd,lat_type,'scf',template_dir)
        if not check_convergence(cmpd,lat_type,'scf'):
            return None, True
        return compare_mag_mom(cmpd,lat_type,template_dir), False

def form_cmpd_dict(cmpd_list):
    """
    Constructs empty dictionary of correct length for testing
    """
    cmpd_diff_dict = {}
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        lat_type = cmpd['lattice_type']
        property_list = [property for property in cmpd if property not in ['formula','lattice_type']]
        cmpd_diff_dict[formula] = {}
        cmpd_diff_dict[formula][lat_type] = {}
        for property in property_list:
            cmpd_diff_dict[formula][lat_type][property] = {}
    return cmpd_diff_dict

def test_cmpd_list(cmpd_list,cmpd_diff_dict,cmpd_template_dir):
    """
    Perform property tests for all compounds
    """
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        lat_type = cmpd['lattice_type']
        property_list = [property for property in cmpd if property not in ['formula','lattice_type']]
        for property in property_list:
            ae_value = cmpd[property]
            cmpd_diff_dict[formula][lat_type][property], error_check = test_property(formula,lat_type,property,ae_value,cmpd_template_dir)
            if error_check:
                bad_run(cmpd_diff_dict)
                return cmpd_diff_dict, True
    return cmpd_diff_dict, False
