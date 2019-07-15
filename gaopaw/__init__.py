import collections
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
import glob


def parse_num_objs(working_dir):
    """
    Read in number of objective functions defined in the dakota.in file
    """
    with open(working_dir+'/../dakota.in') as dakota_input:
        for line in dakota_input:
            if 'num_objective_functions' in line:
                num_objs = line.split()[2]
    return int(num_objs)

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

def unique(value_list):
    """
    Get list of unique elements to be tested
    """
    try:
        value_list = [round(float(value),3) for value in value_list]
    except ValueError:
        list = [str(value) for value in value_list]
    unique_list = []
    for value in value_list:
        if value not in unique_list:
            unique_list.append(value)
    return unique_list

def get_num_objs(cmpd_list,element_list):
    """
    Parse input.json and return total number of objective functions
    """
    num_elems = len(element_list)
    if num_elems == 1 and len(cmpd_list[0].__dict__.keys()) == 1:
        if element_list[0] in ['N','P']:
            return 2
        else:
            return 3
    cmpd_diff_dict = form_cmpd_dict(cmpd_list)
    num_properties = 0
    for cmpd in cmpd_diff_dict.keys():
        for lat_type in cmpd_diff_dict[cmpd].keys():
            for property in cmpd_diff_dict[cmpd][lat_type].keys():
                num_properties += 1
    num_obj_fns = 0
    for elem in element_list:
        if elem in ['N','P']:
            num_obj_fns += 2
        else:
            num_obj_fns += 3
    num_obj_fns += num_properties
    return int(num_obj_fns)

def form_cmpd_dict(cmpd_list):
    """
    Constructs empty dictionary of correct length for testing
    """
    cmpd_diff_dict = {}
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        lat_type = cmpd['lattice_type']
        property_list = [property for property in \
            cmpd if property not in ['formula','lattice_type']]
        cmpd_diff_dict[formula] = {}
        cmpd_diff_dict[formula][lat_type] = {}
        for property in property_list:
            cmpd_diff_dict[formula][lat_type][property] = {}
    return cmpd_diff_dict

def get_element_info(template_dir):
    """
    Read in AE data for elemental lattice constants
    """
    elemental_data = {}
    df_FCC = pd.read_table(template_dir+'/WIEN2k_FCC',sep='\s+',header=None)
    for (elem,lat_const) in zip(df_FCC[0],df_FCC[1]):
        elemental_data[elem] = {}
        elemental_data[elem]['FCC'] = lat_const
    df_BCC = pd.read_table(template_dir+'/WIEN2k_BCC',sep='\s+',header=None)
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
        copyfile('params.in',elem+'/params.in')
        with fileutils.chdir(elem):
            write_atompaw_input(elem,template_dir)
            run_atompaw(elem)
            if not check_upf():
                return elem_diff_dict, True
            copyfile(elem+'.GGA-PBE-paw.UPF','../'+elem+'.GGA-PBE-paw.UPF')
            elem_diff_dict[elem]['elemental']['log'] = compare_log()
            if elem in ['N','P']:
                if elem == 'N':
                    run_QE(elem,'SC','relax',template_dir)
                    if not check_convergence(elem,'SC','relax'):
                        return elem_diff_dict, True
                    elem_diff_dict[elem]['SC']['atomic_positions'] = \
                         compare_atoms(elem,'SC',template_dir)
                    return elem_diff_dict, False
                if elem == 'P':
                    run_QE(elem,'ortho','relax',template_dir)
                    if not check_convergence(elem,'ortho','relax'):
                        return elem_diff_dict, True
                    elemental_data = get_element_info(template_dir)
                    ae_lat = elemental_data[elem]['ortho']
                    elem_diff_dict[elem]['ortho']['lattice_constant'] = \
                        compare_lat(ae_lat,elem,'ortho')
            else:
                for lat_type in ['FCC','BCC']:
                    run_QE(elem,lat_type,'relax',template_dir)
                    if not check_convergence(elem,lat_type,'relax'):
                        return elem_diff_dict, True
                    elemental_data = get_element_info(template_dir)
                    ae_lat = elemental_data[elem][lat_type]
                    elem_diff_dict[elem][lat_type]['lattice_constant'] = \
                        compare_lat(ae_lat,elem,lat_type)
    return elem_diff_dict, False

def test_cmpd_list(cmpd_list,cmpd_diff_dict,cmpd_template_dir):
    """
    Perform property tests for all compounds
    """
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        lat_type = cmpd['lattice_type']
        property_list = [property for property in \
            cmpd if property not in ['formula','lattice_type']]
        for property in property_list:
            ae_value = cmpd[property]
            cmpd_diff_dict[formula][lat_type][property], error_check = \
                test_property(formula,lat_type,property,ae_value,cmpd_template_dir)
            if error_check:
                bad_run(cmpd_diff_dict)
                return cmpd_diff_dict, True
    return cmpd_diff_dict, False

def merge_dicts(dct, merge_dct):
    """
    Recursively merge two dictionaries
    """
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict) and \
         isinstance(merge_dct[k], collections.Mapping)):
            merge_dicts(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]
    return dct

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
        V0, QE_bulk, B_prime = get_bulk(cmpd,lat_type)
        if property == 'eos':
            qe_data = np.loadtxt('QE_EOS.txt')
            qe_eos = {'element': [cmpd], 'V0': [qe_data[0]],
                'B0': [qe_data[1]], 'BP': [qe_data[2]]}
            ae_eos = {'element': [cmpd], 'V0': [ae_data[0]],
                'B0': [ae_data[1]], 'BP': [ae_data[2]]}
            delta_factor = calcDelta(qe_eos,ae_eos,[cmpd])
            return delta_factor, False
        if property == 'bulk_modulus':
            bulk_diff = abs(QE_bulk - ae_data)/ae_data
    if property == 'phonon_frequency':
        run_QE(cmpd,lat_type,'scf',template_dir)
        if not check_convergence(cmpd,lat_type,'scf'):
            return None, True
        run_phonon(cmpd,lat_type,template_dir)
        phonon_diff = compare_phonon(cmpd,lat_type,ae_data,template_dir)
        return phonon_diff, False
    if property == 'atomic_positions':
        return compare_atoms(cmpd,lat_type,template_dir), False
    if property == 'band_gap':
        run_QE(cmpd,lat_type,'scf',template_dir)
        if not check_convergence(cmpd,lat_type,'scf'):
            return None, True
        qe_gap = get_gap(cmpd,lat_type)
        return abs(ae_data-qe_gap), False
    if property == 'magnetization':
        run_QE(cmpd,lat_type,'scf',template_dir)
        if not check_convergence(cmpd,lat_type,'scf'):
            return None, True
        qe_mag = get_mag(cmpd,lat_type)
        return abs(ae_data-qe_mag)/ae_data, False
    if property == 'magnetic_moment':
        run_QE(cmpd,lat_type,'scf',template_dir)
        if not check_convergence(cmpd,lat_type,'scf'):
            return None, True
        return compare_mag_mom(cmpd,lat_type,ae_data,template_dir), False
    raise ValueError('Your property, '+property+', is not defined')

def check_upf():
    """
    Check if a .UPF file was succesfully created by AtomPAW
    """
    if len([fname for fname in glob.iglob('*UPF')]) != 0:
        return True
    return False

def bad_run(num_obj_fns):
    """
    If something went wrong with the run, e.g., no .UPF file created or
    if running QE raised an error, set the objective function to 100
    """
    params, results = di.read_parameters_file('params.in','results.out')
    for num in range(1,num_obj_fns+1):
        label = 'obj_fn_'+str(num)
        results[label].function = 100.0
    results.write()

def compare_lat(ae_lat,cmpd,lat_type):
    """
    Compute difference between AE and PAW lattice constants.
    For cubic systems, a primitive cell is assumed.
    For all other lattice types, a conventional cell is assumed.
    """
    qe_lat = get_lattice_constant(cmpd,lat_type)
    if isinstance(ae_lat,list) == False:
        return abs(qe_lat-ae_lat)/ae_lat
    lat_diff = 0
    for (qe_val,ae_val) in zip(qe_lat,ae_lat):
        lat_diff += abs(qe_val-ae_val)/ae_val
    return lat_diff/len(ae_lat)

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
    template_file = os.path.join(template_dir, elem+'.atompaw.template')
    new_input_file = elem+'.atompaw.in'
    subprocess.check_call(['run', 'dprepro.py', 'params.in', template_file, 
        new_input_file], env=os.environ.copy())

def run_atompaw(elem):
    """
    Run AtomPAW using (elem).atompaw.in
    """
    with open(elem+'.atompaw.in','r') as input_fin, open('log_atompaw', 'w') as log_fout: 
        subprocess.call(['atompaw'], stdin=input_fin, stdout=log_fout, 
            env=os.environ.copy())

def run_QE(cmpd,lat_type,calc_type,template_dir):
    """
    Write and run QE using cmpd.lat_type.calc_type in template_dir.
    """
    if os.path.exists(cmpd+'.'+lat_type+'.'+calc_type+'.out'):
        return
    template_file = os.path.join(template_dir, cmpd+'.'+lat_type+'.'+calc_type+'.template')
    new_input_file = cmpd+'.'+lat_type+'.'+calc_type+'.in'
    shutil.copy(template_file,new_input_file)
    if calc_type == 'scf':
        update_structure(cmpd,lat_type,'scf')
    qe_input = cmpd+'.'+lat_type+'.'+calc_type+'.in'
    subprocess.call(['run','periodic_dft_gui_dir/runner.py','pw.x',qe_input,
        '-MPICORES','4'], env=os.environ.copy())

def get_lattice_constant(cmpd,lat_type):
    """
    Get relaxed lattice constant from QE run.
    Note some tolerance is allowed.
    """
    qe_reader_path = os.path.join(fileutils.get_mmshare_scripts_dir(),
        'periodic_dft_gui_dir', 'qe2mae.py')
    qe_reader_mod = imputils.import_module_from_file(qe_reader_path)
    qe_reader = qe_reader_mod.QEOutputReader(cmpd+'.'+lat_type+'.relax.out')
    struct = qe_reader.structs[qe_reader.final_struct_id]
    cparams = xtal.get_chorus_properties(struct)
    params = xtal.get_params_from_chorus(cparams)
    if lat_type in ['FCC','ZB','RS','diamond','HH']:
        return math.sqrt(2)*params[0]
    if lat_type == 'BCC':
        return (2./3.)*math.sqrt(3)*params[0]
    if lat_type in ['per','SC','CsCl']:
        return params[0]
    if lat_type in ['tetrag','hex','WZ']:
        unique_lat_list = sorted(unique(params[:3]))
        return unique_lat_list[0], unique_lat_list[1]
    if lat_type == 'ortho':
        unique_lat_list = sorted(params[:3])
        return unique_lat_list[0], unique_lat_list[1], unique_lat_list[2]
    if lat_type == 'rhomb':
        lat = params[0]
        angle = params[4]
        return lat, angle
    if lat_type == 'monoclin':
        unique_lat_list = sorted(unique(params[:3]))
        for value in params[3:]:
            if (value > 90.01) or (value < 89.99):
                angle = value
        return unique_lat_list[0], unique_lat_list[1], unique_lat_list[2], angle
    if lat_type == 'triclin':
        unique_lat_list = sorted(unique(params[:3]))
        unique_angle_list = sorted(unique(params[3:]))
        return unique_lat_list + unique_angle_list

def check_convergence(cmpd,lat_type,calc_type):
    """
    Check if the QE run converged
    """
    with open(cmpd+'.'+lat_type+'.'+calc_type+'.out') as qe_output:
        for line in qe_output:
            if any(error in line for error in \
            ['convergence NOT','S matrix not positive definite','stopping ...']):
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
    log_derivs = [fname for fname in glob.iglob('logd*')]
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

def compare_atoms(cmpd,lat_type,template_dir):
    """
    Compare atomic positions of QE-relaxed structure
    and those of the AE-relaxed structure...
    """
    df_ae = pd.read_table(template_dir+'/AE_Struct.'+cmpd+'.'+lat_type,
        sep='\s+',header=None)
    df_ae = df_ae.drop(0,1).transpose()
    qe_reader_path = os.path.join(fileutils.get_mmshare_scripts_dir(),
        'periodic_dft_gui_dir', 'qe2mae.py')
    qe_reader_mod = imputils.import_module_from_file(qe_reader_path)
    qe_reader = qe_reader_mod.QEOutputReader(cmpd+'.'+lat_type+'.relax.out')
    struct = qe_reader.structs[qe_reader.final_struct_id]
    df_qe = struct.getXYZ()
    distance_list = []
    for index in range(len(df_ae.keys())):
        ae_position = np.array(df_ae[index])
        qe_position = np.array(df_qe[index])
        distance_list.append(np.linalg.norm(ae_position-qe_position))
    return sum(distance_list)

def get_mag(cmpd,lat_type):
    """
    Parse QE output (scf run) to obtain total magnetization
    """
    mag = []
    with open(cmpd+'.'+lat_type+'.scf.out') as qe_output:
        for line in qe_output:
            if 'absolute' in line.split():
                mag.append(line.split()[3])
    return float(mag[-1])

def compare_mag_mom(cmpd,lat_type,ae_mag_mom,template_dir):
    """
    Parse QE output (scf run) to obtain individual
    magnetic moments of atoms in given structure.
    Compare these with AE magnetic moments in AE_mag.
    """
    qe_mag_mom = []
    with open(cmpd+'.'+lat_type+'.scf.out') as qe_output:
        for line in qe_output:
            if 'magn:' in line.split():
                qe_mag_mom.append(line.split()[5])
    qe_mag_mom = [float(value) for value in qe_mag_mom]
    rel_diff = []
    for (qe_val,ae_val) in zip(qe_mag_mom,ae_mag_mom):
        if float(ae_val) != 0.0:
            rel_diff.append(abs((qe_val-ae_val)/ae_val))
        else:
            pass
    net_diff = sum(rel_diff)/len(rel_diff)
    return float(net_diff)

def get_gap(cmpd,lat_type):
    """
    Parse QE output (scf run) to obtain band gap.
    Note that unoccupied bands need to be included
    and occupations need to be fixed in the scf run.
    """
    with open(cmpd+'.'+lat_type+'.scf.out') as qe_output:
        for line in qe_output:
            if 'highest' and 'lowest' in line.split():
                band_gap = (float(line.split()[7]) - float(line.split()[6]))
    return float(band_gap)

def birch_murnaghan(vol, vol_equil, bulk, bulk_prime, energy_equil):
    """
    3rd order Birch-Murnaghan equation of state, in the energy-volume form
    """
    vol = np.array(vol)
    return energy_equil + 9 * vol_equil * bulk / 16. * (
        ((vol_equil / vol) ** (2 / 3.) - 1) ** 3 * bulk_prime +
        ((vol_equil / vol) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (vol_equil / vol) ** (2 / 3.)))

def get_bulk(cmpd,lat_type):
    """
    Reads in energy-volume data from E_V.txt and calculates:
    equilibrium volume, bulk modulus, and dB/dP...
    i.e., fit to Birch Murnaghan equation of state
    """
    df = pd.read_table('E_V.txt',sep='\s+',header=None)
    energy = list(df[0])
    energy = np.array([13.6056980659*value for value in energy]) ## Ry to eV
    volume = list(df[1])
    volume = np.array([0.14818453429566825*value for value in volume]) ## bohr^3 to A^3
    initial_parameters = [volume.mean(), 2.5, 4, energy.mean()]
    fit_eqn = eval('birch_murnaghan')
    popt, pcov = cf(fit_eqn, volume, energy, initial_parameters)
    with open(cmpd+'.'+lat_type+'.relax.in') as qe_output:
        for line in qe_output:
            if 'nat=' in line:
                num_atoms = float(line.split('=')[1][:-1])
    volume = popt[0]/num_atoms ## A^3/atom
    bulk = popt[1]*160.2 ## GPa
    bulk_prime = popt[2] ## Dimensionless
    with open('QE_EOS.txt','w+') as eos_file:
        eos_file.write(str(volume)+' '+str(bulk)+' '+str(bulk_prime))
    return volume, bulk, bulk_prime

def run_scale_lat(cmpd,lat_type,template_dir):
    """
    Read in relaxed cell parameters from cmpd.lat_type.relax.out,
    scale this lattice constant from 94% to 106% (7 values created),
    write input and run QE at each value, write corresponding volumes
    and energies into E_V.txt (units of Bohr^3 and Ry^3)
    """
    scale_num = [0.94,0.96,0.98,1.0,1.02,1.04,1.06]
    relax_file = cmpd+'.'+lat_type+'.relax.in'
    UPF_files = [fname for fname in glob.iglob('*UPF')]
    energies = []
    volumes = []
    folder_index = 1
    for value in scale_num:
        folder = cmpd+'_'+str(folder_index)
        new_cell_params = scale_cell(cmpd,lat_type,value)
        new_cell_matrix = np.matrix(new_cell_params)
        volumes.append(np.linalg.det(new_cell_matrix))
        os.mkdir(folder)
        for file in UPF_files:
            copyfile(file,folder+'/'+file)
        copyfile(relax_file,folder+'/'+relax_file)
        copyfile(relax_file[:-2]+'out',folder+'/'+relax_file[:-2]+'out')
        with fileutils.chdir(folder):
            update_structure(cmpd,lat_type,'relax')
            write_cell(cmpd,lat_type,new_cell_params)
            subprocess.call(['run','periodic_dft_gui_dir/runner.py','pw.x',relax_file,
                '-MPICORES','4'], env=os.environ.copy())
            all_energies = []
            with open(relax_file[:-2]+'out') as qe_output:
                for line in qe_output:
                    if '!    total energy              =' in line:
                        all_energies.append(line.split()[4])
            energies.append(all_energies[-1])
        folder_index += 1
    with open('E_V.txt','w+') as ev_file:
        for (e,v) in zip(energies,volumes):
            ev_file.write(str(e)+' '+str(v)+'\n')

def calcDelta(data_f, data_w, eloverlap):
    """
    Calculate the Delta using the data in data_f, data_w on
    element in eloverlap...taken from Delta Package...
    data_f: QE data (dict) calculated with PAWs
    data_w: AE data (dict) calculated with WIEN2k
    eloverlap: names (list) of elements/compounds
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
    delta = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi))
    return delta

def compare_phonon(cmpd,lat_type,ae_freq,template_dir):
    """
    Parse optical phonon frequencies from QE run
    and compare with AE frequencies
    """
    with open(cmpd+'.'+lat_type+'.scf.in') as qe_input:
        for line in qe_input:
            if 'nat=' in line:
                num_atoms = int(line.split('=')[1][:-1])
    with open('dynmat.out') as dyn_output:
        lines = dyn_output.readlines()
    index = 0
    for line in lines:
        if 'mode' in line:
            mode_line = index
        index += 1
    num_freq = 3*num_atoms
    freq_index = mode_line + 1
    freq = []
    for i in range(num_freq):
        freq.append(lines[freq_index].split()[2])
        freq_index += 1
    qe_freq = sorted([float(value) for value in freq])
    rel_diff = []
    for (qe_val,ae_val) in zip(qe_freq,ae_freq):
        rel_diff.append(abs(ae_val-qe_val))
    net_diff = sum(rel_diff)/len(rel_diff)
    return net_diff

def update_structure(cmpd,lat_type,calc_type):
    """
    Parse equilibrium structure from completed relaxation
    and update the corresponding calculation input file.
    """
    with open(cmpd+'.'+lat_type+'.relax.out') as f:
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
            if 'alat' in line: ## assumes Bohr
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
        v = [alat*value for value in v]
        vectors.append(v)
    cell_header = 'CELL_PARAMETERS bohr\n'
    v1 = str(vectors[0][0])+' '+str(vectors[0][1])+' '+str(vectors[0][2])+'\n'
    v2 = str(vectors[1][0])+' '+str(vectors[1][1])+' '+str(vectors[1][2])+'\n'
    v3 = str(vectors[2][0])+' '+str(vectors[2][1])+' '+str(vectors[2][2])+'\n'
    with open(cmpd+'.'+lat_type+'.'+calc_type+'.in') as qe_input:
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
    with open(cmpd+'.'+lat_type+'.'+calc_type+'.in','w+') as qe_file:
        for line in orig_struct:
            qe_file.write(line)
        qe_file.write(coords_header)
        for line in coords:
            qe_file.write(line)
        qe_file.write(cell_header+v1+v2+v3)

def scale_cell(cmpd,lat_type,scale_factor):
    """
    Scale cell volume according to scale_factor
    """
    with open(cmpd+'.'+lat_type+'.relax.out') as qe_output:
        lines = qe_output.readlines()
    index = 0
    for line in lines:
        if 'CELL_PARAMETERS' in line:
            cell_index = [index+1,index+2,index+3]
            split_line = line.split('=')
            if 'alat' in line: ## assumes Bohr
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

def write_cell(cmpd,lat_type,cell):
    """
    Write given cell to QE relaxation input
    """
    vectors = np.array(cell)
    v1 = str(vectors[0][0])+' '+str(vectors[0][1])+' '+str(vectors[0][2])+'\n'
    v2 = str(vectors[1][0])+' '+str(vectors[1][1])+' '+str(vectors[1][2])+'\n'
    v3 = str(vectors[2][0])+' '+str(vectors[2][1])+' '+str(vectors[2][2])+'\n'
    cell_header = 'CELL_PARAMETERS bohr\n'
    with open(cmpd+'.'+lat_type+'.relax.in') as qe_input:
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
    with open(cmpd+'.'+lat_type+'.relax.in','w+') as qe_file:
        for line in orig_struct:
            qe_file.write(line)
        qe_file.write(cell_header+v1+v2+v3)

def run_phonon(cmpd,cmpd_lat_type,template_dir):
    """
    Run QE phonon calculations using ph.x and dynmat.x.
    """
    copyfile(template_dir+'/phonon.in','./phonon.in')
    if os.path.exists('phonon.save'):
        os.remove('phonon.out')
        shutil.rmtree('phonon.save')
        os.remove('phonon.save.qegz')
    scf_savefile = cmpd+'.'+cmpd_lat_type+'.scf.save.qegz'
    subprocess.call(['run','periodic_dft_gui_dir/runner.py','ph.x','phonon.in',
        '-input_save',scf_savefile,'-MPICORES','4'], env=os.environ.copy())
    copyfile(template_dir+'/dynmat.in','./dynmat.in')
    if os.path.exists('dynmat.save'):
        os.remove('dynmat.out')
        shutil.rmtree('dynmat.save')
        os.remove('dynmat.save.qegz')
    subprocess.call(['run','periodic_dft_gui_dir/runner.py','dynmat.x','dynmat.in',
        '-input_save','phonon.save.qegz','-MPICORES','4'], env=os.environ.copy())

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
    with open('OBJ_FN','w+') as obj_file:
        for (value,label) in zip(obj_fn_list,obj_fn_labels):
            value = round(float(value),6)
            if ('lattice_constant' in label) or ('magnetization' in label):
                value = value*100
                obj_file.write(label+':  '+str(value)+'%\n')
            if 'log' in label:
                obj_file.write(label+':  '+str(value)+'\n')
            if 'band_gap' in label:
                obj_file.write(label+':  '+str(value)+' eV\n')
            if 'eos' in label:
                obj_file.write(label[:-3]+'delta_factor:  '+str(value)+'\n')
            if 'phonon_frequency' in label:
                obj_file.write(label+':  '+str(value)+' THz\n')
            if 'atomic_positions' in label:
                obj_file.write(label+':  '+str(value)+' angstroms\n')
    upf_files = [fname for fname in glob.iglob('*UPF')]
    atompaw_files = [fname for fname in glob.iglob('*atompaw*')]
    if not os.path.isdir('Best_Solution'):
        os.mkdir('../Best_Solution')
    results_df = pd.read_table('results.out',sep='\s+',header=None)
    obj_fn_list = [float(value) for value in list(results_df[0])]
    if os.path.exists('../Best_Solution/results.out'):
        last_results_df = pd.read_table('../Best_Solution/results.out',
            sep='\s+',header=None)
        last_obj_fn_list = [float(value) for value in list(last_results_df[0])]
        index = 1
        for obj_fn in obj_fn_list:
            last_max = float(np.loadtxt('../Best_Solution/Max_Error_'+str(index)))
            if obj_fn > last_max:
                os.remove('../Best_Solution/Max_Error_'+str(index))
                with open('../Best_Solution/Max_Error_'+str(index),'w+') as error_file:
                    error_file.write(str(obj_fn))
            index += 1
    else:
        index = 1
        for obj_fn in obj_fn_list:
            with open('../Best_Solution/Max_Error_'+str(index),'w+') as obj_file:
                obj_file.write(str(obj_fn))
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
    if os.path.exists('../Best_Solution/results.out'):
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
        for fname in os.listdir('../Best_Solution/'):
            if 'Max_Error' not in fname:
                os.remove('../Best_Solution/'+fname)
        copyfile('OBJ_FN','../Best_Solution/results.out')
        for (file_atom,file_upf) in zip(atompaw_files,upf_files):
            copyfile(file_atom,'../Best_Solution/'+file_atom)
            copyfile(file_atom,'../Best_Solution/'+file_upf)
        with open('../Best_Solution/rms_error','w+') as rms_file:
            rms_file.write(str(rms_error))
