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
from gaopaw.calcDelta import *


def parse_num_objs(working_dir):
    """
    Read in number of objective functions defined in the dakota.in file.
    """
    with open(os.path.join(working_dir, os.pardir, 'dakota.in')) as dakota_input:
        for line in dakota_input:
            if 'num_objective_functions' in line:
                num_objs = line.split()[2]
    return int(num_objs)

def parse_elems(formula):
    """
    Parse compound formula to obtain constituent elements.
    """
    assert formula[0].isupper(), 'First letter of cmpd formula should be capitalized'
    letters_only = ''.join([letter for letter in formula if not letter.isdigit()])
    index = -1
    elems = []
    for letter in letters_only:
        if letter.isupper():
            elems.append(letter)
            index += 1
        else:
            elems[index] += letter
    return list(set(elems))

def get_num_objs(cmpd_list, element_list):
    """
    Parse input.json and return total number of objective functions.
    """
    num_elems = len(element_list)
    if num_elems == 1 and len(vars(cmpd_list[0])) == 1:
        if element_list[0] in ['N', 'P']:
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
        if elem in ['N', 'P']:
            num_obj_fns += 2
        else:
            num_obj_fns += 3
    num_obj_fns += num_properties
    return int(num_obj_fns)

def form_cmpd_dict(cmpd_list):
    """
    Constructs empty dictionary of correct length for compound testing.
    """
    cmpd_diff_dict = {}
    for cmpd in cmpd_list:
        formula = cmpd.formula
        lat_type = cmpd.lattice_type
        property_list = [property for property in vars(cmpd) 
            if property not in ['formula', 'lattice_type']]
        try:
            cmpd_diff_dict[formula]
        except KeyError:
            cmpd_diff_dict[formula] = {}
        cmpd_diff_dict[formula][lat_type] = {}
        for property in property_list:
            cmpd_diff_dict[formula][lat_type][property] = {}
    return cmpd_diff_dict

def get_element_info(template_dir):
    """
    Read in AE data for elemental lattice constants from specified 
    template directory. This includes FCC/BCC lattice constants for 
    all non f-block elements excluding N and P, for which atomic positions 
    and lattice constants are tested in simple cubic and orthorhombic 
    structures respectively. For f-block, RS oxides are considered.
    """
    elemental_data = {}
    df_fcc = pd.read_table(os.path.join(template_dir, 'WIEN2k_FCC'), 
        sep='\s+', header=None)
    for (elem, lat_const) in zip(df_fcc[0], df_fcc[1]):
        elemental_data[elem] = {}
        elemental_data[elem]['FCC'] = lat_const
    df_bcc = pd.read_table(os.path.join(template_dir, 'WIEN2k_BCC'), 
        sep='\s+', header=None)
    for (elem, lat_const) in zip(df_bcc[0], df_bcc[1]):
        elemental_data[elem]['BCC'] = lat_const
    df_ren = pd.read_table(os.path.join(template_dir, 'F_BLOCK_RS_NITRIDES'),
        sep='\s+', header=None)
    for (elem, lat_const, mag) in zip(df_ren[0], df_ren[1], df_ren[2]):
        elemental_data['%sN' % elem] = {}
        elemental_data['%sN' % elem]['RS'] = lat_const
        elemental_data['%sN' % elem]['magnetization'] = mag
    elemental_data['N'] = {}
    elemental_data['N']['SC'] = 6.1902
    elemental_data['P'] = {}
    elemental_data['P']['ortho'] = [3.304659, 4.573268, 11.316935]
    return elemental_data

def test_element_list(elem_list, template_dir):
    """
    Perform and check UPF generation with atompaw, compare pseudized
    log derivatives with corresponding AE log derivatives for each orbital,
    and compare QE with AE lattice constants for elemental states.
    """
    elem_diff_dict = {}
    elemental_data = get_element_info(template_dir)
    f_block = ['La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',
        'Er','Tm','Yb','Lu','Ac','Th','Pa','U','Np','Pu','Am']
    for elem in elem_list:
        assert elem in elemental_data.keys(), \
            'No AE data available for your element: %s' % elem
        elem_diff_dict[elem] = {}
        elem_diff_dict[elem]['elemental'] = {}
        elem_diff_dict[elem]['elemental']['log'] = {}
        if elem in ['N', 'P'] + f_block:
            if elem == 'N':
                elem_diff_dict[elem]['SC'] = {}
                elem_diff_dict[elem]['SC']['atomic_positions'] = {}
            if elem == 'P':
                elem_diff_dict[elem]['ortho'] = {}
                elem_diff_dict[elem]['ortho']['lattice_constant'] = {}
            if elem in f_block:
                elem_diff_dict['%sN' % elem] = {}
                elem_diff_dict['%sN' % elem]['RS'] = {}
                elem_diff_dict['%sN' % elem]['RS']['lattice_constant'] = {}
                elem_diff_dict['%sN' % elem]['RS']['magnetization'] = {}
        else:
            for lat_type in ['FCC', 'BCC']:
                elem_diff_dict[elem][lat_type] = {}
                elem_diff_dict[elem][lat_type]['lattice_constant'] = {}
    for elem in elem_list:
        os.mkdir(elem)
        copyfile('params.in', os.path.join(elem, 'params.in'))
        with fileutils.chdir(elem):
            write_atompaw_input(elem, template_dir)
            copyfile('%s.atompaw.in' % elem, os.path.join(os.pardir, '%s.atompaw.in' % elem))
            run_atompaw(elem)
            if not check_upf():
                return elem_diff_dict, True
            copyfile('%s.GGA-PBE-paw.UPF' % elem, os.path.join(os.pardir, '%s.GGA-PBE-paw.UPF' % elem))
            elem_diff_dict[elem]['elemental']['log'] = compare_log()
            if elem in ['N', 'P'] + f_block:
                if elem == 'N':
                    run_qe(elem, 'SC', 'relax', template_dir)
                    if not check_convergence(elem, 'SC', 'relax'):
                        return elem_diff_dict, True
                    elem_diff_dict[elem]['SC']['atomic_positions'] = \
                         compare_atoms(elem, 'SC', template_dir)
                if elem == 'P':
                    run_qe(elem, 'ortho', 'relax', template_dir)
                    if not check_convergence(elem, 'ortho', 'relax'):
                        return elem_diff_dict, True
                    ae_lat = elemental_data[elem]['ortho']
                    elem_diff_dict[elem]['ortho']['lattice_constant'] = \
                        compare_lat(ae_lat, elem, 'ortho')
                if elem in f_block:
                    copyfile(os.path.join(os.pardir,'N.GGA-PBE-paw.UPF'),'./N.GGA-PBE-paw.UPF')
                    run_qe('%sN' % elem, 'RS', 'relax', template_dir)
                    if not check_convergence('%sN' % elem, 'RS', 'relax'):
                        return elem_diff_dict, True
                    ae_lat = elemental_data['%sN' % elem]['RS']
                    elem_diff_dict['%sN' % elem]['RS']['lattice_constant'] = \
                        compare_lat(ae_lat, '%sN' % elem, 'RS')
                    run_qe('%sN' % elem, 'RS', 'scf', template_dir)
                    if not check_convergence('%sN' % elem, 'RS', 'scf'):
                        return elem_diff_dict, True
                    qe_mag = get_mag('%sN' % elem, 'RS')
                    ae_mag = elemental_data['%sN' % elem]['magnetization']
                    elem_diff_dict['%sN' % elem]['RS']['magnetization'] = \
                        abs(ae_mag-qe_mag)
            else:
                for lat_type in ['FCC', 'BCC']:
                    run_qe(elem, lat_type, 'relax', template_dir)
                    if not check_convergence(elem, lat_type, 'relax'):
                        return elem_diff_dict, True
                    ae_lat = elemental_data[elem][lat_type]
                    elem_diff_dict[elem][lat_type]['lattice_constant'] = \
                        compare_lat(ae_lat, elem, lat_type)
    return elem_diff_dict, False

def test_cmpd_list(cmpd_list, cmpd_diff_dict, cmpd_template_dir, elem_template_dir):
    """
    For each compound specified in input.json, perform the given property tests
    and compare with AE values. Differences correspond to objective functions.
    """
    elemental_data = get_element_info(elem_template_dir)
    for cmpd in cmpd_list:
        formula = cmpd.formula
        lat_type = cmpd.lattice_type
        property_list = [property for property in vars(cmpd)
            if property not in ['formula', 'lattice_type']]
        for property in property_list:
            elem_err_mssg = """%s not needed for %s, this is tested
                automatically, therefore property should be
                removed from input.json""" % (property, formula)
            if formula in elemental_data.keys():
                if (formula not in ['N','P']) and (lat_type in ['FCC','BCC']):
                    assert property != 'lattice_constant', elem_err_mssg
                if (cmpd == 'N') and (lat_type == 'SC'):
                    assert property != 'atomic_positions', elem_err_mssg
                if (cmpd == 'P') and (lat_type == 'ortho'):
                    assert property != 'lattice_constant', elem_err_mssg
            f_block = ['La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',
                'Er','Tm','Yb','Lu','Ac','Th','Pa','U','Np','Pu','Am']
            if formula in ['%sN' % f_elem for f_elem in f_block]:
                if lat_type == 'RS':
                    assert property != 'lattice_constant', elem_err_mssg
                    assert property != 'magnetization', elem_err_mssg
            ae_value = getattr(cmpd, property)
            cmpd_diff_dict[formula][lat_type][property], error_check = \
                test_property(formula, lat_type, property, ae_value, cmpd_template_dir)
            if error_check:
                return cmpd_diff_dict, True
    return cmpd_diff_dict, False

def merge_dicts(dct, merge_dct):
    """
    *Recursively* merge two dictionaries.
    """
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict) and \
         isinstance(merge_dct[k], collections.Mapping)):
            merge_dicts(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]
    return dct

def test_property(cmpd, lat_type, property, ae_data, template_dir):
    """
    For a given compound and property, perform the required QE
    calculations to obtain data which will be compared with an
    specified AE value (difference ~ objective function).
    """
    run_qe(cmpd, lat_type, 'relax', template_dir)
    if not check_convergence(cmpd, lat_type, 'relax'):
        return None, True
    if property == 'lattice_constant':
        return compare_lat(ae_data, cmpd, lat_type), False
    if property in ['eos', 'bulk_modulus']:
        run_scale_lat(cmpd, lat_type, template_dir)
        V0, QE_bulk, B_prime = get_bulk(cmpd, lat_type)
        if property == 'eos':
            assert len(ae_data) == 3, \
                'Three parameters required for EOS'
            qe_data = np.loadtxt('QE_EOS.txt')
            qe_eos = {'element': [cmpd], 'V0': [qe_data[0]], 
                'B0': [qe_data[1]], 'BP': [qe_data[2]]}
            ae_eos = {'element': [cmpd], 'V0': [ae_data[0]], 
                'B0': [ae_data[1]], 'BP': [ae_data[2]]}
            delta_factor = calcDelta(qe_eos, ae_eos, [cmpd])[0]
            with open('%s.%s.%s.in' % (cmpd, lat_type, 'relax')) as qe_input:
                for line in qe_input:
                    if 'nat=' in line:
                        natoms = float(line.split('=')[1][:-1])
            delta_factor = delta_factor / natoms
            return delta_factor, False
        if property == 'bulk_modulus':
            bulk_diff = abs(QE_bulk - ae_data)/ae_data
            return bulk_diff, False
    if property == 'phonon_frequency':
        run_qe(cmpd, lat_type, 'scf', template_dir)
        if not check_convergence(cmpd, lat_type, 'scf'):
            return None, True
        run_phonon(cmpd, lat_type, template_dir)
        phonon_diff = compare_phonon(cmpd, lat_type, ae_data, template_dir)
        return phonon_diff, False
    if property == 'atomic_positions':
        return compare_atoms(cmpd, lat_type, template_dir), False
    if property == 'band_gap':
        run_qe(cmpd, lat_type, 'scf', template_dir)
        if not check_convergence(cmpd, lat_type, 'scf'):
            return None, True
        qe_gap = get_gap(cmpd, lat_type)
        return abs(ae_data-qe_gap), False
    if property == 'magnetization':
        run_qe(cmpd, lat_type, 'scf', template_dir)
        if not check_convergence(cmpd, lat_type, 'scf'):
            return None, True
        qe_mag = get_mag(cmpd, lat_type)
        return abs(ae_data-qe_mag), False
    if property == 'magnetic_moment':
        run_qe(cmpd, lat_type, 'scf', template_dir)
        if not check_convergence(cmpd, lat_type, 'scf'):
            return None, True
        return compare_mag_mom(cmpd, lat_type, ae_data, template_dir), False
    raise ValueError('Your property, %s , is not defined' % property)

def check_upf():
    """
    Check if a .UPF file was succesfully created by AtomPAW.
    """
    if len(glob.glob('*UPF')) != 0:
        return True
    return False

def bad_run(num_obj_fns):
    """
    If something went wrong with the run, e.g., no .UPF file created or
    if running QE raised an error, set all objective functions to 100.
    """
    params, results = di.read_parameters_file('params.in', 'results.out')
    for num in range(1, num_obj_fns+1):
        label = 'obj_fn_%s' % num
        results[label].function = 100.0
    results.write()

def compare_lat(ae_lat, cmpd, lat_type):
    """
    Calculate the average difference between QE and QE lattice parameters.
    AE values must be given in terms of a conventional unit cell, however,
    a primitive or conventional cell may be used in the QE calculation.
    get_lattice_constant() will convert into conventional units.
    """
    qe_lat = get_lattice_constant(cmpd, lat_type)
    if isinstance(ae_lat, list) == False:
        return abs(qe_lat-ae_lat)/ae_lat
    assert len(ae_lat) == len(qe_lat), \
        'Wrong number of lattice parameters given for specified lattice type'
    lat_diff = 0
    for (qe_val, ae_val) in zip(qe_lat, ae_lat):
        lat_diff += abs(qe_val-ae_val)/ae_val
    return lat_diff/len(ae_lat)

def update_dakota(diff_dict):
    """
    Set the parameters and results files to be used by Dakota.
    Objective functions are updated according to differences between
    pseudized and AE properties calculated with atompaw and QE.
    """
    params, results = di.read_parameters_file('params.in', 'results.out')
    label_index = 1
    for elem in diff_dict.keys():
        for lat_type in diff_dict[elem].keys():
            for property in diff_dict[elem][lat_type].keys():
                label = 'obj_fn_'+str(label_index)
                results[label].function = diff_dict[elem][lat_type][property]
                label_index += 1
    results.write()

def write_atompaw_input(elem, template_dir):
    """
    Write atompaw input file for elem based on existing
    elem.atompaw.templtae in specified template_dir.
    """
    template_file = os.path.join(template_dir, '%s.atompaw.template' % elem)
    new_input_file = '%s.atompaw.in' % elem
    subprocess.check_call(['run', 'dprepro.py', 'params.in', template_file, 
        new_input_file], env=os.environ.copy())

def run_atompaw(elem):
    """
    Run AtomPAW using elem.atompaw.in.
    """
    with open('%s.atompaw.in' % elem, 'r') as input_fin, open('log_atompaw', 'w') as log_fout: 
        subprocess.call(['atompaw'], stdin=input_fin, stdout=log_fout, 
            env=os.environ.copy())

def run_qe(cmpd, lat_type, calc_type, template_dir):
    """
    Write and run QE using cmpd.lat_type.calc_type from template_dir.
    """
    if os.path.exists('%s.%s.%s.out' % (cmpd, lat_type, calc_type)):
        return
    template_file = os.path.join(template_dir, '%s.%s.%s.template' % (cmpd, lat_type, calc_type))
    qe_input = '%s.%s.%s.in' % (cmpd, lat_type, calc_type)
    shutil.copy(template_file, qe_input)
    if calc_type == 'scf':
        update_structure(cmpd, lat_type, 'scf')
    subprocess.call(['run', 'periodic_dft_gui_dir/runner.py', 'pw.x', qe_input, 
        '-MPICORES', '4'], env=os.environ.copy())

def get_lattice_constant(cmpd, lat_type, tol=3):
    """
    Parse relaxed lattice parameters from QE relaxation run.
    Allowed tolerance is given by tol ~ number of decimal places.
    Values are converted to those of a conventional unit cell.
    * Standard oriention of primitive cells (see QE docs) assumed. *
    * Must use unique axis c for base-centered lattices. *
    """
    qe_reader_path = os.path.join(fileutils.get_mmshare_scripts_dir(), 
        'periodic_dft_gui_dir', 'qe2mae.py')
    qe_reader_mod = imputils.import_module_from_file(qe_reader_path)
    qe_reader = qe_reader_mod.QEOutputReader('%s.%s.relax.out' % (cmpd, lat_type))
    struct = qe_reader.structs[qe_reader.final_struct_id]
    cparams = xtal.get_chorus_properties(struct)
    params = np.array(xtal.get_params_from_chorus(cparams)).round(tol)
    unique_lat = np.array(sorted(list(set(params[:3]))))
    unique_angles = np.array(sorted(list(set(params[3:]))))
    err_mssg = 'Input lattice is incorrect, does not match %s' % lat_type
    if lat_type in ['FCC', 'ZB', 'RS', 'diamond', 'HH']: ## face-centered (F)
        assert len(unique_lat) == 1, err_mssg
        if np.array_equal(unique_angles, [60.0]):
            return math.sqrt(2)*params[0]
        if np.array_equal(unique_angles, [90.0]):
            return params[0]
        raise ValueError(err_mssg)
    if lat_type == 'BCC': ## body-centered (I)
        assert len(unique_lat) == 1, err_mssg
        if len(unique_angles) == 2:
            return (2./3.)*math.sqrt(3)*params[0]
        if np.array_equal(unique_angles, [90.0]):
            return params[0]
        raise ValueError(err_mssg)
    if lat_type in ['per', 'SC', 'CsCl']: ## conv (P)
        assert np.array_equal(unique_angles, [90.0]) and len(unique_lat) == 1, err_mssg
        return params[0]
    if lat_type in ['hex', 'WZ']: ## conv (P)
        assert len(unique_lat) == 2 \
            and np.array_equal(unique_angles, [90.0, 120.0]), err_mssg
        return unique_lat[0], unique_lat[1]
    if lat_type == 'rhomb': ## trig (R)
        assert len(unique_lat) == 1 and len(unique_angles) == 1 \
            and not np.array_equal(unique_angles, [90.0]), err_mssg
        return unique_lat[0], unique_angles[0]
    if lat_type == 'tetrag':
        if len(unique_lat) == 1 and len(unique_angles) == 2: ## body-centered (I)
            cell_vecs = get_cell(cmpd, lat_type, 'relax')
            abs_vec = [abs(value) for value in cell_vecs[0]]
            prim_lengths = sorted(set(abs_vec))
            conv_lengths = np.array([2*value for value in prim_lengths]).round(tol)
            return conv_lengths[0], conv_lengths[1]
        if len(unique_lat) == 2 and np.array_equal(unique_angles, [90.0]): ## conv (P)
            return unique_lat[0], unique_lat[1]
        raise ValueError(err_mssg)
    if lat_type == 'ortho':
        if len(unique_lat) == 3 and np.array_equal(unique_angles, [90.0]): ## conv (P)
            return unique_lat[0], unique_lat[1], unique_lat[2]
        if len(unique_lat) == 2 and len(unique_angles) == 2: ## base-centered (C)
            assert list(params[3:]).count(90.0) == 2, err_mssg
            cell_vecs = get_cell(cmpd, lat_type, 'relax')
            a_lat = cell_vecs[0][0]*2.
            b_lat = cell_vecs[0][1]*2.
            c_lat = cell_vecs[2][2]
            conv_lengths = np.array(sorted([a_lat, b_lat, c_lat])).round(tol)
            return conv_lengths[0], conv_lengths[1], conv_lengths[2]
        if len(unique_lat) == 3 and len(unique_angles) == 3: ## face-centered (F)
            cell_vecs = np.array(get_cell(cmpd, lat_type, 'relax'))
            components = [abs(value) for value in cell_vecs.flatten()]
            prim_lengths = sorted(set(components))[1:]
            conv_lengths = np.array([2*value for value in prim_lengths]).round(tol)
            return conv_lengths[0], conv_lengths[1], conv_lengths[2]
        if len(unique_lat) == 1 and len(unique_angles) == 3: ## body-centered (I)
            cell_vecs = np.array(get_cell(cmpd, lat_type, 'relax'))
            components = [abs(value) for value in cell_vecs.flatten()]
            prim_lengths = sorted(set(components))
            conv_lengths = np.array([2*value for value in prim_lengths]).round(tol)
            return conv_lengths[0], conv_lengths[1], conv_lengths[2]
        raise ValueError(err_mssg)
    if lat_type == 'monoclin':
        if list(params[3:]).count(90.0) == 2: ## conv (P)
            for value in params[3:]:
                if round(value, 2) != 90.0:
                    angle = value
            return unique_lat[0], unique_lat[1], unique_lat[2], angle
        if len(unique_lat) == 2 and len(unique_angles) == 2: ## base-centered (C)
            cell_vecs = np.array(get_cell(cmpd, lat_type, 'relax'))
            a_lat = cell_vecs[0][0]*2
            c_lat = cell_vecs[2][2]*2
            compon_1 = cell_vecs[1][0]*2
            compon_2 = cell_vecs[1][1]*2
            b_lat = math.sqrt(compon_1**2 + compon_2**2)
            angle = round(math.degrees(math.acos(compon_1/b_lat)), 3)
            b_lat = b_lat/2. ## Not sure why yet
            conv_lengths = np.array(sorted([a_lat, b_lat, c_lat])).round(tol)
            return conv_lengths[0], conv_lengths[1], conv_lengths[2], angle
        raise ValueError(err_mssg)
    if lat_type == 'triclin': ## conv (P)
        assert len(unique_lat) == 3 and list(params[3:]).count(90.0) < 2, err_mssg
        return params

def check_convergence(cmpd, lat_type, calc_type):
    """
    Check if the QE calculation ran succesfully
    (i.e., w/o error and convergence acheived).
    """
    qe_error_signs = \
        ['convergence NOT', 'S matrix not positive definite', 'stopping ...']
    with open('%s.%s.%s.out' % (cmpd, lat_type, calc_type)) as qe_output:
        for line in qe_output:
            if any(error in line for error in qe_error_signs):
                return False
    return True

def compare_log():
    """"
    Compare arctan of the logarithmic derivatives of the pseudized and 
    exact wavefunctions produced by atompaw (goal is to minimize difference).
    Comparing arctan works better than using log derivs explicitly.
    Note that columns in logderiv.l correspond to:
    energy, exact logderiv, pseudized logderv, exact arctan, pseudized arctan.
    We exclude the unbound state (highest l) from consideration.
    """
    log_derivs = glob.glob('logd*')
    sum_log = 0
    total_diff = 0
    for fname in log_derivs[:-1]: ## Exclude unbound state
        log_data = np.loadtxt(fname).transpose()
        log_exact = np.array(log_data[3])
        log_pseudo = np.array(log_data[4])
        sum_log += sum([abs(value) for value in log_exact])
        total_diff += sum(abs(log_pseudo - log_exact))
    return total_diff/sum_log

def compare_atoms(cmpd, lat_type, template_dir):
    """
    Compare atomic positions of QE-relaxed structure and those 
    of the AE-relaxed structure. Return sum of the distances.
    """
    df_ae = pd.read_table(os.path.join(template_dir, 
        'AE_Struct.%s.%s' % (cmpd, lat_type)), sep='\s+', header=None)
    df_ae = df_ae.drop(0, 1).transpose()
    qe_reader_path = os.path.join(fileutils.get_mmshare_scripts_dir(), 
        'periodic_dft_gui_dir', 'qe2mae.py')
    qe_reader_mod = imputils.import_module_from_file(qe_reader_path)
    qe_reader = qe_reader_mod.QEOutputReader('%s.%s.relax.out' % (cmpd, lat_type))
    struct = qe_reader.structs[qe_reader.final_struct_id]
    df_qe = struct.getXYZ()
    distance_list = []
    for index in range(len(df_ae.keys())):
        ae_position = np.array(df_ae[index])
        qe_position = np.array(df_qe[index])
        distance_list.append(np.linalg.norm(ae_position-qe_position))
    return sum(distance_list)

def get_mag(cmpd, lat_type):
    """
    Parse QE output (scf run) to obtain total magnetization.
    """
    mag = []
    with open('%s.%s.scf.out' % (cmpd, lat_type)) as qe_output:
        for line in qe_output:
            if 'absolute' in line.split():
                mag.append(line.split()[3])
    assert len(mag) != 0, """No magnetization found, 
        spin-polarization must be considered in SCF calculation"""
    return float(mag[-1])

def compare_mag_mom(cmpd, lat_type, ae_mag_mom, template_dir):
    """
    Parse QE output (scf run) to obtain individual
    magnetic moments of atoms in given structure.
    Compare these with AE magnetic moments.
    """
    qe_mag_mom = []
    with open('%s.%s.scf.out' % (cmpd, lat_type)) as qe_output:
        for line in qe_output:
            if 'magn:' in line.split():
                qe_mag_mom.append(line.split()[5])
    assert len(qe_mag_mom) != 0, """No magnetization found, 
        spin-polarization must be considered in SCF calculation"""
    qe_mag_mom = [float(value) for value in qe_mag_mom]
    rel_diff = []
    for (qe_val, ae_val) in zip(qe_mag_mom, ae_mag_mom):
        rel_diff.append(abs(qe_val-ae_val))
    net_diff = sum(rel_diff)/len(rel_diff)
    return net_diff

def get_gap(cmpd, lat_type):
    """
    Parse QE output (scf run) to obtain band gap.
    Note that sufficient unoccupied bands need to be included
    and occupations set to fixed in the scf run.
    """
    band_gap = None
    with open('%s.%s.scf.out' % (cmpd, lat_type)) as qe_output:
        for line in qe_output:
            if 'highest' and 'lowest' in line.split():
                band_gap = (float(line.split()[7]) - float(line.split()[6]))
    if not band_gap:
        err = """Energies of highest occupied and lowest unoccupied 
            states could not be found; ensure occupation is fixed"""
        raise NameError(err)
    return band_gap

def birch_murnaghan(vol, vol_equil, bulk, bulk_prime, energy_equil):
    """
    3rd order Birch-Murnaghan equation of state, in the energy-volume form.
    """
    vol = np.array(vol)
    return energy_equil + 9 * vol_equil * bulk / 16. * (
        ((vol_equil / vol) ** (2 / 3.) - 1) ** 3 * bulk_prime +
        ((vol_equil / vol) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (vol_equil / vol) ** (2 / 3.)))

def get_bulk(cmpd, lat_type, ev_fname='E_V.txt'):
    """
    Reads in energy-volume data from E_V.txt and calculates:
    equilibrium volume, bulk modulus, and dB/dP,
    i.e., fit to Birch Murnaghan equation of state.
    """
    ev_data = np.loadtxt(ev_fname).transpose()
    energy = list(ev_data[0])
    energy = np.array([13.6056980659*value for value in energy]) ## Ry to eV
    volume = list(ev_data[1])
    volume = np.array([0.14818453429566825*value for value in volume]) ## bohr^3 to A^3
    initial_parameters = [volume.mean(), 2.5, 4, energy.mean()]
    fit_eqn = eval('birch_murnaghan')
    popt, pcov = cf(fit_eqn, volume, energy, initial_parameters)
    with open('%s.%s.relax.in' % (cmpd, lat_type)) as qe_output:
        for line in qe_output:
            if 'nat=' in line:
                num_atoms = float(line.split('=')[1][:-1])
    volume = popt[0]/num_atoms ## A^3/atom
    bulk = popt[1]*160.2 ## GPa
    bulk_prime = popt[2] ## Dimensionless
    with open('QE_EOS.txt', 'w+') as eos_file:
        eos_file.write('%s %s %s' % (volume, bulk, bulk_prime))
    return volume, bulk, bulk_prime

def run_qe_as_is(relax_in):
    """
    Run QE w/o needing to fetch input, to be used in EOS calculations.
    """
    subprocess.call(['run', 'periodic_dft_gui_dir/runner.py', 'pw.x', relax_in,
        '-MPICORES', '4'], env=os.environ.copy())

def run_scale_lat(cmpd, lat_type, template_dir, ev_fname='E_V.txt'):
    """
    Read in relaxed cell parameters from cmpd.lat_type.relax.out, 
    scale this lattice constant from 94% to 106% (7 values created), 
    write input and run QE at each value, write corresponding volumes
    and energies into E_V.txt (units of Bohr^3 and Ry^3)
    Note that crystal coords are preferred in QE input, otherwise
    atomic positions will not be scaled with the cell.
    """
    scale_num = [0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06]
    relax_in = '%s.%s.relax.in' % (cmpd, lat_type)
    relax_out = '%s.%s.relax.out' % (cmpd, lat_type)
    energies = []
    volumes = []
    folder_index = 1
    for value in scale_num:
        folder = '%s_%s' % (cmpd, folder_index)
        new_cell_params = scale_cell(cmpd, lat_type, value)
        new_cell_matrix = np.matrix(new_cell_params)
        volumes.append(np.linalg.det(new_cell_matrix))
        os.mkdir(folder)
        for fname in glob.iglob('*UPF'):
            copyfile(fname, os.path.join(folder, fname))
        copyfile(relax_in, os.path.join(folder, relax_in))
        copyfile(relax_out, os.path.join(folder, relax_out))
        with fileutils.chdir(folder):
            update_structure(cmpd, lat_type, 'relax')
            write_cell(cmpd, lat_type, new_cell_params)
            run_qe_as_is(relax_in)
            all_energies = []
            with open(relax_out) as qe_output:
                for line in qe_output:
                    if '!    total energy              =' in line:
                        all_energies.append(line.split()[4])
            energies.append(all_energies[-1])
        folder_index += 1
    with open(ev_fname, 'w+') as ev_file:
        for (e, v) in zip(energies, volumes):
            ev_file.write('%s %s \n' % (e, v))

def compare_phonon(cmpd, lat_type, ae_freq, template_dir):
    """
    Parse acoustic and optical phonon frequencies (usually at gamma)
    from QE run and compare with given AE frequencies.
    """
    with open('%s.%s.scf.in' % (cmpd, lat_type)) as qe_input:
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
    assert len(ae_freq) == len(qe_freq), \
        'Wrong number of phonon frequencies given'
    rel_diff = []
    for (qe_val, ae_val) in zip(qe_freq, ae_freq):
        rel_diff.append(abs(ae_val-qe_val))
    net_diff = sum(rel_diff)/len(rel_diff)
    return net_diff

def get_cell(cmpd, lat_type, calc_type):
    """
    Parse cell from QE output and write to 3x3 array
    consisting of lattice vectors in angstroms or bohrs.
    """
    with open('%s.%s.relax.out' % (cmpd, lat_type)) as f:
        lines = f.readlines()
    index = 0
    for line in lines:
        if 'CELL_PARAMETERS' in line:
            cell_index = [index+1, index+2, index+3]
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
    return vectors

def update_structure(cmpd, lat_type, calc_type):
    """
    Parse equilibrium structure from completed relaxation
    and update the corresponding calculation input file.
    """
    with open('%s.%s.relax.out' % (cmpd, lat_type)) as f:
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
        if ('End' in line) or (line == '\n'):
            break
        else:
            coords.append(line)
    coords_header = 'ATOMIC_POSITIONS %s \n' % coord_type
    index = 0
    for line in lines:
        if 'CELL_PARAMETERS' in line:
            cell_index = [index+1, index+2, index+3]
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
    v1 = '%s %s %s \n' % tuple(vectors[0])
    v2 = '%s %s %s \n' % tuple(vectors[1])
    v3 = '%s %s %s \n' % tuple(vectors[2])
    with open('%s.%s.%s.in' % (cmpd, lat_type, calc_type)) as qe_input:
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
    with open('%s.%s.%s.in' % (cmpd, lat_type, calc_type), 'w+') as qe_file:
        for line in orig_struct:
            qe_file.write(line)
        qe_file.write(coords_header)
        for line in coords:
            qe_file.write(line)
        qe_file.write('%s%s%s%s' % (cell_header, v1, v2, v3))

def scale_cell(cmpd, lat_type, scale_factor):
    """
    Scale cell volume according to scale_factor.
    """
    with open('%s.%s.relax.out' % (cmpd, lat_type)) as qe_output:
        lines = qe_output.readlines()
    index = 0
    for line in lines:
        if 'CELL_PARAMETERS' in line:
            cell_index = [index+1, index+2, index+3]
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
    scaled_cell = np.array(vectors)*(scale_factor**(1./3.))
    return scaled_cell

def write_cell(cmpd, lat_type, cell):
    """
    Write given cell to QE relaxation input.
    Also keep cell volume fixed during relaxation;
    to be used during equation of state calculations.
    """
    vectors = np.array(cell)
    v1 = '%s %s %s \n' % tuple(vectors[0])
    v2 = '%s %s %s \n' % tuple(vectors[1])
    v3 = '%s %s %s \n' % tuple(vectors[2])
    cell_header = 'CELL_PARAMETERS bohr\n'
    with open('%s.%s.relax.in' % (cmpd, lat_type)) as qe_input:
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
                    orig_struct.append("  calculation='relax'\n")
            else:
                orig_struct.append(line)
            line_index += 1
        else:
            line_index += 4
    with open('%s.%s.relax.in' % (cmpd, lat_type), 'w+') as qe_file:
        for line in orig_struct:
            qe_file.write(line)
        qe_file.write('%s%s%s%s' % (cell_header, v1, v2, v3))

def run_phonon(cmpd, lat_type, template_dir):
    """
    Run QE phonon calculations using ph.x and dynmat.x.
    """
    copyfile(os.path.join(template_dir, 'phonon.in'), 'phonon.in')
    if os.path.exists('phonon.save'):
        os.remove('phonon.out')
        shutil.rmtree('phonon.save')
        os.remove('phonon.save.qegz')
    scf_savefile = '%s.%s.scf.save.qegz' % (cmpd, lat_type)
    subprocess.call(['run', 'periodic_dft_gui_dir/runner.py', 'ph.x', 'phonon.in', 
        '-input_save', scf_savefile, '-MPICORES', '4'], env=os.environ.copy())
    copyfile(os.path.join(template_dir, 'dynmat.in'), 'dynmat.in')
    if os.path.exists('dynmat.save'):
        os.remove('dynmat.out')
        shutil.rmtree('dynmat.save')
        os.remove('dynmat.save.qegz')
    subprocess.call(['run', 'periodic_dft_gui_dir/runner.py', 'dynmat.x', 'dynmat.in', 
        '-input_save', 'phonon.save.qegz', '-MPICORES', '4'], env=os.environ.copy())

def dict_to_list(diff_dict):
    """
    Convert a detailed dictionary of cmpds, lattice types, and properties
    into lists of objective function values and labels.
    """
    obj_fn_list = []
    obj_fn_labels = []
    for formula in diff_dict.keys():
        for lat_type in diff_dict[formula].keys():
            for property in diff_dict[formula][lat_type].keys():
                obj_fn_labels.append('%s_%s_%s' % (formula, lat_type, property))
                obj_fn_list.append(diff_dict[formula][lat_type][property])
    return obj_fn_list, obj_fn_labels

def update_obj_file(diff_dict):
    """
    Write objective functions to two files:
    Detailed_Results is placed in the working directory
    and contains labels and values of obj fns.
    Obj_Fn_Data is placed in parent directory and
    contains only numeric data to be parsed later on.
    """
    obj_fn_list, obj_fn_labels = dict_to_list(diff_dict)
    with open('Detailed_Results', 'w+') as obj_file:
        for (value, label) in zip(obj_fn_list, obj_fn_labels):
            value = round(float(value), 6)
            if ('lattice_constant' in label) or ('bulk_modulus' in label):
                value = value*100
                obj_file.write('%s: %s%%\n' % (label, value))
            if ('magnetization' in label) or ('magnetic_moment' in label):
                obj_file.write('%s: %s bohr mag\n' % (label, value))
            if 'log' in label:
                obj_file.write('%s: %s\n' % (label, value))
            if 'band_gap' in label:
                obj_file.write('%s: %s eV\n' % (label, value))
            if 'eos' in label:
                obj_file.write('%sdelta_factor: %s meV/atom\n' % (label[:-3], value))
            if 'phonon_frequency' in label:
                obj_file.write('%s: %s THz\n' % (label, value))
            if 'atomic_positions' in label:
                obj_file.write('%s: %s angstroms\n' % (label, value))
    if not os.path.exists(os.path.join(os.pardir, 'Obj_Fn_Data')):
        obj_file = open(os.path.join(os.pardir, 'Obj_Fn_Data'), 'w+')
        obj_file.close()
    with open(os.path.join(os.pardir, 'Obj_Fn_Data'), 'a') as obj_file:
        obj_file.write('\n')
        for obj_fn in obj_fn_list:
            obj_file.write('%s ' % obj_fn)

def calc_obj_fn(obj_fn_list):
    """
    Calculate the composite objective function which is to be
    given to dakota for minimization.
    The weighted sum approach, for which results are compared
    to the utopian and nadir points, is employed to normalize
    each individual objective function.
    The final/composite objective function is defined
    as the mean absolute error among all normalized
    objective functions.
    """
    obj_fn_matrix = np.loadtxt(os.path.join(os.pardir, 'Obj_Fn_Data')).transpose()
    min_values = []
    max_values = []
    for all_objs in obj_fn_matrix:
        if isinstance(all_objs, np.ndarray):
            min_values.append(min(all_objs))
            max_values.append(max(all_objs))
        else:
            return 0.0
    norm_objs = []
    for (obj, utopia, nadir) in zip(obj_fn_list, min_values, max_values):
        if nadir - utopia != 0:
            norm_objs.append( abs( (obj - utopia)/(nadir - utopia) ) )
        else:
            norm_objs.append(0.0)
    return sum(norm_objs)/len(norm_objs)

def update_best_result(obj_fn_list):
    """
    If current solution is better than all previous solutions, 
    update the Best_Result folder accordingly.
    """
    current_obj_fn = calc_obj_fn(obj_fn_list)
    upf_files = glob.glob('*UPF')
    atompaw_files = glob.glob('*atompaw*')
    if not os.path.isdir(os.path.join(os.pardir, 'Best_Solution')):
        os.mkdir(os.path.join(os.pardir, 'Best_Solution'))
        for fname in upf_files + atompaw_files:
            copyfile(fname, os.path.join(os.pardir, 'Best_Solution', fname))
        with open(os.path.join(os.pardir, 'Best_Solution', 'Obj_Fn'), 'w+') as obj_file:
            obj_file.write(str(current_obj_fn))
        with open(os.path.join(os.pardir, 'Best_Solution', 'data'), 'w+') as data_file:
            for obj in obj_fn_list:
                data_file.write('%s ' % obj)
        copyfile('Detailed_Results', os.path.join(os.pardir, 'Best_Solution', 'Detailed_Results'))
        return
    last_data = np.loadtxt(os.path.join(os.pardir, 'Best_Solution', 'data'))
    last_obj_fn = calc_obj_fn(last_data)
    if current_obj_fn < last_obj_fn:
        for fname in upf_files + atompaw_files:
            copyfile(fname, os.path.join(os.pardir, 'Best_Solution', fname))
        with open(os.path.join(os.pardir, 'Best_Solution', 'Obj_Fn'), 'w+') as obj_file:
            obj_file.write(str(current_obj_fn))
        with open(os.path.join(os.pardir, 'Best_Solution', 'data'), 'w+') as data_file:
            for obj in obj_fn_list:
                data_file.write('%s ' % obj)
        copyfile('Detailed_Results', os.path.join(os.pardir, 'Best_Solution', 'Detailed_Results'))
    else:
        with open(os.path.join(os.pardir, 'Best_Solution', 'Obj_Fn'), 'w+') as obj_file:
            obj_file.write(str(last_obj_fn))


