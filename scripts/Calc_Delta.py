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


## Usage: python Calc_Delta.py elem_name lat_type
## Make sure all template and AE files are in current directory

def main():
    """
    Calculate EOS with QE and compare with AE data
    """
    cmpd = sys.argv[-2]
    cmpd_lat_type = sys.argv[-1]
    template_dir = '.'
    write_QE_input(cmpd,cmpd_lat_type,'relax',template_dir)
    run_QE(cmpd,cmpd_lat_type,'relax')
    run_scale_lat(cmpd,cmpd_lat_type,template_dir)
    V0, QE_bulk, B_prime = get_bulk(cmpd,cmpd_lat_type)
    QE_EOS_data, AE_EOS_data = read_eos(cmpd,cmpd_lat_type,template_dir)
    delta_factor = calcDelta(QE_EOS_data,AE_EOS_data,[cmpd],False)
    f = open('DELTA-FACTOR','w+')
    f.write(str(delta_factor)+'\n')
    f.close()


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


if __name__=='__main__':
    main()

