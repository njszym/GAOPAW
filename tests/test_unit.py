import mock
import glob
import shutil
from collections import Counter
import pandas as pd
import numpy as np
from schrodinger.utils import fileutils
import os
from unittest.mock import Mock
from types import SimpleNamespace
import gaopaw
import sys

top_dir = os.getcwd()

with fileutils.chdir('ex_templates/Empty'):
    gp_run = gaopaw.Runner()

def test_dir():
    assert gp_run.working_dir == os.path.join(top_dir,'ex_templates/Empty')

def test_readInput():
    assert hasattr(gp_run, 'input_settings')
    assert np.array_equal(list(gp_run.input_settings.__dict__.keys()), 
        ['directories', 'compounds'])
    assert np.array_equal(list(gp_run.input_settings.directories.__dict__.keys()),
        ['elem_template_dir', 'cmpd_template_dir', 'include_paw'])
    assert len(gp_run.input_settings.compounds) == 3

def test_getPaws():
    paw_list = ['%s.GGA-PBE-paw.UPF' % elem for elem in ['Si', 'O']]
    with fileutils.chdir('ex_templates/Empty'):
        for fname in paw_list:
            if os.path.exists(fname):
                os.remove(fname)
        gp_run.getPaws()
    for fname in paw_list:
        assert os.path.exists(os.path.join('ex_templates/Empty',fname))

def test_numDakotaObjs():
    with fileutils.chdir('ex_templates/Empty'):
        gp_run = gaopaw.Runner()
        assert gp_run.num_obj_fns == 15
    with fileutils.chdir('Be_workdir/Be'):
        gp_run = gaopaw.Runner()
        assert gp_run.num_obj_fns == 20
    with fileutils.chdir('FeO_workdir/Fe'):
        gp_run = gaopaw.Runner()
        assert gp_run.num_obj_fns == 8
    with fileutils.chdir('La_workdir/La'):
        gp_run = gaopaw.Runner()
        assert gp_run.num_obj_fns == 3


def test_numUserObjs():
    with fileutils.chdir('ex_templates/Empty'):
        gp_run = gaopaw.Runner()
        assert gp_run.numUserObjs() == 9
    with fileutils.chdir('Be_workdir/Be'):
        gp_run = gaopaw.Runner()
        assert gp_run.numUserObjs() == 20
    with fileutils.chdir('La_workdir/La'):
        gp_run = gaopaw.Runner()
        assert gp_run.numUserObjs() == 3

def test_formCmpdDict():
    with fileutils.chdir('ex_templates/Empty'):
        gp_run = gaopaw.Runner()
    assert gp_run.formCmpdDict() == \
        {'Si': {'SC': {'lattice_constant': {}, 'phonon_frequency': {}}}, 
        'SiC': {'ZB': {'band_gap': {}, 'eos': {}, 'lattice_constant': {}}}, 
        'SiO': {'RS': {'lattice_constant': {}}}}

def test_getElementList():
    with fileutils.chdir('ex_templates/Empty'):
        gp_run = gaopaw.Runner()
        assert Counter(gp_run.element_list) == Counter(['C'])
    with fileutils.chdir('Be_workdir/Be'):
        gp_run = gaopaw.Runner()
        assert Counter(gp_run.element_list) == Counter(['Be', 'O', 'S', 'Al', 'B'])

def test_parseElems():
    assert Counter(gp_run.parseElems('H2O2Na2FeCo34Ag5OSH2')) == \
        Counter(['H', 'O', 'Na', 'Fe', 'Co', 'Ag', 'S'])
    assert Counter(gp_run.parseElems('O')) == Counter(['O'])
    assert Counter(gp_run.parseElems('HOSFe')) == \
        Counter(['H', 'O', 'S', 'Fe'])

def run_atompaw_mock(elem):
    elem_path = os.path.join(os.pardir, os.pardir, os.pardir, 'ex_templates', elem)
    for fname in gaopaw.os.listdir(elem_path):
        shutil.copyfile(os.path.join(elem_path, fname), fname)

def test_testElementList():
    with fileutils.chdir(os.path.join('Si_workdir', 'Si')):
        gp_run = gaopaw.Runner()
        for fname in os.listdir():
            if os.path.isdir(fname):
                shutil.rmtree(fname)
            else:
                os.remove(fname)
        shutil.copyfile(os.path.join(os.pardir, 'params.in'), 'params.in')
        with mock.patch.object(gp_run, 'runAtompaw') as run_mock:
            run_mock.side_effect = run_atompaw_mock
            assert gp_run.testElementList()[1] == False

def test_formElementDict():
    with fileutils.chdir(os.path.join('Si_workdir', 'Si')):
        gp_run = gaopaw.Runner()
        elem_dict = gp_run.formElementDict()
        assert Counter(elem_dict.keys()) == Counter(['C', 'Si', 'O'])
        assert Counter(elem_dict['Si'].keys()) == Counter(['elemental','FCC','BCC'])
        assert list(elem_dict['Si']['FCC'].keys()) == ['lattice_constant']

def test_writeAtompawInput():
    with fileutils.chdir(os.path.join('ex_templates', 'Empty')):
        gp_run = gaopaw.Runner()
        if os.path.exists('Si.atompaw.in'):
            os.remove('Si.atompaw.in')
        gp_run.writeAtompawInput('Si')
        assert os.path.exists('Si.atompaw.in')
        os.remove('Si.atompaw.in')

def test_checkUpf():
    with fileutils.chdir(os.path.join('ex_templates', 'Empty')):
        gp_run = gaopaw.Runner()
        assert gp_run.checkUpf() == True
    assert gp_run.checkUpf() == False

def test_compareLog():
    with fileutils.chdir(os.path.join('Log_tests', 'workingdir')):
        gp_run = gaopaw.Runner()
        assert np.isclose(gp_run.compareLog(), 0.03425, rtol=1e-3)

def test_updateStructure():
    with fileutils.chdir(os.path.join('Si_workdir', 'struct_test')):
        gp_run = gaopaw.Runner()
        gp_run.updateStructure('Si', 'diamond', 'scf')
        assert os.path.exists('Si.diamond.scf.in')

def test_checkConvergence():
    with fileutils.chdir(os.path.join('Convergence_tests', 'good')):
        gp_run = gaopaw.Runner()
        assert gp_run.checkConvergence('Si', 'FCC', 'relax') == True
    with fileutils.chdir(os.path.join('Convergence_tests', 'bad')):
        gp_run = gaopaw.Runner()
        assert gp_run.checkConvergence('Si', 'FCC', 'relax') == False

def test_compareAtoms():
    with fileutils.chdir(os.path.join('N_workdir', 'N')):
        gp_run = gaopaw.Runner()
        elem_template_dir = gp_run.input_settings.directories.elem_template_dir
        atom_diff = gp_run.compareAtoms('N', 'SC', elem_template_dir)
        assert np.isclose(atom_diff, 0.042629, rtol=1e-3)

def test_getLatticeConstant():
    with fileutils.chdir(os.path.join('Conv_Tests', 'Si_BCC_2')):
        gp_run = gaopaw.Runner()
        assert np.isclose(
            gp_run.getLatticeConstant('Si', 'BCC'), 3.085, rtol=1e-3)
    with fileutils.chdir(os.path.join('Prim_tests', 'Tetrag_I')):
        assert np.array_equal(
            np.array(gp_run.getLatticeConstant('Si', 'tetrag')).round(3), 
            [1.588, 1.746])
    with fileutils.chdir(os.path.join('Prim_tests', 'Monoclin_C')):
        assert np.array_equal(
            np.array(gp_run.getLatticeConstant('Si', 'monoclin')).round(3), 
            [1.588, 1.746, 1.905, 60.0])
    with fileutils.chdir(os.path.join('Prim_tests', 'Ortho_C')):
        assert np.array_equal(
            np.array(gp_run.getLatticeConstant('Si', 'ortho')).round(3),
            [3.175, 3.493, 3.810])
    with fileutils.chdir(os.path.join('Prim_tests', 'Ortho_F')):
        assert np.array_equal(
            np.array(gp_run.getLatticeConstant('Si', 'ortho')).round(3),
            [3.175, 3.493, 3.810])
    with fileutils.chdir(os.path.join('Prim_tests', 'Ortho_I')):
        assert np.array_equal(
            np.array(gp_run.getLatticeConstant('Si', 'ortho')).round(3),
            [3.175, 3.493, 3.810])

def test_getCell():
    with fileutils.chdir(os.path.join('Prim_tests', 'Tetrag_I')):
        gp_run = gaopaw.Runner()
        assert np.array_equal(
        np.array(gp_run.getCell('Si', 'tetrag', 'relax')).round(3), 
        [[1.50, -1.50, 1.65], [1.50, 1.50, 1.65], [-1.50, -1.50, 1.65]])

def test_getMag():
    with fileutils.chdir(os.path.join('FeO_workdir', 'Fe')):
        gp_run = gaopaw.Runner()
        assert np.isclose(gp_run.getMag('FeO', 'RS'), 3.83, rtol=1e-3)

def test_checkIfElem():
    with fileutils.chdir(os.path.join('FeO_workdir', 'Fe')):
        gp_run = gaopaw.Runner()
        assert gp_run.is_elem == False
    with fileutils.chdir(os.path.join('N_workdir', 'N')):
        gp_run = gaopaw.Runner()
        assert gp_run.is_elem == True

def test_mergeDicts():
    with fileutils.chdir(os.path.join('FeO_workdir', 'Fe')):
        gp_run = gaopaw.Runner()
    assert gp_run.mergeDicts({'A': 1}, {'B': 2}) == {'A': 1, 'B': 2}
    assert gp_run.mergeDicts({'A': {'B': 1}}, {'A': {'C': 2}}) == {'A': {'B': 1, 'C': 2}}

def test_testCmpdList():
    with fileutils.chdir(os.path.join('Si_workdir', 'Si')):
        gp_run = gaopaw.Runner()
    with mock.patch.object(gp_run, 'testProperty') as run_mock:
        run_mock.return_value = [0, False]
        assert gp_run.testCmpdList()[1] == False

def test_testPhaseStability():
    with fileutils.chdir('Polymorphs'):
        gp_run = gaopaw.Runner('current')
        assert gp_run.testPhaseStability('Be', ['FCC','BCC'], None) == 0.0
        assert gp_run.testPhaseStability('Be', ['BCC','FCC'], None) == 100.0

def test_testProperty():
    with fileutils.chdir(os.path.join('BeO_workdir', 'BeO')):
        gp_run = gaopaw.Runner()
        assert np.isclose(
            gp_run.testProperty('BeO', 'RS', 'lattice_constant', 1.0, '.')[0], 
            2.625, rtol=1e-3)
        assert np.isclose(
            gp_run.testProperty('BeO', 'RS', 'band_gap', 1.0, '.')[0],
            7.364, rtol=1e-3)
    with fileutils.chdir(os.path.join('BeS_workdir', 'BeS')):
        gp_run = gaopaw.Runner()
        with mock.patch.object(gp_run, 'runScaleLat') as run_mock:
            assert np.isclose(
                gp_run.testProperty('BeS', 'ZB', 'eos', [1.0, 1.0, 1.0], '.')[0],
                1285.706, rtol=1e-1)
            assert np.isclose(
                gp_run.testProperty('BeS', 'ZB', 'bulk_modulus', 1.0, '.')[0],
                91.903, rtol=1e-2)
    with fileutils.chdir(os.path.join('Be.SC_workdir', 'Be')):
        gp_run = gaopaw.Runner()
        with mock.patch.object(gp_run, 'runPhonon') as run_mock:
            assert np.isclose(
                gp_run.testProperty('Be', 'SC', 'phonon_frequency', 
                [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], '.')[0],
                12.666, rtol=1e-2)
    with fileutils.chdir(os.path.join('N_workdir', 'N')):
        gp_run = gaopaw.Runner()
        elem_template_dir = gp_run.input_settings.directories.elem_template_dir
        assert np.isclose(
            gp_run.testProperty('N', 'SC', 'atomic_positions', 1.0, elem_template_dir)[0],
            0.0426, rtol=1e-3)
    with fileutils.chdir(os.path.join('FeO_workdir', 'Fe')):
        gp_run = gaopaw.Runner()
        assert np.isclose(
            gp_run.testProperty('FeO', 'RS', 'magnetization', 1.0, elem_template_dir)[0],
            2.83, rtol=1e-3)
        assert np.isclose(
            gp_run.testProperty('FeO', 'RS', 'magnetic_moment', [1.0, 1.0], elem_template_dir)[0],
            1.467, rtol=1e-3)

def qe_mock(placeholder):
    shutil.copyfile(os.path.join(os.pardir, os.pardir, 'ex_templates', 'EOS',
        'output_%s' % qe_mock.counter), './ZnO.ZB.relax.out')
    qe_mock.counter += 1

def test_runScaleLat():
    qe_mock.counter = 1
    with fileutils.chdir('EOS_Test'):
        gp_run = gaopaw.Runner()
        for fname in glob.iglob('ZnO_*'):
            shutil.rmtree(fname)
        if os.path.exists('E_V.txt'):
            os.remove('E_V.txt')
        with mock.patch.object(gp_run, 'runCurrentQE') as run_mock:
            run_mock.side_effect = qe_mock
            gp_run.runScaleLat('ZnO', 'ZB')
        assert os.path.exists('E_V.txt')
        assert np.isclose(
            np.loadtxt('E_V.txt').transpose()[1][0], 
            0.94*np.loadtxt('E_V.txt').transpose()[1][3], rtol=1e-3)
        assert np.isclose(
            np.loadtxt('E_V.txt').transpose()[1][-1],
            1.06*np.loadtxt('E_V.txt').transpose()[1][3], rtol=1e-3)

def test_updateObjFile():
    diff_dict = \
        {'Si': {'SC': {'lattice_constant': 0.01, 'phonon_frequency': 0.2}},
        'SiC': {'ZB': {'band_gap': 0.03}}}
    with fileutils.chdir('Si_workdir'):
        gp_run = gaopaw.Runner()
        if os.path.exists('Detailed_Results'):
            os.remove('Detailed_Results')
        gp_run.updateObjFile(diff_dict)
        assert os.path.exists('Detailed_Results')
        df = pd.read_table('Detailed_Results', sep=':', header=None)
        assert np.array_equal(list(df[0]),
            ['Si_SC_lattice_constant', 'Si_SC_phonon_frequency', 'SiC_ZB_band_gap'])

def test_dictToList():
    with fileutils.chdir('Si_workdir'):
        gp_run = gaopaw.Runner()
    diff_dict = \
        {'Si': {'SC': {'lattice_constant': 0.01, 'phonon_frequency': 0.2}},
        'SiC': {'ZB': {'band_gap': 0.03}}}
    assert np.array_equal(gp_run.dictToList(diff_dict),
        [[0.01, 0.2, 0.03], ['Si_SC_lattice_constant', 'Si_SC_phonon_frequency', 
        'SiC_ZB_band_gap']])

def test_updateDakota():
    with fileutils.chdir(os.path.join('N_workdir', 'N')):
        if os.path.exists('results.out'):
            os.remove('results.out')
        gp_run = gaopaw.Runner()
        diff_dict = \
            {'N': {'SC': {'atomic_positions': 1.0}, 'elemental': {'log': 1.0}}}
        gp_run.updateDakota(diff_dict)
        assert os.path.exists('results.out')

def test_badRun():
    with fileutils.chdir(os.path.join('Si_workdir', 'Si')):
        if os.path.exists('results.out'):
            os.remove('results.out')
        gp_run = gaopaw.Runner()
        gp_run.badRun()
        assert os.path.exists('results.out')

def test_writeDakota():
    with fileutils.chdir('dakota_updates'):
        gp_run = gaopaw.Runner(input_dir='current')
        gp_run.updateNumObjs()
        gp_run.updateVars()
        gp_run.updateLabels()
        assert gp_run.numDakotaObjs() == gp_run.numUserObjs()

def test_getBestSoln():
    with fileutils.chdir('Best_Soln'):
        if os.path.isdir('Best_Solution'):
            shutil.rmtree('Best_Solution')
        gp_run = gaopaw.Runner(input_dir='current')
        with mock.patch.object(gp_run, 'runAtompaw') as run_mock:
            gp_run.getBestSoln()
        assert os.path.isdir('Best_Solution')
        with fileutils.chdir('Best_Solution'):
            assert os.path.exists('Be.atompaw.in')
            with open('Detailed_Results') as results:
                assert len(results.readlines()) == \
                gp_run.numUserObjs()

def test_getAtompawEnergies():
    with fileutils.chdir('ap_energies'):
        gp_run = gaopaw.Runner(input_dir='current')
        pseud, ae = gp_run.getAtompawEnergies('Be')
        assert np.isclose(pseud, -29.2657434130951, rtol=1e-12)
        assert np.isclose(ae, -29.2657430162565, rtol=1e-12)

