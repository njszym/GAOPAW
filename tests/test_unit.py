import mock
import shutil
from collections import Counter
import numpy as np
from schrodinger.utils import fileutils
import os
from unittest.mock import Mock
from types import SimpleNamespace
import gaopaw
import sys

top_dir = os.getcwd()

with gaopaw.fileutils.chdir('ex_templates/Empty'):
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
    with gaopaw.fileutils.chdir('ex_templates/Empty'):
        for fname in paw_list:
            if os.path.exists(fname):
                os.remove(fname)
        gp_run.getPaws()
    for fname in paw_list:
        assert os.path.exists(os.path.join('ex_templates/Empty',fname))

def test_numDakotaObjs():
    with gaopaw.fileutils.chdir('ex_templates/Empty'):
        gp_run = gaopaw.Runner()
        assert gp_run.num_obj_fns == 15
    with gaopaw.fileutils.chdir('Be_workdir/Be'):
        gp_run = gaopaw.Runner()
        assert gp_run.num_obj_fns == 20
    with gaopaw.fileutils.chdir('FeO_workdir/Fe'):
        gp_run = gaopaw.Runner()
        assert gp_run.num_obj_fns == 8
    with gaopaw.fileutils.chdir('La_workdir/La'):
        gp_run = gaopaw.Runner()
        assert gp_run.num_obj_fns == 3


def test_numUserObjs():
    with gaopaw.fileutils.chdir('ex_templates/Empty'):
        gp_run = gaopaw.Runner()
        assert gp_run.numUserObjs() == 15
    with gaopaw.fileutils.chdir('Be_workdir/Be'):
        gp_run = gaopaw.Runner()
        assert gp_run.numUserObjs() == 20
    with gaopaw.fileutils.chdir('La_workdir/La'):
        gp_run = gaopaw.Runner()
        assert gp_run.numUserObjs() == 3

def test_formCmpdDict():
    with gaopaw.fileutils.chdir('ex_templates/Empty'):
        gp_run = gaopaw.Runner()
    assert gp_run.formCmpdDict() == \
        {'Si': {'SC': {'lattice_constant': {}, 'phonon_frequency': {}}}, 
        'SiC': {'ZB': {'band_gap': {}, 'eos': {}, 'lattice_constant': {}}}, 
        'SiO': {'RS': {'lattice_constant': {}}}}

def test_getElementList():
    with gaopaw.fileutils.chdir('ex_templates/Empty'):
        gp_run = gaopaw.Runner()
        assert Counter(gp_run.element_list) == Counter(['C', 'Si', 'O'])
    with gaopaw.fileutils.chdir('Be_workdir/Be'):
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

