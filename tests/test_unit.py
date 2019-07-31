import mock
import numpy as np
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
    with gaopaw.fileutils.chdir('ex_templates/Empty'):
        gp_run.getPaws()
    for fname in ['%s.GGA-PBE-paw.UPF' % elem for elem in ['Si', 'O']]:
        assert os.path.exists(os.path.join('ex_templates/Empty',fname))
