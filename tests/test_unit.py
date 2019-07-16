import mock
from types import SimpleNamespace
import gaopaw

working_dir = '/scr/szymansk/gaopaw/tests/sample_workdir'
elem_template_dir = '/scr/szymansk/gaopaw/Elem_Templates'

def test_parse_num_objs():
    assert gaopaw.parse_num_objs(working_dir) == 15

def test_parse_elems():
    assert gaopaw.parse_elems('O') == ['O']
    assert gaopaw.parse_elems('HOSFe') == ['H','O','S','Fe']
    assert gaopaw.parse_elems('Ag2S') == ['Ag','S']
    assert gaopaw.parse_elems('H2O2Na2FeCo34Ag5OSH2') == ['H','O','Na','Fe','Co','Ag','S']

def test_unique():
    assert gaopaw.unique(['Be','Hf','Hf','O','S','Be','Hf','Ga','O']) == ['Be','Hf','O','S','Ga']
    assert gaopaw.unique([1,4,1,2,2,3,5,1]) == [1,4,2,3,5]
    assert gaopaw.unique([12.000001,12.0]) == [12.0]

def test_get_num_objs():
    cmpd_list = \
    [SimpleNamespace(formula='Si', lattice_constant=3.08, lattice_type='SC', phonon_frequency=[0.0, 0.0, 0.0, 6.15, 6.15, 6.15]),
    SimpleNamespace(formula='SiO', lattice_constant=4.616, lattice_type='RS'),
    SimpleNamespace(band_gap=1.24, eos=[10.50197, 212.71678, 3.692], formula='SiC', lattice_constant=4.38, lattice_type='ZB')]
    assert gaopaw.get_num_objs(cmpd_list,['Si','O','C']) == 15
    cmpd_list = [SimpleNamespace(formula='Si')]
    assert gaopaw.get_num_objs(cmpd_list,['Si']) == 3
    cmpd_list = [SimpleNamespace(formula='N')]
    assert gaopaw.get_num_objs(cmpd_list,['N']) == 2
    cmpd_list = [SimpleNamespace(formula='P')]
    assert gaopaw.get_num_objs(cmpd_list,['P']) == 2
    cmpd_list = [SimpleNamespace(formula='N',lattice_type='triclinic',magnetization=3.0)]
    assert gaopaw.get_num_objs(cmpd_list,['N']) == 3

def test_form_cmpd_dict():
    cmpd_list = \
    [SimpleNamespace(formula='Si', lattice_constant=3.08, lattice_type='SC', phonon_frequency=[0.0, 0.0, 0.0, 6.15, 6.15, 6.15]),
    SimpleNamespace(formula='SiO', lattice_constant=4.616, lattice_type='RS'),
    SimpleNamespace(band_gap=1.24, eos=[10.50197, 212.71678, 3.692], formula='SiC', lattice_constant=4.38, lattice_type='ZB')]
    assert gaopaw.form_cmpd_dict(cmpd_list) == \
    {'Si': {'SC': {'lattice_constant': {}, 'phonon_frequency': {}}},
    'SiC': {'ZB': {'band_gap': {}, 'eos': {}, 'lattice_constant': {}}},
    'SiO': {'RS': {'lattice_constant': {}}}}
    cmpd_list = [SimpleNamespace(formula='N',lattice_type='triclinic',magnetization=3.0)]
    assert gaopaw.form_cmpd_dict(cmpd_list) == {'N': {'triclinic': {'magnetization': {}}}}

def test_get_element_info():
    assert round(gaopaw.get_element_info(elem_template_dir)['Ag']['FCC'],3) == 4.158
    assert round(gaopaw.get_element_info(elem_template_dir)['Ag']['BCC'],3) == 3.299
    assert round(gaopaw.get_element_info(elem_template_dir)['Zr']['FCC'],3) == 4.522
    assert round(gaopaw.get_element_info(elem_template_dir)['Zr']['BCC'],3) == 3.569

@mock.patch.object(gaopaw,'run_atompaw')
def test_test_element_list(run_atompaw_mock):
    with gaopaw.fileutils.chdir(working_dir):
        if gaopaw.os.path.isdir('Si'):
            gaopaw.shutil.rmtree('Si')
        if gaopaw.os.path.exists('Si.atompaw.in'):
            gaopaw.os.remove('Si.atompaw.in')
        assert gaopaw.test_element_list(['Si'], elem_template_dir)[1] == False

#@mock.patch.object(gaopaw.os, 'mkdir')
#def test_element_list(mkdir_mock):
#    assert gaopaw.test_element_list(['Si'], elem_template_dir)

#def test_mock():
#    with mock.patch.object(gaopaw, 'run_QE') as run_QE_mock:
#        run_QE_mock.return_value = False
#        print(gaopaw.test_element_mock())
