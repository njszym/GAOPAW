import mock
from unittest.mock import Mock
from types import SimpleNamespace
import gaopaw


def working_dir(cmpd):
    return gaopaw.os.getcwd()+'/'+cmpd+'_workdir'

elem_template_dir = gaopaw.os.getcwd()+'/../Elem_Templates'
ex_template_dir = gaopaw.os.getcwd()+'/ex_templates'

def test_parse_num_objs():
    assert gaopaw.parse_num_objs(ex_template_dir+'/Empty') == 15

def test_parse_elems():
    assert set(gaopaw.parse_elems('O')) == {'O'}
    assert set(gaopaw.parse_elems('HOSFe')) == {'H', 'O', 'Fe', 'S'}
    assert set(gaopaw.parse_elems('Ag2S')) == {'Ag', 'S'}
    assert set(gaopaw.parse_elems('H2O2Na2FeCo34Ag5OSH2')) == \
        {'H', 'O', 'Na', 'Fe', 'Co', 'Ag', 'S'}

def test_get_num_objs():
    cmpd_list = \
    [SimpleNamespace(formula='Si', lattice_constant=3.08, 
        lattice_type='SC', phonon_frequency=[0.0, 0.0, 0.0, 6.15, 6.15, 6.15]), 
    SimpleNamespace(formula='SiO', lattice_constant=4.616, lattice_type='RS'), 
    SimpleNamespace(band_gap=1.24, eos=[10.50197, 212.71678, 3.692], 
        formula='SiC', lattice_constant=4.38, lattice_type='ZB')]
    assert gaopaw.get_num_objs(cmpd_list, ['Si', 'O', 'C']) == 15
    cmpd_list = [SimpleNamespace(formula='Si')]
    assert gaopaw.get_num_objs(cmpd_list, ['Si']) == 3
    cmpd_list = [SimpleNamespace(formula='N')]
    assert gaopaw.get_num_objs(cmpd_list, ['N']) == 2
    cmpd_list = [SimpleNamespace(formula='P')]
    assert gaopaw.get_num_objs(cmpd_list, ['P']) == 2
    cmpd_list = [SimpleNamespace(formula='N', lattice_type='triclinic', 
        magnetization=3.0)]
    assert gaopaw.get_num_objs(cmpd_list, ['N']) == 3

def test_form_cmpd_dict():
    cmpd_list = \
    [SimpleNamespace(formula='Si', lattice_constant=3.08, 
        lattice_type='SC', phonon_frequency=[0.0, 0.0, 0.0, 6.15, 6.15, 6.15]), 
    SimpleNamespace(formula='SiO', lattice_constant=4.616, lattice_type='RS'), 
    SimpleNamespace(band_gap=1.24, eos=[10.50197, 212.71678, 3.692], 
        formula='SiC', lattice_constant=4.38, lattice_type='ZB')]
    assert gaopaw.form_cmpd_dict(cmpd_list) == \
    {'Si': {'SC': {'lattice_constant': {}, 'phonon_frequency': {}}}, 
    'SiC': {'ZB': {'band_gap': {}, 'eos': {}, 'lattice_constant': {}}}, 
    'SiO': {'RS': {'lattice_constant': {}}}}
    cmpd_list = [SimpleNamespace(formula='N', 
        lattice_type='triclinic', magnetization=3.0)]
    assert gaopaw.form_cmpd_dict(cmpd_list) == {'N': 
        {'triclinic': {'magnetization': {}}}}

def test_get_element_info():
    assert gaopaw.np.allclose(
        gaopaw.get_element_info(elem_template_dir)['Ag']['FCC'], 4.158, atol=0.001)
    assert gaopaw.np.allclose(
        gaopaw.get_element_info(elem_template_dir)['Ag']['BCC'], 3.299, atol=0.001)
    assert gaopaw.np.allclose(
        gaopaw.get_element_info(elem_template_dir)['Zr']['FCC'], 4.522, atol=0.001)
    assert gaopaw.np.allclose(
        gaopaw.get_element_info(elem_template_dir)['Zr']['BCC'], 3.569, atol=0.001)

def test_compare_log():
    with gaopaw.fileutils.chdir('Log_tests'):
        assert round(gaopaw.compare_log(), 5) == 0.00334

def run_atompaw_mock(elem):
    for fname in gaopaw.os.listdir(ex_template_dir+'/'+elem+'/'):
        gaopaw.shutil.copyfile(ex_template_dir+'/'+elem+'/'+fname, './'+fname)

def test_test_element_list():
    for elem in ['Si', 'N', 'Be']:
        with gaopaw.fileutils.chdir(working_dir(elem)):
            if gaopaw.os.path.isdir(elem):
                gaopaw.shutil.rmtree(elem)
            with mock.patch.object(gaopaw, 'run_atompaw') as run_mock:
                run_mock.side_effect = run_atompaw_mock
                if elem == 'Be':
                    assert gaopaw.test_element_list([elem], elem_template_dir)[1] == True
                else:
                    assert gaopaw.test_element_list([elem], elem_template_dir)[1] == False

def test_test_cmpd_list():
    cmpd_list = \
    [SimpleNamespace(formula='Si', lattice_constant=3.08, 
        lattice_type='SC', phonon_frequency=[0.0, 0.0, 0.0, 6.15, 6.15, 6.15]), 
    SimpleNamespace(formula='SiO', lattice_constant=4.616, lattice_type='RS'), 
    SimpleNamespace(band_gap=1.24, eos=[10.50197, 212.71678, 3.692], 
        formula='SiC', lattice_constant=4.38, lattice_type='ZB')]
    cmpd_diff_dict = \
    {'Si': {'SC': {'lattice_constant': {}, 'phonon_frequency': {}}}, 
    'SiC': {'ZB': {'band_gap': {}, 'eos': {}, 'lattice_constant': {}}}, 
    'SiO': {'RS': {'lattice_constant': {}}}}
    with mock.patch.object(gaopaw, 'test_property') as run_mock:
        run_mock.return_value = [0, False]
        assert gaopaw.test_cmpd_list(cmpd_list, cmpd_diff_dict, None, elem_template_dir)[1] == False

def test_merge_dicts():
    assert gaopaw.merge_dicts({'A': 1}, {'B': 2}) == {'A': 1, 'B': 2}
    assert gaopaw.merge_dicts({'A': {'B': 1}}, {'A': {'C': 2}}) == {'A': {'B': 1, 'C': 2}}

def test_test_property():
    with gaopaw.fileutils.chdir(working_dir('BeO')):
        assert round(gaopaw.test_property('BeO', 'RS', 'lattice_constant', 1.0, 
            gaopaw.os.getcwd())[0], 3) == 2.625
        assert gaopaw.test_property('BeO', 'RS', 'lattice_constant', 1.0, 
            gaopaw.os.getcwd())[1] == False
        assert round(gaopaw.test_property('BeO', 'RS', 'band_gap', 1.0, 
            gaopaw.os.getcwd())[0], 3) == 7.364
        assert gaopaw.test_property('BeO', 'RS', 'band_gap', 1.0, 
            gaopaw.os.getcwd())[1] == False
    with gaopaw.fileutils.chdir(working_dir('BeS')):
        with mock.patch.object(gaopaw, 'run_scale_lat') as run_mock:
            assert round(gaopaw.test_property('BeS', 'ZB', 'eos', [1.0, 1.0, 1.0], 
                gaopaw.os.getcwd())[0], 3) == 1285.706
            assert gaopaw.test_property('BeS', 'ZB', 'eos', [1.0, 1.0, 1.0], 
                gaopaw.os.getcwd())[1] == False
            assert round(gaopaw.test_property('BeS', 'ZB', 'bulk_modulus', 1.0, 
                gaopaw.os.getcwd())[0], 3) == 91.903
            assert gaopaw.test_property('BeS', 'ZB', 'bulk_modulus', 1.0, 
                gaopaw.os.getcwd())[1] == False
    with gaopaw.fileutils.chdir(working_dir('Be.SC')):
        with mock.patch.object(gaopaw, 'run_phonon') as run_mock:
            assert round(gaopaw.test_property('Be', 'SC', 'phonon_frequency', 
                [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], gaopaw.os.getcwd())[0], 3) == 12.666
            assert gaopaw.test_property('Be', 'SC', 'phonon_frequency', 
                [0.0, 0.0, 0.0, 1.0, 1.0, 1.0], gaopaw.os.getcwd())[1] == False
    with gaopaw.fileutils.chdir(working_dir('N')):
        assert round(gaopaw.test_property('N', 'SC', 'atomic_positions', 1.0, 
            elem_template_dir)[0], 3) == 0.035
        assert gaopaw.test_property('N', 'SC', 'atomic_positions', 1.0, 
            elem_template_dir)[1] == False
    with gaopaw.fileutils.chdir(working_dir('FeO')):
        assert round(gaopaw.test_property('FeO', 'RS', 'magnetization', 1.0, 
            gaopaw.os.getcwd())[0], 3) == 2.83
        assert gaopaw.test_property('FeO', 'RS', 'magnetization', 1.0, 
            gaopaw.os.getcwd())[1] == False
        assert round(gaopaw.test_property('FeO', 'RS', 'magnetic_moment', [1.0, 1.0], 
            gaopaw.os.getcwd())[0], 3) == 1.467
        assert gaopaw.test_property('FeO', 'RS', 'magnetic_moment', [1.0, 1.0], 
            gaopaw.os.getcwd())[1] == False

def test_check_upf():
    with gaopaw.fileutils.chdir(ex_template_dir+'/Empty'):
        assert gaopaw.check_upf() == False
    with gaopaw.fileutils.chdir(ex_template_dir+'/Si'):
        assert gaopaw.check_upf() == True

def test_bad_run():
    with gaopaw.fileutils.chdir(ex_template_dir+'/bad_run_ex'):
        gaopaw.bad_run()
        results = gaopaw.pd.read_table('results.out', sep='\s+', header=None)[0]
        assert len(results) == 3
        for value in results:
            assert value == 100.0

def test_compare_lat():
    with gaopaw.fileutils.chdir(working_dir('Lattice')):
        assert round(gaopaw.compare_lat([1.0, 1.0], 'H', 'hex'), 3) == 3.494
        assert round(gaopaw.compare_lat(1.0, 'Si', 'FCC'), 3) == 2.858
        assert round(gaopaw.compare_lat(1.0, 'Si', 'BCC'), 3) == 2.078
        assert round(gaopaw.compare_lat([1.0, 1.0], 'Si', 'tetrag'), 3) == 2.694
        assert round(gaopaw.compare_lat([1.0, 1.0, 1.0], 'Si', 'ortho'), 3) == 2.510
        assert round(gaopaw.compare_lat([1.0, 1.0], 'Si', 'rhomb'), 3) == 30.241
        assert round(gaopaw.compare_lat([1.0, 1.0, 1.0, 1.0], 'Si', 'monoclin'), 3) == 14.650
        assert round(gaopaw.compare_lat([1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 'Si', 'triclin'), 3) == 22.142

def test_get_lattice_constant():
    with gaopaw.fileutils.chdir('Conv_Tests/Si_BCC_2'):
        assert round(gaopaw.get_lattice_constant('Si', 'BCC'), 3) == 3.085
    with gaopaw.fileutils.chdir('Prim_tests/Tetrag_I'):
        assert gaopaw.get_lattice_constant('Si', 'tetrag')[0] == 1.588
        assert gaopaw.get_lattice_constant('Si', 'tetrag')[1] == 1.746
    with gaopaw.fileutils.chdir('Prim_tests/Monoclin_C'):
        assert gaopaw.get_lattice_constant('Si', 'monoclin')[0] == 1.588
        assert gaopaw.get_lattice_constant('Si', 'monoclin')[1] == 1.746
        assert gaopaw.get_lattice_constant('Si', 'monoclin')[2] == 1.905
        assert gaopaw.get_lattice_constant('Si', 'monoclin')[3] == 60.0
    with gaopaw.fileutils.chdir('Prim_tests/Ortho_C'):
        assert gaopaw.get_lattice_constant('Si', 'ortho')[0] == 3.175
        assert gaopaw.get_lattice_constant('Si', 'ortho')[1] == 3.493
        assert gaopaw.get_lattice_constant('Si', 'ortho')[2] == 3.810
    with gaopaw.fileutils.chdir('Prim_tests/Ortho_F'):
        assert gaopaw.get_lattice_constant('Si', 'ortho')[0] == 3.175
        assert gaopaw.get_lattice_constant('Si', 'ortho')[1] == 3.493
        assert gaopaw.get_lattice_constant('Si', 'ortho')[2] == 3.810
    with gaopaw.fileutils.chdir('Prim_tests/Ortho_I'):
        assert gaopaw.get_lattice_constant('Si', 'ortho')[0] == 3.175
        assert gaopaw.get_lattice_constant('Si', 'ortho')[1] == 3.493
        assert gaopaw.get_lattice_constant('Si', 'ortho')[2] == 3.810

def test_update_dakota():
    diff_dict = \
    {'Si': {'BCC': {'lattice_constant': 0.01, 'phonon_frequency': 0.2}},
    'SiC': {'ZB': {'band_gap': 0.03}}}
    with gaopaw.fileutils.chdir(working_dir('Si')):
        gaopaw.update_dakota(diff_dict)
        assert gaopaw.os.path.exists('results.out')

def test_check_convergence():
    with gaopaw.fileutils.chdir(working_dir('Conv')):
        assert gaopaw.check_convergence('Si', 'FCC', 'relax') == False
    with gaopaw.fileutils.chdir(working_dir('Lattice')):
        assert gaopaw.check_convergence('Si', 'FCC', 'relax') == True

def test_update_structure():
    with gaopaw.fileutils.chdir(working_dir('Si')):
        gaopaw.update_structure('Si', 'diamond', 'scf')

def test_write_cell():
    with gaopaw.fileutils.chdir(working_dir('Si')):
        gaopaw.write_cell('Si', 'diamond', [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        with open('Si.diamond.relax.in') as qe_input:
            assert qe_input.readlines()[-1] == '0.0 0.0 1.0 \n'

def test_dict_to_list():
    diff_dict = \
    {'Si': {'SC': {'lattice_constant': 0.01, 'phonon_frequency': 0.2}}, 
    'SiC': {'ZB': {'band_gap': 0.03}}}
    assert gaopaw.dict_to_list(diff_dict)[0] == [0.01, 0.2, 0.03]
    assert gaopaw.dict_to_list(diff_dict)[1] == \
    ['Si_SC_lattice_constant', 'Si_SC_phonon_frequency', 'SiC_ZB_band_gap']

def test_scale_cell():
    with gaopaw.fileutils.chdir(working_dir('Si')):
        scaled_cell = gaopaw.scale_cell('Si', 'diamond', 0.9)
        assert len(scaled_cell) == 3
        for vec in scaled_cell:
            assert len(vec) == 3
        scaled_matrix = gaopaw.np.matrix(scaled_cell)
        assert round(gaopaw.np.linalg.det(scaled_matrix), 3) == 243.096

def qe_mock(placeholder):
    gaopaw.shutil.copyfile(gaopaw.os.path.join(ex_template_dir,'EOS',
        'output_%s' % qe_mock.counter),'./ZnO.ZB.relax.out')
    qe_mock.counter += 1

def test_run_scale_lat():
    qe_mock.counter = 1
    with gaopaw.fileutils.chdir('EOS_Test'):
        for fname in gaopaw.glob.glob('ZnO_*'):
            gaopaw.shutil.rmtree(fname)
        if gaopaw.os.path.exists('E_V.txt'):
            gaopaw.os.remove('E_V.txt')
        with mock.patch.object(gaopaw, 'run_qe_as_is') as run_mock:
            run_mock.side_effect = qe_mock
            gaopaw.run_scale_lat('ZnO','ZB',None)
        assert gaopaw.os.path.exists('E_V.txt')
        assert min(gaopaw.np.loadtxt('E_V.txt').transpose()[0]) \
            == gaopaw.np.loadtxt('E_V.txt').transpose()[0][3]
        assert round(gaopaw.np.loadtxt('E_V.txt').transpose()[1][0],3) == \
            round(0.94*gaopaw.np.loadtxt('E_V.txt').transpose()[1][3],3)
        assert round(gaopaw.np.loadtxt('E_V.txt').transpose()[1][-1],3) == \
            round(1.06*gaopaw.np.loadtxt('E_V.txt').transpose()[1][3],3)

def test_update_obj_file():
    diff_dict = \
    {'Si': {'SC': {'lattice_constant': 0.01, 'phonon_frequency': 0.2}}, 
    'SiC': {'ZB': {'band_gap': 0.03}}}
    with gaopaw.fileutils.chdir(working_dir('Si')):
        if gaopaw.os.path.exists('Detailed_Results'):
            gaopaw.os.remove('Detailed_Results')
        gaopaw.update_obj_file(diff_dict)
        assert gaopaw.os.path.exists('Detailed_Results')
        df = gaopaw.pd.read_table('Detailed_Results', sep=':', header=None)
        assert gaopaw.np.array_equal(list(df[0]), 
            ['Si_SC_lattice_constant', 'Si_SC_phonon_frequency', 'SiC_ZB_band_gap'])
        assert gaopaw.np.array_equal(list(df[1]),
            [' 1.0%', ' 0.2 THz', ' 0.03 eV'])


def test_calc_obj_fn():
    with open('Obj_Fn_Data', 'w+') as obj_file:
        obj_file.write('1.0 2.0 3.0\n')
        obj_file.write('1.2 1.8 4.0')
    with gaopaw.fileutils.chdir(working_dir('Si')):
        assert round(gaopaw.calc_obj_fn([1.2, 1.8, 4.0]), 3) == 0.667
    gaopaw.os.remove('Obj_Fn_Data')
    with open('Obj_Fn_Data', 'w+') as obj_file:
        obj_file.write('1.0 2.0 3.0\n')
    with gaopaw.fileutils.chdir(working_dir('Si')):
        assert round(gaopaw.calc_obj_fn([1.2, 1.8, 4.0]), 3) == 0.0
    with open('Obj_Fn_Data', 'w+') as obj_file:
        obj_file.write('1.0 2.0 0.0\n')
        obj_file.write('1.2 1.8 0.0')
    with gaopaw.fileutils.chdir(working_dir('Si')):
        assert round(gaopaw.calc_obj_fn([1.2, 1.8, 4.0]), 3) == 0.333
    gaopaw.os.remove('Obj_Fn_Data')
