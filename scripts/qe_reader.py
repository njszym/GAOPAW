from schrodinger.utils import imputils
from schrodinger.utils import fileutils
from schrodinger.application.matsci import property_names as pnames
from schrodinger.application.matsci.nano import xtal
import sys
import os

qe_reader_path = os.path.join(fileutils.get_mmshare_scripts_dir(),
    'periodic_dft_gui_dir', 'qe2mae.py')
qe_reader_mod = imputils.import_module_from_file(qe_reader_path)

qe_reader = qe_reader_mod.QEOutputReader(sys.argv[1])

print('Last structure id = ', qe_reader.final_struct_id)
struct = qe_reader.structs[qe_reader.final_struct_id]

cparams = xtal.get_chorus_properties(struct)
params = xtal.get_params_from_chorus(cparams)
print('Cell dimensions = ', params)
