# Genetic Algorithm Optimizer for PAW datasets

Workflow designed to utilize AtomPAW, Quantum Espresso (QE), and Dakota to generate optimized PAW potentials through implementation of genetic algorithms.

Objective functions, defined as differences between calculated (with QE) and known all-electron (with WIEN2k) data, are minimized for a given element (as an elemental state and/or as a constituent element of a compositie system).

Several elements or compounds of any periodic, ordered structure may be tested.

Currently supported properties include lattice constants, band gaps, Delta-factor (see https://molmod.ugent.be/deltacodesdft), bulk modulus, atomic positions, total magnetization, individual magnetic moments, and phonon frequencies.

### General usage notes:

(i) Create a working directory under which all calculations will take place.

(ii) Write the input.json file containing information on the compounds and properties to be tested. Several sample input files are provided in the Examples/ folder. The general structure is as follows:
- "directories":
  - "elem_template_dir": [path to elemental template files]
  - "cmpd_template_dir": [path to compound templates files]
- "compounds":
  - For each compound, provide a formula, lattice type, and any properties to be tested/optimized.
  - For each property, provide the corresponding all-electron data.

(iii) If necessary, create a template directory (path specified in input.json) containing all QE input files required to calculate properties of compounds specified in the input.json file. Filenames follow the general format of [compound formula].[lattice type].[calculation type].template. Note that elemental properties (log derivs, FCC/BCC lattice constants) are considered automatically and therefore input files for these runs need not be explicitly provided.

(iv) Create a template dakota.in file (see Examples/Dakota_template) in the working directory. Use scripts/write_dakota.py to parse the input.json file and write the following information to dakota.in:
- Elemental variable bounds (optimal values taken from Elem_Templates/BOUNDS; path should be specified in input.json)
- Total number of objective functions
- Objective function labels
For details on the dakota parameters, see https://dakota.sandia.gov/documentation.html.

(v) Execute Dakota optimization through "dakota dakota.in"; generally submit using job script. Note that the evaluation_concurrency variable, which controls the number of parallel jobs allowed to run at once, should be chosen according to the number of processors available. Currently the code is set up to employ 4 processors per job; hence the evaluation_concurrency should equal the total number of processors available divided by 4.

(vi) Once the optimization is complete, use scripts/get_best_soln.py to to retrieve the best solution obtained throughout all generations considered. The weighted sum approach is currently utilized to normalize all objective functions, for which the mean absolute error (MAE) is calculated. The universal minimum in MAE is chosen as the best solution; corresponding objective functions and atompaw input files are placed in the Best_Solution/ folder.

(vii) Upon obtaining a set of optimized pseudopotentials, a final test may be carried out using scripts/test_PP.py. The procedure follows same as usual according to the input.json file, however, no dakota settings are necessary. Additionally, it is recommended that a new directory containing the .UPF files be created, for which the path is specified as "paw_dir" is the input.json (under "directories").
