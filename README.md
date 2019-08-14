# Genetic Algorithm Optimizer for PAW potentials

Workflow designed to utilize AtomPAW, Quantum Espresso (QE), and Dakota to generate optimized PAW potentials through implementation of a multi-objective genetic algorithm.

Objective functions, defined as differences between calculated data (using pseudopotentials in QE) and known benchmark data (generally from all-electron calculations performed with WIEN2k, but may also be from well-converged pseudopotential calculations), are minimized for a given element (as an elemental state and/or as a constituent element of a compositie system). Several elements or compounds of any periodic, ordered structure may be tested in a single run.

### General usage notes:

(1) Create a working directory under which all calculations will take place.

(2) Write the input.json file containing information on the compounds and properties to be tested. Several sample input files are provided in the Examples/ folder. The general structure is as follows:
- "directories":
  - "elem_template_dir": "*absolute path to elemental template files*"
  - "cmpd_template_dir": "*absolute path to compound templates files*"
- "compounds":
  - For each compound, provide a formula, lattice type, and any properties to be tested/optimized.
  - For each property, provide the corresponding all-electron data.

(3) If necessary, create a template directory (i.e., "cmpd_template_dir" specified in input.json) containing all QE input files required to calculate properties of compounds specified in the input.json file. Filenames follow the format of *compound*.*formula*.*lattice type*.*calculation type*.*template*. Note that elemental properties (log derivatives and FCC/BCC lattice constants) are considered automatically and therefore input files (contained in "elem_template_dir") for these runs need not be explicitly provided. 

(4) Create a template dakota.in file (see Examples/Dakota_template) in the working directory. Use scripts/write_dakota.py to parse the input.json file and write the following information to dakota.in:
- Elemental variable bounds (taken from "elem_template_dir"/BOUNDS)
- Total number of objective functions
- Objective function labels
For details on the dakota parameters, see https://dakota.sandia.gov/documentation.html.

(5) Execute Dakota optimization through *dakota dakota.in*; generally submit using job script (see Examples/job_template). Note that the evaluation_concurrency variable, which controls the number of parallel jobs allowed to run at once, should be chosen according to the number of processors available. Currently the code is set up to employ 4 processors per job; hence the evaluation_concurrency should equal the total number of processors available divided by 4.

(6) Once the optimization is complete, use scripts/get_best_soln.py to to retrieve the best solution obtained throughout all generations considered. The weighted sum approach is currently utilized to normalize all objective functions, for which the mean absolute error (MAE) is calculated. By default, all objectives are weighted equally. Alternatively, user specified weights may be used; to do so, add --weights to the end of the python execution statement. In comparing all result sets, the universal minimum in MAE is chosen as the best solution; corresponding objective functions and atompaw input files are placed in the Best_Solution/ folder.

(7) Once a set of optimized pseudopotentials is obtained, a final test may be carried out using scripts/test_PP.py. The procedure follows same as usual according to the input.json file, however, no dakota settings are necessary. Additionally, it is recommended that a new directory containing the .UPF files be created, for which the path is specified as "paw_dir" in the input.json (under "directories").

### Optional arguments:

The following may be specified under the "directories" section of the input.json file:

- "optimize_log_grid": True or False
    - Whether to include number of logarithmic grid points (see AtomPAW docs) in optimization.

- "include_paw": [list of elements]
    - Which elements would you like to pull from a directory ("paw_dir" if specified, otherwise "cmpd_template_dir") and keep fixed throughout optimization
    
### Currently supported properties with required formats

- Logarithmic derivatives of the pseudized wavefunctions
    - Need not be specified explicitly, automatically considered for each element
    
- Lattice parameter(s)
    - For cubic symmetry: "lattice_constant": *a* in angstroms
    - For all lower-symmetry structures: "lattice_constant": [unique lattice constants (angstroms) in order of increasing magnitude, unique non-90 lattice angles (degrees) in order of increasing magntitude]
    
- Delta-factor (based on equation of state)
    - "eos": [volume (cubic angstroms), bulk modulus (GPa), derivative of bulk modulus with respect to pressure]
    
- Bulk modulus
    - "bulk_modulus": *B* in GPa
    
- Phonon frequencies
    - "phonon_frequency": [individual frequencies in order of increasing value (including sign)]
    
- Atomic positions
    - "atomic_positions": True (data will be taken from file (in "cmpd_template_dir") with positions written in QE format)
    
- Electronic band gap
    - "band_gap": gap in eV
    
- Magnetization
    - "magnetization": net magnetization in Bohr magnetons
    
- Individual magnetic moments
    - "magnetic_moment": [moments (in Bohr magnetons) in order consistent with corresponding atoms in QE input]
