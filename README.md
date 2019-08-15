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

(3) If necessary, create a template directory (i.e., "cmpd_template_dir" specified in input.json) containing all QE input files required to calculate properties of compounds specified in the input.json file. Filenames follow the format of *compound.formula.lattice_type.calculation type.template*. Note that elemental properties (log derivatives and FCC/BCC lattice constants) are considered automatically and therefore input files (contained in "elem_template_dir") for these runs need not be explicitly provided. 

(4) Create a template dakota.in file (see Examples/Dakota_template) in the working directory. Modify the string for analysis_driver to include the absolute path to analysis.py. For all other settings, use scripts/write_dakota.py to parse the input.json file and write the following information to dakota.in:
- Elemental variable bounds (taken from "elem_template_dir"/BOUNDS)
- Total number of objective functions
- Objective function labels (For details on the dakota parameters, see https://dakota.sandia.gov/documentation.html)

(5) Execute Dakota optimization through *dakota dakota.in*; generally submit using job script (see Examples/job_template). Note that the evaluation_concurrency variable, which controls the number of parallel jobs allowed to run at once, should be chosen according to the number of processors available. Currently the code is set up to employ 4 processors per job; hence the evaluation_concurrency should equal the total number of processors available divided by 4.

(6) Once the optimization is complete, use scripts/get_best_soln.py to to retrieve the best solution obtained throughout all generations considered. The weighted sum approach is currently utilized to normalize all objective functions, for which the mean absolute error (MAE) is calculated. By default, all objectives are weighted equally. Alternatively, user specified weights may be used; to do so, add --weights to the end of the python execution statement. In comparing all result sets, the universal minimum in MAE is chosen as the best solution; corresponding objective functions and atompaw input files are placed in the Best_Solution/ folder.

(7) Once a set of optimized pseudopotentials is obtained, a final test may be carried out using scripts/test_PP.py. The procedure follows same as usual according to the input.json file, however, no dakota settings are necessary. Additionally, it is recommended that a new directory containing the .UPF files be created, for which the path is specified as "paw_dir" in the input.json (under "directories").

#### Optional arguments:

The following may be specified under the "directories" section of the input.json file:

- "optimize_log_grid": True or False
    - Whether to include number of logarithmic grid points (see AtomPAW docs) in optimization.

- "include_paw": [list of elements]
    - Which elements would you like to pull from a directory ("paw_dir" if specified, otherwise "cmpd_template_dir") and keep fixed throughout optimization
    
#### Currently supported properties with required formats:

- Logarithmic derivatives of the pseudized wavefunctions
    - Need not be specified explicitly, automatically considered for each element
    
- Lattice parameter(s)
    - For cubic symmetry: "lattice_constant": *a* in angstroms
    - For all lower-symmetry structures: "lattice_constant": [unique lattice constants (angstroms) in order of increasing magnitude, unique non-90 lattice angles (degrees) in order of increasing magntitude]
    
- Delta-factor (based on equation of state, see: https://molmod.ugent.be/deltacodesdft)
    - "eos": [volume (cubic angstroms), bulk modulus (GPa), derivative of bulk modulus with respect to pressure]
    
- Bulk modulus
    - "bulk_modulus": *B* in GPa
    
- Phonon frequencies
    - "phonon_frequency": [individual frequencies in order of increasing value (including sign)]
    
- Atomic positions
    - "atomic_positions": True (data will be taken from a file, in "cmpd_template_dir", with positions written in QE format, specifically in cartesian coordinates, units of angstroms)
    
- Electronic band gap
    - "band_gap": gap in eV
    
- Magnetization
    - "magnetization": net magnetization in Bohr magnetons
    
- Individual magnetic moments
    - "magnetic_moment": [moments (in Bohr magnetons) in order consistent with corresponding atoms in QE input]
    
- Polymorph energetics
    - "lattice_type": [lattice type labels for each polymorph, listed in order of increasing energy]

#### Currently supported lattice type labels:

Unique labels, which may share the same Bravais lattice, are useful to distringuish between polymorphs. While the following are currently implemented, more may be added by simply appending to the lists present in the getLatticeConstant method.

- Cubic: "FCC", "ZB", "RS", "diamond", "HH", "BCC", "per", "SC", "CsCl"

- Hexagonal: "hex", "WZ"

- Rhombohedral: "rhomb"

- Orthorhombic: "ortho"

- Tetragonal: "tetrag"

- Monoclinic: "monoclin"

- Triclinic: "triclin"

#### Some caveats:

- Input structures may be in the conventional or primitive setting, however, if the latter is used, the cell orientation must be consistent with the default QE orientation according to ibrav (see: https://www.quantum-espresso.org/Doc/INPUT_PW.html). To avoid such issues, it is best to explicitly set ibrav (> 0) as opposed to specifying your own cell (with ibrav = 0). This also ensures your lattice parameters are nearly exact (e.g., 90 or 120), which serves to avoid any errors in the getLatticeConstant method (checks are made to ensure correct cell shape, with some tolerance).

- Roughly optimized on AtomPAW varialbes have been obtained for the majority of elements (see Elem_Templates/BOUNDS) typically by considering logarithmic derivatives and FCC/BCC lattice constants. The write_dakota.py script will generate bounds according to these values; 5% above and below (default) seems to work well, however, this range may be changed in the updateVars method.

- Generally, FCC/BCC lattice constants are automatically tested for each element throughout an optimization. However, exceptions include:
    - N: dimer separation is tested by considering atomic positions
    - P: lattice constants are tested with respect to the orthorhombic ground state
    - Hg: lattice constants are tested with respect to the tetragonal ground state
    - f-block (where available): lattice constant and magnetization tested for rocksalt nitrides. Note that variable bounds have not yet been obtained for these elements; once fully-optimized nitrogen PAW has been obtained, these optimizations may be performed with fixed N potential (see Examples/La/).
    
- Further properties (e.g., phonon frequencies) may be tested with respect to elemental FCC/BCC states, however, lattice constant should (typically) be ommitted from the properties, considering this is already tested for most elements.

- GAOPAW will determine what calculations are necessary for a given property. By default, a relaxation is always carried as a first step, followed by a scf calculation if appropriate (e.g., to obtain the band gap). However, if you wish to structure fixed (perhaps you'd like to match with experimental or AE values), this can be acheived by modifying the QE relaxation input accordingly: 'none' for cell_dynamics and/or 0.0 scaling constants for ionic force factors.

- Similarily, in some cases it may be useful to keep the symmetry of the cell fixed (note that if symmetry changes to a different Bravais lattice, GAOPAW will give an error), in which case the cell_dofree tag (under &CELL) may be set to 'ibrav' in the corresponding QE input.
