# Genetic Algorithm Optimizer for PAW datasets

Workflow designed to utilize AtomPAW, Quantum Espresso (QE), and Dakota to generate optimized PAW potentials through implementation of genetic algorithms.

Objective functions, defined as differences between calculated (with QE) and known all-electron (with WIEN2k) data, are minimized for a given element (as an elemental state and/or as a constituent element of a compositie system).

Several elements or compounds of any periodic, ordered structure may be tested.

Currently supported properties include lattice constants, band gaps, Delta-factor (see https://molmod.ugent.be/deltacodesdft), bulk modulus, atomic positions, total magnetization, individual magnetic moments, and phonon frequencies.

### General usage notes:

(i) Create a working directory under which all calculations will take place.

(ii) Write the input.json file containing information on the compounds and properties to be tested. See the Examples/ folder for information regarding parameter names and formatting guidelines.

(iii) Create a template directory (path specified in input.json) containing all AtomPAW and QE input files for the compounds specified in the input.json file. Filenames follow the general format of (compound formula).(lattice type).(calculation type).template. Note that elemental properties (log derivs, FCC/BCC lattice constants) are considered automatically and therefore input files for these runs need not be explicitly provided.

(iv) Create a template dakota.in file (may take from Examples/ folder) in current directory. Use the write_dakota.py code from the scripts/ folder to parse input.json file and write information regarding elemental variable bounds (optimal values for each element taken from Elem_Templates folder; path should be specified in the input.json file) and total number of objective functions, along with corresponding labels, to the dakota.in file. For details on the dakota parameters, see https://dakota.sandia.gov/documentation.html.

(v) Execute Dakota optimization through "dakota dakota.in" (or submit using job script).

(vi) Once the optimization is complete, use the get_best_soln.py to retrieve the best solution (must be executed in same folder as the dakota.in and input.json files). The weighted sum approach is currently utilized to normalize all objective functions. The universal minimum in the mean absolute error of these noramlize objective functions is then chosen as the best solution, for which detailed results and atompaw input files are placed in the Best_Solution/ folder of the current working directory. 

(vii) Upon obtaining a set of optimized pseudopotentials, a final test may be carried out using the test_PP.py script. The procedure follows same as usual according to the input.json file, however, no dakota settings are necessary. Additionally, a new directory containing the .UPF files must be created and should be specified in the header of the test_PP.py code.
