# Genetic Algorithm Optimizer for PAW datasets

Workflow designed to utilize AtomPAW, Quantum Espresso (QE), and Dakota to generate optimized PAW potentials through implementation of genetic algorithms.

Objective functions, defined as differences between calculated (with QE) and known all-electron (with WIEN2k) data, are minimized for a given element.

Elements or compounds of any periodic, ordered structure maybe tested.

Currently supported properties include lattice constants, band gaps, Delta-factor (see https://molmod.ugent.be/deltacodesdft), bulk modulus, atomic positions, total magnetization, individual magnetic moments, and phonon frequencies.

### General usage notes:

(i) Create a working directory under which all calculations will take place.

(ii) Write a gaopaw.yaml file containing information on the elements, compounds, and properties to be tested.

(iii) Create a template directory (path specified in gaopaw.yaml) containing all AtomPAW and QE input files.

(iv) Write a dakota.in file containing all control parameters for the genetic algorithm as implemented in Dakota. An example dakota.in file with (roughly) optimized parameters may be found in the "Dakota_Template" folder. For details, see https://dakota.sandia.gov/documentation.html.

(v) Execute the Dakota through the following command: "dakota dakota.in".

(vi) As the calculations take place, the best result will be continuously updated and placed in the "Best_Result" folder along with the relevant input parameters.
