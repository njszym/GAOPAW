from gaopaw import *

## Usage: python optimize_loggrid.py [elem name]

## May pass energy_tol and lat_tol args to optimizeLogGrid()
## to change toterances; defaults are 1e-6 eV and 1e-3 angstroms
## for energy and lattice constants respectively.

def main():
    """
    Iteratively decrease loggrid density until
    tolerances to changes in energy and lattice
    constants are reached.
    """
    elem = sys.argv[-1]
    gp_run = Runner(input_dir='current')
    gp_run.optimizeLogGrid(elem)

if __name__=='__main__':
    main()

