import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
from gaopaw import *


def main():
    """
    Get best solution from dakota_tabular.dat
    """
    gp_run = Runner(input_dir='current')
    if sys.argv[-1] == '--weights':
        gp_run.getBestSoln(use_weights=True)
    else:
        gp_run.getBestSoln()

if __name__=='__main__':
    main()
