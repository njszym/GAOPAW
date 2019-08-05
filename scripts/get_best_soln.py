import os
import sys
sys.path.append(os.path.join(os.getcwd(), os.pardir))
from gaopaw import *


def main():
    """
    Get best solution from dakota_tabular.dat
    """
    gp_run = Runner(input_dir='current')
    gp_run.getBestSoln()

if __name__=='__main__':
    main()
