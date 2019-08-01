from gaopaw import *


def main():
    """
    Get best solution from dakota_tabular.dat
    """
    gp_run = Runner(input_dir='current')
    Runner.getBestSoln()

if __name__=='__main__':
    main()
