import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
from gaopaw import *


def main():
    """
    Parse information from input.json and update dakota.in accordingly.
    """
    gp_run = Runner(input_dir='current', writing_dakota=True)
    gp_run.updateNumObjs()
    gp_run.updateVars()
    gp_run.updateLabels()

if __name__=='__main__':
    main()

