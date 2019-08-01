from gaopaw import *


def main():
    """
    Parse information from input.json and update dakota.in accordingly.
    """
    gp_run = Runner(input_dir='current')
    gp_run.updateNumObjs()
    gp_run.updateVars()
    gp_run.updateLabels()

if __name__=='__main__':
    main()

