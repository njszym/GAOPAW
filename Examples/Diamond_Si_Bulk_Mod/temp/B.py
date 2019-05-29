import pandas as pd
import numpy as np
from scipy.optimize import curve_fit as cf


def birch_murnaghan(V, V0, B0, B0_prime, E0):
    """
    3rd order Birch-Murnaghan equation of state, in the energy-volume form
    """
    V = np.array(V)
    return E0 + 9 * V0 * B0 / 16. * (
        ((V0 / V) ** (2 / 3.) - 1) ** 3 * B0_prime +
        ((V0 / V) ** (2 / 3.) - 1) ** 2 * (6 - 4 * (V0 / V) ** (2 / 3.)))

def get_bulk(num_atoms):
    """
    Reads in energy-volume data from E_V.txt and calculates bulk modulus
    """
    df = pd.read_table('E_V.txt',sep='\s+',header=None)
    E = list(df[0])
    V = list(df[1])
    V = np.array([0.14818453429566825*value for value in V]) ## bohr^3 to A^3
    Y = np.array([13.6056980659*value for value in E]) ## Ry to eV
    initial_parameters = [V.mean(), 2.5, 4, Y.mean()]
    fit_eqn = eval('birch_murnaghan')
    popt, pcov = cf(fit_eqn, V, Y, initial_parameters)
    volume = popt[0]/num_atoms
    bulk = popt[1]*160.2
    B_prime = popt[2]
    return float(volume), float(bulk), float(B_prime)

