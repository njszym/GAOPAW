import os
from schrodinger.utils import fileutils
import numpy as np
import glob
import shutil
import math


## DEPRACATED

def scale_factor(ibrav,scale_vol,celldm):
    """
    Convert scale_vol into a value suitable to
    mutliply celldm(1) such the the resulting cell
    is scaled by the desired amount. Depends in the
    bravais lattice, hence ibrav is needed.
    """
    if ibrav == 1: ## Cubic (P)
        return scale_vol**(1./3.)
    if ibrav == 2: ## Cubic (F)
        return (4*scale_vol)**(1./3.)
    if ibrav == 3: ## Cubic (I)
        return (2*scale_vol)**(1./3.)
    if ibrav == 4: ## Hexagonal (P)
        return (scale_vol*((2./math.sqrt(3))/celldm[3]))**(1./3.)
    if ibrav == 5: ## Rhombohedral (R)
        return (scale_vol / (1-3*celldm[4]**2+2*celldm[4]**3)**(1./2.))**(1./3.)
    if ibrav == 6: ## Tetragonal (P)
        return (scale_vol / celldm[3])**(1./3.)
    if ibrav == 7: ## Tetragonal (I)
        return (2*scale_vol / celldm[3])**(1./3.)
    if ibrav == 8: ## Orthorhombic (P)
        return (scale_vol / (celldm[2]*celldm[3]))**(1./3.)
    if ibrav in [9, 11]: ## Orthorhombic (C, I)
        return (2*scale_vol / (celldm[2]*celldm[3]))**(1./3.)
    if ibrav == 10: ## Orthorhombic (P)
        return (4*scale_vol / (celldm[2]*celldm[3]))**(1./3.)
    if ibrav == 12: ## Monoclinic (P)
        return  (scale_vol / (celldm[2]*celldm[3]*math.sqrt(1-celldm[4]**2)))**(1./3.)
    if ibrav == 13: ## Monoclinic (C)
        return (2*scale_vol / (celldm[2]*celldm[3]*math.sqrt(1-celldm[4]**2)))**(1./3.)
    if ibrav == 14: ## Triclinic (P)
        return (scale_vol / ( celldm[2]*celldm[3]*math.sqrt(
            1 - celldm[4]**2 - celldm[6]**2 + 2*celldm[4]*celldm[5]*celldm[6])))

