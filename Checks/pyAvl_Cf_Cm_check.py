# Check file for pyAVL_Cf_Cm


import numpy as np
import sys
import avlwrapper as avl
sys.path.append('../AVL_Automation')
from pyAvl_Cf_Cm import pyAvl_Cf_Cm


# x = np.random.random((12, 1))
x = np.zeros((12,1))
x[0] = 42
x[1] = 4*np.pi/180
x[2] = 3*np.pi/180

print(x)

# Path to aircraft geometry file for AVL
# acftpath = '/Users/thomasgreenhill/SUAVE/ASW22BL/ASW_22_Demo/Geometry/ASW22BL_Wing_With_Tips_Tail_Surfaces_v1.4.avl'
acftpath = '/Users/thomasgreenhill/Documents/AVL/C182RG/C182RG.avl'

global rho, T, heng
rho = 0.002377

controlvec = np.array([
    [10],
    [0],
    [0]
])

global Ixx, Iyy, Izz, Ixz, cbar, S, b, g, m, X_cg

def acft_props():
    # Properties for F-4 Phantom (all in IMP units)
    Ixx = 2.5e4
    Iyy = 1.22e5
    Izz = 1.4e5
    Ixz = 9.75e2

    cbar = 16
    S = 530
    bspan = 38.67
    g = 32.2
    m = 38924/32.2
    X_cg = 0.289

    return Ixx, Iyy, Izz, Ixz, cbar, S, bspan, g, m, X_cg

Ixx, Iyy, Izz, Ixz, cbar, S, bspan, g, m, X_cg = acft_props()


CX, CY, CZ, Cm, Cn, Cl = pyAvl_Cf_Cm(acftpath, x, controlvec, Ixx, Iyy, Izz, Ixz, m, X_cg, rho, g, bspan, cbar, runfileexport=False)

print(CX)
