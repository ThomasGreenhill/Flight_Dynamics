# Check file for pyAVL_Cf_Cm


import numpy as np
import sys
import avlwrapper as avl
sys.path.append('../AVL_Automation')
from pyAvl_Cf_Cm_FAER import pyAvl_Cf_Cm_FAER


# x = np.random.random((12, 1))
x = np.zeros((12,1))
x[0] = 42
x[1] = 4*np.pi/180
x[2] = 3*np.pi/180


# Path to aircraft geometry file for AVL
# acftpath = '/Users/thomasgreenhill/SUAVE/ASW22BL/ASW_22_Demo/Geometry/ASW22BL_Wing_With_Tips_Tail_Surfaces_v1.4.avl'
acftpath = '/Users/thomasgreenhill/Documents/AVL/C182RG/C182RG_SI.avl'

aircraft = avl.FileWrapper(acftpath)

# Show the aicraft geometry
# session = avl.Session(geometry=aircraft)
# if 'gs_bin' in session.config.settings:
#     img = session.save_geometry_plot()[0]
#     avl.show_image(img)
# else:
#     session.show_geometry()
    


global rho, T, heng
rho = 1.225

controlvec = np.array([
    [0],
    [0],
    [0],
    [0]
])

global Ixx, Iyy, Izz, Ixz, cbar, S, b, g, m, X_cg

def acft_props():
   # Properties for Cessna 172 Skyhawk http://jsbsim.sourceforge.net/MassProps.html
    Ixx = 948*1.35581795
    Iyy = 1346*1.35581795
    Izz = 1967*1.35581795
    Ixz = 0

    cbar = 4.9*0.3048
    S = 174*0.3048*0.3048
    bspan = 35.8*0.3048
    g = 9.81
    m = 1700*0.453592
    X_cg = 41/12*0.3048

    return Ixx, Iyy, Izz, Ixz, cbar, S, bspan, g, m, X_cg

Ixx, Iyy, Izz, Ixz, cbar, S, bspan, g, m, X_cg = acft_props()


CX, CY, CZ, Cm, Cn, Cl = pyAvl_Cf_Cm_FAER(acftpath, x, controlvec, Ixx, Iyy, Izz, Ixz, m, X_cg, rho, g, bspan, cbar, runfileexport=False)

print(CX)
