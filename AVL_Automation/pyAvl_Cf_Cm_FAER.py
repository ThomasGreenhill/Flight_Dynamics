# Calls AVL to produce the force and moment coefficients about the body axis frame, given the current states 
#
#   Function calls:
#       {none}
#
#   References:
#       1. Garza, Morelli: *A collection of Nonlinear Aircraft Simulations in MATLAB*
#
#   Notes:
#       1. Current control inputs are mapped for those of B737
#
#   History:
#       1. Created and debugged, TVG 01.08.2020

import sys
sys.path.append('../avlwrapper')
import avlwrapper as avl
import numpy as np

global r2d
r2d = 180/np.pi


def pyAvl_Cf_Cm_FAER(acftpath, x, controlvec, Ixx, Iyy, Izz, Ixz, m, X_cg, rho, g, bspan, cbar, runfileexport):
    # print(x)
    # print(controlvec)

    # Path to aircraft geometry
    # NOTE: Airfoil file references in the aircraft geometry .avl file must be updated according to your computer's local file system
    aircraft = avl.FileWrapper(acftpath)

    # Show the aicraft geometry
    # session = avl.Session(geometry=aircraft)
    # if 'gs_bin' in session.config.settings:
    #     img = session.save_geometry_plot()[0]
    #     avl.show_image(img)
    # else:
    #     session.show_geometry()

    # Unpack the state variables
    Vt = float(x[0])
    a = float(x[1])
    b = float(x[2])
    Phi = float(x[3])
    The = float(x[4])
    Psi = float(x[5])
    p = float(x[6])
    q = float(x[7])
    r = float(x[8])
    xE = float(x[9])
    yE = float(x[10])
    h = float(x[11])

    # Unpack the control inputs
    df = float(controlvec[0])
    da = float(controlvec[1])
    de = float(controlvec[2])
    dr = float(controlvec[3])



    current_params = avl.Case(name='Get_Current_Coefficients',
                                   alpha=a*r2d, beta=b*r2d,
                                   roll_rate=p*bspan/(2*Vt), pitch_rate=q*cbar/(2*Vt), yaw_rate=r*bspan/(2*Vt),
                                   velocity=Vt, density=rho, gravity=g, 
                                   mass=m, Ixx=Ixx, Iyy=Iyy, Izz=Izz, Izx=Ixz,
                                   bank=Phi, elevation=The, heading=Psi,
                                   flap=df, aileron=da, elevator=de, rudder=dr
                                   )

    

    session = avl.Session(geometry=aircraft, cases=[current_params])

    
    result = session.run_all_cases()

    # UNCOMMENT to print out the force and moment coefficients
    # print("CX = {}".format(
    #     result['Get_Current_Coefficients']['Totals']['CXtot']))
    # print("CY = {}".format(
    #     result['Get_Current_Coefficients']['Totals']['CYtot']))
    # print("CZ = {}".format(
    #     result['Get_Current_Coefficients']['Totals']['CZtot']))
    # print("Cm = {}".format(
    #     result['Get_Current_Coefficients']['Totals']['Cmtot']))
    # print("Cn = {}".format(
    #     result['Get_Current_Coefficients']['Totals']['Cntot']))
    # print("Cl = {}".format(
    #     result['Get_Current_Coefficients']['Totals']['Cltot']))

    ## Export the run file (optional)
    if runfileexport == True:
        session.export_run_files()

    CX = result['Get_Current_Coefficients']['Totals']['CXtot']
    CY = result['Get_Current_Coefficients']['Totals']['CYtot']
    CZ = result['Get_Current_Coefficients']['Totals']['CZtot']
    Cm = result['Get_Current_Coefficients']['Totals']['Cmtot']
    Cn = result['Get_Current_Coefficients']['Totals']['Cntot']
    Cl = result['Get_Current_Coefficients']['Totals']['Cltot']

    return CX, CY, CZ, Cm, Cn, Cl
