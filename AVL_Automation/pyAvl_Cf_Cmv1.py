# Calls AVL to produce the
#
#   Function calls:
#       ../Utilities/rk4.py (Runge-Kutta 4th Order Solver)
#
#   References:
#       1. Garza, Morelli: *A collection of Nonlinear Aircraft Simulations in MATLAB*
#
#
#   TVG 01.08.2020
#
#

import avlwrapper as avl
import numpy as np

global r2d
r2d = 180/np.pi


def pyAvl_Cf_Cm(acftpath, x, controlvec, Ixx, Iyy, Izz, Ixz, m, X_cg, rho, g, bspan, cbar, runfileexport):
    aircraft = avl.FileWrapper(acftpath)
    
    # Unpack the state variables
    # Vt = x[0]
    # print('VT', Vt)
    # a = x[1]
    # b = x[2]
    # Phi = x[3]
    # The = x[4]
    # Psi = x[5]
    # p = x[6]
    # q = x[7]
    # r = x[8]
    # xE = x[9]
    # yE = x[10]
    # h = x[11]

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
    de = float(controlvec[0])
    da = float(controlvec[1])
    dr = float(controlvec[2])
    df = float(0)
    dw = float(0)

    

    # Show the aicraft geometry
    # session = avl.Session(geometry=aircraft)
    # if 'gs_bin' in session.config.settings:
    #     img = session.save_geometry_plot()[0]
    #     avl.show_image(img)
    # else:
    #     session.show_geometry()

    # Set the control deflections
    # da_param = avl.Parameter(name='d1', setting='d1', value=da)
    # df_param = avl.Parameter(name='d2', setting='d2', value=df)
    # dw_param = avl.Parameter(name='d3', setting='d3', value=dw)
    # de_param = avl.Parameter(name='d4', setting='d4', value=de)
    # dr_param = avl.Parameter(name='d5', setting='d5', value=dr)

    # Input the current states

    # alpha_const = avl.State(name='alpha', value=a*r2d, unit='deg')

    # alpha_param = avl.Case(name='alpha', setting='alpha', value=a*r2d)

    # beta_param = avl.Parameter(name='beta', setting='beta', value=b*r2d)
    # # print('p{}'.format(p*bspan/(2*Vt)))
    # p_param = avl.Parameter(name='pb/2V', setting='pb/2V', value=p*bspan/(2*Vt))

    # # print('q{}'.format(q*cbar/(2*Vt)))
    # q_param = avl.Parameter(name='pitch_rate', setting='qc/2V', value=q*cbar/(2*Vt))

    # # print('r{}'.format(r*bspan/(2*Vt)))
    # r_param = avl.Parameter(name='yaw_rate', setting='rb/2V', value=r*bspan/(2*Vt))

    current_params = avl.Case(name='Get_Current_Coefficients',
                                   alpha=a*r2d, beta=b*r2d,
                                   roll_rate=p*bspan/(2*Vt), pitch_rate=q*cbar/(2*Vt), yaw_rate=r*bspan/(2*Vt),
                                   velocity=Vt, density=rho, gravity=g, 
                                   mass=m, Ixx=Ixx, Iyy=Iyy, Izz=Izz, Izx=Ixz,
                                   bank=Phi, elevation=The, heading=Psi,
                                   aileron=da, flaps=df, winglet=dw, elevator=de, rudder=dr
                                   )

    # current_states = avl.Case(name='Get_Current_Coefficients', alpha=alpha_const)

    session = avl.Session(geometry=aircraft, cases=[current_params])
    # print(session)
    
    result = session.run_all_cases()

    # print(result)
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
    # print("x = {}".format(x))
    # print("u = {}".format(controlvec))

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
