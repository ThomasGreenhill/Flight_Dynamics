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


def pyAvl_Cf_Cm(acftpath, x, controlvec, Ixx, Iyy, Izz, Ixz, m, X_cg, rho, bspan, cbar):
    # Unpack the state variables
    Vt = x[0]
    a = x[1]
    b = x[2]
    Phi = x[3]
    The = x[4]
    Psi = x[5]
    p = x[6]
    q = x[7]
    r = x[8]
    xE = x[9]
    yE = x[10]
    h = x[11]

    # Unpack the control inputs
    de = controlvec[0]
    da = controlvec[1]
    dr = controlvec[2]
    df = 0
    dw = 0

    aircraft = avl.FileWrapper(acftpath)

    # Show the aicraft geometry
    # session = avl.Session(geometry=aircraft)
    # if 'gs_bin' in session.config.constraints:
    #     img = session.save_geometry_plot()[0]
    #     avl.show_image(img)
    # else:
    #     session.show_geometry()

    # Set the control deflections
    da_param = avl.Parameter(name='aileron', constraint='d1', value=da)
    df_param = avl.Parameter(name='flaps', constraint='d2', value=df)
    dw_param = avl.Parameter(name='winglet', constraint='d3', value=dw)
    de_param = avl.Parameter(name='elevator', constraint='d4', value=de)
    dr_param = avl.Parameter(name='rudder', constraint='d5', value=dr)

    # Input the current states

    alpha_const = avl.State(name='alpha', value=a*r2d, unit='deg')

    alpha_param = avl.Parameter(name='alpha', constraint='alpha', value=a*r2d)
    beta_param = avl.Parameter(name='beta', constraint='beta', value=b*r2d)
    # print('p{}'.format(p*bspan/(2*Vt)))
    p_param = avl.Parameter(name='p', constraint='pb/2V', value=p*bspan/(2*Vt))

    # print('q{}'.format(q*cbar/(2*Vt)))
    q_param = avl.Parameter(name='q', constraint='qc/2V', value=q*cbar/(2*Vt))

    # print('r{}'.format(r*bspan/(2*Vt)))
    r_param = avl.Parameter(name='r', constraint='rb/2V', value=r*bspan/(2*Vt))

    current_params = avl.Case(name='Get_Current_Coefficients',
                                   alpha=alpha_param, beta=beta_param,
                                   roll_rate=p_param, pitch_rate=q_param, yaw_rate=r_param,
                                   aileron=da_param, flaps=df_param, winglet=dw_param, elevator=de_param, rudder=dr_param)

    current_states = avl.Case(name='Get_Current_Coefficients', alpha=alpha_const)

    # create session with the geometry object and the cases
    session = avl.Session(geometry=aircraft, cases=[current_states, current_params])

    # get results and write the resulting dict to a JSON-file
    # session.show_geometry()
    results = session.get_results()
    with open('out.json', 'w') as f:
        f.write(json.dumps(results))

    # print(result)
    print("CX = {}".format(
        result['Get_Current_Coefficients']['Totals']['CXtot']))
    print("CY = {}".format(
        result['Get_Current_Coefficients']['Totals']['CYtot']))
    print("CZ = {}".format(
        result['Get_Current_Coefficients']['Totals']['CZtot']))
    print("Cm = {}".format(
        result['Get_Current_Coefficients']['Totals']['Cmtot']))
    print("Cn = {}".format(
        result['Get_Current_Coefficients']['Totals']['Cntot']))
    print("Cl = {}".format(
        result['Get_Current_Coefficients']['Totals']['Cltot']))
    print("x = {}".format(x))
    print("u = {}".format(controlvec))

    session.export_run_files()
