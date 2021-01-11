# Simulation code for aircraft with linear kinematics (constant air properties)
#
#   Function calls:
#       ../Utilities/rk4.py (Runge-Kutta 4th Order Solver)
#       ../AVL_Automation/pyAvl_Cf_Cm.py (Calls AVL to come up with force and moment coefficients at current time step)
#
#   References:
#       1. Garza, Morelli: *A collection of Nonlinear Aircraft Simulations in MATLAB* 
#
#
#   TVG 01.08.2020
#
#   

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../Utilities')
from rk4 import rk4
sys.path.append('../AVL_Automation')
from pyAvl_Cf_Cm import pyAvl_Cf_Cm
sys.path.append("../Figure_Formatting")
import formatfigures

# Figure formatting preamble
returnval = formatfigures.formatfigures()
colors = returnval[0]
params = returnval[1]

# Path to aircraft geometry file for AVL
acftpath = '/Users/thomasgreenhill/SUAVE/ASW22BL/ASW_22_Demo/Geometry/ASW22BL_Wing_With_Tips_Tail_Surfaces_v1.4.avl'

global rho, T, heng
rho = 1.225
# Neglect thrust and engine angular momentum for now 
T = 0
heng = 0

global Ixx, Iyy, Izz, Ixz, cbar, S, bspan, g, m, X_cg

def acft_props():
    # Bogus properties for ASW-22 BL
    Ixx = 100000
    Iyy = 200000
    Izz = 300000
    Ixz = 0.1

    cbar = 0.659
    S = 17.36575
    bspan = 29.35
    g = 9.81
    m = 850
    X_cg = 0.3

    return Ixx, Iyy, Izz, Ixz, cbar, S, bspan, g, m, X_cg

Ixx, Iyy, Izz, Ixz, cbar, S, bspan, g, m, X_cg = acft_props()

def control(t):
    '''Control input depending on current time
    
    Control input vector is organized as follows:
        controlvec = np.array([
            [de],
            [da],
            [dr]
        ])
    
    Note: 
        AVL expects units in degrees
    '''
    tstart = 5
    tend = 20

    if t > tstart and t < tend:
        controlvec = np.array([
            [10],
            [0],
            [0]
        ])
    else:
        controlvec = np.array([
            [0],
            [0],
            [0]
        ])
    return controlvec

def nonlin_eom(t, x):
    '''
    The states are organized as follows:
        x = np.array([
            [Vt],
            [a],
            [b],
            [Phi],
            [The],
            [Psi],
            [p],
            [q],
            [r],
            [xE],
            [yE],
            [h]
        ])
    '''
    
    # Unpack the state variables
    # print(x)
    # print('Next \n')
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

    qbar = 0.5*rho*Vt**2

    controlvec = control(t)

    # Inertial constants
    c1 = ((Iyy-Izz)*Izz-Ixz**2)/(Ixx*Izz-Ixz**2)
    c2 = ((Ixx-Iyy+Izz)*Ixz)/(Ixx*Izz-Ixz**2)
    c3 = Izz/(Ixx*Izz-Ixz**2)
    c4 = Ixz/(Ixx*Izz-Ixz**2)
    c5 = (Izz-Ixx)/Iyy
    c6 = Ixz/Iyy
    c7 = 1/Iyy
    c8 = ((Ixx-Iyy)*Ixx-Ixz**2)/(Ixx*Izz-Ixz**2)
    c9 = Ixx/(Ixx*Izz-Ixz**2)

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

    # Translational pseudo-states
    u = Vt*np.cos(a)*np.cos(b)
    v = Vt*np.sin(b)
    w = Vt*np.sin(a)*np.cos(b)

    # Calculate the body-axis force coefficients with AVL
    Cx, Cy, Cz, Cm, Cn, Cl = pyAvl_Cf_Cm(acftpath, x, controlvec, Ixx, Iyy, Izz, Ixz, m, X_cg, rho, g, bspan, cbar, runfileexport=False)

    # Compute the state derivatives

    udot = r*v - q*w - g*np.sin(The) + (qbar*S*Cx+T)/m
    vdot = p*w - r*u + g*np.cos(The)*np.sin(Phi) + (qbar*S*Cy)/m
    wdot = q*u - p*v + g*np.cos(The)*np.cos(Phi) + (qbar*S*Cz)/m

    Vtdot = (u*udot + v*vdot + w*wdot)/Vt
    adot = (u*wdot - w*udot)/(u**2 + w**2)
    bdot = (Vt*vdot - v*Vtdot)/(Vt**2*np.sqrt(1-(v/Vt)**2))
    
    pdot = (c1*r + c2*p + c4*heng)*q + qbar*S*bspan*(c3*Cl + c4*Cn)
    qdot = (c5*p - c7*heng)*r - c6*(p**2 - r**2) + qbar*S*cbar*c7*Cm
    rdot = (c8*p - c2*r + c9*heng)*q + qbar*S*bspan*(c4*Cl + c9*Cn)

    Phidot = p + np.tan(The)*(q*np.sin(Phi) + r*np.cos(Phi))
    Thedot = q*np.cos(Phi) - r*np.sin(Phi)
    Psidot = (q*np.sin(Phi) + r*np.cos(Phi))/np.cos(Phi)

    xEdot = u*np.cos(Phi)*np.cos(Phi) + v*(np.cos(Psi)*np.sin(The)*np.sin(Phi) - np.sin(Psi)*np.cos(Phi)) + w*(np.cos(Psi)*np.sin(The)*np.cos(Phi) + np.sin(Psi)*np.sin(Phi))
    yEdot = u*np.sin(Psi)*np.cos(The) + v*(np.sin(Psi)*np.sin(The)*np.sin(Phi) + np.cos(Psi)*np.cos(Phi)) + w*(np.sin(Psi)*np.sin(The)*np.cos(Phi) - np.cos(Psi)*np.sin(Phi))
    hdot = u*np.sin(The) - v*np.cos(The)*np.sin(Phi) - w*np.cos(The)*np.cos(Phi)

    xdot = np.array([
            Vtdot,
            adot,
            bdot,
            Phidot,
            Thedot,
            Psidot,
            pdot,
            qdot,
            rdot,
            xEdot,
            yEdot,
            hdot
        ])

    return xdot


# Initial conditions and time vector
tspan = np.linspace(0,30,31)
x0 = np.zeros((12,1))
x0[0] = 42
x0[1] = 4*np.pi/180

# Solve using RK4
x = rk4(nonlin_eom, tspan, x0, print_timestamp=True)

# Plotting
labels = ["Total Velocity", "alpha", "beta", "Phi", "Theta", "Psi", "Pitch Rate", "Roll Rate", "Yaw Rate", "x-Earth Axis Location", "y-Earth Axis Location", "Altitude Above Earth"]

plt.figure(figsize=(12, 12))
plt.plot(tspan[:],np.transpose(x[:,:]), label=labels)
plt.legend(loc='best')
plt.show()
