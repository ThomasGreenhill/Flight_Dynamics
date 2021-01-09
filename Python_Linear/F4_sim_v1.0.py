# Simulation code for aircraft with linear kinematics
#
#   Function calls:
#       ../Utilities/rk4.py (Runge-Kutta 4th Order Solver)
#
#
#

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../Utilities')
from rk4 import rk4

global Xu, Xa, Zu, Za, Mu, Ma, Madot, Mq, Xde, Zde, Mde, Yb, Lb, Lp, Lr, Nb, Np, Nr, Ydr, Ldr, Ndr, Yda, Lda, Nda, g, U0, V0, W0, The0, Phi0, Psi0, Ixz, Ixx, Izz, xnb, ynb, znb

xnb = 0
ynb = 0
znb = 0

# Data from F4 Phantom at Sea Level and Ma = 0.8
Xu = -0.0162
Xa = -0.82
Zu = -0.073
Za = -1374.9
Mu = -0.0017
Ma = -17.76
Madot = -0.592
Mq = -1.36
Xde = 0
Zde = -141.0
Mde = -32.30

Yb = -299.2
Lb = -27.21
Lp = -3.034
Lr = 0.876
Nb = 16.03
Np = 0.009
Nr = -0.753
Ydr = 39.47
Ldr = 7.774
Ndr = -7.920
Yda = -6.644
Lda = 22.15
Nda = 0.556

g = 32.2

# Initial Conditions
U0 = 893  # ft/s
V0 = 0
W0 = 0

The0 = 0
Phi0 = 0
Psi0 = 0

# Inertial Properties
Ixz = 9.75*10**2
Ixx = 2.5*10*4
Izz = 1.4*10**5

# Additional unused variables
global Xda, Xb, Xdr, Xv, Xw, Xp, Xq, Xr, Yde, Ya, Yu, Yv, Yw, Yp, Yq, Yr, Zda, Zb, Zdr, Zv, Zw, Zp, Zq, Zr, Lde, La, Lu, Lv, Lw, Lq, Mda, Mb, Mdr, Mv, Mw, Mp, Mr, Nde, Na, Nu, Nv, Nw, Nq
Xda = 0
Xb = 0
Xdr = 0
Xv = 0
Xw = 0
Xp = 0
Xq = 0
Xr = 0

Yde = 0
Ya = 0
Yu = 0
Yv = 0
Yw = 0
Yp = 0
Yq = 0
Yr = 0

Zda = 0
Zb = 0
Zdr = 0
Zv = 0
Zw = 0
Zp = 0
Zq = 0
Zr = 0

Lde = 0
La = 0
Lu = 0
Lv = 0
Lw = 0
Lq = 0

Mda = 0
Mb = 0
Mdr = 0
Mv = 0
Mw = 0
Mp = 0
Mr = 0

Nde = 0
Na = 0
Nu = 0
Nv = 0
Nw = 0
Nq = 0

global Xvdot, Xwdot, Xpdot, Xqdot, Xrdot, Yudot, Ywdot, Ypdot, Yqdot, Yrdot, Zudot, Zvdot, Zpdot, Zqdot, Zrdot, Ludot, Lvdot, Lwdot, Lqdot, Lrdot, Mudot, Mvdot, Mwdot, Mpdot, Mrdot, Nudot, Nvdot, Nwdot, Npdot, Nqdot, Nrdot
Xvdot = 0
Xwdot = 0
Xpdot = 0
Xqdot = 0
Xrdot = 0

Yudot = 0
Ywdot = 0
Ypdot = 0
Yqdot = 0
Yrdot = 0

Zudot = 0
Zvdot = 0
Zpdot = 0
Zqdot = 0
Zrdot = 0

Ludot = 0
Lvdot = 0
Lwdot = 0
Lqdot = 0
Lrdot = 0

Mudot = 0
Mvdot = 0
Mwdot = 0
Mpdot = 0
Mrdot = 0

Nudot = 0
Nvdot = 0
Nwdot = 0
Npdot = 0
Nqdot = 0

# Input time series
tspan = np.linspace(0, 50, 5001)

# No perturbations initially
x0 = np.zeros((12, 1))


def kinematics(t, x):
    '''
    x is arranged with the perturbed speeds and angular rates as follows

        x =    [[x],
                [u],
                [y],
                [v],
                [z],
                [w],
                [phi],
                [p],
                [the],
                [q],
                [psi]
                [r]]
    '''

    C1 = np.sqrt(U0 ^ 2+V0 ^ 2+W0 ^ 2)

    E = np.array(
        [
            [1,    0,       0,   0,          0,        0,
             0,   0,                    0,   0,      0,   0],
            [0,    1,       0,   -Xvdot,     0,        -Xwdot,
             0,   -Xpdot,                0,   -Xqdot, 0,   -Xrdot],
            [0,    0,       1,   0,          0,        0,
             0,   0,                     0,   0,      0,   0],
            [0,    -Yudot,  0,   1-Yb/C1,    0,        -Ywdot,
             0,   -Ypdot,                0,   -Yqdot, 0,   -Yrdot],
            [0,    0,       0,   0,          1,        0,
             0,   0,                     0,   0,      0,   0],
            [0,    -Zudot,  0,   -Zvdot,     0,        1,
             0,   -Zpdot,                0,   -Zqdot, 0,   -Zrdot],
            [0,    0,       0,   0,          0,        0,
             1,   0,                     0,   0,      0,   0],
            [0,    -Ludot,  0,   -Lvdot,     0,        -Lwdot,      0,
             1,                     0,   -Lqdot, 0,   -Lrdot-Ixz/Ixx],
            [0,    0,       0,   0,          0,        0,
             0,   0,                     1,   0,      0,   0],
            [0,    -Mudot,  0,   -Mvdot,     0,        -Mwdot,
             0,   -Mpdot,                0,   1,      0,   -Mrdot],
            [0,    0,       0,   0,          0,        0,
             0,   0,                     0,   0,      1,   0],
            [0,   -Nudot,   0,   -Nvdot,     0,        -Nwdot,      0,   -Npdot-Ixz/Izz,        0,   -Nqdot, 0,   1]])

    A = np.array(
        [
            [0,   1,       0,   0,      0,   0,   0,
             0,       0,               0,         0,   0],
            [0,   Xu,      0,   Xv,     0,   Xw,  0,
             Xp,      -g*np.cos(The0),    Xq-W0,	0,   Xr+V0],
            [0,   0,       0,   1,      0,   0,   0,
             0,       0,               0,         0,   0],
            [0,   Yu,      0,   Yv,     0,   Yw,  g *
             np.cos(The0), Yp+W0,   0,               Yq,        0,   Yr-U0],
            [0,   0,       0,   0,      0,   1,   0,
             0,       0,               0,         0,   0],
            [0,   Zu,      0,   Zv,     0,   Zw,  0,
             Zp-V0,   -g*np.sin(The0),    Zq+U0,	0,   Zr],
            [0,   0,       0,   0,      0,   0,   0,
             1,       0,               0,     	0,   0],
            [0,   Lu,      0,   Lv,     0,   Lw,  0,
             Lp,      0,               Lq,        0,   Lr],
            [0,   0,       0,   0,      0,   0,   0,
             0,       0,               1,     	0,   0],
            [0,   Mu,      0,   Mv,     0,   Mw,  0,
             Mp,      0,               Mq,    	0,   Mr],
            [0,   0,       0,   0,      0,   0,   0,
             0,       0,               0,   	0,   1],
            [0,   Nu,      0,   Nv,     0,   Nw,  0,           Np,      0,               Nq,        0,   Nr]])

    B = np.array(
        [
            [0,   0,   0],
            [Xda, Xde, Xdr],
            [0,   0,   0],
            [Yda, Yde, Ydr],
            [0,   0,   0],
            [Zda, Zde, Zdr],
            [0,   0,   0],
            [Lda, Lde, Ldr],
            [0,   0,   0],
            [Mda, Mde, Mdr],
            [0,   0,   0],
            [Nda, Nde, Ndr]])

    # Nonlinear part
    a = np.arctan((x[5, [0]]+x[7, [0]]*ynb-x[9, [0]]*xnb)/U0)
    b = np.arcsin((x[3, [0]]-x[7, [0]]*znb+x[11, [0]]*xnb)/C1)

    abderivs = np.array(
        [
            [0, 0],
            [Xa, Xb],
            [0, 0],
            [Ya, Yb],
            [0, 0],
            [Za, Zb],
            [0, 0],
            [La, Lb],
            [0, 0],
            [Ma, Mb],
            [0, 0],
            [Na, Nb]
        ]
    )

    # Define the input sequence
    ti = 1
    tf = 5
    if t > ti and t < tf:
        u = np.array([[0],
                      [5],
                      [0]])
    else:
        u = np.zeros((3, 1))

    Einv = np.linalg.inv(E)
    xdot = np.matmul(np.matmul(Einv, A), x) + np.matmul(np.matmul(Einv,
                                                                  B), u) + np.matmul(abderivs, np.array([a, b]))
    # print('np.matmul(np.matmul(Einv,A),x)')
    # print(np.matmul(abderivs, np.array([a, b])))
    # print('iter \n')
    return xdot


x = rk4(kinematics, tspan, x0)

plt.figure
plt.plot(tspan,np.transpose(x))
plt.legend
plt.show()
