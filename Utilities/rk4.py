import numpy as np


def rk4(xdot, tspan, x0):
    """
    Runge-Kutta 4th Order Solver

    Inputs:
        xdot: kinematics function vector (ODE)
        tspan: time vector spanning t0:dt:tf
        x0: initial conditions vector

    Outputs:
        x: numerical solution to ODE input 'xdot'

    History: 
        Created and debugged, 01.05.2020

    """

    # Allocate space
    x = np.zeros((len(x0),len(tspan)))
    x[:,0] = x0[:,0]

    for ii in range(0,len(tspan)-1):
        dt = tspan[ii+1]-tspan[ii]
        
        # Come up with the RK4 coefficients
        temp1 = xdot(tspan[ii], x[:,[ii]])
        k1 = dt*temp1

        temp2 = xdot(tspan[ii] + 0.5*dt, x[:,[ii]] + 0.5*k1[:,[0]])
        k2 = dt*temp2
        
        temp3 = xdot(tspan[ii] + 0.5*dt, x[:,[ii]] + 0.5*k2[:,[0]])
        k3 = dt*temp3

        temp4 = xdot(tspan[ii+1], x[:,[ii]] + k3[:,[0]])
        k4 = dt*temp4

        # Solve for the next time step
        
        x[:,[ii+1]] = x[:,[ii]] + (k1[:,[0]] + 2*k2[:,[0]] + 2*k3[:,[0]] + k4[:,[0]])/6
    
    return x



