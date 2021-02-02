import numpy as np


def rk2(xdot, tspan, x0, print_timestamp):
    """
    Runge-Kutta 2nd Order Solver

    Inputs:
        xdot: kinematics function vector (ODE)
        tspan: time vector spanning t0:dt:tf
        x0: initial conditions vector

    Outputs:
        x: numerical solution to ODE input 'xdot'

    History: 
        Created and debugged, 01.11.2020

    """

    # Allocate space
    x = np.zeros((len(x0),len(tspan)))
    x[:,[0]] = x0[:,[0]]

    print("Starting Simulation")
    for ii in range(0,len(tspan)-1):
        dt = tspan[ii+1]-tspan[ii]
        
        # Come up with the RK4 coefficients
        temp1 = xdot(tspan[ii], x[:,[ii]])
        k1 = dt*temp1

        temp2 = xdot(tspan[ii] + 0.5*dt, x[:,[ii]] + 0.5*k1[:,[0]])
        k2 = dt*temp2

        # Solve for the next time step
        
        x[:,[ii+1]] = x[:,[ii]] + k2[:,[0]]

        if print_timestamp == True:
            print('\rRunning t = {:1.5f} (s)'.format(tspan[ii]), end='', flush=True)


    print("\n Simulation Complete")
    return x




