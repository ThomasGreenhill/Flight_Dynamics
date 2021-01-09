import numpy as np
import matplotlib.pyplot as plt
import sys
import time

sys.path.append('../Utilities')

from rk4 import rk4

## Check using duffing's equation

tspan = np.linspace(0,50,5001)
x0 = np.array([[0], [1]])

def duffing(t, x):
    w = 1 
    a = 1
    b = 1
    d = 0.22
    ff = 0.3

    xdot = np.array([[x[1]], [ff*np.cos(w*t)-d*x[1]+b*x[0]-a*(x[0])**3]])

    return xdot

t = time.time()

x = rk4(duffing, tspan, x0)

elapsed = time.time() - t

# print(elapsed)

plt.figure()
plt.plot(tspan,x[1,:])
plt.show()

plt.figure()
plt.plot(tspan,x[0,:])
plt.show()

import cProfile
cProfile.run('rk4(duffing, tspan, x0)')