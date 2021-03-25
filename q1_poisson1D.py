#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW20 - Question 1. Laplace equation and membranes in 1D

Created on Wed Mar 24 23:04:28 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt

def gauss1D(u, P, T, dx, maxiter=25):
    """ Gauss-Seidel solution to Poisson elliptic equation in 1-dimension
    u: array, x-axis coordinates list
    P: float, uniform pressure on the membrane
    T: float, tension of the membrane
    dx: float, increment in the x-axis
    maxiter: int, max iterations
    """
    for j in range(maxiter):
        for i in range(1, len(u)-1):  # do not change the edges
            u[i] = 0.5 * (u[i+1] + u[i-1] - (P/T) * dx**2)
    return u

resolution = 10  # number of points along x-axis
minX = 0
maxX = 1
P = -1  # Uniform pressure on the membrane
T = 1  # Tension of the membrane

x = np.linspace(0, resolution, resolution+1)
u = np.zeros(resolution+1, dtype=float)
dx = (maxX - minX) / resolution  # increment

# Modifications
u[-1] = 2
P = -100
u[5] = 4
print(u)

# Solve using Gauss-Seidel function
u2 = gauss1D(u, P, T, dx)

plt.plot(x, u2)
plt.xlabel('x')
plt.ylabel('u')
plt.savefig('q1_poisson1D.png', dpi=300, bbox_inches='tight')