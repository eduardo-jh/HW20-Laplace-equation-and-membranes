#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BE523 Biosystems Analysis & Design
HW20 - Question 1. Laplace equation and membranes in 2D with dx=dy

Created on Thu Mar 25 00:18:48 2021
@author: eduardo
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def gauss2D(u, P, T, dx, maxiter=25):
    """ Gauss-Seidel solution to Poisson elliptic equation in 2-dimensions
    u: array, x-axis coordinates list
    P: float, uniform pressure on the membrane
    T: float, tension of the membrane
    dx: float, increment in the x-axis
    maxiter: int, max iterations
    """
    for k in range(maxiter):
        for i in range(1, len(u)-1):  # do not change the edges
            for j in range(1, len(u)-1):
                u[i][j] = (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - (P/T) * dx**2)/4
    return u

resolution = 10  # number of points along x-axis
minX = minY = -1
maxX = maxY = 1
P = -1  # Uniform pressure on the membrane
T = 1  # Tension of the membrane

x = np.linspace(0, resolution, resolution+1)
u = np.zeros(shape=(resolution, resolution), dtype=float)
dx = (maxX - minX) / resolution  # increment

# change values at the left and right sides
for i in range(0, 5):
    u[i][0] = u[i][-1] = i
for i in range(5, 10):
    u[i][0] = u[i][-1] = 10 - i
P = -100

# Solve using Gauss-Seidel function
sol = gauss2D(u, P, T, dx)

# Create a 3D-plot
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 1, 1, projection='3d')
dy = dx
X = np.arange(minX, maxX, dx)
Y = np.arange(minY, maxY, dy)
X, Y = np.meshgrid(X, Y)
Z = sol

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
ax.set_zlim3d(np.min(sol), np.max(sol))
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
fig.colorbar(surf, shrink=0.5, aspect=10)
fig.savefig('q2_poisson2D.png', dpi=300, bbox_inches='tight')