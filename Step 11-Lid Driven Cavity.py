# LID DRIVEN CAVITY PROBLEM USING FINITE DIFFERENCE METHOD

import numpy
from matplotlib import pyplot

nx = 51  # Number of mesh points in x direction
ny = 51  # Number of mesh points in y direction
nt = 400  # Number of time steps
nit = 50  # Number of time step for pressure iterations
domain = 1.0  # Domain of the mesh system

dx = domain / (nx - 1)
dy = domain / (ny - 1)

rho = 1  # Density
nu = 0.1  # Kinematic Viscosity
dt = 0.001  # Length of time step

u = numpy.zeros((ny, nx))  # Creating the initial array for values of u
v = numpy.zeros((ny, nx))  # Creating the initial array for values of v
p = numpy.zeros((ny, nx))  # Creating the initial array for values of pressure

un = numpy.empty_like(u)  # Creating a null valued array for storing values of u for nth iteration
vn = numpy.empty_like(v)  # Creating a null valued array for storing values of v for nth iteration
pn = numpy.empty_like(p)  # Creating a null valued array for storing values of p for nth iteration

# LOOPING FOR TIME ITERATIONS
'''
Using the following schemes of Finite Difference Method:

Forward Difference in Time
Backward Difference in Space
Central Difference for Poisson Equation of Pressure
'''
for i in range(1, nt + 1):

    # Copying values of u and v from previous step
    vn = v.copy()
    un = u.copy()

    for n in range(1, nit + 1):  # Inner loop covering the Poisson Pressure term

        # Copying values of Pressure from previous step
        pn = p.copy()

        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy ** 2 + (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx ** 2) /
                         (2 * (dx ** 2 + dy ** 2)) - dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) *
                         (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) /
                                           (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                                 ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)) ** 2 -
                                 2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                                      (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) -
                                 ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) ** 2)))

        # BCs of Pressure
        p[:, -1] = p[:, -2]  # dp/dx = 0 at x = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]  # dp/dx = 0 at x = 0
        p[-1, :] = 0  # p = 0 at y = 2

    # ITERATING FOR VALUES OF u
    u[1:-1, 1:-1] = (un[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                     dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                     nu * (dt / dx ** 2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                           dt / dy ** 2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

    # ITERATING FOR VALUES OF v
    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                     dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                     nu * (dt / dx ** 2 * (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) + dt / dy ** 2 *
                           (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

    # BCs of u
    u[0, :] = 0  # u = 0 at y = 0
    u[-1, :] = 1  # u = 1 at y = 2, i.e on the lid
    u[:, 0] = 0  # u = 0 at x = 0
    u[:, -1] = 0  # u = 0 at x = 2

    # BCs of v
    v[0, :] = 0  # v = 0 at y = 0
    v[-1, :] = 0  # v = 0 at y = 2
    v[:, 0] = 0  # v = 0 at x = 0
    v[:, -1] = 0  # v = 0 at x = 2


# PLOTTING THE PRESSURE CONTOURS AND VELOCITY STREAMLINES TOGETHER

x = numpy.linspace(0, domain, nx)  # Grid points along x-axis
y = numpy.linspace(0, domain, ny)  # Grid points along y-axis
X, Y = numpy.meshgrid(x, y)  # Creating the grid for plotting
V = numpy.sqrt(numpy.square(u) + numpy.square(v))  # Absolute velocity

pyplot.contourf(X, Y, p)  # Plotting the pressure contours
pyplot.colorbar()  # Showing the legend for value comparision
pyplot.streamplot(X, Y, u, v, color='white')  # Showing velocity streamlines in white colour
pyplot.show()  # Executing the data visualization
