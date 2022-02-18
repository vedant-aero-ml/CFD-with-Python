import numpy
from matplotlib import pyplot

nx = 41  # Number of grid points
dx = 2/(nx - 1)  # Distance between two grid points...working range is (0, 2)...dx = range/total panels
nt = 20  # Number of time steps we want
nu = 0.2
sigma = 0.2  # Sigma is a parameter, we'll learn more about it later
dt = sigma * dx**2 / nu


u = numpy.ones(nx)      # Array storing the initial values at every points = 1; its 2 at 0.5<=x<=1..them defined below
u[int(0.5 / dx):int(1 / dx + 1)] = 2  # Setting u = 2 between 0.5 and 1 as per our I.C.s

un = numpy.ones(nx)  # Temporary array to store values of Un

for j in range(nt):  # Running the temporal loop
    un = u.copy()  # Copying values of u into un
    for i in range(1, nx-1):
        u[i] = un[i] + nu * dt / dx ** 2 * (un[i + 1] - 2 * un[i] + un[i - 1])

pyplot.plot(numpy.linspace(0, 2, nx), u)
pyplot.grid()
pyplot.show()
