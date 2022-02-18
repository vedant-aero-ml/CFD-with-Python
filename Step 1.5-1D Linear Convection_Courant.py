import numpy
from matplotlib import pyplot


def linearconv(nx):
    dx = 2 / (nx - 1)
    nt = 20  # nt is the number of timesteps we want to calculate
    c = 1
    sigma = 0.5  # Courant Number

    dt = sigma * dx  # Calculating safe dt value by sigma (Courant Number)

    u = numpy.ones(nx)
    u[int(.5 / dx):int(1 / dx + 1)] = 2

    un = numpy.ones(nx)

    for n in range(nt):  # Temporal Loop
        un = u.copy()  # Copy the existing values of u into un
        for i in range(1, nx):
            u[i] = un[i] - c * dt / dx * (un[i] - un[i - 1])
    return u


nx = 89
u = linearconv(nx)
pyplot.plot(numpy.linspace(0, 2, nx), u)
pyplot.show()