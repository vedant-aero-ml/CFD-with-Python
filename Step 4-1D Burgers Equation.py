import numpy
from matplotlib import pyplot
import sympy
import math

x, nu, t = sympy.symbols('x nu t')

phi = (sympy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) + sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * nu * (t + 1))))

phi_prime = phi.diff(x)  # Differentiating phi wrt x

# Defining  the IC
u = -2 * nu * (phi_prime / phi) + 4

# Now allowing th function to be operated by numerical inputs, i.e. converting SymPy to numerical functions
ufunc = sympy.utilities.lambdify((t, x, nu), u)

nx = 101
nt = 100
dx = 2 * math.pi / (nx - 1)
nu = 0.07
dt = nu * dx

# Forward marching solutions

un = numpy.empty(nx)  # The array to hold values for Un

for i in range(nt):
    un = u.copy()
    for j in range(1, nx-1):
        u[j] = un[j] - un[j] * dt / dx * (un[j] - un[j-1]) + nu * dt / dx**2 * (un[j+1] - 2 * un[j] + un[j-1])

    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * dt / dx ** 2 * (un[1] - 2 * un[0] + un[-2])
    u[-1] = u[0]  # Solution after the grid is exceeded reverts..as its a periodic BC


u_analytical = numpy.asarray([ufunc(nt * dt, xi, nu) for xi in x])  # Analytical solution

pyplot.plot(x, u, marker='o', lw=2, label='Computational')
pyplot.plot(x, u_analytical, label='Analytical')
pyplot.legend()
pyplot.show()
