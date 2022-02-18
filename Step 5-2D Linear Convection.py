from matplotlib import pyplot
import numpy


nx = 101  # Number of grid points
ny = 101
dx = 2/(nx - 1)  # Distance between two grid points...working range is (0, 2)...dx = range/total panels
dy = 2/(ny - 1)
c = 1
nt = 25  # Number of time steps we want
sigma = 0.2
dt = sigma * dx
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)
X, Y = numpy.meshgrid(x, y)

u = numpy.ones((ny, nx))  # Array storing the initial values at every points = 1
u[int(0.5 / dy):int(1 / dy + 1), int(0.5 / dx):int(1 / dx + 1)] = 2  # Setting u = 2 between 0.5 and 1 as per our I.C.s

un = numpy.ones((ny, nx))  # Temporary array to store values of Un

'''
# Plotting the space at t = 0
fig = pyplot.figure(dpi=100)
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, u[:])
'''

# Calculating the values in temporal stepping
for i in range(1, nt + 1):
    un = u.copy()
    # (ny, nx) hai to ulta chalega
    u[1:, 1:] = (un[1:, 1:] - (c * dt / dx * (un[1:, 1:] - un[1:, :-1])) - (c * dt / dy * (un[1:, 1:] - un[:-1, 1:])))

    # Boundary Conditions
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

# Plotting the space at the end of the time
fig = pyplot.figure(dpi=100)
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, u)
pyplot.show()
