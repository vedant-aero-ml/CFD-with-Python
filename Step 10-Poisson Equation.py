import numpy
from matplotlib import pyplot

nx = 51
ny = 51
nt = 50
dx = 2/(nx-1)
dy = 2/(ny-1)

p = numpy.zeros((nx, ny))  # Pressure term array
b = numpy.zeros((nx, ny))  # Source term array

b[int(nx/4), int(ny/40)] = 100  # Source term spikes
b[int(3*nx/4), int(3*ny/4)] = -100

pn = numpy.empty_like(p)  # Creating the array to hold values for Pn


for i in range(1, nt + 1):
    pn = p.copy()

    p[1:-1, 1:-1] = ((pn[2:, 1:-1] + pn[:-2, 1:-1])*dy**2 +
                     (pn[1:-1, 2:] + pn[1:-1, :-2])*dx**2 -
                     b[1:-1, 1:-1]*dx**2*dy**2) / (2*(dx**2 + dy**2))

    # BCs
    p[0, :] = 0
    p[-1, :] = 0
    p[:, 0] = 0
    p[:, -1] = 0

x, y = numpy.linspace(0.0, 2.0, nx), numpy.linspace(0.0, 1.0, ny)
X, Y = numpy.meshgrid(x, y)

# Plotting
pyplot.figure()
ax = pyplot.gca(projection='3d')
ax.plot_surface(X, Y, p)
pyplot.show()
