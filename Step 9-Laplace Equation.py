from matplotlib import pyplot
import numpy

nx = 31
ny = 31
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
l1_norm_target = 0.0001
l1_norm = 1

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 1, ny)
X, Y = numpy.meshgrid(x, y)

p = numpy.zeros((ny, nx))  # Creating the initial arrays of zeros

p[:, 0] = 0  # p = 0 @ x = 0
p[:, -1] = y  # p = y @ x = 2
p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

pn = numpy.empty_like(p)

while l1_norm > l1_norm_target:
    pn = p.copy()

    p[1:-1, 1:-1] = ((dy ** 2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) + dx ** 2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) /
                     (2 * (dx ** 2 + dy ** 2)))

    p[:, 0] = 0  # p = 0 @ x = 0
    p[:, -1] = y  # p = y @ x = 2
    p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
    p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

    l1_norm = (numpy.sum(numpy.abs(p[:]) - numpy.abs(pn[:])) / numpy.sum(numpy.abs(pn[:])))


# Plotting
pyplot.figure()
ax = pyplot.gca(projection='3d')
ax.plot_surface(X, Y, p)
pyplot.show()
