import grids
import numpy as np

from helpers import plot_grid
import pyvista as pv

from hexalattice.hexalattice import *
import time


nx = 12
ny = 4
# plt.show()  # import matplotlib.pyplot as plt


area_bounds = np.array([0, 0, nx, ny])


# start_sir = time.time()
grid = grids.gen_grid(area_bounds, r=1, n_sides=3)

plotter = pv.Plotter()

# plotter.add_mesh(grid, color='blue')
plotter.add_mesh(grid)


plotter.set_background('gray')
plotter.view_xy()
plotter.add_camera_orientation_widget()
plotter.show()

# end_sir = time.time()
# plot_grid.plot_grid_3d(grid)



