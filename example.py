import grids
import numpy as np

from helpers import plot_grid
import pyvista as pv

from hexalattice.hexalattice import *
import time


nx = 20
ny = 20
# plt.show()  # import matplotlib.pyplot as plt


area_bounds = np.array([0, 0, nx, ny])


# start_sir = time.time()
grid = grids.gen_grid(area_bounds, r=1, n_sides=4, vertex='down')

plot_grid.plot_grid_3d(grid, cell_data='cellid')



