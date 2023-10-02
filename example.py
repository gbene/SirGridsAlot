import grids
import Grids
import numpy as np

from vtkmodules.all import *
from helpers import plot_grid
import pyvista as pv

from hexalattice.hexalattice import *
import time


nx = 20
ny = 20
# plt.show()  # import matplotlib.pyplot as plt


area_bounds = np.array([0, 0, nx, ny])


# start_sir = time.time()
# grid = grids.gen_grid(area_bounds, r=1, n_sides=4)

grid = Grids.HexGrid(radius=1, bounds='/home/gabriele/STORAGE/Progetti/github/Py21/data/shp/grid_pesce_grosso.shp')


plotter = pv.Plotter()
plotter.add_mesh(grid.grid())
# grid.boundary_object.plot()
plotter.set_background('gray')
plotter.view_xy()
plotter.add_camera_orientation_widget()
plotter.show()



