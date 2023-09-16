import grids
import numpy as np

from hexalattice.hexalattice import *
import time


nx = 1000
ny = 1000


start_lat = time.time()
hex_centers, _ = create_hex_grid(nx=nx,
                                 ny=ny,
                                 do_plot=True)
end_lat = time.time()

print(end_lat-start_lat)
# plt.show()  # import matplotlib.pyplot as plt


area_bounds = np.array([0, 0, nx, ny])


start_sir = time.time()
grid = grids.gen_grid(area_bounds, r=1, n_sides=6)
end_sir = time.time()

print(end_sir-start_sir)

print((end_lat-start_lat)/(end_sir-start_sir))


# print(grid.cell_centers().points)


