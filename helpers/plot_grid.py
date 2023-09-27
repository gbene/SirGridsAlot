import pyvista as pv
import matplotlib.pyplot as plt


def plot_grid_3d(grid:  pv.PolyData, cell_data=None):

    plotter = pv.Plotter()
    if cell_data:
        plotter.add_mesh(grid, scalars=cell_data)
    else:
        plotter.add_mesh(grid, style='wireframe', color='red')
    plotter.set_background('gray')
    plotter.view_xy()
    plotter.add_camera_orientation_widget()
    plotter.show()