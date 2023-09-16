import pyvista as pv
import numpy as np

from vtkmodules.all import *
from shapely.geometry import Polygon, LineString
import geopandas as gpd

from helpers import shp_helpers
import concurrent.futures
import itertools


def filt_cells(args):
    cells, indexes, b_area = args
    rem_cell = []
    for index in indexes:
        points = cells.get_cell(index).points
        cell = Polygon(points)
        if not b_area.contains(cell):
            rem_cell.append(index)
    return rem_cell


def get_boundary(path: str) -> np.ndarray:

    data = gpd.read_file('data/shp/grid_pesce_grosso.shp')

    center = np.array(data.geometry[0].centroid.xy).flatten()

    int_area = data.translate(-center[0], -center[1]).geometry
    area_bounds = np.array(int_area.bounds).flatten()

    return area_bounds


# grid.cell_data['id'] = np.arange(0, grid.n_cells)
# area = shp_helpers.shp2vtk(int_area)
# area_loop = shp_helpers.shp2vtk(int_area.boundary)
# boundary = int_area.boundary
#
#
# with concurrent.futures.ProcessPoolExecutor() as executor:
#     chunks = np.array_split(np.arange(grid.n_cells), 12)
#
#     args = [[grid, chunk, int_area] for chunk in chunks]
#     results = executor.map(filt_cells, args)
#     indx_list = list(results)
#
# indx_list = list(itertools.chain(*indx_list))