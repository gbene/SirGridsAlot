import pyvista as pv
import numpy as np


def gen_grid(bounds: np.ndarray, r: float = 0.3, n_sides: int = 4, return_centers: bool= False):

    """
    Efficiently create a diamond(90) or hex(60) grid in a rectangle boundary.

    The process works as follows:
        1. Define the angle theta and radius r
        2. Create a ranges of angles spanning 90-theta and 360 with a step theta
        3. Get the n of sides e.g. For theta = 90 n sides will be 4
        4. Create a zeros matrix of n_sides x 4 x 4
        5. Fill for each angle a translation matrix with the translation vector = [r*np.cos(angle), r*np.sin(angle), 0, 1]
        6. Create a grid of center points in range of bounds. The n of points in each direction is function of theta
        7. Offset the x coords every row by a factor function of theta
        8. The dot product between the center and the matrix results in n points = n sides translated by the angle

    :param n_sides: Number of sides of the grid cell
    :param r: Radius of circle enclosing the cell
    :param bounds: Boundary limits (xmin, ymin, xmax, ymax)
    :return:
    """

    theta = 360/n_sides
    theta_rad = np.deg2rad(theta)
    angles = np.deg2rad(np.arange((90-theta), 360, theta))
    n_sides = len(angles)

    trans_matrix = np.zeros((n_sides, 4, 4))

    for i, angle in enumerate(angles):
        trans_matrix[i] = np.array([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [r*np.cos(angle), r*np.sin(angle), 0, 1]
        ])

    nx = np.sin(theta_rad)*2*r
    ny = r + np.cos(theta_rad)*r

    xmin, ymin, xmax, ymax = bounds
    x = np.arange(xmin, xmax, nx)
    y = np.arange(ymin, ymax, ny)

    xv, yv = np.meshgrid(x, y)
    zv = np.zeros_like(xv)
    ones = np.ones_like(xv)
    xv[::2, :] += np.sin(theta_rad)*r

    points = np.hstack((xv.reshape(-1, 1), yv.reshape(-1, 1), zv.reshape(-1, 1), ones.reshape(-1, 1)))
    if return_centers:
        return points
    else:
        hex_grid = points.dot(trans_matrix)[:, :, :-1]
        n_points = np.shape(hex_grid)[0]*np.shape(hex_grid)[1]
        conn = np.insert(np.arange(0, n_points), np.arange(0, n_points, n_sides), n_sides)

        vtk_obj = pv.PolyData(hex_grid.reshape(-1, 3), faces=conn)
        return vtk_obj





#
# grid.remove_cells(indx_list, inplace=True)
# transform_matrix = np.array(
#     [
#         [1, 0, 0, 0],
#         [0, 1, 0, 0],
#         [0, 0, 1, -0.5],
#         [0, 0, 0, 1],
#     ])
#
# area_trans = area.transform(transform_matrix, inplace=False)
# # test_bound = shp2vtk(data.geometry[0].boundary)
# extr = vtkLinearExtrusionFilter()
# extr.SetInputData(area_trans)
# extr.CappingOn()
# extr.SetScaleFactor(1)
# extr.Update()
#
# select = vtkSelectEnclosedPoints()
# select.SetInputData(grid)
# select.SetSurfaceData(extr.GetOutput())
# threshold = vtkMultiThreshold()
# insideId = threshold.AddBandpassIntervalSet(
#         1, 1,
#         vtkDataObject.FIELD_ASSOCIATION_POINTS, 'SelectedPoints',
#         0, 1)
# threshold.SetInputConnection(select.GetOutputPort())
# threshold.OutputSet(insideId)
# threshold.Update()
# print(threshold.GetOutput())

