import pyvista as pv
import numpy as np
from vtkmodules.vtkFiltersCore import vtkAppendPolyData, vtkCleanPolyData


def gen_grid(bounds: np.ndarray, r: float = 0.3, n_sides: int = 4, return_centers: bool= False, vertex: str = 'up'):

    """
    Efficiently create a diamond(90) or hex(60) grid in a rectangle boundary.

    The process works as follows:
        1. Define the angle theta and radius r
        2. Create a ranges of angles spanning 90-theta and 360 with a step theta
        3. Get the n of sides e.g. For theta = 90 n sides will be 4
        4. Create a zeros matrix of n_sides x 4 x 4
        5. Fill for each angle a translation matrix with the translation vector = [r*np.cos(angle), r*np.sin(angle), 0, 1]
        6. Create a grid of center centers in range of bounds. The n of centers in each direction is function of theta
        7. Offset the x coords every row by a factor function of theta
        8. The dot product between the center and the matrix results in n centers = n sides translated by the angle

    :param n_sides: Number of sides of the grid cell
    :param r: Radius of circle enclosing the cell
    :param bounds: Boundary limits (xmin, ymin, xmax, ymax)
    :param return_centers: Return the centers of the cells
    :param vertex: Decide the orientation of the squares and hexagons. For vertex='up' (default), square grids will be
    diamond grids and hexagons will be "pointing" up. For vertex='down', square grids will be rectilinear and hexagons
    will lie flat.

    :return:
    """
    if n_sides == 3:
        r *= np.sqrt(3)
        tri = 1
        n_sides = 6
        bounds[:2] = bounds[:2][::-1]
        bounds[2:] = bounds[2:][::-1]
        vertex = 'up'
    else:
        tri = 0

    theta = 360/n_sides
    theta_rad = np.deg2rad(theta)

    fac = 1
    if vertex == 'up':
        angles = np.deg2rad(np.arange((90-theta), 360, theta))

        nx = np.sin(theta_rad)*2*r
        ny = r + np.absolute(np.cos(theta_rad))*r

    elif vertex == 'down':
        if n_sides == 4:
            angles = np.deg2rad(np.arange(45, 360, theta))
            nx = ny = 2*r/np.sqrt(2)
            fac = 0
        elif n_sides == 6:
            angles = np.deg2rad(np.arange(0, 360, theta))
            nx = r + np.absolute(np.cos(theta_rad))*r
            ny = np.sin(theta_rad)*2*r

    n_angles = len(angles)
    trans_matrix = np.zeros((n_angles, 4, 4))

    for i, angle in enumerate(angles):
        trans_matrix[i] = np.array([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [r*np.cos(angle), r*np.sin(angle), 0, 1]
        ])

    xmin, ymin, xmax, ymax = bounds
    x = np.arange(xmin, xmax, nx)
    y = np.arange(ymin, ymax, ny)

    xv, yv = np.meshgrid(x, y)

    n_col, n_rows = np.shape(yv)

    zv = np.zeros_like(xv)
    ones = np.ones_like(xv)

    if vertex == 'up':
        xv[::2, :] -= np.sin(theta_rad)*r
    elif vertex == 'down':
        yv[:, ::2] -= fac*np.sin(theta_rad)*r


    centers = np.hstack((xv.reshape(-1, 1), yv.reshape(-1, 1), zv.reshape(-1, 1), ones.reshape(-1, 1)))

    if return_centers:
        return centers
    else:
        hex_grid = centers.dot(trans_matrix)[:, :, :-1]

        n_points = np.shape(hex_grid)[0]*np.shape(hex_grid)[1]
        hex_grid = hex_grid.reshape(-1, 3)

        if tri:

            new_hexgrid = hex_grid.reshape(n_col, n_rows, 18)

            new_hexgrid[:, 0][::2][:, 6:12] = -1
            new_hexgrid[:, -1][1::2][:, :3] = -1
            new_hexgrid[:, -1][1::2][:, 15:] = -1
            new_hexgrid = new_hexgrid.reshape(-1, 3)
            rem_index = np.where(np.all(new_hexgrid == -1, axis=1) == 1)

            vert_idx = np.arange(0, n_points).reshape(-1, 6)

            vert_idx = np.roll(np.repeat(vert_idx, 2, axis=1), -1, axis=1).reshape(-1, 2)

            zeros_vert = np.zeros((n_points, 1), dtype=int)

            part1 = np.append(zeros_vert, vert_idx, axis=1).flatten()

            len_index = max(np.delete(vert_idx, rem_index[0], axis=0).flatten())

            center_grid = np.append(new_hexgrid, centers[:, :-1], axis=0)

            centers_idx = np.repeat(np.arange(n_points, n_points + len(centers)), 6).reshape(-1, 1)

            zeros_centers = np.zeros((len(centers_idx), 2), dtype=int)

            part2 = np.append(centers_idx, zeros_centers, axis=1).flatten()

            tri_conn = (part1+part2).reshape(-1, 3)

            tri_conn = np.insert(tri_conn, np.arange(0, len(tri_conn.flatten()), 3), 3)

            vtk_obj = pv.PolyData(center_grid, faces=tri_conn)
            vtk_obj.remove_points(rem_index[0], inplace=True)
            vtk_obj.rotate_z(angle=90, point=vtk_obj.center, inplace=True)
            vtk_obj.flip_x(inplace=True)

        else:
            conn = np.insert(np.arange(0, n_points), np.arange(0, n_points, n_angles), n_angles)
            vtk_obj = pv.PolyData(hex_grid, faces=conn)



        # #
        # #
        # # # if n_sides == 3:
        # # #     appender = vtkAppendPolyData()
        # # #     vtk_copy = vtk_obj.rotate_z(point=vtk_obj.center, angle=180)
        # # #     # appender.AddInputData(vtk_obj)
        # # #     # appender.AddInputData(vtk_copy)
        # #
        # #
        # #
        # #
        # #
        # #
        # #
        vtk_obj.cell_data['cellid'] = np.arange(vtk_obj.n_cells)

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

