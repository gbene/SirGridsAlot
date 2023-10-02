from abc import ABC, abstractmethod
import numpy as np
import pyvista as pv
import geopandas as gpd
from shapely.geometry import MultiLineString

class BaseGrid(ABC):
    """
    Abstract class for grids:

    1. TriGrid
    2. QuadGrid
    3. HexGrid
    4. GridCollection
    """

    def __init__(self, radius: float, bounds, n_sides: int, orientation: str):
        """
        Init the entity.
        """

        self._radius = radius

        if isinstance(bounds, np.ndarray):
            self._boundary = bounds
        elif isinstance(bounds, str):

            boundary_shp = gpd.read_file(bounds)['geometry'][0]

            self._boundary = np.array(boundary_shp.bounds)

        self._n_sides = n_sides
        self._orientation = orientation

        self._centers = None
        self._nx = self._ny = 0

        self._grid_object = self.gen_grid()

    @property
    def radius(self) -> float:
        """
        Get the radius of the circumscribed circle.

        :return: float
        """
        return self._radius

    def side(self) -> float:
        """
        Get the length of the side of a cell
        :return: float
        """
        theta_rad = self.theta(rad_only=True)

        side = 2*self.radius*np.sin(theta_rad/2)

        return side

    def height(self) -> float:
        """
        Get the distance between the center of the cell and a side
        :return: float
        """
        theta_rad = self.theta(rad_only=True)

        height = self.radius * np.cos(theta_rad/2)

        return height

    @property
    def bounds(self) -> np.ndarray:
        """
        Get the boundaries of the grid (xmin, ymin, xmax, ymax).

        :return: numpy array
        """
        return self._boundary

    @property
    def boundary_object(self) -> pv.PolyData:
        """
        Get the representation of the shp boundary.

        :return: Polydata object
        """

        return self._boundary_object

    @property
    def number_of_sides(self) -> int:
        """
        Get the number of sides of the cell of the grid

        :return: integer
        """
        return self._n_sides

    @property
    def cell_orientation(self) -> str:
        """
        Get the orientation for hex and square grids ("up" or "flat")

        :return: string
        """
        return self._orientation

    @property
    def cell_type(self) -> str:
        """
        Get the name of the cell of the grid.

        + 3: Triangle
        + 4:
            + Diamond if orientation is up
            + Square if orientation is flat
        + 6:
            + Hexagon if orientation is up
            + Flaxagon if orientation is down

        :return: string
        """

        if self.number_of_sides == 3:
            return "triangle"
        elif self.number_of_sides == 4:
            if self.cell_orientation == 'up':
                return "diamond"
            else:
                return "square"
        else:
            if self.cell_orientation == 'up':
                return "hexagon"
            else:
                return "flaxagon"

    @property
    def spacing(self) -> list:
        """
        Get the spacing between center points of the grid [nx, ny]

        :return: list
        """
        return [self._nx, self._ny]

    @property
    def centers(self) -> np.ndarray:
        """
        Get the grid cell center coordinates

        :return: numpy array
        """
        centers = self._centers[:, :-1]
        return centers.reshape(-1, 3)

    @abstractmethod
    def gen_grid(self) -> pv.PolyData:
        """
        Generate the grid given the different options. A pyvista polydata object is created with m cells of sides n_sides
        of the grid.

        Each grid has a different way to generate the grid so this needs to be implemented in the different grid classes

        :return: pyvista Polydata
        """
        pass

    @abstractmethod
    def area(self) -> float:
        """
        Calculate the area of the given cell. Each grid will have its own area calculation

        :return: float
        """
        pass

    def grid(self) -> pv.PolyData:
        """
        Get the grid pyvista entity.

        :return: pyvista PolyData
        """

        return self._grid_object

    def theta(self, rad_only: bool = False) -> list:
        """
        Get a list of angular values [degree, radians] of angle between two points in the grid
        :return: list
        """

        theta = 360 / self.number_of_sides
        theta_rad = np.deg2rad(theta)

        if rad_only:
            return theta_rad
        else:
            return [theta, theta_rad]

    def matrix(self, angles: np.ndarray) -> np.ndarray:

        n_angles = len(angles)
        r = self.radius
        trans_matrix = np.zeros((n_angles, 4, 4))

        for i, angle in enumerate(angles):
            trans_matrix[i] = np.array([
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [r * np.cos(angle), r * np.sin(angle), 0, 1]
            ])

        return trans_matrix


class TriGrid(BaseGrid):
    """
    Class for tri grid specific
    """

    def __init__(self, radius: float, bounds, n_sides: int = 6, orientation: str= 'up'):
        super().__init__(radius, bounds, n_sides, orientation)

    def gen_grid(self) -> pv.PolyData:

        theta, theta_rad = self.theta()
        r = self.radius

        self.bounds[:2] = self.bounds[:2][::-1]
        self.bounds[2:] = self.bounds[2:][::-1]

        angles = np.deg2rad(np.arange((90 - theta), 360, theta))

        self._nx = np.sin(theta_rad) * 2 * r
        self._ny = r + np.absolute(np.cos(theta_rad)) * r

        n_angles = len(angles)
        trans_matrix = self.matrix(angles)

        x_min, y_min, x_max, y_max = self.bounds
        x = np.arange(x_min, x_max, self._nx)
        y = np.arange(y_min, y_max, self._ny)

        xv, yv = np.meshgrid(x, y)

        n_col, n_rows = np.shape(yv)

        zv = np.zeros_like(xv)
        ones = np.ones_like(xv)

        xv[::2, :] -= np.sin(theta_rad) * r

        self._centers = np.hstack((xv.reshape(-1, 1), yv.reshape(-1, 1), zv.reshape(-1, 1), ones.reshape(-1, 1)))

        hex_grid = self._centers.dot(trans_matrix)[:, :, :-1]

        n_points = np.shape(hex_grid)[0] * np.shape(hex_grid)[1]
        hex_grid = hex_grid.reshape(-1, 3)

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

        center_grid = np.append(new_hexgrid, self.centers, axis=0)

        centers_idx = np.repeat(np.arange(n_points, n_points + len(self.centers)), 6).reshape(-1, 1)

        zeros_centers = np.zeros((len(centers_idx), 2), dtype=int)

        part2 = np.append(centers_idx, zeros_centers, axis=1).flatten()

        tri_conn = (part1 + part2).reshape(-1, 3)

        tri_conn = np.insert(tri_conn, np.arange(0, len(tri_conn.flatten()), 3), 3)

        vtk_obj = pv.PolyData(center_grid, faces=tri_conn)

        vtk_obj.remove_points(rem_index[0], inplace=True)

        vtk_obj.rotate_z(angle=90, point=vtk_obj.center, inplace=True)
        vtk_obj.flip_x(inplace=True)

        vtk_obj.cell_data['cellid'] = np.arange(vtk_obj.n_cells)

        return vtk_obj

    def area(self) -> float:

        """
        The area of the hexagon

        :return: float
        """
        side = self.side()
        height = self.height()
        area = side*height/2
        return area


class QuadGrid(BaseGrid):
    """
    Class for quad grid specific
    """

    def __init__(self, radius: float, bounds, n_sides: int = 4, orientation: str= 'up'):
        super().__init__(radius, bounds, n_sides, orientation)

    def gen_grid(self) -> pv.PolyData:

        theta, theta_rad = self.theta()
        r = self.radius

        if self.cell_orientation == 'up':
            angles = np.deg2rad(np.arange((90 - theta), 360, theta))

            self._nx = np.sin(theta_rad) * 2 * r
            self._ny = r + np.absolute(np.cos(theta_rad)) * r

        else:
            angles = np.deg2rad(np.arange(45, 360, theta))

            self._nx = self._ny = 2*r/np.sqrt(2)

        n_angles = len(angles)
        trans_matrix = self.matrix(angles)

        x_min, y_min, x_max, y_max = self.bounds
        x = np.arange(x_min, x_max, self._nx)
        y = np.arange(y_min, y_max, self._ny)

        xv, yv = np.meshgrid(x, y)

        n_col, n_rows = np.shape(yv)

        zv = np.zeros_like(xv)
        ones = np.ones_like(xv)

        if self.cell_orientation == 'up':
            xv[::2, :] -= np.sin(theta_rad) * r

        self._centers = np.hstack((xv.reshape(-1, 1), yv.reshape(-1, 1), zv.reshape(-1, 1), ones.reshape(-1, 1)))

        hex_grid = self._centers.dot(trans_matrix)[:, :, :-1]

        n_points = np.shape(hex_grid)[0] * np.shape(hex_grid)[1]
        hex_grid = hex_grid.reshape(-1, 3)

        conn = np.insert(np.arange(0, n_points), np.arange(0, n_points, n_angles), n_angles)
        vtk_obj = pv.PolyData(hex_grid, faces=conn)

        vtk_obj.cell_data['cellid'] = np.arange(vtk_obj.n_cells)
        return vtk_obj

    def area(self) -> float:

        """
        The area of the hexagon

        :return: float
        """
        side = self.side()
        area = np.power(side, 2)
        return area


class HexGrid(BaseGrid):
    """
    Class for hex grid specific
    """

    def __init__(self, radius: float, bounds, n_sides: int = 6, orientation: str= 'up'):
        super().__init__(radius, bounds, n_sides, orientation)

    def gen_grid(self) -> pv.PolyData:

        theta, theta_rad = self.theta()
        r = self.radius

        if self.cell_orientation == 'up':
            angles = np.deg2rad(np.arange((90 - theta), 360, theta))

            self._nx = np.sin(theta_rad) * 2 * r
            self._ny = r + np.absolute(np.cos(theta_rad)) * r

        else:
            angles = np.deg2rad(np.arange(0, 360, theta))

            self._nx = r + np.absolute(np.cos(theta_rad))*r
            self._ny = np.sin(theta_rad)*2*r

        n_angles = len(angles)
        trans_matrix = self.matrix(angles)

        x_min, y_min, x_max, y_max = self.bounds
        x = np.arange(x_min, x_max, self._nx)
        y = np.arange(y_min, y_max, self._ny)

        xv, yv = np.meshgrid(x, y)

        n_col, n_rows = np.shape(yv)

        zv = np.zeros_like(xv)
        ones = np.ones_like(xv)

        if self.cell_orientation == 'up':
            xv[::2, :] -= np.sin(theta_rad) * r
        else:
            yv[:, ::2] -= np.sin(theta_rad) * r

        self._centers = np.hstack((xv.reshape(-1, 1), yv.reshape(-1, 1), zv.reshape(-1, 1), ones.reshape(-1, 1)))

        hex_grid = self._centers.dot(trans_matrix)[:, :, :-1]

        n_points = np.shape(hex_grid)[0] * np.shape(hex_grid)[1]
        hex_grid = hex_grid.reshape(-1, 3)

        conn = np.insert(np.arange(0, n_points), np.arange(0, n_points, n_angles), n_angles)
        vtk_obj = pv.PolyData(hex_grid, faces=conn)

        vtk_obj.cell_data['cellid'] = np.arange(vtk_obj.n_cells)
        return vtk_obj

    def area(self) -> float:

        """
        The area of the hexagon

        :return: float
        """
        side = self.side()
        area = np.power(side, 2)*3*np.sqrt(3)/2
        return area

