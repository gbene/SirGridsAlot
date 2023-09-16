from shapely.geometry import Point, MultiPoint, LineString, MultiLineString
from vtkmodules.all import *
from pyvista import PolyData, lines_from_points
from numpy import array, zeros, hstack


def p2vtk(point: (Point, MultiPoint)) -> PolyData:
    out_xyz = array(point.coords)
    n_points = len(out_xyz)
    if out_xyz.ndim == 2:
        z = zeros((n_points, 1))
        out_xyz = hstack((out_xyz, z))
    vtk_obj = PolyData(out_xyz)
    return vtk_obj


def l2vtk(line: LineString) -> PolyData:

    out_xyz = array(line.coords)
    n_points = len(out_xyz)
    if out_xyz.ndim == 2:
        z = zeros((n_points, 1))
        out_xyz = hstack((out_xyz, z))
    vtk_obj = lines_from_points(out_xyz)
    return vtk_obj


def ml2vtk(multiline: MultiLineString) -> PolyData:

    appender = vtkAppendPolyData()
    for line in multiline:
        vtk_obj = l2vtk(line)
        appender.AddInputData(vtk_obj)
    appender.Update()
    vtk_obj = PolyData(appender.GetOutput())
    return vtk_obj


def shp2vtk(shp_obj):

    if shp_obj.geom_type == 'Point' or shp_obj.geom_type == 'MultiPoint':
        vtk_obj = p2vtk(shp_obj)
    elif shp_obj.geom_type == 'LineString':
        vtk_obj = l2vtk(shp_obj)
    elif shp_obj.geom_type == 'MultiLineString':
        vtk_obj = ml2vtk(shp_obj)
    elif shp_obj.geom_type == 'Polygon':
        if shp_obj.boundary.geom_type == "MultiLineString":
            vtk_obj = ml2vtk(shp_obj.boundary)
        else:
            vtk_obj = l2vtk(shp_obj.boundary)

        tr = vtkContourTriangulator()
        tr.AddInputData(vtk_obj)
        tr.Update()
        vtk_obj = tr.GetOutput()
    else:
        vtk_obj = None
    return PolyData(vtk_obj)
