# core stuff
import pickle as pk
import numpy as np

# meshing
import pyvista
import tetgen
import meshio as mio
import openmesh as om

"""
General info:
This script contains methods to transform/ convert between different kinds of mesh files.
"""


def read_pk_file(filename):
    """
    Reads in a .pk file and returns the vertices and triangles (faces)
    Args:
        filename (str): The filename of the .pk file

    Returns:
        mesh_points, mesh_triangles (tuple): list of mesh points and list of mesh triangles

    Notes:
        Adapted from
        https://github.com/darioizzo/geodesyNets/blob/master/archive/Modelling%20Bennu%20with%20mascons.ipynb

    """
    with open(filename, "rb") as f:
        mesh_points, mesh_triangles = pk.load(f)
    mesh_points = np.array(mesh_points)
    mesh_triangles = np.array(mesh_triangles)
    # Characteristic dimension
    L = max(mesh_points[:, 0]) - min(mesh_points[:, 0])

    # Non dimensional units
    mesh_points = mesh_points / L * 2 * 0.8
    return mesh_points, mesh_triangles


def write_to_node_faces_ele_file(path, filename, nodes, faces, ele):
    """
    Write to a .node, .face and .ele file.
    Args:
        path (str): the path where to write
        filename (str): the filename of the file to create
        nodes (list): list of vertices (cartesian coordinates)
        faces (list): list of faces (triangles consisting of vertices' indices)
        ele (list): list of polyhedral elements

    Returns:
        None

    """
    with open(path + filename + ".node", "w") as f:
        f.write("# Node count, 3 dimensions, no attribute, no boundary marker\n")
        f.write("{} {} {} {}\n".format(int(nodes.size / 3), 3, 0, 0))
        f.write("# Node index, node coordinates\n")
        index = 0
        for n in nodes:
            f.write("{} {} {} {}\n".format(index, n[0], n[1], n[2]))
            index += 1
    with open(path + filename + ".face", "w") as f:
        f.write("# Number of faces, boundary marker off\n")
        f.write("{} {}\n".format(int(faces.size / 3), 0))
        f.write("# Face index, nodes of face\n")
        index = 0
        for fac in faces:
            f.write("{} {} {} {}\n".format(index, fac[0], fac[1], fac[2]))
            index += 1
    with open(path + filename + ".ele", "w") as f:
        f.write("# number of tetrahedra, number of nodes per tet, no region attribute\n")
        f.write("{} {} {}\n".format(int(ele.size / 4), 4, 0))
        f.write("# tetrahedra index, nodes\n")
        index = 0
        for tet in ele:
            f.write("{} {} {} {} {}\n".format(index, tet[0], tet[1], tet[2], tet[3]))
            index += 1


def write__tsoulis_fortran_files(path, nodes, faces):
    """
    Writes to a topout, xyposnew and dataut file - exactly the input file which the
    FORTRAN implementation by Tsoulis requires as input.
    Args:
        path (str):  the path where to write
        nodes (list): list of vertices (cartesian coordinates)
        faces (list): list of faces (triangles consisting of vertices' indices)

    Returns:
        None
    """
    with open(path + "topoaut", "w") as f:
        for fac in faces:
            f.write(" {} {} {}\n".format(fac[0] + 1, fac[1] + 1, fac[2] + 1))
    with open(path + "xyzposnew", "w") as f:
        for n in nodes:
            f.write(" {} {} {}\n".format(n[0], n[1], n[2]))
    with open(path + "dataut", "w") as f:
        for fac in faces:
            f.write(" 3")


def write_other(path, name, nodes, faces):
    """
    Writes other files like .ply files.
    Args:
        path (str):  the path where to write
        name (str): the name of the file to be created, the suffix determines the concrete write operation
        nodes (list): list of vertices (cartesian coordinates)
        faces (list): list of faces (triangles consisting of vertices' indices)

    Returns:
        None
    """
    mesh = om.TriMesh()
    mesh.add_vertices(nodes)
    mesh.add_faces(faces)
    om.write_mesh((path+name), mesh)
    # Alternative Framework - Code
    # cells = [("triangle", faces)]
    # mesh = mio.Mesh(nodes, cells)
    # mesh.write((path + name))


def main():
    # TODO Pass file names via terminal?
    # Read the input .pk file
    print("Reading file...")
    mesh_points, mesh_triangles = read_pk_file("../mesh/input/Eros.pk")

    # Tetrahralize the mesh
    print("Tetrahralize...")
    tgen = tetgen.TetGen(mesh_points, mesh_triangles)
    nodes, elems = tgen.tetrahedralize()

    # Print a vtk file, diplayable with ParaView
    print("Writing vtk files...")
    tgen.write("../mesh/Eros/Eros.vtk")

    # Print the input files for the C++ implementation of the Polyhedral Model
    print("Writing polyhedron-cpp input files..")
    write_to_node_faces_ele_file("../mesh/Eros/", "Eros", nodes, mesh_triangles, elems)

    # Print the input files for the FORTRAN implementation of the Polyhedral Model by Tsoulis
    print("Writing FORTRAN files..")
    write__tsoulis_fortran_files("../mesh/Eros/", mesh_points, mesh_triangles)

    # Write other files
    print("Writing other files...")
    write_other("../mesh/Eros/", "Eros.ply", nodes, mesh_triangles)

    # Plot the tetrahralized mesh
    # print("Showing Plot")
    # # tgen.grid.plot(show_edges=True)
    # pl = pyvista.Plotter()
    # pl.add_mesh(tgen.grid)
    # # pl.add_axes_at_origin(labels_off=True)
    # pl.show()

    print("Finished.")


if __name__ == "__main__":
    main()
