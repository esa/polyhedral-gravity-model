import sys
import time

import open3d as o3d
import numpy as np
import math

class Node:
    idx: int
    x: float
    y: float
    z: float

    def __init__(self, idx: int, x: float, y: float, z: float):
        self.idx = idx
        self.x = x
        self.y = y
        self.z = z

    @staticmethod
    def from_line(line):
        index, *floats = line.split(" ")
        index = int(index)
        floats = list(map(float, floats))
        return Node(index, *floats)

    def __add__(self, other):
        return Node(-1, self.x + other.x, self.y + other.y, self.z + other.z)
    
    def vertex(self):
        return [self.x,self.y,self.z]


class Face:
    idx: int
    v_idx1: int
    v_idx2: int
    v_idx3: int

    def __init__(self, idx: int, v_idx1: int, v_idx2: int, v_idx3: int):
        self.idx = idx
        self.v_idx1 = v_idx1
        self.v_idx2 = v_idx2
        self.v_idx3 = v_idx3

    def mid_point(self, nodes):
        node = nodes[self.v_idx1] + nodes[self.v_idx2] + nodes[self.v_idx3]
        node.idx = len(nodes)
        return node

    def generate_split_faces(self, node_idx, next_face_idx):
        return [
            Face(next_face_idx, self.v_idx1, self.v_idx2, node_idx),
            Face(next_face_idx + 1, self.v_idx2, self.v_idx3, node_idx),
            Face(next_face_idx + 2, self.v_idx3, self.v_idx1, node_idx),
        ]

    @staticmethod
    def from_line(line):
        results = list(map(int, line.split(" ")))
        return Face(*results)
    
    def triangle(self):
        return [self.v_idx1, self.v_idx2, self.v_idx3]


def fetch_data(node_file, face_file):
    with open(node_file, "r") as f:
        node_lines = f.readlines()
    with open(face_file, "r") as f:
        face_lines = f.readlines()
    nodes = []
    faces = []
    skipped_first = False
    for line in node_lines:
        if line.strip().startswith("#"):
            continue
        if not skipped_first:
            skipped_first = True
            continue
        nodes.append(Node.from_line(line))

    skipped_first = False
    for line in face_lines:
        if line.strip().startswith("#"):
            continue
        if not skipped_first:
            skipped_first = True
            continue
        faces.append(Face.from_line(line))
    return nodes, faces

def scale_mesh(nodes, faces, number_of_faces):
    vertices = list(map(lambda n: n.vertex(), nodes))
    triangles = list(map(lambda n: n.triangle(), faces))
    mesh = o3d.geometry.TriangleMesh(vertices=o3d.utility.Vector3dVector(vertices), triangles=o3d.utility.Vector3iVector(triangles))
    iterations = max(0, math.ceil(math.log(number_of_faces / len(triangles)) / math.log(4)))
    if iterations > 0:
        mesh = mesh.subdivide_loop(number_of_iterations=iterations)
    mesh = mesh.simplify_quadric_decimation(target_number_of_triangles=number_of_faces)
    mesh = mesh.remove_unreferenced_vertices()
    vertices = np.asarray(mesh.vertices).tolist()
    triangles = np.asarray(mesh.triangles).tolist()
    nodes = []
    faces = []
    for i in range(len(vertices)):
        nodes.append(Node(i,*vertices[i]))
    for i in range(len(triangles)):
        faces.append(Face(i, *triangles[i]))
    return nodes, faces


def write_to_file(nodes, faces, filename):
    with open(f"{filename}.node", "w") as f:
        f.write("# Adapted from https://github.com/darioizzo/geodesyNets/tree/master/3dmeshes\n")
        f.write("# Node count, 3 dimensions, no attribute, no boundary marker\n")
        f.write(f"{len(nodes)} 3 0 0\n")
        f.write("# Node index, node coordinates\n")

        for node in nodes:
            line = f"{node.idx} {node.x} {node.y} {node.z}\n"
            f.write(line)

    with open(f"{filename}.face", "w") as f:
        f.write("# Adapted from https://github.com/darioizzo/geodesyNets/tree/master/3dmeshes\n")
        f.write("# Number of faces, boundary marker off\n")
        f.write(f"{len(faces)} 0\n")
        f.write("# Face index, nodes of face\n")

        for face in faces:
            line = f"{face.idx} {face.v_idx1} {face.v_idx2} {face.v_idx3}\n"
            f.write(line)


def main():
    if len(sys.argv) != 3:
        print("Usage: python scale_mesh.py <node_file> <face_file>")
        return -1
    nodes, faces = fetch_data(sys.argv[1], sys.argv[2])
    filename = sys.argv[1].split(".")[0]
    face_amounts = [round(1000* math.sqrt(3)**k) for k in range(10)]
    for amount in face_amounts:
        nodes, faces = scale_mesh(nodes, faces, amount)
        write_to_file(nodes, faces, f"{filename}_scaled-{amount}")


if __name__ == "__main__":
    main()
