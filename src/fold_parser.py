import json
import numpy as np
from faces import Face
from graph import Graph

def parse_fold_file(filename):
    """
    parse the fold file and return a list of faces
    :param filename: file name (str) of the .fold file
    :return: list of faces
    """
    with open(filename) as f:
        data = json.load(f)

    vertices_coords = data["vertices_coords"]
    faces_vertices = data["faces_vertices"]
    faces = []
    for f in faces_vertices:
        vertices = [np.array(vertices_coords[i]) for i in f]
        # determine the normal direction of the face
        direction = np.cross(vertices[1] - vertices[0], vertices[3] - vertices[0])
        direction /= np.linalg.norm(direction)
        direction_str = None
        if abs(direction[0]) > .9:
            if direction[0] < 0:
                direction_str = "-x"
            else:
                direction_str = "+x"
        elif abs(direction[1]) > .9:
            if direction[1] < 0:
                direction_str = "-y"
            else:
                direction_str = "+y"
        elif abs(direction[2]) > .9:
            if direction[2] < 0:
                direction_str = "-z"
            else:
                direction_str = "+z"
        faces.append(Face(vertices, direction_str))
    return faces

def create_dual_graph(faces):
    """
    creates a dual graph from a list of faces. Nodes are faces and edges exist between nodes if the two faces
    are adjacent on the polyhedra
    :param faces: list of Face objects
    :return: Graph object
    """
    edges = []
    for u in range(len(faces)):
        adjacent = []
        for v in range(u, len(faces)):
            # check to see how many vertices are shared between face u and face v, 6 unique vertices means adjacent
            if len(set(faces[u].get_vertices(as_tuple=True) + faces[v].get_vertices(as_tuple=True))) == 6:
                adjacent.append(v)
        edges.append(adjacent)

    return Graph(faces, E=edges)


if __name__ == "__main__":
    faces = parse_fold_file("../data/unit_cube.fold")
    # for f in faces:
    #     print f
    print create_dual_graph(faces)
