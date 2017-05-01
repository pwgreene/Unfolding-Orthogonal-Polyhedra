from graph import Graph
from faces import Face
import json
import numpy as np


class Polyhedron(object):

    def __init__(self, vertices=None, faces=None, filename=None):
        """
        initialize a polyhedron either with a list of vertices and faces, OR with a fold file filename
        faces - list of Face objects
        faces_vertices - list of indices into self.vertices defining each face
        vertices - list of coordinates (np array) of each vertex
        """
        if filename is not None:
            self.faces, self.faces_vertices, self.vertices = self.parse_fold_file(filename)
        elif vertices is not None and faces is not None:
            self.faces = faces
            self.vertices = vertices
        else:
            raise Exception("must pass in either both vertices and faces or just filename to constructor")
        self.primal_graph = self.create_primal_graph()
        self.dual_graph = self.create_dual_graph()
        # TODO: initialize Components
        self.components = None


    def parse_fold_file(self, filename):
        """
        parse the fold file and return a list of faces
        :param filename: file name (str) of the .fold file
        :return: None, but update self.faces and self.vertices
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

        return faces, faces_vertices, vertices_coords

    def create_primal_graph(self):
        """
        create primal graph, where nodes are vertices and edges are edges in the polyhedron
        :return: Graph object
        """
        edges = set()
        for face in self.faces_vertices:
            for i in range(len(face)):
                next_v = (i+1) % len(face)
                if (face[i], face[next_v]) not in edges and (face[next_v], face[i]) not in edges:
                    edges.add((face[i], face[next_v]))

        # convert set of tuples to list of lists
        edges = list(edges)
        edges = [list(e) for e in edges]
        print edges

        return Graph(self.vertices, E_long_list=edges)

    def create_dual_graph(self):
        """
        creates a dual graph from a list of faces. Nodes are faces and edges exist between nodes if the two faces
        are adjacent on the polyhedron (share an edge)
        :return: Graph object
        """
        edges = []
        for u in range(len(self.faces)):
            adjacent = []
            for v in range(u, len(self.faces)):
                # check to see how many vertices are shared between face u and face v, 6 unique vertices means adjacent
                if len(set(self.faces[u].get_vertices(as_tuple=True) + self.faces[v].get_vertices(as_tuple=True))) == 6:
                    adjacent.append(v)
            edges.append(adjacent)

        return Graph(self.faces, E=edges)

if __name__ == "__main__":
    p = Polyhedron(filename="../data/unit_cube.fold")
    print p.primal_graph