from graph import Graph
from faces import Face
import polyhedra_generation
from component import Component
import json
import numpy as np


class Polyhedron(object):

    def __init__(self, vertices=None, faces=None, filelist=None):
        """
        initialize a polyhedron either with a list of vertices and faces OR
        a list of FOLD files where each file is a separate component of the polyhedron
        faces - list of Face objects
        faces_vertices - list of indices into self.vertices defining each face
        vertices - list of coordinates (np array) of each vertex
        components - list of lists indices into self.faces where each list of indices defines a component
        """

        if vertices is not None and faces is not None:
            self.faces = faces
            self.vertices = vertices
        elif filelist is not None:
            self.components = []
            total_faces = 0
            for f in range(len(filelist)):
                with open(filelist[f]) as f_o:
                    data = json.load(f_o)
                    faces, faces_vertices, vertices = self.parse_fold_file(filelist[f])
                    self.components.append(range(total_faces, total_faces+len(faces)))
                    total_faces += len(faces)
                if f == 0:
                    polyhedra_generation.create_fold_file("tmp.fold", data)
                else:
                    polyhedra_generation.create_fold_file("tmp.fold", data, append=True)
            self.faces, self.faces_vertices, self.vertices = self.parse_fold_file("tmp.fold")
        else:
            raise Exception("must pass in either both vertices and faces or filename(s) to constructor")

        self.primal_graph = self.create_primal_graph()
        self.dual_graph = self.create_dual_graph()
        for c in self.components:
            print "Component:"
            for f in [self.faces[i] for i in c]:
                print f
        self.layers = self.get_layers()

        # all_faces = self.dual_graph.get_V()
        # component_graph = self.dual_graph.copy()
        # for i in xrange(len(self.layers) - 1):
        #   faces_between = [key for key in all_faces if all_faces[key].between_layers(self.layers[i + 1], self.layers[i])]
        #   subgraph_between = self.dual_graph.subgraph(faces_between)
        #   for face in faces_between:
        #     if face not in component_graph.get_V():
        #       continue
        #     connections = subgraph_between.get_reachable(face)
        #     component_dual_graph = subgraph_between.subgraph(connections)
        #     component = Component(component_dual_graph, self.layers[i + 1], self.layers[i])
        #     component_graph = component_graph.combine_vertices(connections, component)
        # self.components = component_graph


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
            direction = direction.astype(np.float)
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

        return Graph(self.vertices, E_list=edges)

    def create_dual_graph(self):
        """
        creates a dual graph from a list of faces. Nodes are faces and edges exist between nodes if the two faces
        are adjacent on the polyhedron (share an edge)
        :return: Graph object
        """
        edges = []
        for u in range(len(self.faces)):
            for v in range(u, len(self.faces)):
                # check to see how many vertices are shared between face u and face v, 6 unique vertices means adjacent
                if len(set(self.faces[u].get_vertices(as_tuple=True) + self.faces[v].get_vertices(as_tuple=True))) == 6:
                    edges.append([u, v])

        return Graph(self.faces, E_list=edges)

    def write_dual_graph(self, filename):
        """
        write dual graph to FOLD format in linakge form for visualization
        :param filename: string
        :return: None
        """
        data_out = {"vertices_coords": [],
                "edges_vertices": []}
        for u in self.dual_graph.E:
            center = self.dual_graph.get_V()[u].get_center()
            data_out["vertices_coords"].append([x for x in center])
            for v in self.dual_graph.E[u]:
                data_out["edges_vertices"].append([u, v])
        polyhedra_generation.create_fold_file(filename, data_out, frame_class="linkage")
    
    # returns a list of layers as a list of y-values
    def get_layers(self):
        y_values = [vertex[1] for vertex in self.vertices]
        return sorted(list(set(y_values)))

if __name__ == "__main__":
    p = Polyhedron(filelist=["../data/test/unit_cube_open.fold", "../data/test/rect_box.fold"])
