from graph import Graph
from faces import Face
import polyhedra_generation
from component import Component
from component_node import ComponentNode
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
        self.layers = self.get_layers()

        all_faces = self.dual_graph.get_V()
        component_graph = self.dual_graph.copy()
        for i in xrange(len(self.layers) - 1):
            faces_between = [key for key in all_faces if all_faces[key].between_layers(self.layers[i + 1], self.layers[i])]
            subgraph_between = self.dual_graph.subgraph(faces_between)
            for face in faces_between:
                if face not in component_graph.get_V():
                    continue
                connections = subgraph_between.get_reachable(face)
                component_dual_graph = subgraph_between.subgraph(connections)
                component = Component(component_dual_graph, self.layers[i + 1], self.layers[i])
                component_graph = component_graph.combine_vertices(connections, component)
        self.component_graph = component_graph
        self.unfolding_tree = self.create_unfolding_tree()

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
    
    def create_unfolding_tree(self):
      components_dict = self.component_graph.get_V()
      root = components_dict.keys()[0]
      root_component = self.component_graph.get_vertex(root)
      root_node = ComponentNode(root_component)
      remaining = [c for c in components_dict if c != root]
      root_node, _ = self.create_unfolding_subtree(root_node, remaining)
      return root_node
      
    def create_unfolding_subtree(self, root_node, remaining_components):
      remaining = remaining_components
      for c in remaining_components:
        if not remaining:
          break
        if c not in remaining:
          continue
        remaining.remove(c)
        component = self.component_graph.get_vertex(c)
        bridge = self.get_bridge(root_node.component, component)
        if not bridge:
          continue
        c_node = ComponentNode(component)
        child, remaining = self.create_unfolding_subtree(c_node, remaining)
        child.add_parent_bridge(bridge[1])
        root_node.add_child(child)
        root_node.add_child_bridge(bridge[0])
      return root_node, remaining

    # returns the bridge connecting component c1 and component c2
    # returned as a list of length 2 of lists where the first sublist is the portion of the bridge in c1 as a sequence/path of faces and the second sublist is the portion not in c1
    def get_bridge(self, c1, c2):
      if c1.y == c2.y_minus_1:
        y = c1.y
      else:
        y = c1.y_minus_1

      c1_faces = c1.get_faces()
      c2_faces = c2.get_faces()
      c1_z = [face for face in c1_faces if c1.get_face(face).direction == '+z' or c1.get_face(face).direction == '-z']
      c2_z = [face for face in c2_faces if c2.get_face(face).direction == '+z' or c2.get_face(face).direction == '-z']
      all_faces_dict = self.dual_graph.get_V()
      y_faces = [face for face in all_faces_dict if all_faces_dict[face].in_layer(y)]
      face_subgraph = self.dual_graph.subgraph(c1_z + c2_z + y_faces)
      c2_z_set = set(c2_z)
      
      for face in c1_z:
        layers = [[face]]
        all_faces = [face]
        while layers[-1]:
          next_layer = []
          for vertex in layers[-1]:
            connections = face_subgraph.get_connections(vertex)
            connections = [c for c in connections if c not in all_faces]
            all_faces.extend(connections)
            if not c2_z_set.isdisjoint(connections):
              for c in connections:
                if c in c2_z_set:
                  next_layer = [c]
                  break
              break
            next_layer.extend([c for c in connections if face_subgraph.get_V()[c].in_layer(y)])
          layers.append(next_layer)
          if next_layer and next_layer[0] in c2_z_set:
            break
        if layers[-1] and layers[-1][0] in c2_z_set:
          break
      
      if not layers[-1]:
        return []

      path = []
      path.append(layers[-1][0])
      prev = layers[-1][0]
      for i in reversed(xrange(len(layers) - 1)):
        connections = face_subgraph.get_connections(prev)
        for face in layers[i]:
          if face in connections:
            path.insert(0, face)
            prev = face
            break
      c1_bridge = []
      c2_bridge = []
      for face in path:
        if face in c1_faces:
          c1_bridge.append(face)
        else:
          c2_bridge.append(face)
      return [c1_bridge, c2_bridge]

if __name__ == "__main__":
    p = Polyhedron(filelist=["../data/test/unit_cube_open.fold", "../data/test/rect_box.fold"])
    #c1 = p.component_graph.get_V()[0]
    #c2 = p.component_graph.get_V()[5]
    #c3 = Component(p.dual_graph.subgraph([6]), 2, 1)
    #print p.get_bridge(c1, c3)

