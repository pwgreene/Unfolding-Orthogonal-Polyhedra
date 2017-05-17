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
                    faces, _, _, _ = self.parse_fold_file(filelist[f])
                    self.components.append(range(total_faces, total_faces+len(faces)))
                    total_faces += len(faces)
                if f == 0:
                    polyhedra_generation.create_fold_file("tmp.fold", data)
                else:
                    polyhedra_generation.create_fold_file("tmp.fold", data, append=True)
            self.faces, self.faces_vertices, self.vertices, self.edges = self.parse_fold_file("tmp.fold")
        else:
            raise Exception("must pass in either both vertices and faces or filename(s) to constructor")
        
        faces_vertices, edges_vertices, vertices_coords = self.grid_divide2(self.faces_vertices, self.edges, self.vertices)
        polyhedra_generation.create_fold_file("tmp2.fold", {'faces_vertices': faces_vertices, 'edges_vertices': edges_vertices, 'vertices_coords': vertices_coords})
        self.faces, self.faces_vertices, self.vertices, self.edges = self.parse_fold_file("tmp2.fold")

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
        edges_vertices = data["edges_vertices"]
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

        return faces, faces_vertices, vertices_coords, edges_vertices
    '''
    def grid_divide(self, faces, edges, vertices):
      new_vertices = [tuple(v) for v in vertices]
      new_faces = [tuple(f) for f in faces]
      new_edges = [tuple(e) for e in edges]
      
      e_f_map = {}
      f_e_map = {}
      for e in new_edges:
        e_f_map[e] = []
        for f in new_faces:
          if e[0] in f and e[1] in f:
            e_f_map[e].append(f)
            if f in f_e_map:
              f_e_map[f].append(e)
            else:
              f_e_map[f] = [e]
      cuts = [set(), set(), set()]
      for v in new_vertices:
        cuts[0].add(v[0])
        cuts[1].add(v[1])
        cuts[2].add(v[2])
      
      for i in cuts:
        for cut in cuts[i]:
          removal_edges = []
          adding_edges = []
          for e in new_edges:
            if new_vertices[e[0]][i] == new_vertices[e[1]][i]:
              continue
            v_1 = min(new_vertices[e[0]][i], new_vertices[e[1]][i])
            v_2 = max(new_vertices[e[0]][i], new_vertices[e[1]][i])
            if v_1 >= cut or v_2 <= cut:
              continue
            if e in removal_edges:
              continue
            new_vertex = list(new_vertices[e[0]])
            new_vertex[i] = cut
            new_vertex = tuple(new_vertex)
            new_vertex_idx = len(new_vertices)
            removal_edges.add(e)
            adding_edges.add((e[0], new_vertex_idx))
            adding_edges.add((e[1], new_vertex_idx))
            new_vertices.add(new_vertex)
            f_e_update = {}
            e_f_update = {}
            f_e_remove = []
            e_f_remove = [e]
            for f in e_f_map[e]:
              for f_e in f_e_map[f]:
                if e[0] not in f_e and e[1] not in f_e:
                  # TODO: check if endpoint has cut value
                  new_vertex2 = list(new_vertices[f_e[0]])
                  new_vertex2[i] = cut
                  new_vertex2 = tuple(new_vertex2)
                  new_vertex2_idx = len(new_vertices)
                  if new_vertex2 not in new_vertices:
                    new_vertices.add(new_vertex2)
                    removal_edges.add(f_e)
                    adding_edges.add((f_e[0], new_vertex2_idx))
                    adding_edges.add((f_e[1], new_vertex2_idx))
                    e_f_remove.append(f_e)
                  adding_edges.add((new_vertex_idx, new_vertex2_idx))
                  new_faces.remove(f) 
                  f_e_remove.append(f)
                  for j in xrange(len(f) - 1):
                    if f[j] in e and f[(j + 1) % len(f)] in e:
                      new_face1 = (f[j], new_vertex_idx, new_vertex2_idx, f[j - 1])
                      new_face2 = (new_vertex_idx, f[(j + 1) % len(f)], f[(j + 2) % len(f)], new_vertex2_idx)
                      for face_e in f_e_map[f]:
                        if f[j] in face_e and f[j - 1] in face_e:
                          last_edge1 = face_e
                        elif f[j + 1) % len(f)] in face_e and f[(j + 2) % len(f)] in face_e:
                          last_edge2 = face_e
                      f_e_update[new_face1] = [(f[j], new_vertex_idx), (new_vertex_idx, new_vertex2_idx), (f[j - 1], new_vertex2_idx), last_edge1]
                      f_e_update[new_face2] = [(f[(j + 1) % len(f)], new_vertex_idx), last_edge2, (f[(j + 2) % len(f)], new_vertex2_idx), (new_vertex_idx, new_vertex2_idx)]
                      alt_faces_e = list(e_f_map[e])
                      alt_faces_e.remove(f)
                      for alt_face in alt_faces_e:
                        connected_edges = list(f_e_map[alt_face])
                        connected_edges.remove(e)
                        connected_edges.extend([(f[j], new_vertex_idx), (f[(j + 1) % len(f)], new_vertex_idx)])
                        f_e_update[alt_face] = connected_edges
                      alt_faces_f_e = list(e_f_map[f_e])
                      alt_faces_f_e.remove(f)
                      for alt_face in alt_faces_f_e:
                        connected_edges = list(f_e_map[alt_face])
                        connected_edges.remove(f_e)
                        connected_edges.extend([(f[j - 1], new_vertex2_idx), (f[(j + 2) % len(f)], new_vertex2_idx)])
                        f_e_update[alt_face] = connected_edges
                      alt_faces_last_edge1 = list(e_f_map[last_edge1])
                      alt_face_last_edge1.remove(f)
                      alt_faces_last_edge2 = list(e_f_map[last_edge2])
                      alt_faces_last_edge2.remove(f)
                      e_f_update[last_edge1] = alt_faces_last_edge1 + [new_face1]
                      e_f_update[(f[j], new_vertex_idx)] = alt_faces_e + [new_face1]
                      e_f_update[(f[j - 1], new_vertex2_idx)] = alt_faces_f_e + [new_face1]
                      e_f_update[last_edge2] = alt_faces_last_edge2 + [new_face2]
                      e_f_update[(f[(j + 1) % len(f)], new_vertex_idx)] = alt_faces_e + [new_face2]
                      e_f_update[(f[(j + 2) % len(f)], new_vertex2_idx)] = alt_faces_f_e + [new_face2]
                      e_f_update[(new_vertex_idx, new_vertex2_idx)] = [new_face1, new_face2]
                      break
                    elif f[j] in f_e and f[(j + 1) % len(f)] in f_e:
                      new_face1 = (f[j], new_vertex2_idx, new_vertex_idx, f[j - 1])
                      new_face2 = (new_vertex2_idx, f[(j + 1) % len(f)], f[(j + 2) % len(f)], new_vertex_idx)
                      for face_e in f_e_map[f]:
                        if f[j] in face_e and f[j - 1] in face_e:
                          last_edge1 = face_e
                        elif f[j + 1) % len(f)] in face_e and f[(j + 2) % len(f)] in face_e:
                          last_edge2 = face_e
                      f_e_update[new_face1] = [(f[j], new_vertex2_idx), (new_vertex_idx, new_vertex2_idx), (f[j - 1], new_vertex_idx), last_edge1]
                      f_e_update[new_face2] = [(f[(j + 1) % len(f)], new_vertex2_idx), last_edge2, (f[(j + 2) % len(f)], new_vertex_idx), (new_vertex_idx, new_vertex2_idx)]
                      alt_faces_e = list(e_f_map[e])
                      alt_faces_e.remove(f)
                      for alt_face in alt_faces_e:
                        connected_edges = list(f_e_map[alt_face])
                        connected_edges.remove(e)
                        connected_edges.extend([(f[j - 1], new_vertex_idx), (f[(j + 2) % len(f)], new_vertex_idx)])
                        f_e_update[alt_face] = connected_edges
                      alt_faces_f_e = list(e_f_map[f_e])
                      alt_faces_f_e.remove(f)
                      for alt_face in alt_faces_f_e:
                        connected_edges = list(f_e_map[alt_face])
                        connected_edges.remove(f_e)
                        connected_edges.extend([(f[j], new_vertex2_idx), (f[(j + 1) % len(f)], new_vertex2_idx)])
                        f_e_update[alt_face] = connected_edges
                      alt_faces_last_edge1 = list(e_f_map[last_edge1])
                      alt_face_last_edge1.remove(f)
                      alt_faces_last_edge2 = list(e_f_map[last_edge2])
                      alt_faces_last_edge2.remove(f)
                      e_f_update[last_edge1] = alt_faces_last_edge1 + [new_face1]
                      e_f_update[(f[j], new_vertex2_idx)] = alt_faces_f_e + [new_face1]
                      e_f_update[(f[j - 1], new_vertex_idx)] = alt_faces_e + [new_face1]
                      e_f_update[last_edge2] = alt_faces_last_edge2 + [new_face2]
                      e_f_update[(f[(j + 1) % len(f)], new_vertex2_idx)] = alt_faces_f_e + [new_face2]
                      e_f_update[(f[(j + 2) % len(f)], new_vertex_idx)] = alt_faces_e + [new_face2]
                      e_f_update[(new_vertex_idx, new_vertex2_idx)] = [new_face1, new_face2]
                      break
                  new_faces.append(new_face1)
                  new_faces.append(new_face2)
                  break
            f_e_map.update(f_e_update)
            e_f_map.update(e_f_update)
            for remove_f in f_e_remove:
              del f_e_map[remove_f]
            for remove_e in e_f_remove:
              del e_f_map[remove_e]
          for e in removal_edge:
            new_edges.remove(e)
          new_edges.extend(adding_edges)
      return new_faces, new_edges, new_vertices
    '''

    def grid_divide2(self, faces, edges, vertices):
      all_vertices = vertices
      all_edges = edges
      all_faces = faces

      cuts = [set(), set(), set()]
      for v in all_vertices:
        cuts[0].add(v[0])
        cuts[1].add(v[1])
        cuts[2].add(v[2])
      
      for d in xrange(len(cuts)):
        for cut in cuts[d]:
          new_vertices = []
          new_edges = []
          new_faces = []
          cut_edges = []
          old_faces = []
          for e in all_edges:
            if all_vertices[e[0]][d] == all_vertices[e[1]][d]:
              continue
            v_1 = min(all_vertices[e[0]][d], all_vertices[e[1]][d])
            v_2 = max(all_vertices[e[0]][d], all_vertices[e[1]][d])
            if v_1 < cut and v_2 > cut:
              new_vertex = list(all_vertices[e[0]])
              new_vertex[d] = cut
              new_vertex_id = len(all_vertices) + len(new_vertices)
              new_vertices.append(new_vertex)
              new_edges.extend([[e[0], new_vertex_id], [e[1], new_vertex_id]])
              cut_edges.append(e)
          for i in xrange(len(new_vertices)):
            for j in xrange(len(new_vertices)):
              if i >= j:
                continue
              for f in all_faces:
                if self.e_in_f_FOLD(cut_edges[i], f) and self.e_in_f_FOLD(cut_edges[j], f):
                  new_edges.append([len(all_vertices) + i, len(all_vertices) + j])
                  old_faces.append(f)
                  new_face1, new_face2 = self.split_FOLD_face(f, len(all_vertices) + i, len(all_vertices) + j, cut_edges[i], cut_edges[j])
                  new_faces.extend([new_face1, new_face2])
                  # break
          all_vertices.extend(new_vertices)
          for ce in cut_edges:
            all_edges.remove(ce)
          all_edges.extend(new_edges)
          for of in old_faces:
            all_faces.remove(of)
          all_faces.extend(new_faces)
      return all_faces, all_edges, all_vertices

    def e_in_f_FOLD(self, e, f):
      if e[0] in f and e[1] in f:
        return True
      return False

    # v1 and v2 are vertex ids (not actual vertex)
    def split_FOLD_face(self, f, v1, v2, e1, e2):
      for i in xrange(len(f)):
        if f[i] in e1 and f[(i + 1) % len(f)] in e1:
          new_face1 = [f[i], v1, v2, f[i - 1]]
          new_face2 = [v1, f[(i + 1) % len(f)], f[(i + 2) % len(f)], v2]
          break
        elif f[i] in e2 and f[(i + 1) % len(f)] in e2:
          new_face1 = [f[i], v2, v1, f[i - 1]]
          new_face2 = [v2, f[(i + 1) % len(f)], f[(i + 2) % len(f)], v1]
          break
      return new_face1, new_face2
    
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

def write_to_off(self, out_filename):
    """
    output this polyhedron as an .off file
    :param out_filename: (str) filename of output .off file
    :return: None
    """
    nv, nt = len(self.vertices), len(self.faces)*2  # will split each face into 2 triangles
    tri = []  # list of triangles (indices into vertex array)
    for v1, v2, v3, v4 in self.faces_vertices:
      tri.append([v1, v2, v3])
      tri.append([v1, v3, v4])
    with open(out_filename, 'w') as out:
      out.write("OFF\n")
      out.write("%s %s 0\n" % (nv, nt))
      for v in self.vertices:
        out.write("%s %s %s\n" % (v[0], v[1], v[2]))
      for f in tri:
        out.write("3   %s %s %s\n" % (f[0], f[1], f[2]))
    print "wrote %s vertices and %s faces to %s" % (nv, nt, out_filename)

if __name__ == "__main__":
    p = Polyhedron(filelist=["../data/test/unit_cube_open.fold", "../data/test/rect_box.fold"])
    c = p.component_graph.get_V()[5]
    # c.unfold_strip_leaf(11, "-y")
    # c.unfold_strip_intermediate([13, 15, 14], ["+y", "+y", "+y"], 11, "+y", [4, 4, 2])
    c.unfold_strip_root(12, "-y", 2)
    p.write_to_off("../out/poly.off")
    
    #p = Polyhedron(filelist=["../data/boxes2.fold"])

