from graph import Graph
import numpy as np
import polyhedra_generation

class Component:
  
  # face_graph is the graph with faces as vertices of the component
  def __init__(self, face_graph, y, y_minus_1, is_leaf=False):
    self.is_leaf = is_leaf
    self.y_minus_1 = y_minus_1
    self.full_graph = face_graph
    faces = self.full_graph.get_V()

    prot_faces = [key for key in faces if not faces[key].in_layer(y) and not faces[key].in_layer(y_minus_1)]
    self.prot = self.full_graph.subgraph(prot_faces)
    # TODO: make self.prot a tree and choose the root using self.unfolding_path

    self.f_0 = 11  #TODO: find f_0 (bridge z+ direction)

    back_rim_faces = [key for key in faces if faces[key].in_layer(y)]
    self.back_rim = self.full_graph.subgraph(back_rim_faces)

    front_rim_faces = [key for key in faces if faces[key].in_layer(y_minus_1)]
    self.front_rim = self.full_graph.subgraph(front_rim_faces)

    self.depth = abs(y - y_minus_1)

    # if is_leaf:
    try:
      # self.connector_path is list of indices into full graph, where f_i is first element and f_j is last
      self.f_i = 9
      self.protrusion_path = self.unfolding_path()[:-1]
      self.connector_path = self.find_connector()
      self.f_j = self.connector_path[-1]
      self.is_leaf = True
    except Exception as e:
      print "error:",e
    if self.is_leaf:
      self.strip_paths_leaf()

  def unfolding_path(self):
    """
    create a tree that navigates the protrusion ccw
    :param root: root face (should be +z or -z direction)
    :param is_leaf: bool
    :return: None
    """
    root = self.f_0
    p = root

    vertices = self.prot.get_V()
    edges = {v: [] for v in vertices}
    visited = {}
    path = [root]
    back_at_root = False
    while not back_at_root:
      children = self.prot.E[p]
      for c in children:
        p_face = self.prot.get_V()[p]
        c_face = self.prot.get_V()[c]
        if p_face.direction == c_face.direction:
          if ((p_face.direction == "+x" and p_face.center[2] > c_face.center[2]) or
                (p_face.direction == "-z" and p_face.center[0] > c_face.center[0]) or
                (p_face.direction == "-x" and p_face.center[2] < c_face.center[2]) or
                (p_face.direction == "+z" and p_face.center[0] < c_face.center[0])):
            if c == root:
              back_at_root = True
            edges[p].append(c)
            path.append(c)
            p = c
            visited[c] = 1
            break
        elif ((p_face.direction == "+x" and c_face.direction == "-z") or
                (p_face.direction == "-z" and c_face.direction == "-x") or
                (p_face.direction == "-x" and c_face.direction == "+z") or
                (p_face.direction == "+z" and c_face.direction == "+x")):
          if c == root:
            back_at_root = True

          edges[p].append(c)
          path.append(c)
          p = c
          visited[c] = 1
          break
    return path
    # return Graph(vertices, E_dict=edges)


  def find_connector(self):
    """
    find the connector of the component if starting from f_i
    :param f_i: index of starting face
    :return: list of Faces
    """
    back_rim = set(self.back_rim.get_V())
    face_path = [self.f_i]
    next_face = back_rim.intersection(set(self.full_graph.get_connections(face_path[-1])))

    while next_face:
      face_path.append(next_face.pop())
      next_face = back_rim.intersection(set(self.full_graph.get_connections(face_path[-1])))
      next_face = {f for f in next_face if self.full_graph.get_V()[f].center[2] < self.full_graph.get_V()[face_path[-1]].center[2]}

    next_face = self.full_graph.get_connections(face_path[-1])
    for f in next_face:
      if self.full_graph.get_V()[f].direction == "-z":
        face_path.append(f)
        break
    return face_path

  def strip_paths_leaf(self):
    # cut_1 is the left-most cut

    cut_1 = []
    cut_2 = []
    cut_3 = []

    f_0_vertices = self.full_graph.get_V()[self.f_0].vertices
    # first 3 points of s and t, respectively
    cut_1.append(f_0_vertices[0])
    cut_2.append((f_0_vertices[0]+f_0_vertices[1])/2.0)

    # go along protrusion until hitting f_j
    # FIRST TURN
    path_index = 0
    cur_face_index = self.protrusion_path[path_index]
    while cur_face_index != self.f_j:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + self.depth/2.0, cur_face.vertices[0][2]])
      cut_1.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index]

    # 1st side of f_j
    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + self.depth/2.0, cur_face.vertices[0][2]])
    cut_1.append(next_p)
    #
    path_index += 1
    cur_face_index  = self.protrusion_path[path_index]
    # # 2nd side of f_j
    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + self.depth/2.0, cur_face.vertices[0][2]])
    cut_1.append(next_p)

    # SECOND/THIRD TURN
    cur_face_index = self.protrusion_path[path_index]
    while cur_face_index != self.f_i:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[3][0], cur_face.vertices[3][1], cur_face.vertices[3][2]])
      cut_1.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]

    # FOURTH TURN - go along connector
    for f in self.connector_path[1:]:
      cur_face = self.full_graph.get_V()[f]
      next_p = np.array([cur_face.vertices[2][0], cur_face.vertices[2][1], cur_face.vertices[2][2]])
      cut_1.append(next_p)


    # second cut
    cut_2.append(np.array([cut_2[0][0], cut_2[0][1] + .25*self.depth, cut_2[0][2]]))

    # FIRST TURN, go until hitting f_k (face before f_0)
    path_index = 1
    cur_face_index = self.protrusion_path[path_index]
    while cur_face_index != self.f_0:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + .25*self.depth, cur_face.vertices[0][2]])
      cut_2.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]

    # SECOND TURN
    next_p = np.array([(next_p[0]+f_0_vertices[0][0])/2.0, next_p[1], (next_p[2]+f_0_vertices[0][2])/2.0])
    cut_2.append(next_p)

    next_p = np.array([next_p[0], next_p[1]+self.depth/2.0, next_p[2]])
    cut_2.append(next_p)

    # THIRD TURN - go from f_0 to f_i
    path_index = 0
    cur_face_index = self.protrusion_path[path_index]
    while cur_face_index != self.f_i:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[0][0], next_p[1], cur_face.vertices[0][2]])
      cut_2.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]

    # cross f_i
    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[0][0], next_p[1], cur_face.vertices[0][2]])
    cut_2.append(next_p)
    path_index += 1

    cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[0][0], next_p[1], cur_face.vertices[0][2]])
    cut_2.append(next_p)

    # FORTH TURN - go down along connector
    for f in self.connector_path[1:]:
      cur_face = self.full_graph.get_V()[f]
      next_p = np.array([cur_face.vertices[3][0], cur_face.vertices[3][1], cur_face.vertices[3][2]])
      cut_2.append(next_p)

    # FIFTH TURN - back to f_i
    # TODO: this is an edge of the connector when the connector touches the protrusion on 3 sides. Not implemented for
    # general case
    path_index = 1
    cur_face_index = self.protrusion_path[path_index]
    while cur_face_index != self.f_0:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[0][0], self.y_minus_1, cur_face.vertices[0][2]])
      cut_3.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]

    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[0][0], self.y_minus_1, cur_face.vertices[0][2]])
    cut_3.append(next_p)
    # cut 3 - go around entire protrusion starting from f_0


    # self.write_cut_path([cut_1, cut_2])
    self.write_strip(cut_3, [])

  def write_cut_path(self, paths):
    vertices = []
    edges = []
    faces = []
    for path in paths:
      for i in range(len(path)-1):
        u, v = path[i].astype(float), path[i+1].astype(float)
        # u1 = u.copy()
        u1 = np.array([u[0]+.01, u[1]+.2, u[2]+.01])
        # v1 = v.copy()
        v1 = np.array([v[0]+.01, v[1]+.2, v[2]+.01])
        faces.append([len(vertices), len(vertices)+1, len(vertices)+2, len(vertices)+3])
        edges.append([len(vertices), len(vertices)+1])
        vertices.extend([list(x) for x in [u, v, v1, u1]])
    polyhedra_generation.create_fold_file("cuts.fold", {"vertices_coords": vertices, "faces_vertices":faces,
                                                          "edges_vertices": edges})
  def write_strip(self, path1, path2):
    vertices = [list(x) for x in path1]
    vertices.extend([list(x) for x in path2])
    faces = []

    for i in range(len(path1) - 1):
      faces.append([i, i+1, len(path1)+i+1, len(path1)+i])
    polyhedra_generation.create_fold_file("path.fold", {"vertices_coords": vertices, "faces_vertices":faces,
                                                          "edges_vertices": []}, append=True)



  def __str__(self):
    return self.full_graph.__str__()
