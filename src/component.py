from graph import Graph
import numpy as np
import polyhedra_generation

class Component:
  
  # face_graph is the graph with faces as vertices of the component
  def __init__(self, face_graph, y, y_minus_1, is_root=False):
    self.num_children = 1
    self.is_root = is_root
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
      # TODO: define the children bridge faces
      self.num_children = 1
      self.children = []
    except Exception as e:
      print "error:",e
    if self.num_children == 0:
      self.strip_paths_leaf()
    elif self.is_root:
      pass  #TODO implement unfolding of root
    else:
      self.unfold_strip_intermediate()

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
    last_face = None  # to avoid adding unnecessary vertices to cut if staying on same plane
    while cur_face_index != self.f_j:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + self.depth/2.0, cur_face.vertices[0][2]])
      if last_face is None or cur_face.direction != last_face.direction:
        cut_1.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index]
      last_face = cur_face

    # first side of f_j
    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + self.depth/2.0, cur_face.vertices[0][2]])
    cut_1.append(next_p)

    path_index += 1
    cur_face_index  = self.protrusion_path[path_index]
    # # 2nd side of f_j
    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + self.depth/2.0, cur_face.vertices[0][2]])
    cut_1.append(next_p)

    # SECOND/THIRD TURN
    cur_face_index = self.protrusion_path[path_index]
    last_face = None
    while cur_face_index != self.f_i:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[3][0], cur_face.vertices[3][1], cur_face.vertices[3][2]])
      if last_face is None or cur_face.direction != last_face.direction:
        cut_1.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
      last_face = cur_face

    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[3][0], cur_face.vertices[3][1], cur_face.vertices[3][2]])
    cut_1.append(next_p)

    # FOURTH TURN - go along connector
    f = self.connector_path[-1]
    cur_face = self.full_graph.get_V()[f]
    next_p = np.array([cur_face.vertices[2][0], cur_face.vertices[2][1], cur_face.vertices[2][2]])
    cut_1.append(next_p)

    # second cut
    cut_2.append(np.array([cut_2[0][0], cut_2[0][1] + .25*self.depth, cut_2[0][2]]))

    # FIRST TURN, go until hitting f_k (face before f_0)
    path_index = 1
    cur_face_index = self.protrusion_path[path_index]
    last_face = self.full_graph.get_V()[cur_face_index]
    while cur_face_index != self.f_0:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + .25*self.depth, cur_face.vertices[0][2]])
      if last_face is None or cur_face.direction != last_face.direction:
        cut_2.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
      last_face = cur_face

    # SECOND TURN
    next_p = np.array([(next_p[0]+f_0_vertices[0][0])/2.0, next_p[1], (next_p[2]+f_0_vertices[0][2])/2.0])
    cut_2.append(next_p)

    next_p = np.array([next_p[0], next_p[1]+self.depth/2.0, next_p[2]])
    cut_2.append(next_p)

    # THIRD TURN - go from f_0 to f_i
    path_index = 0
    cur_face_index = self.protrusion_path[path_index]
    last_face = self.full_graph.get_V()[cur_face_index]
    while cur_face_index != self.f_i:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[0][0], next_p[1], cur_face.vertices[0][2]])
      if last_face is None or cur_face.direction != last_face.direction:
        cut_2.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
      last_face = cur_face

    # cross f_i
    path_index += 1
    cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[0][0], next_p[1], cur_face.vertices[0][2]])
    cut_2.append(next_p)

    # FOURTH TURN - go down along connector
    f = self.connector_path[1]
    cur_face = self.full_graph.get_V()[f]
    next_p = np.array([cur_face.vertices[3][0], cur_face.vertices[3][1], cur_face.vertices[3][2]])
    cut_2.append(next_p)


    f = self.connector_path[-1]
    cur_face = self.full_graph.get_V()[f]
    next_p = np.array([cur_face.vertices[3][0], cur_face.vertices[3][1], cur_face.vertices[3][2]])
    cut_2.append(next_p)

    # FIFTH TURN - back to f_i from f_j
    # TODO: this is an edge of the connector when the connector touches the protrusion on 3 sides. Not implemented for
    # general case
    # find index of f_j in prot_path
    general_case = False
    if general_case:
      path_index = [i for i in range(len(self.protrusion_path)) if self.protrusion_path[i] == self.f_j].pop()
      cur_face_index = self.protrusion_path[path_index]
      last_face = None
      while cur_face_index != self.f_i:
        cur_face = self.full_graph.get_V()[cur_face_index]
        next_p = np.array([cur_face.vertices[0][0], next_p[1], cur_face.vertices[0][2]])
        if last_face is None or cur_face.direction != last_face.direction:
          cut_2.append(next_p)
        path_index -= 1
        cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
        last_face = cur_face

    # cut 3 - go around entire protrusion starting from f_0
    path_index = 1
    cur_face_index = self.protrusion_path[path_index]
    last_face = None
    while cur_face_index != self.f_0:
      cur_face = self.full_graph.get_V()[cur_face_index]
      next_p = np.array([cur_face.vertices[0][0], self.y_minus_1, cur_face.vertices[0][2]])
      if last_face is None or last_face.direction != cur_face.direction:
        cut_3.append(next_p)
      path_index += 1
      cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
      last_face = cur_face

    cur_face = self.full_graph.get_V()[cur_face_index]
    next_p = np.array([cur_face.vertices[0][0], self.y_minus_1, cur_face.vertices[0][2]])
    cut_3.append(next_p)

    # print cut_2
    # self.write_cut_path([cut_1, cut_2, cut_3])
    # self.write_strip(cut_3, [])

  def unfold_strip_intermediate(self):
    B1 = []
    B2 = []
    all_cuts = []
    try:
      f_0_face = self.full_graph.get_V()[self.f_0]
    except:
      return

    for child in self.children:
      bridge_index = child.bridge_face() # TODO: define this - should be bridge between this component and its children
      bridge_face = self.full_graph.get_V()[bridge_index]
      # is child bridge in B1 or B2?
      if bridge_face.in_layer(self.y_minus_1):
        B1.append(child)
      else:
        B2.append(child)

    B1 = [13, 12]
    B2 = [14]
    layer_width = float(self.depth) / (len(B1) + 2*len(B2))
    f_0_strip_width = abs(f_0_face.vertices[0][0] - f_0_face.vertices[0][1])/float(len(B1) + len(B2))

    # sort B1, B2 by distance from f_k, respective to ccw direction around protrusion
    # B1 = sorted(B1, key=lambda c: [i for i in range(len(self.protrusion_path)) if self.protrusion_path[i] == c].pop())
    # B2 = sorted(B2, key=lambda c: [i for i in range(len(self.protrusion_path)) if self.protrusion_path[i] == c].pop())
    for i in range(len(B1)):
      child = B1[i]
      # bridge_index = child.bridge_face() # TODO: define this - should be bridge between this component and its children
      bridge_index = B1[i]
      bridge_face = self.full_graph.get_V()[bridge_index]

      width = abs(bridge_face.vertices[0][0] - bridge_face.vertices[1][0])
      # num_cuts = child.num_children * 2 + 1  # number of strips in child + 1
      num_cuts = 5
      cut_width = float(width) / (num_cuts - 1)
      cut_paths = [[] for _ in range(num_cuts)]

      strip_depth = layer_width * i  # how far the closest cut goes into the protrusion
      for cut_num in range(num_cuts):
        # first and second point
        cut_depth = strip_depth + (layer_width/(num_cuts-1))*cut_num
        if not (i == 0 and cut_num == 0):
          print bridge_face.direction
          if bridge_face.direction == "+z":
            cut_point = np.array([bridge_face.vertices[0][0] + cut_width * cut_num,
                                  bridge_face.vertices[0][1], bridge_face.vertices[0][2]])
          elif bridge_face.direction == "-z":
            cut_point = np.array([bridge_face.vertices[0][0] - cut_width * cut_num,
                                  bridge_face.vertices[0][1], bridge_face.vertices[0][2]])
          cut_paths[cut_num].append(cut_point)
          cut_paths[cut_num].append(np.array([cut_point[0], cut_point[1] + cut_depth, cut_point[2]]))

        # second point - then turn cw to f_0

        starting_face_index = [j for j in range(len(self.protrusion_path)) if self.protrusion_path[j] == bridge_index].pop()

        path_index = starting_face_index
        cur_face_index = self.protrusion_path[path_index]
        # last_face = None
        while cur_face_index != self.f_0:
          cur_face = self.full_graph.get_V()[cur_face_index]
          next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + cut_depth, cur_face.vertices[0][2]])
          # if last_face is None or last_face.direction != cur_face.direction:
          cut_paths[cut_num].append(next_p)
          path_index -= 1
          cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
          # last_face = cur_face

        # go to middle of f_0, then down to parent bridge
        if not (i == 0 and cut_num == 0):
          if f_0_face.direction == "+z":
            x = f_0_face.vertices[1][0] - (f_0_strip_width * i) - f_0_strip_width/(num_cuts-1) * cut_num
          else:
            x = f_0_face.vertices[1][0] + (f_0_strip_width * i) + f_0_strip_width/(num_cuts-1) * cut_num
          y = next_p[1]
          z = next_p[2]

          cut_paths[cut_num].append(np.array([x, y, z]))
          y = f_0_face.vertices[0][1]
          cut_paths[cut_num].append(np.array([x, y, z]))

      all_cuts.extend(cut_paths)

    # strips in B2
    for i in range(len(B2)):
      child = B2[i]
      # bridge_index = child.bridge_face() # TODO: define this - should be bridge between this component and its children
      bridge_index = B2[i]
      bridge_face = self.full_graph.get_V()[bridge_index]

      width = abs(bridge_face.vertices[0][0] - bridge_face.vertices[1][0])
      # num_cuts = child.num_children * 2 + 1  # number of strips in child + 1
      num_cuts = 3  # TODO: get number of cuts for children
      cut_width = float(width) / (num_cuts - 1)
      cut_paths = [[] for _ in range(num_cuts)]

      strip_depth = layer_width * i # how far the closest cut goes into the protrusion
      for cut_num in range(num_cuts):
        cut_depth = strip_depth + (layer_width/(num_cuts-1))*cut_num
        cut_depth2 = layer_width + cut_depth
        if not (i == 0 and cut_num == 0):
          if bridge_face.direction == "+z":
            cut_point = np.array([bridge_face.vertices[3][0] + cut_width * cut_num,
                                  bridge_face.vertices[3][1], bridge_face.vertices[3][2]])
          elif bridge_face.direction == "-z":
            cut_point = np.array([bridge_face.vertices[3][0] - cut_width * cut_num,
                                  bridge_face.vertices[3][1], bridge_face.vertices[3][2]])
          cut_paths[cut_num].append(cut_point)
          cut_paths[cut_num].append(np.array([cut_point[0], cut_point[1] - cut_depth, cut_point[2]]))

        # continue to f_k
        f_k = self.protrusion_path[-1]
        f_k_face = self.full_graph.get_V()[f_k]
        if f_k_face.direction == "+z":
          f_k_strip_width = abs(f_k_face.vertices[0][0] - f_k_face.vertices[1][0]) / float(len(B2))
        else:  # -x direction
          f_k_strip_width = abs(f_k_face.vertices[0][2] - f_k_face.vertices[1][2]) / float(len(B2))
        starting_face_index = [j for j in range(len(self.protrusion_path)) if self.protrusion_path[j] == bridge_index].pop()

        path_index = starting_face_index
        cur_face_index = self.protrusion_path[path_index]
        # last_face = None
        while cur_face_index != f_k:
          cur_face = self.full_graph.get_V()[cur_face_index]
          next_p = np.array([cur_face.vertices[3][0], cur_face.vertices[3][1] - cut_depth, cur_face.vertices[3][2]])
          # if last_face is None or last_face.direction != cur_face.direction:
          cut_paths[cut_num].append(next_p)
          path_index -= 1
          cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
          # last_face = cur_face

        # go down f_k
        if f_k_face.direction == "+z":
          x = f_k_face.vertices[3][0] + (f_k_strip_width * i) + f_k_strip_width/(num_cuts - 1) * cut_num
        else:  # -z direction
          x = f_k_face.vertices[3][0] - (f_k_strip_width * i) - f_k_strip_width/(num_cuts - 1) * cut_num
        y = next_p[1]
        z = next_p[2]

        cut_paths[cut_num].append(np.array([x, y, z]))
        y = f_k_face.vertices[3][1] - cut_depth2  # move down into protrusion along f_k
        if (cut_num != 0):
          cut_paths[cut_num].append(np.array([x, y, z]))

        cur_face_index = self.protrusion_path[-1]
        # last_face = None
        while cur_face_index != self.f_0:
          cur_face = self.full_graph.get_V()[cur_face_index]
          next_p = np.array([cur_face.vertices[3][0], y, cur_face.vertices[3][2]])
          # if last_face is None or last_face.direction != cur_face.direction:
          cut_paths[cut_num].append(next_p)
          path_index -= 1
          cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
          # last_face = cur_face

        # go down to parent bridge along f_0
        if f_0_face.direction == "+z":
          x = f_0_face.vertices[0][0] + (f_0_strip_width * i) + f_0_strip_width/(num_cuts-1) * cut_num
        else:
          x = f_0_face.vertices[0][0] - (f_0_strip_width * i) - f_0_strip_width/(num_cuts-1) * cut_num
        z = next_p[2]

        cut_paths[cut_num].append(np.array([x, y, z]))
        y = f_0_face.vertices[0][1]
        cut_paths[cut_num].append(np.array([x, y, z]))
        all_cuts.extend(cut_paths)

    # self.write_cut_path_as_fold(all_cuts)
    self.write_cut_path(all_cuts, "../out/cuts.txt")
    # for p in all_cuts:
    #   print p

  def write_cut_path(self, paths, filename):
    n_paths = len(paths)
    with open(filename, 'w') as f:
      f.write("%s\n" % n_paths)
      for path in paths:
        f.write("%s\n" % len(path))
        for point in path:
          f.write("%s %s %s\n" % (point[0], point[1], point[2]))

  def write_cut_path_as_fold(self, paths):
    vertices = []
    edges = []
    faces = []
    for path in paths:
      last_u, last_v = None, None
      for i in range(len(path)-1):
        u, v = path[i].astype(float), path[i+1].astype(float)

        u1 = np.array([u[0]+.01, u[1]+.02, u[2]+.01])
        v1 = np.array([v[0]+.01, v[1]+.02, v[2]+.01])

        faces.append([len(vertices), len(vertices)+1, len(vertices)+2, len(vertices)+3])
        edges.append([len(vertices), len(vertices)+1])
        vertices.extend([list(x) for x in [u, v, v1, u1]])
        last_u, last_v = u, v
    polyhedra_generation.create_fold_file("../out/cuts.fold", {"vertices_coords": vertices, "faces_vertices":faces,
                                                          "edges_vertices": edges})
  def write_strip(self, path1, path2):
    vertices = [list(x) for x in path1]
    vertices.extend([list(x) for x in path2])
    faces = []

    for i in range(len(path1) - 1):
      faces.append([i, i+1, len(path1)+i+1, len(path1)+i])
    polyhedra_generation.create_fold_file("../out/path.fold", {"vertices_coords": vertices, "faces_vertices":faces,
                                                          "edges_vertices": []}, append=True)



  def __str__(self):
    return self.full_graph.__str__()
