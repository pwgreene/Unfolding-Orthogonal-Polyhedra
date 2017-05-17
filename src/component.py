from graph import Graph
import numpy as np
import polyhedra_generation

class Component:
  
  # face_graph is the graph with faces as vertices of the component
  def __init__(self, face_graph, y, y_minus_1, is_root=False):
    self.num_children = 1
    self.is_root = is_root
    self.y_minus_1 = y_minus_1
    self.y = y
    self.full_graph = face_graph
    faces = self.full_graph.get_V()

    prot_faces = [key for key in faces if not faces[key].in_layer(y) and not faces[key].in_layer(y_minus_1)]
    self.prot = self.full_graph.subgraph(prot_faces)
    # TODO: make self.prot a tree and choose the root using self.unfolding_path
    # TODO: initialize the total number of strips in this component
    # self.num_strips = ?

    # self.f_0 = 11  #TODO: find f_0 (bridge z+ direction)

    back_rim_faces = [key for key in faces if faces[key].in_layer(y)]
    self.back_rim = self.full_graph.subgraph(back_rim_faces)

    front_rim_faces = [key for key in faces if faces[key].in_layer(y_minus_1)]
    self.front_rim = self.full_graph.subgraph(front_rim_faces)

    self.depth = abs(y - y_minus_1)

    '''
    try:
      # self.connector_path is list of indices into full graph, where f_i is first element and f_j is last
      self.f_i = 9
      self.protrusion_path = self.unfolding_path()[:-1]
      self.connector_path = self.find_connector()
      self.f_j = self.connector_path[-1]
      # TODO: define the children bridge faces
      self.num_children = 1
      self.children = []
      self.is_root = True
    except Exception as e:
      print "error:",e
    if self.num_children == 0:
      self.unfold_strip_leaf()
    elif self.is_root:
      #TODO implement unfolding of root
      self.unfold_strip_root()
    else:
      self.unfold_strip_intermediate()
    '''


  def compute_protrusion_path(self, f_0):
    """
    create a tree that navigates the protrusion ccw
    :return: None
    """
    root = f_0
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

    self.protrusion_path = path[:-1]
    # return Graph(vertices, E_dict=edges)


  def find_connector(self, f_i):
    """
    find the connector of the component if starting from f_i
    :param f_i: index of starting face
    :return: list of Faces
    """
    back_rim = set(self.back_rim.get_V())
    face_path = [f_i]
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

  def unfold_strip_leaf(self, f_0):
    # cut_1 is the left-most cut

    cut_1 = []
    cut_2 = []
    cut_3 = []

    self.compute_protrusion_path(f_0)

    f_0_vertices = self.full_graph.get_V()[f_0].vertices
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
    while cur_face_index != f_0:
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
    while cur_face_index != f_0:
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

    # TODO: uncomment this after defining self.parent face or equivalent
    if True:#not self.parent_face.in_layer(self.y_minus_1):
      # flip orientation
      new_cut1, new_cut2, new_cut3 = [], [], []
      for v in cut_1:
        v_new = np.array([v[0], self.depth - v[1], v[2]])
        new_cut1.append(v_new)
      for v in cut_2:
        v_new = np.array([v[0], self.depth - v[1], v[2]])
        new_cut2.append(v_new)
      for v in cut_3:
        v_new = np.array([v[0], self.depth - v[1], v[2]])
        new_cut3.append(v_new)
      cut_1, cut_2, cut_3 = new_cut1, new_cut2, new_cut3


    # print cut_2
    self.write_cut_path([cut_1, cut_2, cut_3], "../out/cuts.txt")
    # self.write_strip(cut_3, [])

  def unfold_strip_intermediate(self, child_faces, child_face_directions, f_0, parent_face_direction,
                                num_leaves_children):
    """
    unfolds the strip given locations of faces that connect to rest of tree. Directions are +y or -y
    :param child_faces: list of indices of face on bridges to children
    :param child_face_directions: direction of children (list of strings)
    :param f_0: face on bridge to parent (defined in paper)
    :param parent_face_direction: direction of bridge to parent (string)
    :param num_leaves_children: number of leaves of the tree rooted at each child (list of ints)
    :return: None
    """

    self.compute_protrusion_path(f_0)
    print self.protrusion_path

    B1 = []
    B2 = []
    all_cuts = []
    f_0_face = self.full_graph.get_V()[f_0]
    # map each face index of child to corresponding direction/number of leaves
    strips_on_face = dict([(child_faces[i], (child_face_directions[i], num_leaves_children[i]))
                         for i in range(len(child_faces))])

    for c in range(len(child_faces)):
      # is child bridge in B1 or B2?
      if child_face_directions[c] == "-y":
        B1.append(child_faces[c])
      else:
        B2.append(child_faces[c])

    # B1 = [13, 15]
    # B2 = [14]
    # B1_cuts = 5 # for testing
    # B2_cuts = 3
    assert (set(B1).isdisjoint(set(B2)))

    # sort B1, B2 by distance from f_k, respective to ccw direction around protrusion
    B1 = sorted(B1, key=lambda c: [i for i in range(len(self.protrusion_path)) if self.protrusion_path[i] == c].pop())
    B2 = sorted(B2, key=lambda c: [i for i in range(len(self.protrusion_path)) if self.protrusion_path[i] == c].pop())

    flip_sides = parent_face_direction == "+y"
    if flip_sides: # reverse orientation with respect to paper's description
      B1, B2 = B2[::-1], B1[::-1]
    cuts_so_far = 0  #keeps track of number of cuts going back to parent
    print B1, B2, flip_sides

    f_0_strip_width = abs(f_0_face.vertices[0][0] - f_0_face.vertices[0][1])/float(len(B1) + len(B2))

    # total_strips = (B1_cuts-1)*len(B1) + (B2_cuts-1)*len(B2)
    total_strips = sum(num_leaves_children)
    # strip_width_to_parent = abs(f_0_face.vertices[0][0] - f_0_face.vertices[0][1])/float(self.total_strips)

    strip_width_to_parent = abs(f_0_face.vertices[0][0] - f_0_face.vertices[0][1])/float(total_strips)
    print strip_width_to_parent, total_strips

    layer_width = float(self.depth) / (len(B1) + 2*len(B2))
    for i in range(len(B1)):
      # child = B1[i]
      # bridge_index = child.bridge_face() # TODO: define this - should be bridge between this component and its children
      bridge_index = B1[i]
      bridge_face = self.full_graph.get_V()[bridge_index]

      width = abs(bridge_face.vertices[0][0] - bridge_face.vertices[1][0])
      num_cuts = strips_on_face[bridge_index][1] + 1
      # num_cuts = B1_cuts
      cut_width = float(width) / (num_cuts - 1)
      cut_paths = [[] for _ in range(num_cuts)]

      side_index = 0

      strip_depth = layer_width * i  # how far the closest cut goes into the protrusion
      for cut_num in range(num_cuts):
        # first and second point
        # if flip_sides:
        #   cut_depth = strip_depth + (layer_width/(num_cuts-1))*(num_cuts-cut_num)
        # else:
        cut_depth = strip_depth + (layer_width/(num_cuts-1))*cut_num
        if not (i == 0 and cut_num == 0):
          print bridge_face.direction
          if bridge_face.direction == "+z":
            if flip_sides:
              cut_point = np.array([bridge_face.vertices[2][0] - cut_width * cut_num,
                                    bridge_face.vertices[2][1], bridge_face.vertices[2][2]])
            else:
              cut_point = np.array([bridge_face.vertices[0][0] + cut_width * cut_num,
                                    bridge_face.vertices[0][1], bridge_face.vertices[0][2]])
          elif bridge_face.direction == "-z":
            if flip_sides:
              cut_point = np.array([bridge_face.vertices[2][0] + cut_width * cut_num,
                                    bridge_face.vertices[2][1], bridge_face.vertices[2][2]])
            else:
              cut_point = np.array([bridge_face.vertices[0][0] - cut_width * cut_num,
                                    bridge_face.vertices[0][1], bridge_face.vertices[0][2]])
          cut_paths[cut_num].append(cut_point)
          if flip_sides:
            cut_paths[cut_num].append(np.array([cut_point[0], cut_point[1] - cut_depth, cut_point[2]]))
          else:
            cut_paths[cut_num].append(np.array([cut_point[0], cut_point[1] + cut_depth, cut_point[2]]))

        # second point - then turn to f_0

        starting_face_index = [j for j in range(len(self.protrusion_path)) if self.protrusion_path[j] == bridge_index].pop()

        path_index = starting_face_index
        cur_face_index = self.protrusion_path[path_index]
        # last_face = None
        while cur_face_index != f_0:
          cur_face = self.full_graph.get_V()[cur_face_index]
          if flip_sides:
            next_p = np.array([cur_face.vertices[2][0], cur_face.vertices[2][1] - cut_depth,
                              cur_face.vertices[2][2]])
          else:
            next_p = np.array([cur_face.vertices[0][0], cur_face.vertices[0][1] + cut_depth,
                              cur_face.vertices[0][2]])
          # if last_face is None or last_face.direction != cur_face.direction:
          cut_paths[cut_num].append(next_p)
          if flip_sides:
            path_index += 1
          else:
            path_index -= 1
          cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
          # last_face = cur_face

        # go to middle of f_0, then down to parent bridge
        if not (i == 0 and cut_num == 0):
          if f_0_face.direction == "+z":
            if flip_sides:
              x = f_0_face.vertices[3][0] + strip_width_to_parent*(cuts_so_far + cut_num)#(f_0_strip_width * i) + f_0_strip_width/(num_cuts-1) * cut_num
            else:
              x = f_0_face.vertices[1][0] - strip_width_to_parent*(cuts_so_far + cut_num)#(f_0_strip_width * i) - f_0_strip_width/(num_cuts-1) * cut_num
          else:
            if flip_sides:
              x = f_0_face.vertices[3][0] - strip_width_to_parent*(cuts_so_far + cut_num)#(f_0_strip_width * i) - f_0_strip_width/(num_cuts-1) * cut_num
            else:
              x = f_0_face.vertices[1][0] + strip_width_to_parent*(cuts_so_far + cut_num)#(f_0_strip_width * i) + f_0_strip_width/(num_cuts-1) * cut_num
          y = next_p[1]
          z = next_p[2]

          cut_paths[cut_num].append(np.array([x, y, z]))
          if flip_sides:
            y = f_0_face.vertices[3][1]
          else:
            y = f_0_face.vertices[0][1]
          cut_paths[cut_num].append(np.array([x, y, z]))

      all_cuts.extend(cut_paths)
      cuts_so_far += num_cuts-1

    # strips in B2
    cuts_so_far = 0
    for i in range(len(B2)):
      child = B2[i]

      bridge_index = B2[i]
      bridge_face = self.full_graph.get_V()[bridge_index]
      print bridge_face

      width = abs(bridge_face.vertices[0][0] - bridge_face.vertices[1][0])
      print bridge_index
      print "width:", width
      # num_cuts = child.num_children * 2 + 1  # number of strips in child + 1
      num_cuts = strips_on_face[bridge_index][1] + 1
      cut_width = float(width) / (num_cuts - 1)
      cut_paths = [[] for _ in range(num_cuts)]

      strip_depth = layer_width * i # how far the closest cut goes into the protrusion
      for cut_num in range(num_cuts):
        cut_depth = strip_depth + (layer_width/(num_cuts-1))*cut_num
        if flip_sides and len(B2) > 1:  # haven't figured out why, but this is necessary
          x = strip_depth + (layer_width/(num_cuts-1))*(cut_num+num_cuts-1)
          cut_depth2 = layer_width + x
        else:
          cut_depth2 = layer_width + cut_depth
        print cut_depth, cut_width
        if not (i == 0 and cut_num == 0):
          if bridge_face.direction == "+z":
            if flip_sides:
              cut_point = np.array([bridge_face.vertices[1][0] - cut_width * cut_num,
                                    bridge_face.vertices[1][1], bridge_face.vertices[1][2]])
            else:
              cut_point = np.array([bridge_face.vertices[3][0] + cut_width * cut_num,
                                    bridge_face.vertices[3][1], bridge_face.vertices[3][2]])
          elif bridge_face.direction == "-z":
            if flip_sides:
              cut_point = np.array([bridge_face.vertices[1][0] + cut_width * cut_num,
                                    bridge_face.vertices[1][1], bridge_face.vertices[1][2]])
            else:
              cut_point = np.array([bridge_face.vertices[3][0] - cut_width * cut_num,
                                    bridge_face.vertices[3][1], bridge_face.vertices[3][2]])
          cut_paths[cut_num].append(cut_point)
          if flip_sides:
            cut_paths[cut_num].append(np.array([cut_point[0], cut_point[1] + cut_depth, cut_point[2]]))
          else:
            cut_paths[cut_num].append(np.array([cut_point[0], cut_point[1] - cut_depth, cut_point[2]]))

        # continue to f_k
        if flip_sides:
          f_k = self.protrusion_path[1]
        else:
          f_k = self.protrusion_path[-1]
        print f_k, self.protrusion_path
        f_k_face = self.full_graph.get_V()[f_k]
        if f_k_face.direction == "+z":
          f_k_strip_width = abs(f_k_face.vertices[0][0] - f_k_face.vertices[1][0]) / float(len(B2))
        else:  # -x direction
          f_k_strip_width = abs(f_k_face.vertices[0][2] - f_k_face.vertices[1][2]) / float(len(B2))
        starting_face_index = [j for j in range(len(self.protrusion_path)) if self.protrusion_path[j] == bridge_index].pop()

        path_index = starting_face_index
        cur_face_index = self.protrusion_path[path_index]

        while cur_face_index != f_k:
          cur_face = self.full_graph.get_V()[cur_face_index]
          if flip_sides:
            next_p = np.array([cur_face.vertices[1][0], cur_face.vertices[1][1] + cut_depth, cur_face.vertices[1][2]])
          else:
            next_p = np.array([cur_face.vertices[3][0], cur_face.vertices[3][1] - cut_depth, cur_face.vertices[3][2]])
          cut_paths[cut_num].append(next_p)
          if flip_sides:
            path_index += 1
          else:
            path_index -= 1
          cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]
        '''
        # go down f_k
        if f_k_face.direction == "+z":
          if flip_sides:
            x = f_k_face.vertices[1][0] - (f_k_strip_width * i) - f_k_strip_width/(num_cuts - 1) * cut_num
          else:
            x = f_k_face.vertices[3][0] + (f_k_strip_width * i) + f_k_strip_width/(num_cuts - 1) * cut_num
        else:  # -z direction
          if flip_sides:
            x = f_k_face.vertices[1][0] + (f_k_strip_width * i) + f_k_strip_width/(num_cuts - 1) * cut_num
          else:
            x = f_k_face.vertices[3][0] - (f_k_strip_width * i) - f_k_strip_width/(num_cuts - 1) * cut_num
        y = next_p[1]
        z = next_p[2]

        cut_paths[cut_num].append(np.array([x, y, z]))

        # move down into protrusion along f_k
        if flip_sides:
          y = f_k_face.vertices[1][1] + cut_depth2
          cur_face_index = self.protrusion_path[1]
        else:
          y = f_k_face.vertices[3][1] - cut_depth2
          cur_face_index = self.protrusion_path[-1]
        cut_paths[cut_num].append(np.array([x, y, z]))

        # go along until reaching f_0
        while cur_face_index != f_0:
          cur_face = self.full_graph.get_V()[cur_face_index]
          if flip_sides:
            next_p = np.array([cur_face.vertices[1][0], y, cur_face.vertices[1][2]])
            path_index += 1
          else:
            next_p = np.array([cur_face.vertices[3][0], y, cur_face.vertices[3][2]])
            path_index -= 1
          cut_paths[cut_num].append(next_p)

          cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]

        # go down to parent bridge along f_0
        if f_0_face.direction == "+z":
          if flip_sides:
            x = f_0_face.vertices[2][0] - strip_width_to_parent*(cut_num + cuts_so_far)#(f_0_strip_width * i) - f_0_strip_width/(num_cuts-1) * cut_num
          else:
            x = f_0_face.vertices[0][0] + strip_width_to_parent*(cut_num + cuts_so_far)#(f_0_strip_width * i) + f_0_strip_width/(num_cuts-1) * cut_num
        else:
          if flip_sides:
            x = f_0_face.vertices[2][0] + strip_width_to_parent*(cut_num + cuts_so_far)#(f_0_strip_width * i) + f_0_strip_width/(num_cuts-1) * cut_num
          else:
            x = f_0_face.vertices[0][0] - strip_width_to_parent*(cut_num + cuts_so_far)#(f_0_strip_width * i) - f_0_strip_width/(num_cuts-1) * cut_num
        z = next_p[2]

        cut_paths[cut_num].append(np.array([x, y, z]))
        if flip_sides:
          y = f_0_face.vertices[3][1]
        else:
          y = f_0_face.vertices[0][1]
        cut_paths[cut_num].append(np.array([x, y, z]))
        '''

      all_cuts.extend(cut_paths)

      cuts_so_far += num_cuts-1
    for c in all_cuts:
        print c
    # self.write_cut_path_as_fold(all_cuts)
    self.write_cut_path(all_cuts, "../out/cuts.txt")
    # for p in all_cuts:
    #   print p

  def unfold_strip_root(self, f_0):

    self.compute_protrusion_path(f_0)
    f_0_vertices = self.full_graph.get_V()[f_0].vertices

    #num_cuts = self.num_leaves*2 + 1 #TODO: define this correctly - there are 6 strips but cuts come together, so only 4 cuts (3 pairs + L(Q_1))
    num_cuts = 3
    # num_cuts_on_connector = self.num_leaves + 1
    num_cuts_on_connector = 2

    cut_paths = [[] for cut_num in range(num_cuts)]
    # bridge_index = self.child_bridge  # TODO
    bridge_index = f_0
    bridge_face = self.full_graph.get_V()[bridge_index]

    # flip_sides = self.parent.in_layer(self.y_minus_1)
    flip_sides = True  # define direction to move

    if flip_sides:  # connector is flipped
      self.f_i, self.f_j = self.f_j, self.f_i

    f_i_face = self.full_graph.get_V()[self.f_i]
    # strip width only depends on number of cuts on connector (num leaves + 1)
    f_i_strip_width = abs(f_i_face.vertices[0][0] - f_i_face.vertices[1][0])/float(num_cuts_on_connector - 2)
    print f_i_strip_width

    for cut_num in range(num_cuts):
      cut_width = abs(bridge_face.vertices[0][0] - bridge_face.vertices[1][0]) / float(num_cuts - 1)
      cut_depth = ((self.depth/2.)/(num_cuts-1))*cut_num

      # start at f_0 and go down to depth
      if not (cut_num == 0):
        if bridge_face.direction == "+z":
          if flip_sides:
            cut_point = np.array([bridge_face.vertices[1][0] - cut_width * cut_num,
                                  bridge_face.vertices[1][1], bridge_face.vertices[1][2]])
          else:
            cut_point = np.array([bridge_face.vertices[3][0] + cut_width * cut_num,
                                  bridge_face.vertices[3][1], bridge_face.vertices[3][2]])
        elif bridge_face.direction == "-z":
          if flip_sides:
            cut_point = np.array([bridge_face.vertices[1][0] + cut_width * cut_num,
                                  bridge_face.vertices[1][1], bridge_face.vertices[1][2]])
          else:
            cut_point = np.array([bridge_face.vertices[3][0] - cut_width * cut_num,
                                  bridge_face.vertices[3][1], bridge_face.vertices[3][2]])

        cut_paths[cut_num].append(cut_point)
        if flip_sides:
          y = cut_point[1] + cut_depth
        else:
          y = cut_point[1] - cut_depth
        cut_paths[cut_num].append(np.array([cut_point[0], y, cut_point[2]]))

      if flip_sides:
        f_1 = self.protrusion_path[-1]
      else:
        f_1 = self.protrusion_path[1]

      starting_face_index = [j for j in range(len(self.protrusion_path)) if self.protrusion_path[j] == bridge_index].pop()

      path_index = starting_face_index
      cur_face_index = self.protrusion_path[path_index]

      while cur_face_index != f_1:
        cur_face = self.full_graph.get_V()[cur_face_index]
        if flip_sides:
          next_p = np.array([cur_face.vertices[1][0], cur_face.vertices[1][1] + cut_depth, cur_face.vertices[1][2]])
        else:
          next_p = np.array([cur_face.vertices[3][0], cur_face.vertices[3][1] - cut_depth, cur_face.vertices[3][2]])
        cut_paths[cut_num].append(next_p)
        if flip_sides:
          path_index += 1
        else:
          path_index -= 1
        cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]

      f_1_face = self.full_graph.get_V()[f_1]
      if f_1_face.direction == "+z":
        f_1_strip_width = abs(f_1_face.vertices[0][0] - f_1_face.vertices[1][0]) / float(num_cuts-1)
      else:  # -x direction
        f_1_strip_width = abs(f_1_face.vertices[0][2] - f_1_face.vertices[1][2]) / float(num_cuts-1)
      print f_1, f_1_strip_width, num_cuts

      # go down f_1
      if f_1_face.direction == "+z":
        if flip_sides:
          x = f_1_face.vertices[1][0] - f_1_strip_width * cut_num#(num_cuts-cut_num-1)
        else:
          x = f_1_face.vertices[3][0] + f_1_strip_width * cut_num
        z = next_p[2]


      elif f_1_face.direction == "-z":
        if flip_sides:
          x = f_1_face.vertices[1][0] + f_1_strip_width * cut_num#(num_cuts-cut_num-1)
        else:
          x = f_1_face.vertices[3][0] - f_1_strip_width * cut_num
        z = next_p[2]
      elif f_1_face.direction == "+x":
        if flip_sides:
          z = f_1_face.vertices[0][2] - f_1_strip_width * (num_cuts-cut_num-1)
        else:
          z = f_1_face.vertices[2][2] + f_1_strip_width * (num_cuts-cut_num-1)
        x = next_p[0]
      else:
        if flip_sides:
          z = f_1_face.vertices[0][2] + f_1_strip_width * (num_cuts-cut_num-1)
        else:
          z = f_1_face.vertices[2][2] - f_1_strip_width * (num_cuts-cut_num-1)
        x = next_p[0]

      y = next_p[1]

      cut_paths[cut_num].append(np.array([x, y, z]))

      if flip_sides:
        y += self.depth/2.
      else:
        y -= self.depth/2.

      cut_paths[cut_num].append(np.array([x, y, z]))

      #continue to f_j, then stop half the cuts
      if not f_1_face.direction in ["+x", "-x"]:
        if flip_sides:
          path_index += 1
        else:
          path_index -= 1
      cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]

      while cur_face_index != self.f_j:
        cur_face = self.full_graph.get_V()[cur_face_index]
        if flip_sides:
          next_p = np.array([cur_face.vertices[1][0], y, cur_face.vertices[1][2]])
        else:
          next_p = np.array([cur_face.vertices[3][0], y, cur_face.vertices[3][2]])
        cut_paths[cut_num].append(next_p)
        if flip_sides:
          path_index += 1
        else:
          path_index -= 1
        cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]

      if 0 < cut_num < num_cuts_on_connector:
        # follow ~half of cuts along to f_i then connector and piece them together with others
        while cur_face_index != self.f_i:
          cur_face = self.full_graph.get_V()[cur_face_index]
          if flip_sides:
            next_p = np.array([cur_face.vertices[1][0], y, cur_face.vertices[1][2]])
          else:
            next_p = np.array([cur_face.vertices[3][0], y, cur_face.vertices[3][2]])
          cut_paths[cut_num].append(next_p)
          if flip_sides:
            path_index += 1
          else:
            path_index -= 1
          cur_face_index = self.protrusion_path[path_index % len(self.protrusion_path)]

        cur_face = self.full_graph.get_V()[cur_face_index]
        # go to correct width, then down connector
        if cur_face.direction == "+z":
          if flip_sides:
            x = cur_face.vertices[1][0] - f_i_strip_width * (cut_num-1)#(num_cuts-cut_num-1)
          else:
            x = cur_face.vertices[3][0] + f_i_strip_width * (cut_num-1)

        elif cur_face.direction == "-z":
          if flip_sides:
            x = cur_face.vertices[1][0] + f_i_strip_width * (cut_num-1)#(num_cuts-cut_num-1)
          else:
            x = cur_face.vertices[3][0] - f_i_strip_width * (cut_num-1)
        z = cur_face.vertices[0][2]
        cut_paths[cut_num].append(np.array([x, y, z]))

        # to connector edge
        if flip_sides:
          y = f_i_face.vertices[3][1]
        else:
          y = f_i_face.vertices[0][1]
        cut_paths[cut_num].append(np.array([x, y, z]))

        # down connector
        f_j_face = self.full_graph.get_V()[self.f_j]
        if flip_sides:
          z = f_j_face.vertices[3][2]
        else:
          z = f_j_face.vertices[1][2]
        cut_paths[cut_num].append(np.array([x, y, z]))

        # up f_j, then connect to other cuts
        if flip_sides:
          y = f_j_face.vertices[3][1] - (self.depth/2.)/(num_cuts-1)*(cut_num-1)
        else:
          y = f_j_face.vertices[1][1] + (self.depth/2.)/(num_cuts-1)*(cut_num-1)
        cut_paths[cut_num].append(np.array([x, y, z]))

        if flip_sides:
          x = f_j_face.vertices[3][0]
        else:
          x = f_j_face.vertices[1][0]
        cut_paths[cut_num].append(np.array([x, y, z]))

    # finally, check to see we connected the cuts
    self.write_cut_path(cut_paths, "../out/cuts.txt")
    print "total_cuts:", len(cut_paths)
    # if self.num_leaves != 1:
    for cut_num in range(1, num_cuts_on_connector):
      if not np.allclose(cut_paths[cut_num][-1], cut_paths[num_cuts-cut_num][-1]):
        raise Exception("%s in cut %s is not equal to %s in cut %s" %
                        (cut_paths[cut_num][-1], cut_num+1, cut_paths[num_cuts-cut_num][-1], num_cuts-cut_num+1))



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
  def get_faces(self):
    return [face for face in self.full_graph.get_V()]

  def get_face(self, face):
    return self.full_graph.get_V()[face]


  def __str__(self):
    return self.full_graph.__str__()

