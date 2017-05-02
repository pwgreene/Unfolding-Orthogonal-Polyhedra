from graph import Graph

class Component:
  
  # face_graph is the graph with faces as vertices of the component
  def __init__(self, face_graph, y, y_minus_1):
    self.full_graph = face_graph
    faces = self.full_graph.get_V()
    
    prot_faces = [key for key in faces if not faces[key].in_layer(y) and not faces[key].in_layer(y_minus_1)]
    self.prot = self.full_graph.subgraph(prot_faces)
    # TODO: make self.prot a tree and choose the root using self.protrusion_tree
    
    back_rim_faces = [key for key in faces if faces[key].in_layer(y)]
    self.back_rim = self.full_graph.subgraph(back_rim_faces)
    
    front_rim_faces = [key for key in faces if faces[key].in_layer(y_minus_1)]
    self.front_rim = self.full_graph.subgraph(front_rim_faces)


  def protrusion_tree(self, root):
    """
    create a tree that navigates the protrusion ccw
    :param root: root face (should be +z or -z direction)
    :return: None
    """
    p = root
    vertices = self.prot.get_V()
    edges = {v: [] for v in vertices}
    visited = {}
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
            p = c
            visited[c] = 1
            break
    g = Graph(vertices, E_dict=edges)

   
  def __str__(self):
    return self.full_graph.__str__()
