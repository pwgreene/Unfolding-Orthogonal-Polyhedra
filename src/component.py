class Component:
  
  # face_graph is the graph with faces as vertices of the component
  def __init__(self, face_graph, y, y_minus_1):
    self.full_graph = face_graph
    faces = self.full_graph.get_V()
    
    prot_faces = [key for key in faces if not faces[key].in_layer(y) and not faces[key].in_layer(y_minus_1)]
    self.prot = self.full_graph.subgraph(prot_faces)
    
    back_rim_faces = [key for key in faces if faces[key].in_layer(y)]
    self.back_rim = self.full_graph.subgraph(back_rim_faces)
    
    front_rim_faces = [key for key in faces if faces[key].in_layer(y_minus_1)]
    self.front_rim = self.full_graph.subgraph(front_rim_faces)
   
  def __str__(self):
    return self.full_graph.__str__()
