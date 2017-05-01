class Graph:
  
  # Input V as a list of vertices
  # Input exactly 1 of E or E_long_list (not both or none)
  # E is edges as a list of lists where list i is the list of indices of vertices vertex i is connected to
  # E_long_list is a list of all edges where an edge is a list of length 2 of the indices of the vertices it connects, edges are assumed bidirectional and should not be input twice (eg. if [0, 1] is given [1, 0] should not be)
  def __init__(self, V, E=None, E_long_list=None):
    if not E and not E_long_list:
      raise Exception("No edges inputted")
    elif E and E_long_list:
      raise Exception("Can't input edges in multiple formats")
    self.V = V
    if E:
      self.E = E
    else:
      self.E = []
      for _ in xrange(len(V)):
        self.E.append([])
      for edge in E_long_list:
        # assumes bi-directional
        self.E[edge[0]].append(edge[1])
        self.E[edge[1]].append(edge[0])

  def __str__(self):
    rep = ""
    for i in range(len(self.E)):
      for j in range(len(self.E[i])):
        rep += "(%s)->(%s)\n" % (self.V[i], self.V[self.E[i][j]])
    return rep
