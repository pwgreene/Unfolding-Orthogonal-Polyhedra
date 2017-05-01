import copy

class Graph:
  
  # Input V as a list of vertices
  # Input exactly 1 of E or E_long_list (not both or none)
  # E is edges as a list of lists where list i is the list of indices of vertices vertex i is connected to
  # E_long_list is a list of all edges where an edge is a list of length 2 of the indices of the vertices it connects, edges are assumed bidirectional and should not be input twice (eg. if [0, 1] is given [1, 0] should not be)
  def __init__(self, V, E=None, E_list=None):
    if not E and not E_list:
      raise Exception("No edges inputted")
    elif E and E_list:
      raise Exception("Can't input edges in multiple formats")
    self.V = V
    if E:
      self.E = E
    else:
      self.E = []
      for _ in xrange(len(V)):
        self.E.append([])
      for edge in E_list:
        # assumes bi-directional
        self.E[edge[0]].append(edge[1])
        self.E[edge[1]].append(edge[0])

  def get_V(self):
    return self.V

  # returns list of indicies of vertices connected to vertex i
  def get_connections(self, i):
    return self.E[i]
  
  # returns list of vertices with indices from index_list
  def get_vertices(self, index_list):
    return [self.V[i] for i in index_list]

  def copy(self):
    return copy.deepcopy(self)

  # deletes vertex i and all associated connections
  # reindexes edges to reflect new vertex indexing
  def del_vertex(self, i):
    self.V.pop(i)
    self.E.pop(i)
    for j in xrange(len(self.E)):
      new_edge_list = []
      for e in self.E[j]:
        if e > i:
          new_edge_list.append(e - 1)
        elif e < i:
          new_edge_list.append(e)
      self.E[j] = new_edge_list

  def __str__(self):
    rep = ""
    for i in range(len(self.E)):
      for j in range(len(self.E[i])):
        rep += "{%s} -> {%s}\n" % (self.V[i], self.V[self.E[i][j]])
    return rep
