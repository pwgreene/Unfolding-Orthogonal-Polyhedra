import copy

class Graph:
  
  # Input V as a list of vertices or a dictionary of label:vertex pairs
  # Input exactly 1 of E, E_list or E_dict
  # E is edges as a list of lists where list i is the list of indices of vertices vertex i is connected to (adjacency list)
  # E_list is a list of all edges where an edge is a list of length 2 of the indices of the vertices it connects, edges are assumed bidirectional and should not be input twice (eg. if [0, 1] is given [1, 0] should not be)
  # E_dict is dictionary where key is the label of the first vertex in the edge and the value is a list of vertex labels for vertices it is connected to
  def __init__(self, V, E=None, E_list=None, E_dict=None):
    if sum(map(bool, [E, E_list, E_dict])) > 1:
      raise Exception("Can't input edges in multiple formats")
    
    if type(V) is list:
      self.V = {i:V[i] for i in xrange(len(V))}
    elif type(V) is dict:
      self.V = V
    else:
      raise Exception("Invalid value for V")
    
    if E:
      self.E = {i:E[i] for i in xrange(len(E))}
    elif E_list:
      self.E = {i:[] for i in xrange(len(V))}
      for edge in E_list:
        # assumes bi-directional
        self.E[edge[0]].append(edge[1])
        self.E[edge[1]].append(edge[0])
    elif E_dict:
      self.E = E_dict
    else:
      self.E = {i:[] for i in xrange(len(V))}

  def get_V(self):
    return self.V

  # returns list of indicies of vertices connected to vertex i
  def get_connections(self, i):
    return list(self.E[i])
  
  # returns list of vertices with indices from index_list
  def get_vertices(self, index_list):
    return [self.V[i] for i in index_list]
  
  # returns vertex with inputted index
  def get_vertex(self, index):
    return self.V[index]
  
  def copy(self):
    return copy.deepcopy(self)

  # deletes vertex i and all associated connections
  def del_vertex(self, i):
    self.V.pop(i)
    self.E.pop(i)
    for j in self.E:
      self.E[j] = [e for e in self.E[j] if e != i]

  # returns the subgraph containing only vertices with labels in vertex_list
  # returns a new graph object
  def subgraph(self, vertex_list):
    vertices = {key:self.V[key] for key in self.V if key in vertex_list}
    edges = {key:self.E[key] for key in self.E if key in vertex_list}
    for j in edges:
      edges[j] = [e for e in edges[j] if e in vertex_list]
    return Graph(vertices, E_dict=edges)
  
  # returns a list of vertices reachable when starting from start (start is included in this dict)
  def get_reachable(self, start):
    visited = [start]
    queue = self.get_connections(start)
    while queue:
      current = queue.pop()
      if current in visited:
        continue
      else:
        visited.append(current)
        queue.extend(self.get_connections(current))
    return visited

  # combines vertices with label in vertex_list and represents combo as new_vertex with vertex_list[0] as the label
  # returns a new graph object
  def combine_vertices(self, vertex_list, new_vertex):
    vertices = {key:self.V[key] for key in self.V if key not in vertex_list}
    vertices[vertex_list[0]] = new_vertex
    new_vertex_edges = [self.E[vertex] for vertex in vertex_list]
    new_vertex_edges_flattened = list(set([e for sublist in new_vertex_edges for e in sublist if e not in vertex_list]))
    edges = {key:self.E[key] for key in self.E if key not in vertex_list}
    for j in edges:
      edges[j] = [e for e in edges[j] if e not in vertex_list]
    for i in new_vertex_edges_flattened:
      edges[i].append(vertex_list[0])
    edges[vertex_list[0]] = new_vertex_edges_flattened
    return Graph(vertices, E_dict=edges)

  def __str__(self):
    rep = ""
    for i in self.E:
      for j in self.E[i]:
        rep += "{%s} -> {%s}\n" % (self.V[i], self.V[j])
    return rep
