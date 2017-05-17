import numpy as np

class Face(object):

    def __init__(self, vertices, direction, sort_vertices=True):
        """
        initialize the face, vertices v1, v2, v3, v4 (np arrays) should be defined in ccw order
        edges of the face are (v1, v2), (v2, v3), (v3, v4), and (v4, v1)
        vertices: list of vertex coordinates
        type: defines the type of face
        sort_vertices: (bool) sort list of vertices by self.sorted_vertices or not. default=True
        direction: +- x, +- y, or +- z (string)
        """
        self.vertices = vertices
        self.principal_directions = ["x", "y", "z"]
        if len(direction) == 1:
            self.direction = "+"+direction
        else:
            self.direction = direction
        self.center = self.compute_center()
        if sort_vertices:
            self.vertices = self.sorted_vertices()

    def get_direction(self):
        return self.direction

    def get_vertices(self, as_tuple=False):
        if as_tuple:
            return [tuple(self.vertices[i]) for i in range(len(self.vertices))]
        else:  # as numpy arrays
            return [self.vertices[i].copy() for i in range(len(self.vertices))]

    def get_direction_as_vector(self):
        sign = self.direction[0]
        axis = self.direction[1]
        direction = np.zeros((3, ))
        if axis == "x":
            direction[0] = 1
        elif axis == "y":
            direction[1] = 1
        elif axis == "z":
            direction[2] = 1
        else:
            raise ValueError

        if sign == "-":
            return -direction
        else:
            return direction

    def __str__(self):
        return ("(%s, %s, %s, %s) direction: "+self.direction) % tuple(self.vertices)

    # returns true if this face is in layer y=y, returns false otherwise
    def in_layer(self, y):
      for vertex in self.vertices:
        if vertex[1] != y:
          return False
      return True
    
    # returns true if this face is between layers y=y and y=y_minus_1 (inclusive)
    # if the face is entirely in layer y (resp. y_minus_1), its normal must be +y (resp. -y)
    # assumes y_minus_1 < y
    def between_layers(self, y, y_minus_1):
      if self.in_layer(y):
        if self.direction == '+y':
          return True
        else:
          return False
      elif self.in_layer(y_minus_1):
        if self.direction == '-y':
          return True
        else:
          return False
      for vertex in self.vertices:
        if vertex[1] < y_minus_1 or vertex[1] > y:
          return False
      return True

    def divide_face(self, axis, ratio):
        """
        divide the face into two pieces on axis according to ratio
        axis: (string) "x", "y", "z"
        ratio: len 2 tup of float or string
        :return: two new Faces
        """
        assert axis in self.principal_directions and len(ratio) == 2
        total = sum(ratio)
        w1 = ratio[0]/float(total)
        w2 = ratio[1]/float(total)
        # check the axis and find midpoint accordingly
        for i in range(len(self.principal_directions)):
            if axis == self.principal_directions[i]:
                # first/third edges parallel to axis
                if (self.vertices[0]-self.vertices[1])[i] != 0:
                    mid = w1*self.vertices[0][i] + w2*self.vertices[1][i]
                    new_v1, new_v2 = self.vertices[1].copy(), self.vertices[2].copy()
                    new_v1[i], new_v2[i] = mid, mid
                    return (Face([self.vertices[0], new_v1, new_v2, self.vertices[3]], self.direction),
                            Face([new_v1, self.vertices[1], self.vertices[2], new_v2], self.direction))
                else:  # second/fourth edges parallel to axis
                    mid = w1*self.vertices[1][i] + w2*self.vertices[2][i]
                    new_v2, new_v3 = self.vertices[2].copy(), self.vertices[3].copy()
                    # print new_v2, new_v3
                    new_v2[i], new_v3[i] = mid, mid
                    print new_v2, new_v3
                    return (Face([self.vertices[0], self.vertices[1], new_v2, new_v3], self.direction),
                            Face([self.vertices[3], new_v3, new_v2, self.vertices[2]], self.direction))


    def compute_center(self):
        """
        :return: average of all the vertices
        """
        avg = np.zeros((3,))
        for v in self.vertices:
            avg += .25*v
        return avg

    def sorted_vertices(self):
        """
        Return new list of vertices such that first vertex is the closet (resp. ccw) vertex
        and smallest y-value (still listed in ccw order).
        :return: new list of vertices. Doesn't modify self.vertices!
        """
        if self.direction == "+z":
            first = sorted(self.vertices, key=lambda z: (z[0], z[1]))[0]
        elif self.direction == "-z":
            first = sorted(self.vertices, key=lambda z: (-z[0], z[1]))[0]
        elif self.direction == "+x":
            first = sorted(self.vertices, key=lambda x: (-x[2], x[1]))[0]
        elif self.direction == "-x":
            first = sorted(self.vertices, key=lambda x: (x[2], x[1]))[0]
        elif self.direction == "+y":
            first = sorted(self.vertices, key=lambda y: (-y[0], y[2]))[0]
        elif self.direction == "-y":
            first = sorted(self.vertices, key=lambda y: (y[0], y[2]))[0]
        else:
            raise TypeError("no direction inputted")

        for i in range(len(self.vertices)):
            if np.array_equal(first, self.vertices[i]):
                return self.vertices[i:] + self.vertices[:i]


    def is_on_face(self, vertex):
        """
        check if the vertex is on the face
        :param vertex: np array, tuple, or list of coordinates
        :return: (bool) True if vertex on face, False o.w.
        """
        # TODO: implement this
        return False


    def __eq__(self, other):
        """Override the default Equals behavior"""
        if not isinstance(other, self.__class__):
            return False

        this_vertices = {tuple(self.vertices[i]) for i in range(len(self.vertices))}
        other_vertices = {tuple(other.vertices[i]) for i in range(len(other.vertices))}

        if len(this_vertices.union(other_vertices)) == len(this_vertices) and self.direction == other.direction:
            return True
        else:
            return False

if __name__ == "__main__":
    # TODO: test this with rectangular faces
    v_1 = np.array([0., 0., 0.])
    v_2 = np.array([1., 0., 0.])
    v_3 = np.array([1., 1., 0.])
    v_4 = np.array([0., 1., 0.])
    f1 = Face([v_1, v_2, v_3, v_4], "+z")
    f2 = Face([v_4, v_3, v_2, v_1], "-z")
    f3 = Face([v_4, v_3, v_2, v_1], "+z")
    assert f1 != f2
    assert f1 == f3
