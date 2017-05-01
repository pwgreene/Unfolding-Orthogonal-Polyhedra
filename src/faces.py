import numpy as np


class Face(object):

    def __init__(self, vertices, direction):
        """
        initialize the face, vertices v1, v2, v3, v4 (np arrays) should be defined in ccw order
        edges of the face are (v1, v2), (v2, v3), (v3, v4), and (v4, v1)
        vertices: list of vertex coordinates
        type: defines the type of face
        direction: +- x, +- y, or +- z (string)
        """
        self.vertices = vertices
        self.principal_directions = ["x", "y", "z"]
        if len(direction) == 1:
            self.direction = "+"+direction
        else:
            self.direction = direction

    def get_direction(self):
        return self.direction

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


if __name__ == "__main__":
    # TODO: test this with rectangular faces
    v_1 = np.array([0., 0., 0.])
    v_2 = np.array([1., 0., 0.])
    v_3 = np.array([1., 1., 0.])
    v_4 = np.array([0., 1., 0.])
    f = Face([v_1, v_2, v_3, v_4], "+z")
    print f
    for face in f.divide_face("y", (1, 1)):
        print face