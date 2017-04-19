import json
import sys

def create_box(corner, length, width, height):
    """
    create a rectangular box
    :param top_corner: tup (x, y, z) coordinates of bottom left front corner
    :param length: float
    :param width: float
    :param height: float
    :return: a dictionary to mapping vertices, edges, and faces (list of lists) as specified in the FOLD format
    """
    vertices = []
    faces = []
    edges = []

    assert length > 0 and width > 0 and height > 0
    corner_x, corner_y, corner_z = corner
    for x in (corner_x, length):
        for y in (corner_y, width):
            for z in (corner_z, height):
                vertices.append([x, y, z])

    edges.extend([[0, 1], [1, 2], [2, 3], [3, 1]])
    edges.extend([[4, 5], [5, 6], [6, 7], [7, 4]])
    edges.extend([[3, 6], [7, 2], [1, 4], [5, 0]])

    faces.append([0, 1, 2, 3])
    faces.append([4, 5, 6, 7])
    faces.append([2, 3, 6, 7])
    faces.append([0, 1, 4, 5])
    faces.append([0, 2, 4, 6])
    faces.append([1, 3, 5, 7])

    return {"vertices_coords": vertices,
            "faces_vertices": faces,
            "edges_vertices": edges}


if __name__ == "__main__":
    # FOR TESTING
    create_box([0, 0, 0], 1, 1, 1)
    # corner, l, w, h, filename = sys.argv[1:6]

