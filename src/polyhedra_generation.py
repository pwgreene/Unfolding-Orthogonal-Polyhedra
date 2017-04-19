import json
import sys

def create_box(corner, length, width, height):
    """
    create a rectangular box
    :param corner: tup (x, y, z) coordinates of bottom left front corner
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
    for x in (0, length):
        for y in (0, width):
            for z in (0, height):
                vertices.append([corner_x+x, corner_y+y, corner_z+z])

    edges.extend([[0, 1], [1, 3], [3, 2], [2, 0]])
    edges.extend([[4, 5], [5, 7], [7, 6], [6, 4]])
    edges.extend([[6, 2], [4, 0], [5, 1], [7, 3]])
    # faces oriented with normals pointed outwards
    faces.append([0, 1, 3, 2])
    faces.append([6, 7, 5, 4])
    faces.append([2, 3, 7, 6])
    faces.append([4, 5, 1, 0])
    faces.append([0, 2, 6, 4])
    faces.append([5, 7, 3, 1])

    return {"vertices_coords": vertices,
            "faces_vertices": faces,
            "edges_vertices": edges}



def create_fold_file(filename, d, append=False):
    """
    create a json object out of a given python dictionary
    :param filename: string
    :param d: dictionary
    :return: json object
    """
    if append:
        with open(filename) as f:
            data = json.load(f)
        # need to offset the indices of the new vertices
        num_existing_vertices = len(data["vertices_coords"])
        num_existing_edges = len(data["edges_vertices"])
        for e in range(len(d["edges_vertices"])):
            for i in range(2):
                d["edges_vertices"][e][i] += num_existing_vertices
        for f in range(len(d["faces_vertices"])):
            for i in range(4):
                d["faces_vertices"][f][i] += num_existing_vertices
        data["vertices_coords"].extend(d["vertices_coords"])
        data["edges_vertices"].extend(d["edges_vertices"])
        data["faces_vertices"].extend(d["faces_vertices"])
        data["edges_assignment"].extend("B" for i in range(num_existing_edges))
        print d

    else:
        # start with default stuff
        data = {"file_spec":1,
                "file_creator": "Python script",
                "file_author": "pwgreene",
                "frame_title": filename,
                "file_classes": ["singleModel"],
                "frame_classes": ["foldedForm"],
                "frame_attributes": ["3D"],
                }
        data.update(d)
        # we're not dealing with foldings, so edge assignments don't matter
        num_edges = len(d["edges_vertices"])
        edges_assignment = {"edges_assignment": ["B" for x in range(num_edges)]}
        data.update(edges_assignment)

    with open(filename, 'w') as f:
        json.dump(data, f, ensure_ascii=True, indent=2)


if __name__ == "__main__":
    # FOR TESTING

    new_box = create_box([-1, 0, 0], 1, 1, 1)
    create_fold_file("../data/boxes.fold", new_box)
    create_fold_file("../data/boxes.fold", create_box([0, 0, 0], 1, 3, 5), append=True)
    # corner, l, w, h, filename = sys.argv[1:6]

