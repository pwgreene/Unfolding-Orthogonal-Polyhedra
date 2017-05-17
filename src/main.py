from polyhedron import Polyhedron
import numpy as np

def unfold_polyhedron(polyhedron, root, bridge_cuts=[]):
    """
    unfolds each component of the polyhedron
    :param polyhedron: Polyhedron object
    :param root: ComponentNode object representing node rooted here
    :param bridge_cuts: list of tuples of the format ((f_0 for child, parent), number of leaves)
    :return:
    """
    if len(root.children) == 0:  # leaf
        print "processing leaf"
        f_0 = root.parent_bridge[-1]
        print "f_0", f_0
        print "bridge:", root.parent_bridge
        bridge_normal = polyhedron.faces[root.parent_bridge[0]].direction
        # check direction of parent component
        if bridge_normal == "+y":
            parent_direction = "+y"
        else:
            parent_direction = "-y"
        root.component.unfold_strip_leaf(f_0, parent_direction)
        # bridge_cuts.append(((f_0, root.parent), 1))
        root.f_0 = f_0
        root.num_leaves = 1
        return 1

    elif root.parent is None:  # root
        bridge_normal = polyhedron.faces[root.children[0].parent_bridge[0]].direction
        # check direction of child component
        if bridge_normal == "+y":
            child_direction = "-y"
        else:
            child_direction = "+y"
        num_leaves = unfold_polyhedron(polyhedron, root.children[0])
        print "processing root"
        root.component.unfold_strip_root(root.children_bridges[0][0], child_direction, num_leaves)
        root.num_leaves = num_leaves
        return num_leaves

    else:  # intermediate, recursive case
        num_leaves = 0
        child_faces = []
        child_face_directions = []
        f_0 = root.parent_bridge[-1]
        print "f_0", f_0
        print "parent:", root.parent_bridge
        # print "chlidren:", root.children[0].parent_bridge
        bridge_normal = polyhedron.faces[f_0].direction
        if bridge_normal == "+y":
            parent_direction = "-y"
        else:
            parent_direction = "+y"
        num_leaves_children = []
        # iterate through children and unfold them recursively
        for i in range(len(root.children)):
            child = root.children[i]
            child_bridge = root.children_bridges[i] + root.children[i].parent_bridge
            print "child_bridge", child_bridge
            bridge_normal = polyhedron.faces[child_bridge[-1]].direction
            if bridge_normal == "+y":
                child_direction = "-y"
            else:
                child_direction = "+y"

            child_faces.append(child_bridge[0])
            child_face_directions.append(child_direction)
            child_leaves = unfold_polyhedron(polyhedron, child)
            num_leaves += child_leaves
            num_leaves_children.append(child_leaves)
            print "-------"
            print child_faces
            print child_face_directions
            print child_leaves
            print num_leaves
            print num_leaves_children
            print "------"
        # bridge_cuts.append(((f_0, root.parent), num_leaves))
        root.f_0 = f_0
        root.num_leaves = num_leaves
        print "processing intermediate"
        root.component.unfold_strip_intermediate(child_faces, child_face_directions, f_0,
                                                 parent_direction, num_leaves_children)
        return num_leaves

def gather_cuts(polyhedron, root, cuts=[]):

    cuts.extend(root.component.cuts)
    if len(root.children) > 0:
        for i in range(len(root.children)):
            gather_cuts(polyhedron, root.children[i], cuts)
            # bridge cuts
            bridge_start_face = polyhedron.faces[root.children_bridges[i][0]]
            bridge_end_face = polyhedron.faces[root.children[i].f_0]
            bridge_middle_face = polyhedron.faces[root.children[i].parent_bridge[0]]

            x_min, x_max = sorted([bridge_start_face.vertices[0][0], bridge_start_face.vertices[1][0]])
            z_min = bridge_start_face.vertices[0][2]
            z_max = bridge_end_face.vertices[0][2]
            y = bridge_middle_face.vertices[0][1]
            num_cuts = root.children[i].num_leaves * 2 + 1
            cut_width = abs(x_max - x_min) / float(num_cuts - 1)
            for cut in range(num_cuts):
                v1 = np.array([x_min, y, z_min])
                v2 = np.array([x_min, y, z_max])
                cuts.append([v1, v2])
                x_min += cut_width
            print bridge_start_face.vertices, bridge_end_face.vertices
    return cuts

def write_cut_path(paths, filename):
    n_paths = len(paths)
    with open(filename, 'w') as f:
      f.write("%s\n" % n_paths)
      for path in paths:
        f.write("%s\n" % len(path))
        for point in path:
          f.write("%s %s %s\n" % (point[0], point[1], point[2]))

if __name__ == "__main__":
    # p1 = Polyhedron(filelist=["../data/test/unit_cube_open.fold", "../data/test/rect_box.fold"])
    p = Polyhedron(filelist=["../data/the_box.fold"])
    # print unfold_polyhedron(p1, p1.unfolding_tree)
    print p.unfolding_tree.children_bridges
    unfold_polyhedron(p, p.unfolding_tree)
    print p.unfolding_tree.children[0].children
    p.write_to_off("../out/poly.off")
    write_cut_path(gather_cuts(p, p.unfolding_tree), "../out/cuts.txt")
