from polyhedron import Polyhedron

def unfold_polyhedron(polyhedron, root, bridge_cuts=[]):
    """
    unfolds each component of the polyhedron
    :param polyhedron: Polyhedron object
    :param root: ComponentNode object representing node rooted here
    :param bridge_cuts: list of tuples of the format ((f_0 for child, parent), number of leaves)
    :return:
    """
    if len(root.children) == 0:  # leaf
        f_0 = root.parent_bridge[-1]
        bridge_normal = polyhedron.faces[f_0].direction
        # check direction of parent component
        if bridge_normal == "+y":
            parent_direction = "+y"
        else:
            parent_direction = "-y"
        root.component.unfold_strip_leaf(f_0, parent_direction)
        bridge_cuts.append(((f_0, root.parent), 1))
        return 1

    elif root.parent is None:  # root
        print "processing root"
        bridge_normal = polyhedron.faces[root.children[0].parent_bridge[0]].direction
        # check direction of child component
        if bridge_normal == "+y":
            child_direction = "-y"
        else:
            child_direction = "+y"
        num_leaves = unfold_polyhedron(polyhedron, root.children[0])
        root.component.unfold_strip_root(root.children_bridges[0][0], child_direction, num_leaves)
        return num_leaves

    else:  # intermediate, recursive case
        num_leaves = 0
        child_faces = []
        child_face_directions = []
        f_0 = root.parent_bridge[0]
        bridge_normal = polyhedron.faces[f_0].direction
        if bridge_normal == "+y":
            parent_direction = "-y"
        else:
            parent_direction = "+y"
        num_leaves_children = []
        # iterate through children and unfold them recursively
        for i in range(len(root.children)):
            child = root.children[i]
            child_bridge = root.children_bridges[i]
            bridge_normal = polyhedron.faces(child_bridge[1]).direction
            if bridge_normal == "+y":
                child_direction = "-y"
            else:
                child_direction = "+y"

            child_faces.append(child_bridge[0])
            child_face_directions.append(child_direction)
            child_leaves = unfold_polyhedron(polyhedron, child)
            num_leaves += child_leaves
            num_leaves_children.append(child_leaves)
        bridge_cuts.append(((f_0, root.parent), num_leaves))
        root.component.unfold_strip_intermediate(child_faces, child_face_directions, f_0,
                                                 parent_direction, num_leaves_children)
        return num_leaves

def gather_cuts(polyhedron, root, cuts=[]):

    cuts.extend(root.component.cuts)
    if len(root.children) > 0:
        for child in root.children:
            gather_cuts(polyhedron, child, cuts)
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
    p1 = Polyhedron(filelist=["../data/test/unit_cube_open.fold", "../data/test/rect_box.fold"])
    print unfold_polyhedron(p1, p1.unfolding_tree)
    write_cut_path(gather_cuts(p1, p1.unfolding_tree), "../out/cuts.txt")