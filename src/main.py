from polyhedron import Polyhedron

def unfold_polyhedron(p):
    """
    unfolds each component of the polyhedron
    :param p: Polyhedron object
    :return:
    """
    root = p.unfolding_tree
    print root.children_bridges
    if len(root.children) == 0:  # leaf
        f_0 = root.parent_bridge[0]
        bridge_normal = p.faces(f_0).direction
        # check direction of parent component
        if bridge_normal == "+y":
            parent_direction = "-y"
        else:
            parent_direction = "+y"
        root.component.unfold_strip_leaf(f_0, parent_direction)
        return 1

    elif root.parent is None:  # root
        bridge_normal = p.faces(root.children_bridges[0][1]).direction
        # check direction of child component
        if bridge_normal == "+y":
            child_direction = "-y"
        else:
            child_direction = "+y"
        num_leaves = unfold_polyhedron(root.children[0])
        root.component.unfold_strip_root(root.children_bridges[0], child_direction, num_leaves)
        return num_leaves

    else:  # intermediate, recursive case
        num_leaves = 0
        child_faces = []
        child_face_directions = []
        f_0 = root.parent_bridge[0]
        bridge_normal = p.faces(f_0).direction
        if bridge_normal == "+y":
            parent_direction = "-y"
        else:
            parent_direction = "+y"
        num_leaves_children = []
        # iterate through children and unfold them recursively
        for i in range(len(root.children)):
            child = root.children[i]
            child_bridge = p.children_bridges[i]
            bridge_normal = p.faces(child_bridge[1]).direction
            if bridge_normal == "+y":
                child_direction = "-y"
            else:
                child_direction = "+y"

            child_faces.append(child_bridge[0])
            child_face_directions.append(child_direction)
            child_leaves = unfold_polyhedron(child)
            num_leaves += child_leaves
            num_leaves_children.append(child_leaves)
        root.component.unfold_strip_intermediate(child_faces, child_face_directions, f_0,
                                                 parent_direction, num_leaves_children)
        return num_leaves

if __name__ == "__main__":
    p = Polyhedron(filelist=["../data/test/unit_cube_open.fold", "../data/test/rect_box.fold"])
    print unfold_polyhedron(p)