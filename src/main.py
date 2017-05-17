from polyhedron import Polyhedron

def unfold_polyhedron(p, num_leaves=0):
    """
    unfolds each component of the polyhedron
    :param p: Polyhedron object
    :return:
    """
    root = p.unfolding_tree
    if len(root.children) == 0:


    elif len(root.parent) == 0:
        bridge_normal = p.faces(root.children_bridges[1])
        # check direction of child component
        if bridge_normal == "+y":
            child_direction = "-y"
        else:
            child_direction = "+y"
        num_leaves = unfold_polyhedron(root.children[0])
        root.component.unfold_strip_root(root.children_bridges[0], child_direction, num_leaves)
        return num_leaves