class ComponentNode:
  
  def __init__(self, component):
    self.component = component
    self.children = []
    self.parent = None
    self.parent_bridge = None
    self.children_bridges = []

  def add_child(self, child):
    self.children.append(child)
    child.add_parent(self)

  def add_parent(self, parent):
    self.parent = parent
  
  # bridge is the list of face keys in component that are part of the bridge to the most recently added child
  def add_child_bridge(self, bridge):
    self.children_bridges.append(bridge)

  # bridge is the list of face keys in component that are part of the bridge to parent
  def add_parent_bridge(self, bridge):
    self.parent_bridge = bridge
