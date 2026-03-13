import numpy as np

def from_navis(navis_skel) -> "Skeleton":
  """
  Convert navis skeletons to osteoid.Skeleton. This
  should be more efficient than the SWC interchange method.
  """
  from .. import Skeleton

  vertex_types = None
  if len(navis_skel.nodes.type):
    # 'root', 'slab', 'branch', 'end', first letter of each
    vertex_types = [ ord(t[0]) for t in navis_skel.nodes.type ] 
    mapping = np.zeros(ord('t') + 1, dtype=np.uint8)
    mapping[ord('s')] = 7 # custom
    mapping[ord('r')] = 8 # custom
    mapping[ord('b')] = 5 # fork point
    mapping[ord('t')] = 6 # end point
    vertex_types = mapping[vertex_types]

  return Skeleton(
    vertices=navis_skel.vertices, 
    edges=navis_skel.edges - 1,
    radii=navis_skel.nodes.radius.to_numpy(),
    vertex_types=vertex_types,
  )
