from typing import Union
import datetime
import re

import numpy as np

def to_swc(
  skeleton:"Skeleton", 
  contributors:str = "", 
  soma_threshold:float = np.inf,
) -> str:
  """
  SWC file generator. 

  The SWC format was first defined in 
  
  R.C Cannona, D.A Turner, G.K Pyapali, H.V Wheal. 
  "An on-line archive of reconstructed hippocampal neurons".
  Journal of Neuroscience Methods
  Volume 84, Issues 1-2, 1 October 1998, Pages 49-54
  doi: 10.1016/S0165-0270(98)00091-0

  This website is also helpful for understanding the format:

  https://web.archive.org/web/20180423163403/http://research.mssm.edu/cnic/swc.html

  soma_threshold: use radius information to determine if there 
    is a soma in this data and if so, use it as the root of its 
    component

  Returns: swc as a string
  """
  sx, sy, sz = np.diag(skeleton.transform)[:3]

  swc_header = f"""# ORIGINAL_SOURCE Osteoid {__VERSION__}
# CREATURE 
# REGION
# FIELD/LAYER
# TYPE
# CONTRIBUTOR {contributors}
# REFERENCE
# RAW 
# EXTRAS 
# SOMA_AREA
# SHINKAGE_CORRECTION 
# VERSION_NUMBER {__VERSION__}
# VERSION_DATE {datetime.datetime.now(datetime.timezone.utc).isoformat()}
# SCALE {sx:.6f} {sy:.6f} {sz:.6f}
"""

  def generate_swc(skel, offset):
    if skel.edges.size == 0:
      return ""

    index = defaultdict(set)
    visited = defaultdict(bool)
    for e1, e2 in skel.edges:
      index[e1].add(e2)
      index[e2].add(e1)

    root = skel.edges[0,0]
    if soma_threshold < np.inf and hasattr(skel.radii):
      if np.any(skel.radii >= soma_threshold):
        root = np.argmax(skel.radii)

    stack = [ root ]
    parents = [ -1 ]

    pairs = []

    while stack:
      node = stack.pop()
      parent = parents.pop()

      if visited[node]:
        continue

      pairs.append([node,parent])
      visited[node] = True
      
      for child in index[node]:
        stack.append(child)
        parents.append(node)

    return pairs

  def create_row(node, parent, offset):
    return [
      (node + offset),
      skel.vertex_types[node],
      skel.vertices[node][0],
      skel.vertices[node][1],
      skel.vertices[node][2],
      skel.radii[node],
      (parent if parent == -1 else (parent + offset)),        
    ]

  def render_row(row):
    return "{n} {T} {x:0.6f} {y:0.6f} {z:0.6f} {R:0.6f} {P}".format(
      n=row[0],
      T=row[1],
      x=row[2],
      y=row[3],
      z=row[4],
      R=row[5],
      P=row[6],
    )

  def renumber(rows):
    mapping = { -1: -1 }
    N = 1
    for row in rows:
      node = row[0]
      if node in mapping:
        row[0] = mapping[node]
        continue
      else:
        row[0] = N
        mapping[node] = N
        N += 1

    for row in rows:
      row[-1] = mapping[row[-1]]

    return rows

  skels = skeleton.components()
  swc = swc_header + "\n"
  offset = 0
  all_rows = []
  for skel in skels: 
    pairs = generate_swc(skel, offset)
    rows = [ 
      create_row(node, parent, offset)
      for node, parent in pairs
    ]
    del pairs
    all_rows.extend(rows)
    offset += skel.vertices.shape[0]

  all_rows = renumber(all_rows)
  swc += "\n".join((
    render_row(row)
    for row in all_rows
  ))

  return swc

def from_swc(skeleton:"Skeleton", swcstr:str) -> "Skeleton":
  """
  The SWC format was first defined in 
  
  R.C Cannona, D.A Turner, G.K Pyapali, H.V Wheal. 
  "An on-line archive of reconstructed hippocampal neurons".
  Journal of Neuroscience Methods
  Volume 84, Issues 1-2, 1 October 1998, Pages 49-54
  doi: 10.1016/S0165-0270(98)00091-0

  This website is also helpful for understanding the format:

  https://web.archive.org/web/20180423163403/http://research.mssm.edu/cnic/swc.html

  Returns: Skeleton
  """
  from .. import Skeleton
  
  lines = swcstr.split("\n")
  while len(lines) and (lines[0] == '' or re.match(r'[#\s]', lines[0][0])):
    l = lines.pop(0)

  if len(lines) == 0:
    return Skeleton()

  vertices = []
  edges = []
  radii = []
  vertex_types = []

  label_index = {}
  
  N = 0

  for line in lines:
    if line.replace(r"\s", '') == '':
      continue
    (vid, vtype, x, y, z, radius, parent_id) = line.split(" ")
    
    coord = ( float(x), float(y), float(z) )
    vid = int(vid)
    parent_id = int(parent_id)

    label_index[vid] = N

    if parent_id >= 0:
      if vid < parent_id:
        edge = [vid, parent_id]
      else:
        edge = [parent_id, vid]

      edges.append(edge)

    vertices.append(coord)
    vertex_types.append(int(vtype))

    try:
      radius = float(radius)
    except ValueError:
      radius = -1 # e.g. radius = NA or N/A

    radii.append(radius)

    N += 1

  for edge in edges:
    edge[0] = label_index[edge[0]]
    edge[1] = label_index[edge[1]]

  return Skeleton(vertices, edges, radii, vertex_types)
