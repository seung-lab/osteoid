from typing import Optional, Any

import numpy as np

def to_precomputed(skeleton:"Skeleton") -> bytes:
  edges = skeleton.edges.astype(np.uint32)
  vertices = skeleton.vertices.astype(np.float32)
  
  result = BytesIO()

  # Write number of positions and edges as first two uint32s
  result.write(struct.pack('<II', vertices.size // 3, edges.size // 2))
  result.write(vertices.tobytes('C'))
  result.write(edges.tobytes('C'))

  def writeattr(attr, dtype, text):
    if attr is None:
      return

    attr = attr.astype(dtype)

    if attr.shape[0] != vertices.shape[0]:
      raise SkeletonEncodeError("Number of {} {} ({}) must match the number of vertices ({}).".format(
        dtype, text, attr.shape[0], vertices.shape[0]
      ))
    
    result.write(attr.tobytes('C'))

  for attr in skeleton.extra_attributes:
    arr = getattr(skeleton, attr['id'])
    writeattr(arr, np.dtype(attr['data_type']), attr['id'])

  return result.getvalue()

def from_precomputed(
  skelbuf:bytes, 
  segid:Optional[int] = None, 
  vertex_attributes:Optional[dict[str, Any]] = None,
) -> "Skeleton":
  """
  Convert a buffer into a Skeleton object.

  Format:
  num vertices (Nv) (uint32)
  num edges (Ne) (uint32)
  XYZ x Nv (float32)
  edge x Ne (2x uint32)

  Default Vertex Attributes:

    radii x Nv (optional, float32)
    vertex_type x Nv (optional, req radii, uint8) (SWC definition)

  Specify your own:

  vertex_attributes = [
    {
      'id': name of attribute,
      'num_components': int,
      'data_type': dtype,
    },
  ]

  More documentation: 
  https://github.com/seung-lab/cloud-volume/wiki/Advanced-Topic:-Skeletons-and-Point-Clouds
  """
  from .. import Skeleton
  
  if len(skelbuf) < 8:
    raise SkeletonDecodeError("{} bytes is fewer than needed to specify the number of verices and edges.".format(len(skelbuf)))

  num_vertices, num_edges = struct.unpack('<II', skelbuf[:8])
  min_format_length = 8 + 12 * num_vertices + 8 * num_edges

  if len(skelbuf) < min_format_length:
    raise SkeletonDecodeError("The input skeleton was {} bytes but the format requires {} bytes.".format(
      len(skelbuf), min_format_length
    ))

  vstart = 2 * 4 # two uint32s in
  vend = vstart + num_vertices * 3 * 4 # float32s
  vertbuf = skelbuf[ vstart : vend ]

  estart = vend
  eend = estart + num_edges * 4 * 2 # 2x uint32s

  edgebuf = skelbuf[ estart : eend ]

  vertices = np.frombuffer(vertbuf, dtype='<f4').reshape( (num_vertices, 3) )
  edges = np.frombuffer(edgebuf, dtype='<u4').reshape( (num_edges, 2) )

  skeleton = Skeleton(vertices, edges, segid=segid)

  if len(skelbuf) == min_format_length:
    return skeleton

  if vertex_attributes is None:
    vertex_attributes = kls._default_attributes()

  start = eend
  end = -1
  for attr in vertex_attributes:
    num_components = int(attr['num_components'])
    data_type = np.dtype(attr['data_type'])
    end = start + num_vertices * num_components * data_type.itemsize
    attrbuf = np.frombuffer(skelbuf[start : end], dtype=data_type)

    if num_components > 1:
      attrbuf = attrbuf.reshape( (num_vertices, num_components) )

    setattr(skeleton, attr['id'], attrbuf)
    start = end

  skeleton.extra_attributes = vertex_attributes

  return skeleton