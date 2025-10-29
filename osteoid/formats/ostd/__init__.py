import numpy as np

from .header import (
	OstdHeader, OstdAttribute,
	TO_DATATYPE, FROM_DATATYPE,
	EdgeRepresentation, GraphType, 
	CompressionType, PhysicalUnit,
	SpaceType,
)

def to_ostd(
  skel:"Skeleton", 
  physical_unit:str = "nm",
) -> bytes:

  transform = np.eye(4)
  transform[:3,:] = skel.transform
  transform = transform.astype(np.float32, copy=False)

  header = OstdHeader(
    Nv=int(skel.vertices.shape[0]),
    Ne=int(skel.vertices.shape[0]),
    Nattr=len(skel.extra_attributes),
    attr_name_width=max([ len(attr["name"]) for attr in skel.extra_attributes ]),
    physical_unit=physical_unit,
    id=skel.id,
    vertex_dtype=TO_DATATYPE[skel.vertices.dtype],
    edge_representation=EdgeRepresentation.PAIR,
    graph_type=GraphType.GRAPH,
    compression=CompressionType.NONE,
    num_components=len(skel.components()),
    transform=transform,
  )

  vertices = skel.vertices.tobytes("F")
  vcrc = lib.crc32c(vertices)

  edges = skel.edges.astype(header.edge_dtype, copy=False).tobytes("F")
  ecrc = lib.crc32c(edges)

  skeleton_binary = b''.join([
    header.to_bytes(),
    vertices, vcrc,
    edges, ecrc,
  ])

  attr_binaries = []
  attr_components = []
  offset = len(skeleton_binary)
  for attr in skel.extra_attributes:
    attr_binary = getattr(skel, attr["id"]).tobytes("F")
    attr_binary_crc = lib.crc32c(attr_binary)
    attr_binary = b''.join([ attr_binary, attr_binary_crc ])
    attr_binaries.append(attr_binary)

    name = bytearray(attr_name_width)
    name[:len(attr["id"])] = attr["id"]
    dtype = np.uint8(TO_DATATYPE[np.dtype(attr["data_type"]).type])
    type_field = dtype | (np.uint8(attr["num_components"]) << 4)
    attr_components.append(b''.join([ name, type_field, offset.to_bytes(8, 'little') ]))

    offset += len(attr_binary)

  attr_header = b''.join(attr_components)

  final_binary = b''.join([
    skeleton_binary,
    *attr_binaries,
    attr_header,
    lib.crc16(attr_header),
  ])
  whole_file_crc = lib.crc32c(final_binary)
  final_binary += whole_file_crc
  return final_binary

def from_ostd(binary:bytes) -> "Skeleton":
  from .. import Skeleton
  
  header = OstdHeader.from_bytes(binary)

  dtype = FROM_DATATYPE[header.vertex_dtype]
  Nvb = header.Nvb
  hb = OstdHeader.HEADER_BYTES

  offset = hb + Nvb
  vert_stored_crc32c = int.from_bytes(binary[offset:offset+4], 'little')
  vert_computed_crc32 = lib.crc32c(binary[hb:hb+Nvb])

  if vert_stored_crc32c != vert_computed_crc32:
    raise ValueError(f"vertex corruption detected. stored crc32c: {vert_stored_crc32c} computed crc32c: {vert_computed_crc32}")

  vertices = np.frombuffer(
    binary, 
    dtype=dtype, 
    offset=HEADER_LENGTH,
    count=Nv * 3,
  )
  vertices = vertices.reshape((Nv,3), order="C")

  edge_dtype = find_edge_dtype(Nv)
  num_edge_bytes = Ne * 2 * np.dtype(edge_dtype).itemsize

  computed_edge_crc32c = lib.crc32c(binary[offset:offset + num_edge_bytes ])

  if flags["edge_repr"] == EdgeRepresentation.PAIR:
    edges = np.frombuffer(
      binary,
      dtype=,
      offset=offset,
      count=(Ne * 2)
    )
    edges = edges.reshape((Ne, 2), order='C')
  else: # PARENT
    raise ValueError("PARENT format is not yet supported.")

  offset += num_edge_bytes
  stored_edge_crc32c = read_int(4)

  if stored_edge_crc32c != computed_edge_crc32c:
    raise ValueError(f"edge corruption detected. stored crc32c: {stored_edge_crc32c} computed crc32c: {computed_edge_crc32c}")

  if flags["space"] == SpaceType.VOXEL:
    space = 'voxel'
  else:
    space = 'physical'

  attr_width = attr_name_width + 1 + 8

  attr_bytes = num_attr * attr_width
  attributes_binary = binary[-(attr_bytes + 2):-2]
  stored_attr_crc16 = int.from_bytes(binary[-2:], 'little')
  computed_attr_crc16 = lib.crc16(attributes_binary)

  if stored_attr_crc16 != computed_attr_crc16:
    raise ValueError(f"attribute corruption detected. stored crc16: {stored_attr_crc16} computed crc16: {computed_attr_crc16}")

  skel = Skeleton(
    vertices=vertices,
    edges=edges,
    segid=skel_id,
    space=space,
    transform=transform[:3,:],
  )

  extra_attributes = []

  for i in range(num_attr):
    single_attr_binary = attributes_binary[i * attr_width: (i+1) * attr_width]
    name = single_attr_binary[:attr_name_width]
    name = name.split(b'\x00', 1)[0].decode('utf-8')
    field = single_attr_binary[attr_name_width]
    dtype = DataType(field & 0b1111)
    dtype = np.dtype(FROM_DATATYPE[dtype])
    num_components = int((field >> 4) & 0b1111)

    offset = int.from_bytes(single_attr_binary[-8:], 'little')
    attr = np.frombuffer(
      binary, 
      dtype=dtype, 
      offset=offset, 
      count=(num_components * Nv)
    )

    offset += attr.nbytes
    stored_attr_crc32c = int.from_bytes(binary[offset:offset+4], 'little')
    computed_attr_crc32c = lib.crc32c(attr)

    if stored_attr_crc32c != computed_attr_crc32c:
      raise ValueError(f"attribute {name} was corrupted. stored crc32c: {stored_attr_crc32c} computed crc32c: {computed_attr_crc32c}")

    attr = attr.reshape((Nv, num_components), order="C")
    setattr(skel, name, attr)
    extra_attributes.append({
      "id": name,
      "data_type": dtype.name,
      "num_components": num_components,
    })

  return skel

