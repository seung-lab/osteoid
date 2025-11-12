from enum import IntEnum
import uuid 

import numpy as np

class PhysicalUnit(IntEnum):
  VOXEL = 0

  ANGSTROM = 1
  FEMTOMETER = 2
  PICOMETER = 3
  NANOMETER = 4
  MICROMETER = 5
  MILLIMETER = 6
  CENTIMETER = 7
  METER = 8
  KILOMETER = 9
  MEGAMETER = 10
  LIGHTYEAR = 11
  PARSEC = 12

  MIL = 13
  INCH = 14
  FOOT = 15
  YARD = 16
  STATUTE_MILE = 17
  NAUTICAL_MILE = 18

TO_UNIT = {
  "vx": PhysicalUnit.VOXEL,
  "voxel": PhysicalUnit.VOXEL,

  "A": PhysicalUnit.ANGSTROM,
  "angstrom": PhysicalUnit.ANGSTROM,
  
  "fm": PhysicalUnit.FEMTOMETER,
  "femtometer": PhysicalUnit.FEMTOMETER,
  
  "pm": PhysicalUnit.PICOMETER,
  "picometer": PhysicalUnit.PICOMETER,
  
  "nm": PhysicalUnit.NANOMETER,
  "nanometer": PhysicalUnit.NANOMETER,
  
  "um": PhysicalUnit.MICROMETER,
  "micrometer": PhysicalUnit.MICROMETER,
  "micron": PhysicalUnit.MICROMETER,
  
  "mm": PhysicalUnit.MILLIMETER,
  "millimeter": PhysicalUnit.MILLIMETER,

  "cm": PhysicalUnit.CENTIMETER,
  "centimeter": PhysicalUnit.CENTIMETER,
  
  "m": PhysicalUnit.METER,
  "meter": PhysicalUnit.METER,
  
  "km": PhysicalUnit.KILOMETER,
  "kilometer": PhysicalUnit.KILOMETER,

  "Mm": PhysicalUnit.MEGAMETER,
  
  "ly": PhysicalUnit.LIGHTYEAR,
  "lightyear": PhysicalUnit.LIGHTYEAR,

  "pc": PhysicalUnit.PARSEC,
  "parsec": PhysicalUnit.PARSEC,

  "mil": PhysicalUnit.MIL,
  
  "in": PhysicalUnit.INCH,
  "inch": PhysicalUnit.INCH,
  "inches": PhysicalUnit.INCH,
  
  "ft": PhysicalUnit.FOOT,
  "foot": PhysicalUnit.FOOT,
  "feet": PhysicalUnit.FOOT,

  "yd": PhysicalUnit.YARD,
  "yard": PhysicalUnit.YARD,
  
  "mi": PhysicalUnit.STATUTE_MILE,
  "mile": PhysicalUnit.STATUTE_MILE,
  
  "nmi": PhysicalUnit.NAUTICAL_MILE,
}

class CompressionType(IntEnum):
  NONE = 0
  GZIP = 1
  BZIP2 = 2
  ZSTD = 3

class DataType(IntEnum):
  F8 = 0
  F16 = 1
  F32 = 2
  F64 = 3
  U8 = 4
  U16 = 5
  U32 = 6
  U64 = 7
  I8 = 8
  I16 = 9
  I32 = 10
  I64 = 11

TO_DATATYPE = {
  np.float16: DataType.F16,
  np.float32: DataType.F32,
  np.float64: DataType.F64,
  np.uint8: DataType.U8,
  np.uint16: DataType.U16,
  np.uint32: DataType.U32,
  np.uint64: DataType.U64,
  np.int8: DataType.I8,
  np.int16: DataType.I16,
  np.int32: DataType.I32,
  np.int64: DataType.I64,
}

FROM_DATATYPE = { v,k for k,v in TO_DATATYPE.items() }

class EdgeRepresentation(IntEnum):
  PAIR = 0
  PARENT = 1

class GraphType(IntEnum):
  UNKNOWN = 0
  GRAPH = 1
  TREE = 2

class SpaceType(IntEnum):
  VOXEL = 0
  PHYSICAL = 1

FORMAT_VERSION = 0

MAGIC_NUMBER = b'ostd'

HEADER_BYTE_SCHEMA = [
  ['magic', len(MAGIC_NUMBER) ],
  ['format_version', 1 ],
  ['id', 16 ],
  ['flags', 4 ],
  ['Nv', 8 ], # num vertices
  ['Ne', 8 ], # num edges
  ['Nattr', 2 ],
  ['attr_name_width', 1 ],
  ['num_components', 4 ],
  ['transform', 4*4*4 ],
  ['crc16', 2 ],
]

HEADER_OFFSETS = {}
HEADER_LENGTH = 0
for name, amt in HEADER_BYTE_SCHEMA:
  HEADER_OFFSETS[name] = (HEADER_LENGTH, HEADER_LENGTH + amt)
  HEADER_LENGTH += amt

def parse_flags(flags:int) -> dict[str,IntEnum]:
  return {
    'vertex_dtype': DataType(flags & 0b1111),
    'edge_repr': EdgeRepresentation((flags & 0b10000) >> 4),
    'compression': CompressionType((flags & 0b111100000) >> 5),
    'graph_type': GraphType((flags & 0b1110000000000) >> 9),
    'unit': PhysicalUnit((flags & 0b11111000000000000) >> 12),
    'space': SpaceType((flags & 0b1100000000000000000) >> 17),
  }

def find_edge_dtype(num_verts:int) -> int:
  if num_verts < np.dtype(np.uint8).max:
    return np.uint8
  elif num_verts < np.dtype(np.uint16).max:
    return np.uint16
  elif num_verts < np.dtype(np.uint32).max:
    return np.uint32
  else:
    return np.uint64

def to_ostd(
  skel:"Skeleton", 
  physical_unit:str = "nm",
) -> bytes:
  physical_unit = physical_unit or "nm"
  if physical_unit[-1] == "s":
    physical_unit = physical_unit[:-1]
  physical_unit = TO_UNIT[physical_unit]

  skeleton_id = int(skel.id).to_bytes(16, 'little')
  if skeleton_id is None:
    skeleton_id = uuid.uuid4().bytes

  space = SpaceType.VOXEL
  if skel.space == "physical":
    space = SpaceType.PHYSICAL

  flags = np.uint32(0)
  flags |= np.uint32(TO_DATATYPE[skel.vertices.dtype]) # 4 bits
  flags |= np.uint32(EdgeRepresentation.PAIR << 4) # 1 bit
  flags |= np.uint32(CompressionType.NONE << 5) # 4 bits
  flags |= np.uint32(GraphType.UNKNOWN << 9) # 3 bits
  flags |= np.uint32(physical_unit << 12) # 5 bits
  flags |= np.uint32(space << 17) # 2 bits

  Nv = int(skel.vertices.shape[0])
  Ne = int(skel.vertices.shape[0])
  # edge_data_width = find_edge_data_width(Nv)

  num_attributes = len(skel.extra_attributes)
  num_components = len(skel.components())

  transform = np.eye(4)
  transform[:3,:] = skel.transform
  transform = transform.astype(np.float32, copy=False)

  attr_name_width = 0
  if len(skel.extra_attributes):
    attr_name_width = max([ len(attr["name"]) for attr in skel.extra_attributes ])

  if attr_name_width > 255:
    raise ValueError("longest attribute name is longer than 255 chars")

  header = b''.join([
    FORMAT_VERSION.to_bytes(1),
    skeleton_id,
    flags.tobytes(),
    int(Nv).to_bytes(8, 'little'),
    int(Ne).to_bytes(8, 'little'),
    int(num_attributes).to_bytes(2, 'little'),
    int(attr_name_width).to_bytes(1, 'little'),
    int(num_components).to_bytes(4, 'little'),
    transform.tobytes("F"),
  ])
  header_crc16 = lib.crc16(header)

  vertices = skel.vertices.tobytes("F")
  vcrc = lib.crc32c(vertices)

  edges = skel.edges.tobytes("F")
  ecrc = lib.crc32c(edges)

  skeleton_binary = b''.join([
    MAGIC_NUMBER,
    header,
    header_crc16,
    vertices,
    vcrc,
    edges,
    ecrc,
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

def validate_ostd(binary:bytes):
  if len(binary) < HEADER_LENGTH:
    raise ValueError(f"buffer is too short to be an osteoid file. buffer: {len(binary)} bytes, needed at least: {HEADER_LENGTH} bytes")

  magic = binary[:4]
  if magic != MAGIC_NUMBER:
    raise ValueError(f"magic numbers did not match. Expected: {MAGIC_NUMBER}, Got: {magic}")

  format_version = binary[4]
  if format_version > FORMAT_VERSION:
    raise ValueError(f"format version in buffer exceeded maximum supported version. Got: {format_version} Maximum Supported Version: {FORMAT_VERSION}")

  crc_start, crc_end = HEADER_OFFSETS['crc16']
  stored_crc16 = int.from_bytes(binary[crc_start:crc_end])
  computed_crc16 = lib.crc16(binary[4:crc_start])

  if stored_crc16 != computed_crc16:
    raise ValueError(f"Header corruption detected. Stored CRC16: {stored_crc16}, Computed CRC16: {computed_crc16}")

def from_ostd(binary:bytes) -> "Skeleton":
  from .. import Skeleton
  
  validate_ostd(binary)

  offset = 5

  def read_int(N):
    nonlocal offset
    x = int.from_bytes(binary[offset:offset+N], 'little')
    offset += N
    return x

  skel_id = read_int(16)
  flags = read_int(2)
  Nv = read_int(8)
  Ne = read_int(8)
  num_attr = read_int(2)
  attr_name_width = read_int(1)
  num_components = read_int(4)

  flags = parse_flags(flags)

  transform = np.frombuffer(binary, dtype=np.float32, offset=offset, count=4 * 4)
  transform = transform.reshape((4,4), order="F")
  offset += transform.nbytes

  dtype = flags["vertex_dtype"]
  dtype = FROM_DATATYPE[dtype]

  num_vertex_bytes = Nv * 3 * np.dtype(dtype).itemsize

  offset = HEADER_LENGTH + num_vertex_bytes
  vert_stored_crc32c = read_int(4)
  vert_computed_crc32 = lib.crc32c(binary[HEADER_LENGTH:HEADER_LENGTH+num_vertex_bytes])

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

