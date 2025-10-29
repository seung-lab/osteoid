from typing import Optional, Any

from enum import IntEnum
import uuid 

import numpy as np
import numpy.typing as npt

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
  BOOLEAN = 12
  PACKED_BOOLEAN = 13

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
  GRAPH = 0
  TREE = 1

class SpaceType(IntEnum):
  VOXEL = 0
  PHYSICAL = 1

def find_edge_dtype(num_verts:int) -> np.dtype:
  if num_verts < np.dtype(np.uint8).max:
    return np.uint8
  elif num_verts < np.dtype(np.uint16).max:
    return np.uint16
  elif num_verts < np.dtype(np.uint32).max:
    return np.uint32
  else:
    return np.uint64

class OstdHeader:
  MAGIC = b'ostd'
  FORMAT_VERSION = 0
  HEADER_BYTES = 114
  # SCHEMA = [
  #   ['magic', 4 ],
  #   ['format_version', 1 ],
  #   ['id', 16 ],
  #   ['flags', 4 ],
  #   ['Nv', 8 ], # num vertices
  #   ['Ne', 8 ], # num edges
  #   ['Nattr', 2 ],
  #   ['attr_name_width', 1 ],
  #   ['num_components', 4 ],
  #   ['transform', 4*4*4 ],
  #   ['crc16', 2 ],
  # ]

  def __init__(
    self,
    Nv:int,
    Ne:int,
    space:Union[str, SpaceType] = SpaceType.VOXEL,
    vertex_datatype:DataType = DataType.F32,
    edge_representation:EdgeRepresentation = EdgeRepresentation.PAIR,
    graph_type:GraphType = GraphType.GRAPH,
    physical_unit:Union[str, PhysicalUnit] = PhysicalUnit.NANOMETER,
    compression:CompressionType = CompressionType.NONE,
    crc16:Optional[int] = None,
    transform:npt.NDArray[np.float32] = np.eye(4),
    Nattr:int = 0,
    attr_name_width:int = 1,
    num_components:Optional[int] = None,
    format_version:int = 0,
    id:Optional[int] = None,
  ):
    self.Nv = int(Nv)
    self.Ne = int(Ne)
    self.space = space
    if isinstance(space, str) and space == "physical":
      self.space = SpaceType.PHYSICAL

    self.vertex_datatype = vertex_datatype
    self.edge_dtype = find_edge_dtype(Nv)
    self.edge_representation = edge_representation
    self.graph_type = graph_type

    if isinstance(physical_unit, str):
      self.physical_unit = TO_UNIT[physical_unit]
    else:
      self.physical_unit = physical_unit

    self.compression = compression
    self.transform = np.asfortranarray(transform)
    self.num_components = int(num_components)

    self.format_version = format_version

    assert 0 < attr_name_width <= 255
    assert Nattr < (1 << 15)

    self.attr_name_width = attr_name_width
    self.Nattr = Nattr

    if id is None:
      self.id = uuid.uuid4().bytes
      self.id = int.from_bytes(id, 'little')
    elif isinstance(id, (bytes, bytearray)):
      self.id = int.from_bytes(id, 'little')
    else:
      self.id = int(id)

    self.crc16 = crc16

  @property
  def vertex_dtype(self):
    return np.dtype(FROM_DATATYPE[self.vertex_datatype])

  @property
  def edge_dtype(self):
    if self.Nv < np.dtype(np.uint8).max:
      return np.uint8
    elif self.Nv < np.dtype(np.uint16).max:
      return np.uint16
    elif self.Nv < np.dtype(np.uint32).max:
      return np.uint32
    else:
      return np.uint64

  @property
  def Nvb(self):
    return self.Nv * 3 * self.vertex_dtype.itemsize

  def encode_flags(self) -> int:
    flags = np.uint32(0)
    flags |= np.uint32(self.vertex_datatype.value & 0b1111) # 4 bits
    flags |= np.uint32((self.edge_representation.value & 0b1) << 4) # 1 bit
    flags |= np.uint32((self.compression.value & 0b1111) << 5) # 4 bits
    flags |= np.uint32((self.graph_type.value & 0b111) << 9) # 3 bits
    flags |= np.uint32((self.physical_unit.value & 0b11111) << 12) # 5 bits
    flags |= np.uint32((self.space.value & 0b11) << 17) # 2 bits
    return flags

  def decode_flags(self, flags:int):
    self.vertex_datatype = DataType(flags & 0b1111)
    self.edge_repr = EdgeRepresentation((flags & 0b10000) >> 4)
    self.compression = CompressionType((flags & 0b111100000) >> 5)
    self.graph_type = GraphType((flags & 0b1110000000000) >> 9)
    self.unit = PhysicalUnit((flags & 0b11111000000000000) >> 12)
    self.space = SpaceType((flags & 0b1100000000000000000) >> 17)

  @classmethod
  def validate_header(kls, binary:bytes, ):
    if len(binary) < OstdHeader.HEADER_BYTES:
      raise ValueError(f"buffer is too short to be an osteoid file. buffer: {len(binary)} bytes, needed at least: {OstdHeader.HEADER_BYTES} bytes")

    magic = binary[:4]
    if magic != OstdHeader.MAGIC:
      raise ValueError(f"magic numbers did not match. Expected: {OstdHeader.MAGIC}, Got: {magic}")

    format_version = binary[4]
    if format_version > OstdHeader.FORMAT_VERSION:
      raise ValueError(f"format version in buffer exceeded maximum supported version. Got: {format_version} Maximum Supported Version: {OstdHeader.FORMAT_VERSION}")

    crc_start, crc_end = HEADER_OFFSETS['crc16']
    stored_crc16 = int.from_bytes(binary[crc_start:crc_end])
    computed_crc16 = lib.crc16(binary[4:crc_start])

    if stored_crc16 != computed_crc16:
      raise ValueError(f"Header corruption detected. Stored CRC16: {stored_crc16}, Computed CRC16: {computed_crc16}")

  @classmethod
  def from_bytes(kls, binary:bytes, crc_check:bool = True) -> "OstdHeader":
    OstdHeader.validate_header(binary)

    offset = len(OstdHeader.MAGIC) + 1 # magic + format_version

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

    header = OstdHeader(Nv, Ne)
    header.decode_flags(flags)

    transform = np.frombuffer(binary, dtype=np.float32, offset=offset, count=4 * 4)
    header.transform = transform.reshape((4,4), order="F")

  def to_bytes(self) -> bytes:
    header = b''.join([
      self.format_version.to_bytes(1),
      int(self.id).to_bytes(16, 'little'),
      self.encode_flags().tobytes(),
      int(self.Nv).to_bytes(8, 'little'),
      int(self.Ne).to_bytes(8, 'little'),
      int(self.Nattr).to_bytes(2, 'little'),
      int(self.attr_name_width).to_bytes(1, 'little'),
      int(self.num_components).to_bytes(4, 'little'),
      self.transform.tobytes("F"),
    ])
    header_crc16 = lib.crc16(header)

    return b''.join([
      OstdHeader.MAGIC,
      header,
      header_crc16,
    ])

@dataclass
class OstdAttribute:
  id: str
  dtype: str
  num_components: int

  def to_bytes(self, name_width:int):
    if self.num_components >= 16:
      raise ValueError(f"Only 15 components supported per an attribute. Got {self.num_components} for {self.id}.")

    name = bytearray(name_width)
    name[:len(self.id)] = self.id
    dtype = np.uint8(TO_DATATYPE[np.dtype(self.dtype).type])
    type_field = dtype | (np.uint8(self.num_components) << 4)
    return b''.join([ name, type_field ])

  def to_dict(self) -> dict[str,Any]:
    return {
      "id": self.id,
      "data_type": np.dtype(self.dtype).name,
      "num_components": self.num_components,
    }

  @classmethod
  def from_bytes(self, binary:bytes, name_width:int) -> "OstdAttribute":
    name = binary[:name_width]
    name = name.split(b'\x00', 1)[0].decode('utf-8')
    field = binary[name_width]
    dtype = DataType(field & 0b1111)
    dtype = np.dtype(FROM_DATATYPE[dtype])
    num_components = int((field >> 4) & 0b1111)
    return OstdAttribute(id=name, dtype=dtype, num_components=num_components)

