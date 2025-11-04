from typing import Optional, Any, Literal

from .types import (
  CompressionType,
  DataType,
  EdgeRepresentation,
  GraphType,
  PhysicalUnit,
  SpaceType,
  TO_DATATYPE, 
  FROM_DATATYPE,
)

import uuid

import numpy as np
import numpy.typing as npt

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
  HEADER_BYTES = 81

  def __init__(
    self,
    Nv:int,
    Ne:int,
    append_mode:bool = False,
    attribute_header_bytes:int = 0,
    crc16:Optional[int] = None,
    edge_datatype:DataType = DataType.U32,
    edge_compression:CompressionType = CompressionType.NONE,
    edge_representation:EdgeRepresentation = EdgeRepresentation.PAIR,
    edge_bytes:int = 0,
    format_version:int = 0,
    graph_type:GraphType = GraphType.GRAPH,
    has_transform:bool = True,
    id:Optional[int] = None,
    num_axes:Literal[2,3] = 3,
    num_components:Optional[int] = None,
    physical_unit:Union[str, PhysicalUnit] = PhysicalUnit.NANOMETER,
    space:Union[str, SpaceType] = SpaceType.VOXEL,
    spatial_index_bytes:int = 0,
    total_bytes:int = 0,
    vertex_compression:CompressionType = CompressionType.NONE,
    vertex_datatype:DataType = DataType.F32,
    vertex_bytes:int = 0,
  ):
    self.Nv = int(Nv)
    self.Ne = int(Ne)
    self.space = space
    if isinstance(space, str) and space == "physical":
      self.space = SpaceType.PHYSICAL

    self.append_mode = append_mode

    self.vertex_datatype = vertex_datatype
    self.edge_dtype = edge_dtype
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
    return np.dtype(FROM_DATATYPE[self.edge_datatype])

  @property
  def Nvb(self):
    return self.Nv * self.num_axes * self.vertex_dtype.itemsize

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
  compression: CompressionType
  content_length: int

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

