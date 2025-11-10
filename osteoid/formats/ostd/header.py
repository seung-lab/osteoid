from dataclasses import dataclass
from typing import Optional, Any, Literal, Union

from enum import IntEnum

from ... import lib

from .types import (
  AreaType,
  AttributeType,
  AxisPermutationType,
  CompressionType,
  CurrentSpaceType,
  DataType,
  EdgeRepresentationType,
  ElectricalType,
  EnergyType,
  GraphType,
  LengthType,
  MassType,
  SIPrefixType,
  SpaceType,
  TemperatureType,
  TimeType,
  VolumeType,
  TO_LENGTH_UNIT,
  TO_DATATYPE, 
  FROM_DATATYPE,
  TO_AXIS_PERMUTATION,
  QUANTITY_TYPE,
)

import uuid
import struct

import numpy as np
import numpy.typing as npt

@dataclass
class CoordinateFrame:
  sign_x:bool
  sign_y:bool
  sign_z:bool
  permutation:AxisPermutationType

def parse_coordinate_frame_orientation(orientation:str) -> AxisPermutationType:
  if len(orientation) > 6:
    raise ValueError(f"Unable to parse orientation: {orientation[:100]}")

  orientation = orientation.upper()
  normalized = orientation.replace('+', '').replace('-', '')

  if not (2 <= len(normalized) <= 3):
    raise ValueError(f"Unable to parse orientation: {normalized}")

  POSITIVE = 0
  NEGATIVE = 1

  signs = [ POSITIVE, POSITIVE, POSITIVE ]
  mapping = { "X": 0, "Y": 1, "Z": 2 }

  for i in range(len(orientation) - 1):
    if orientation[i] == "-":
      signs[mapping[orientation[i+1]]] = NEGATIVE

  permutation = TO_AXIS_PERMUTATION[normalized]

  return CoordinateFrame(*signs, permutation)

class OstdHeader:
  MAGIC = b'ostd'
  FORMAT_VERSION = 0
  HEADER_BYTES = 88

  def __init__(
    self,
    Nv:int,
    Ne:int,
    append_mode:bool = False,
    attribute_header_bytes:int = 0,
    coordinate_frame_orientation:str = '+X+Y+Z',
    crc16:Optional[int] = None,
    edge_data_type:DataType = DataType.U32,
    edge_compression:CompressionType = CompressionType.NONE,
    edge_representation:EdgeRepresentationType = EdgeRepresentationType.PAIR,
    edge_bytes:int = 0,
    format_version:int = 0,
    graph_type:GraphType = GraphType.GRAPH,
    has_transform:bool = True,
    id:Optional[int] = None,
    num_axes:Literal[2,3] = 3,
    num_components:int = np.iinfo(np.uint32).max,
    length_unit:Union[str, tuple[SIPrefixType, IntEnum]] = (SIPrefixType.NANO, LengthType.METER),
    physical_path_length:float = float('NaN'),
    space:int = 0,
    spatial_index_bytes:int = 0,
    total_bytes:int = 0,
    vertex_compression:CompressionType = CompressionType.NONE,
    vertex_data_type:DataType = DataType.F32,
    vertex_bytes:int = 0,
    voxel_centered:bool = True,
  ):
    self.Nv = int(Nv)
    self.Ne = int(Ne)
    self.space = space
    self.append_mode = append_mode

    self.vertex_data_type = vertex_data_type
    self.edge_data_type = edge_data_type
    self.edge_representation = edge_representation
    self.graph_type = graph_type

    if isinstance(length_unit, str):
      self.length_unit = TO_LENGTH_UNIT[length_unit]
    else:
      self.length_unit = length_unit

    if isinstance(coordinate_frame_orientation, str):
      self.coordinate_frame_orientation = parse_coordinate_frame_orientation(coordinate_frame_orientation)
    else:
      self.coordinate_frame_orientation = coordinate_frame_orientation

    self.voxel_centered = voxel_centered
    self.vertex_compression = vertex_compression
    self.edge_compression = edge_compression
    self.has_transform = bool(has_transform)
    self.num_axes = num_axes
    self.num_components = int(num_components)
    self.physical_path_length = float(physical_path_length)

    self.vertex_bytes = vertex_bytes
    self.edge_bytes = edge_bytes
    self.spatial_index_bytes = spatial_index_bytes
    self.attribute_header_bytes = attribute_header_bytes
    self.total_bytes = total_bytes
    self.format_version = format_version

    if id is None:
      self.id = uuid.uuid4().bytes
      self.id = int.from_bytes(self.id, 'little')
    elif isinstance(id, (bytes, bytearray)):
      self.id = int.from_bytes(id, 'little')
    else:
      self.id = int(id)

    self.crc16 = crc16

  @property
  def vertex_dtype(self):
    return np.dtype(FROM_DATATYPE[self.vertex_data_type])

  @property
  def edge_dtype(self):
    return np.dtype(FROM_DATATYPE[self.edge_data_type])

  @property
  def Nvb(self):
    return self.Nv * self.num_axes * self.vertex_dtype.itemsize

  def encode_flags(self) -> int:
    flags = np.uint64(0)
    offset = 0
    def write_int(value:int, nbits:int):
      nonlocal flags
      nonlocal offset
      mask = (1 << nbits) - 1
      flags |= np.uint64(value & mask) << offset
      offset += nbits

    write_int(self.vertex_data_type.value, 4)
    write_int(self.edge_data_type.value, 4)

    write_int(self.vertex_compression.value, 4)
    write_int(self.edge_compression.value, 4)

    write_int(self.graph_type.value, 3)
    write_int(self.length_unit[0].value, 5)
    write_int(self.length_unit[1].value, 4)

    write_int(self.coordinate_frame_orientation.sign_x, 1)
    write_int(self.coordinate_frame_orientation.sign_y, 1)
    write_int(self.coordinate_frame_orientation.sign_z, 1)
    write_int(self.coordinate_frame_orientation.permutation.value, 3)

    write_int(self.num_axes, 3)
    write_int(bool(self.has_transform), 1)
    write_int(self.voxel_centered, 1)
    write_int(self.edge_representation.value, 2)

    return flags

  def decode_flags(self, flags:int):
    offset = 0
    def read_int(nbits):
      nonlocal offset
      mask = (1 << nbits) - 1
      res = (flags >> offset) & mask
      offset += nbits
      return res

    self.vertex_data_type = DataType(read_int(4))
    self.edge_data_type = DataType(read_int(4))
    self.vertex_compression = CompressionType(read_int(4))
    self.edge_compression = CompressionType(read_int(4))

    self.graph_type = GraphType(read_int(3))
    self.length_unit = (SIPrefixType(read_int(5)), LengthType(read_int(4)))

    self.coordinate_frame_orientation = CoordinateFrame(
      sign_x=bool(read_int(1)),
      sign_y=bool(read_int(1)),
      sign_z=bool(read_int(1)),
      permutation=AxisPermutationType(read_int(3)),
    )

    self.num_axes = read_int(3)

    self.has_transform = bool(read_int(1))
    self.voxel_centered = bool(read_int(1))
    self.append_mode = bool(read_int(1))
    self.edge_representation = EdgeRepresentationType(read_int(2))

  @classmethod
  def validate_header(kls, binary:bytes):
    if len(binary) < OstdHeader.HEADER_BYTES:
      raise ValueError(f"buffer is too short to be an osteoid file. buffer: {len(binary)} bytes, needed at least: {OstdHeader.HEADER_BYTES} bytes")

    magic = binary[:4]
    if magic != OstdHeader.MAGIC:
      raise ValueError(f"magic numbers did not match. Expected: {OstdHeader.MAGIC}, Got: {magic}")

    format_version = binary[4]
    if format_version > OstdHeader.FORMAT_VERSION:
      raise ValueError(f"format version in buffer exceeded maximum supported version. Got: {format_version} Maximum Supported Version: {OstdHeader.FORMAT_VERSION}")

    total_bytes = int.from_bytes(binary[5:13])
    if len(binary) < total_bytes:
      raise ValueError("Stream contained too few bytes.")

    # hb = OstdHeader.HEADER_BYTES
    # stored_crc16 = int.from_bytes(binary[hb-2:hb])
    # computed_crc16 = lib.crc16(binary[len(OstdHeader.MAGIC):hb-2])

    # if stored_crc16 != computed_crc16:
    #   raise ValueError(f"Header corruption detected. Stored CRC16: {stored_crc16}, Computed CRC16: {computed_crc16}")

  @classmethod
  def from_bytes(kls, binary:bytes, crc_check:bool = True) -> "OstdHeader":
    OstdHeader.validate_header(binary)

    offset = len(OstdHeader.MAGIC)

    def read_int(N):
      nonlocal offset
      x = int.from_bytes(binary[offset:offset+N], 'little')
      offset += N
      return x

    def read_float():
      nonlocal offset
      N = 4
      x = struct.unpack('f', binary[offset:offset+N])[0]
      offset += N
      return x 

    format_version = read_int(1)
    total_bytes = read_int(8)
    skel_id = read_int(16)
    flags = read_int(8)
    Nv = read_int(8)
    Ne = read_int(8)
    vertex_bytes = read_int(8)
    edge_bytes = read_int(8)
    spatial_index_bytes = read_int(4)
    attribute_header_bytes = read_int(4)
    num_components = read_int(4)
    physical_path_length = read_float()
    current_space = read_int(1)
    crc16 = read_int(2)

    header = OstdHeader(
      Nv, Ne,
      id=skel_id,
      total_bytes=total_bytes,
      vertex_bytes=vertex_bytes,
      edge_bytes=edge_bytes,
      spatial_index_bytes=spatial_index_bytes,
      attribute_header_bytes=attribute_header_bytes,
      num_components=num_components,
      physical_path_length=physical_path_length,
      space=current_space,
      format_version=format_version,
      crc16=crc16,
    )
    header.decode_flags(flags)

    return header

  def to_bytes(self) -> bytes:
    header = b''.join([
      int(self.format_version).to_bytes(1),
      int(self.total_bytes).to_bytes(8),
      int(self.id).to_bytes(16, 'little'),
      self.encode_flags().tobytes(),
      int(self.Nv).to_bytes(8, 'little'),
      int(self.Ne).to_bytes(8, 'little'),
      int(self.vertex_bytes).to_bytes(8, 'little'),
      int(self.edge_bytes).to_bytes(8, 'little'),
      int(self.spatial_index_bytes).to_bytes(4, 'little'),
      int(self.attribute_header_bytes).to_bytes(4, 'little'),
      int(self.num_components).to_bytes(4, 'little'),
      struct.pack('f', self.physical_path_length),
      int(self.space).to_bytes(1),
    ])
    header_crc16 = lib.crc16(header)
    
    return b''.join([
      OstdHeader.MAGIC,
      header,
      header_crc16.to_bytes(2, 'little'),
    ])

@dataclass
class OstdSpatialIndex:
  minpt:npt.NDArray[np.float32]
  maxpt:npt.NDArray[np.float32]
  chunk_size:npt.NDArray[np.float32]
  index_binary:bytes
  paths_binary:bytes

@dataclass
class OstdAttribute:
  def __init__(self):
    name = ""
    attribute_type = None
    data_type = None
    compression = None
    unit = (None, None)
    num_components = 0
    content_length = 0

  def to_bytes(self, name_width:int):
    if self.num_components >= 16:
      raise ValueError(f"Only 15 components supported per an attribute. Got {self.num_components} for {self.id}.")

    name = bytearray(name_width)
    name[:len(self.id)] = self.id
    dtype = np.uint8(TO_DATATYPE[np.dtype(self.dtype).type])
    type_field = dtype | (np.uint8(self.num_components) << 4)
    return b''.join([ name, type_field ])

  # def to_dict(self) -> dict[str,Any]:
  #   return {
  #     "id": self.id,
  #     "data_type": np.dtype(self.dtype).name,
  #     "num_components": self.num_components,
  #   }

  def decode_flags(self, flags:int):
    offset = 0
    masks = [ 0b0, 0b1, 0b11, 0b111, 0b1111, 0b11111 ]
    def read_int(N):
      nonlocal offset
      res = (flags >> N) & masks[N] 
      offset += N
      return res

    self.data_type = DataType(read_int(4))
    self.compression = CompressionType(read_int(4))
    self.attribute_type = AttributeType(read_int(1))

  @classmethod
  def from_bytes(self, binary:bytes, name_width:int) -> "OstdAttribute":

    attr = OstdAttribute()

    name = binary[:name_width]
    name = name.split(b'\x00', 1)[0].decode('utf-8')
    attr.name = name

    flags = int.from_bytes(binary[name_width:name_width+2])
    attr.decode_flags(flags)

    offset = name_width + 2

    unit_info = int.from_bytes(binary[offset:offset + 2])
    DimensionTypeClass = QUANTITY_TYPE[int(unit_info & 0b111)]
    si_prefix = SIPrefixType((unit_info >> 3) & 0b11111)
    dimension = DimensionTypeClass((unit_info >> (3+5)) & 0b1111)
    attr.unit = (si_prefix, dimension)

    offset += 2
    attr.num_components = int.from_bytes(binary[offset:offset+1])
    offset += 1
    attr.content_length = int.from_bytes(binary[offset:offset+8])

    hb = len(OstdHeader.HEADER_BYTES)
    attr.crc16 = int.from_bytes(binary[hb-2:hb])

    return attr


@dataclass
class OstdAttributeHeader:
  name_width:int
  attributes:list[OstdAttribute]
  crc16:int

  def to_bytes(self) -> bytes:
    header = [
      len(self.attributes).to_bytes(2),
      int(self.name_width).to_bytes(1),
    ]
    header += [ attr.to_bytes() for attr in self.attributes ]
    binary = b''.join(header)
    binary += lib.crc16(binary)
    return binary

  @classmethod
  def from_bytes(kls, binary:bytes) -> "OstdAttributeHeader":
    pass
