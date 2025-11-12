from typing import Optional, Any

from collections import OrderedDict
from dataclasses import dataclass
from functools import partial

from enum import IntEnum

import numpy as np
import numpy.typing as npt

from .header import (
  OstdAttribute,
  OstdAttributeSection,
  OstdHeader, 
  OstdSpatialIndex,
  OstdTransform,
  OstdTransformSection,
)

from .types import (
  AttributeType,
  CompressionType, 
  EdgeRepresentationType,
  LengthType,
  SIPrefixType,
  SpaceType,
  length_conversion_factor,
  TO_DATATYPE,
  SI_PREFIX_VALUE,
  FROM_LENGTH_UNIT,
)

# represents one skeleton section
# in a possibly multipart file
@dataclass
class OstdSkeletonPart:
  header:OstdHeader
  vertices:npt.NDArray[np.generic]
  edges:npt.NDArray[np.unsignedinteger]
  spaces:Optional[OstdTransformSection] = None
  spatial_index:Optional[OstdSpatialIndex] = None
  attributes:Optional[
    OrderedDict[str,tuple[tuple[SIPrefixType, IntEnum], npt.NDArray[np.generic]]]
  ] = None

  def to_bytes(self) -> bytes:
    vertex_binary = self.vertices.tobytes("C")
    vertex_binary += lib.crc32c(vertex_binary)

    edge_binary = self.edges.tobytes("C")
    edge_binary += lib.crc32c(edge_binary)

    self.header.vertex_bytes = len(vertex_binary)
    self.header.edge_bytes = len(edge_binary)
    self.header.total_bytes = (
      OstdHeader.HEADER_BYTES + 
      header.vertex_bytes + 
      header.edge_bytes + 
      header.spatial_index_bytes +
      header.attribute_header_bytes
    )
    return b''.join([
      header.to_bytes(),
      vertex_binary,
      edge_binary,
    ])

  @classmethod
  def from_bytes(kls, binary:bytes, offset:int = 0) -> "OstdSkeletonPart":
    header = OstdHeader.from_bytes(binary, offset=offset)

    off = offset + OstdHeader.HEADER_BYTES

    if header.has_transform:
      spaces = OstdTransformSection.from_bytes(
        binary[off:off + 1 + 65 * 255]
      )
      off += transform.nbytes
    else:
      spaces = OstdTransformSection([ 
        OstdTransform(SpaceType.VOXEL, np.eye(4,4, dtype=np.float32))
      ])
      off += 0

    spatial_index = None
    off += 0
    if header.spatial_index_bytes > 0:
      raise ValueError("Spatial index not implemented.")

    if header.vertex_compression != CompressionType.NONE:
      raise ValueError(f"Compression not yet supported.")
    if header.edge_compression != CompressionType.NONE:
      raise ValueError(f"Compression not yet supported.")

    vertex_bytes = binary[off:off+header.vertex_bytes]
    off += header.vertex_bytes

    vertices = np.frombuffer(
      vertex_bytes, 
      count=(header.Nv * header.num_axes),
      dtype=header.vertex_dtype,
    ).reshape((header.Nv, header.num_axes), order="C")
    del vertex_bytes

    if header.edge_representation != EdgeRepresentationType.PAIR:
      raise ValueError("Only edge lists are currently supported.")

    edge_bytes = binary[off:off+header.edge_bytes]
    off += header.edge_bytes

    edges = np.frombuffer(
      edge_bytes, 
      count=header.Ne, 
      dtype=header.edge_dtype
    ).reshape((header.Ne, 2), order="C")

    attribute_header = None
    attributes = OrderedDict()
    if header.attribute_header_bytes > 0:
      attribute_header = OstdAttributeSection.from_bytes(binary[-header.attribute_header_bytes:])
      for attribute in attribute_header.attributes:

        if attribute.attribute_type == AttributeType.EDGE:
          raise ValueError("Edge vertex attributes are not yet supported.")
        if attribute.compression != CompressionType.NONE:
          raise ValueError("Attribute compression not yet supported.")

        arr = np.frombuffer(
          binary,
          offset=off,
          count=(header.Nv * attribute.num_components),
          dtype=attribute.dtype,
        ).reshape((header.Nv, attribute.num_components), order="C")

        attributes[attribute.name] = (attribute.unit, arr)

    return OstdSkeletonPart(
      header=header,
      spaces=spaces,
      vertices=vertices,
      edges=edges,
      attributes=attributes,
    )

# represents a full skeleton including
# multiple parts
class OstdSkeleton:
  def __init__(self, parts:list[OstdSkeletonPart] = []):
    self.parts = parts

  @property
  def id(self):
    return self.parts[0].header.id

  @property
  def num_vertices(self) -> int:
    return sum(( part.header.Nv for part in self.parts ))

  @property
  def num_edges(self) -> int:
    return sum(( part.header.Ne for part in self.parts ))

  @property
  def num_components(self) -> int:
    """Note: This may be an overestimate for multi-part files."""
    return sum(( part.header.num_components for part in self.parts ))

  @property
  def unit(self) -> tuple[SIPrefixType, LengthType]:
    if len(self.parts[0]) == 0:
      return (SIPrefixType.NONE, LengthType.VOXEL)
    return self.parts[0].header.length_unit

  @property
  def human_readable_unit(self) -> str:
    return FROM_LENGTH_UNIT[self.unit]

  @property
  def physical_length(self) -> float:
    master_si_unit, master_base_unit = self.parts[0].header.length_unit
    master_si_value = SI_PREFIX_VALUE[master_si_unit]

    physical_length = 0
    for part in self.parts:
      if np.isnan(part.header.physical_length):
        # NB: Would be possible to compute this on demand....
        raise ValueError("NaN encountered in physical length reporting.")

      if part.header.unit == master_unit:
        physical_length += part.header.physical_length
        continue

      si_unit, base_unit = part.header.unit
      si_conversion = (master_si_value / SI_PREFIX_VALUE[si_unit])
      base_conversion = length_conversion_factor(base_unit, master_base_unit)
      physical_length += part.header.physical_length * (si_conversion * base_conversion)

    return physical_length

  @property
  def vertices(self) -> npt.NDArray[Any]:
    master_si_prefix, master_unit = self.parts[0].header.length_unit
    master_si_value = SI_PREFIX_VALUE[master_si_prefix]

    verts = []
    for part in self.parts:
      if part.header.length_unit == master_unit:
        verts.append(part.vertices)
        continue

      si_prefix, base_unit = part.header.length_unit
      factor = (master_si_value / SI_PREFIX_VALUE[si_prefix])
      factor *= length_conversion_factor(base_unit, master_unit)

      if isinstance(part.vertices, np.floating):
        verts.append(
          part.vertices * factor
        )
      elif factor >= 0:
        verts.append(
          part.vertices * factor
        )
      else:
        verts.append(
          part.vertices // (1/factor)
        )

    return np.concatenate(verts)

  @property
  def edges(self) -> npt.NDArray[np.uint64]:
    if len(self.parts) <= 1:
      return self.parts[0].edges

    offset = 0
    edges = []
    for part in parts:
      edges.append(part.edges + offset)
      offset += part.header.Nv

    return np.concatenate(edges)

  def is_multipart(self) -> bool:
    return len(self.parts) > 1

  def append(self, skel:"OstdSkeleton"):
    self.parts.append(skel)

  def save(self, filename:str):
    pass

  @classmethod
  def load(kls, filename:str) -> "OstdSkeleton":
    pass

  @classmethod
  def create(kls, 
    vertices:npt.NDArray[np.generic], 
    edges:npt.NDArray[np.unsignedinteger],
    id:Optional[int] = None,
    coordinate_frame_orientation:str = "+X+Y+Z",
    voxel_centered:bool = True,
  ):
    header = OstdHeader(
      Nv = vertices.shape[0],
      Ne = edges.shape[0],
      coordinate_frame_orientation = coordinate_frame_orientation,
      edge_data_type = TO_DATATYPE[np.dtype(edges.dtype).type],
      edge_compression = CompressionType.NONE,
      edge_representation = EdgeRepresentationType.PAIR,
      has_transform = False,
      id = id,
      num_axes = vertices.shape[1],
      vertex_data_type = TO_DATATYPE[np.dtype(vertices.dtype).type],
      voxel_centered = bool(voxel_centered),
    )
    return OstdSkeleton([  
      OstdSkeletonPart(header, vertices=vertices, edges=edges)
    ])

  def to_bytes(self) -> bytes:
    return b''.join([
      part.to_bytes()
      for part in self.parts
    ])

  @classmethod
  def from_bytes(kls, binary:bytes) -> "OstdSkeleton":
    offset = 0
    parts = []
    while offset < len(binary):
      part = OstdSkeletonPart.from_bytes(binary, offset=offset)
      parts.append(part)
      offset += part.header.total_bytes

    skel = OstdSkeleton(parts)

    if len(parts) == 0:
      return skel

    def getattribute(skel:OstdSkeleton, name:str) -> np.ndarray:
      return np.concatenate(( part.attributes[name][1] for part in skel.parts ))

    for name in parts[0].attributes:
      prop = property(fget=partial(getattribute, skel, name))
      setattr(skel, name, prop)

    return skel







