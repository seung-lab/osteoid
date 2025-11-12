from typing import Optional, Any

from collections import OrderedDict
from dataclasses import dataclass
from functools import partial

from enum import IntEnum

import numpy as np
import numpy.typing as npt

from ... import lib

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
  TO_LENGTH_UNIT,
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
    vertex_binary += lib.crc32c(vertex_binary).to_bytes(4, 'little')

    edge_binary = self.edges.tobytes("C")
    edge_binary += lib.crc32c(edge_binary).to_bytes(4, 'little')

    self.header.vertex_bytes = len(vertex_binary)
    self.header.edge_bytes = len(edge_binary)
    self.header.total_bytes = (
      OstdHeader.HEADER_BYTES + 
      self.header.vertex_bytes + 
      self.header.edge_bytes + 
      self.header.spatial_index_bytes +
      self.header.attribute_header_bytes
    )
    return b''.join([
      self.header.to_bytes(),
      vertex_binary,
      edge_binary,
    ])

  def change_space(self, idx:int):
    if idx < 0 or idx >= len(self.spaces.spaces) or idx > 255:
      raise ValueError(f"No transform exists for this space index. {idx} (max: {len(self.spaces.spaces)})")

    if idx == self.header.space:
      return

    verts = self.vertices
    src_space = self.header.space
    identity = np.eye(4, dtype=verts.dtype)
    total_transform = identity

    if src_space != 0:
      # transform from src space to voxels
      src_to_vox = np.linalg.inv(self.spaces[src_space].transform)
      total_transform = total_transform @ src_to_vox

    if idx != 0:
      # transform from voxels to target space
      vox_to_dst = self.spaces[idx - 1].transform
      total_transform = total_transform @ vox_to_dst

    if not np.allclose(total_transform, identity):
      ones = np.ones((verts.shape[0], 1), dtype=verts.dtype)
      verts = np.hstack([verts, ones])
      del ones
      verts = verts @ total_transform.T
      verts[:, :3] /= verts[:, 3:4]
      self.vertices = verts[:, :3]

    self.header.space = idx

  @classmethod
  def from_bytes(kls, binary:bytes, offset:int = 0) -> "OstdSkeletonPart":
    header = OstdHeader.from_bytes(binary, offset=offset)

    off = offset + OstdHeader.HEADER_BYTES

    if header.has_transform:
      spaces = OstdTransformSection.from_bytes(
        binary[off:off + 1 + 65 * 255]
      )
      off += spaces.nbytes
    else:
      spaces = OstdTransformSection([ 
        OstdTransform(SpaceType.GENERIC, np.eye(4,4, dtype=np.float32))
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

    vertices = np.frombuffer(
      binary,
      offset=off,
      count=(header.Nv * header.num_axes),
      dtype=header.vertex_dtype,
    ).reshape((header.Nv, header.num_axes), order="C")

    off += header.vertex_bytes

    if header.edge_representation != EdgeRepresentationType.PAIR:
      raise ValueError("Only edge lists are currently supported.")

    edges = np.frombuffer(
      binary, 
      offset=off,
      count=header.Ne * 2, 
      dtype=header.edge_dtype
    ).reshape((header.Ne, 2), order="C")

    off += header.edge_bytes

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

class OstdSkeletonProperties:
  def __init__(self):
    self.props = {}

  def __getattr__(self, key:str):
    if key not in self.props:
      raise AttributeError(f"skeleton property has no attribute \"{key}\"")

    val = self.props[key]
    if callable(val):
      return val()
    return val

# represents a full skeleton including
# multiple parts
class OstdSkeleton:
  def __init__(self, parts:list[OstdSkeletonPart] = []):
    self.parts = parts
    self.a = OstdSkeletonProperties()

  @property
  def id(self) -> int:
    return self.parts[0].header.id

  @property
  def spaces(self):
    return self.parts[0].spaces

  @property
  def transforms(self):
    return [ space.transform for space in self.parts[0].spaces ]

  @property
  def coordinate_frame_orientation(self):
    return self.parts[0].header.coordinate_frame_orientation

  @property
  def voxel_centered(self):
    return self.parts[0].header.voxel_centered

  def change_space(self, idx:int):
    for part in self.parts:
      if not part.header.has_transform and idx != 0:
        raise ValueError("Skeleton part does not have a transform section in its binary.")
      part.change_space(idx)

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
    if len(self.parts) == 0:
      return (SIPrefixType.NONE, LengthType.VOXEL)
    return self.parts[0].header.length_unit

  @property
  def human_readable_unit(self) -> str:
    return FROM_LENGTH_UNIT[self.unit]

  @property
  def physical_length(self) -> float:
    master_unit = self.parts[0].header.length_unit
    master_si_unit, master_base_unit = master_unit
    master_si_value = SI_PREFIX_VALUE[master_si_unit]

    physical_length = 0
    for part in self.parts:
      if np.isnan(part.header.physical_length):
        # NB: Would be possible to compute this on demand....
        raise ValueError("NaN encountered in physical length reporting.")

      if part.header.length_unit == master_unit:
        physical_length += part.header.physical_length
        continue

      si_unit, base_unit = part.header.length_unit
      si_conversion = (master_si_value / SI_PREFIX_VALUE[si_unit])
      base_conversion = length_conversion_factor(base_unit, master_base_unit)
      physical_length += part.header.physical_length * (si_conversion * base_conversion)

    return physical_length

  @property
  def vertices(self) -> npt.NDArray[Any]:
    master_unit = self.parts[0].header.length_unit
    master_si_prefix, master_base_unit = self.parts[0].header.length_unit
    master_si_value = SI_PREFIX_VALUE[master_si_prefix]

    verts = []
    for part in self.parts:
      if part.header.length_unit == master_unit:
        verts.append(part.vertices)
        continue

      si_prefix, base_unit = part.header.length_unit
      factor = (master_si_value / SI_PREFIX_VALUE[si_prefix])
      factor *= length_conversion_factor(base_unit, master_base_unit)

      if np.issubdtype(part.vertices.dtype, np.floating):
        verts.append(
          part.vertices * factor
        )
      else:
        casted = part.vertices.astype(np.float32)
        casted *= factor
        verts.append(
          casted.astype(part.vertices.dtype, copy=False)
        )

    return np.concatenate(verts, axis=0)

  @property
  def edges(self) -> npt.NDArray[np.uint64]:
    if len(self.parts) <= 1:
      return self.parts[0].edges

    offset = 0
    edges = []
    for part in self.parts:
      edges.append(part.edges + offset)
      offset += part.header.Nv

    return np.concatenate(edges)

  def is_multipart(self) -> bool:
    return len(self.parts) > 1

  def append(self, skel:"OstdSkeleton"):
    self.parts.extend(skel.parts)

  def save(self, filename:str):
    pass

  @classmethod
  def load(kls, filename:str) -> "OstdSkeleton":
    pass

  @classmethod
  def create(kls, 
    vertices:npt.NDArray[np.generic], 
    edges:npt.NDArray[np.unsignedinteger],
    length_unit:str = "nm",
    id:Optional[int] = None,
    spaces:list = [],
    coordinate_frame_orientation:str = "+X+Y+Z",
    voxel_centered:bool = True,
    attributes:dict[str,npt.NDArray[np.generic]] = {},
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
      length_unit = TO_LENGTH_UNIT[length_unit.lower()],
      num_axes = vertices.shape[1],
      vertex_data_type = TO_DATATYPE[np.dtype(vertices.dtype).type],
      voxel_centered = bool(voxel_centered),
    )

    spaces = OstdTransformSection([
        OstdTransform(space=transform[0], transform=transform[1])
        if isinstance(transform, tuple)
        else OstdTransform(space=SpaceType.GENERIC, transform=transform)
      for transform in spaces
    ])

    return OstdSkeleton([  
      OstdSkeletonPart(
        header=header, 
        vertices=vertices, 
        edges=edges,
        spaces=spaces,
        attributes=attributes,
      )
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







