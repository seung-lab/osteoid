from typing import Optional, Any, Literal

from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from functools import partial

from enum import IntEnum

import numpy as np
import numpy.typing as npt

import fastremap
import fastosteoid

from ... import lib

from .header import (
  OstdAttribute,
  OstdAttributeSection,
  OstdHeader, 
  OstdTransform,
  OstdTransformSection,
)

from .types import (
  AttributeType,
  CompressionType, 
  GraphType,
  SIUnit,
  SpaceType,
  TO_DATATYPE,
  TRANSLATE_UNIT,
)

def check_crc32c(binary:bytes, stored_crc:int):
    computed_crc = lib.crc32c(binary)    
    if stored_crc != computed_crc:
      raise ValueError(f"Header corruption detected. Stored crc32c: {stored_crc}, Computed crc32c: {computed_crc}")

# represents one skeleton section
# in a possibly multipart file
@dataclass
class OstdSkeletonPart:
  header:OstdHeader
  vertices:npt.NDArray[np.generic]
  edges:npt.NDArray[np.unsignedinteger]
  spaces:Optional[OstdTransformSection] = None
  attributes:Optional[
    OrderedDict[str,tuple[SIUnit, npt.NDArray[np.generic]]]
  ] = None

  def _encode_vertices(self, vertices:np.ndarray) -> bytes:
    if self.header.vertex_compression == CompressionType.NONE:
      vertex_binary = vertices.tobytes("C")
      vertex_binary += lib.crc32c(vertex_binary).to_bytes(4, 'little')
      return vertex_binary
    elif self.header.vertex_compression == CompressionType.DRACO:
      import DracoPy
      return DracoPy.compress(vertices, preserve_order=True, quantization_bits=12)
    elif self.header.vertex_compression == CompressionType.GZIP:
      import deflate
      return deflate.gzip_compress(vertices.tobytes("C"))
    else:
      raise ValueError(f"Unsupported compression type: {self.header.vertex_compression}")

  def _encode_linked_paths(self) -> tuple[bytes, bytes]:
    (all_paths, all_edges, has_cycle, N) = fastosteoid.linked_paths(self.edges)

    # flatten all_paths
    all_paths = [ 
      item 
      for sublist in all_paths 
        for item in sublist
    ]

    # sort vertices by path locations
    reorder = np.concatenate(all_paths)
    verts = self.vertices[reorder]
    edges = np.concatenate(all_edges).astype(self.header.edge_dtype, copy=False)

    inv = np.empty([ self.header.Nv ], dtype=reorder.dtype, order="C")
    inv[reorder] = np.arange(len(reorder), dtype=reorder.dtype)
    edges = inv[edges]
    del inv

    pdtype = lib.compute_dtype(verts.shape[0])
    path_lengths = np.asarray([
      len(path) for path in all_paths
    ], dtype=pdtype)

    self.header.num_components = N
    if has_cycle:
      self.header.graph_type = GraphType.CYCLIC
    else:
      self.header.graph_type = GraphType.TREE

    edge_binary = b''.join([
      len(path_lengths).to_bytes(8, 'little'),
      path_lengths.tobytes(),
      edges.tobytes(),
    ])
    edge_binary += lib.crc32c(edge_binary).to_bytes(4, 'little')

    vertex_binary = self._encode_vertices(verts)

    return (vertex_binary, edge_binary, reorder)

  def _encode_geometry(self) -> tuple[bytes,bytes,Optional[np.ndarray]]:
      return self._encode_linked_paths()

  @classmethod
  def _decode_edge_representation_linked_paths(
    kls,
    header:OstdHeader,
    binary:bytes,
    offset:int,
  ) -> npt.NDArray[np.unsignedinteger]:
    check_buf = np.frombuffer(
      binary, 
      offset=offset, 
      count=(header.edge_bytes - 4), 
      dtype=np.uint8
    )
    crc_off = offset + header.edge_bytes - 4
    stored_crc32c = int.from_bytes(binary[crc_off:crc_off+4], 'little')
    check_crc32c(check_buf,  stored_crc32c)

    edges_buf = np.frombuffer(
      binary, 
      offset=offset, 
      count=header.edge_bytes, 
      dtype=np.uint8
    )

    edge_width = np.dtype(header.edge_dtype).itemsize
    return fastosteoid.decode_linked_path_edges(edges_buf, header.Nv, edge_width)

  @classmethod
  def _decode_edge_representation(
    kls, 
    header:OstdHeader,
    binary:bytes,
    offset:int,
  ) -> npt.NDArray[np.unsignedinteger]:
    return kls._decode_edge_representation_linked_paths(header, binary, offset)
  
  def to_bytes(self) -> bytes:
    vertex_binary, edge_binary, reorder = self._encode_geometry()
    self.header.cable_length = self.cable_length()

    self.header.has_transform = False
    transform_binary = b''
    if self.spaces is not None and len(self.spaces.spaces):
      transform_binary = self.spaces.to_bytes()
      self.header.has_transform = True

    self.header.attribute_header_bytes = 0
    attributes_header_binary = b''
    if self.attributes is not None and len(self.attributes) > 0:
      attributes_header_binary = self._create_attributes_header().to_bytes()
      self.header.attribute_header_bytes = len(attributes_header_binary)

    attributes_binary = b''
    if self.attributes is not None and len(self.attributes) > 0:
      attributes_binary = []
      for name, (unit, arr) in self.attributes.items():
        assert arr.shape[0] == self.header.Nv
        if reorder is not None:
          arr = arr[reorder]
        attr_binary = arr.tobytes("C")
        attr_binary += lib.crc32c(attr_binary).to_bytes(4, 'little')
        attributes_binary.append(attr_binary)
      attributes_binary.append(attributes_header_binary)
      attributes_binary = b''.join(attributes_binary)

    self.header.vertex_bytes = len(vertex_binary)
    self.header.edge_bytes = len(edge_binary)
    self.header.total_bytes = (
      OstdHeader.HEADER_BYTES +
      len(transform_binary) +  
      self.header.vertex_bytes + 
      self.header.edge_bytes + 
      len(attributes_binary)
    )

    return b''.join([
      self.header.to_bytes(),
      transform_binary,
      vertex_binary,
      edge_binary,
      attributes_binary,
    ])

  def _create_attributes_header(self) -> OstdAttributeSection:
    attrs = []
    for name, (unit, arr) in self.attributes.items():
      attr = OstdAttribute(
        name=name,
        attribute_type=AttributeType.VERTEX,
        data_type=TO_DATATYPE[np.dtype(arr.dtype).type],
        compression=CompressionType.NONE,
        unit=unit,
        num_components=(1 if arr.ndim == 1 else arr.shape[1]),
        content_length=arr.nbytes,
      )
      attrs.append(attr)
    return OstdAttributeSection(attrs)

  def cable_length(self):
    if not np.isnan(self.header.cable_length):
      return self.header.cable_length

    original_space = self.header.space
    self.physical_space()

    verts = self.vertices
    edges = self.edges

    dtype = self.vertices.dtype
    if np.dtype(dtype).itemsize < 4:
      dtype = np.float32

    v1 = verts[edges[:,0]].astype(dtype, copy=False)
    v2 = verts[edges[:,1]].astype(dtype, copy=False)
    dist = np.linalg.norm(v2 - v1, axis=1)

    self.header.cable_length = np.sum(dist.astype(np.float64)).astype(np.float32)

    self.change_space(original_space)

    return self.header.cable_length

  def _graph_type(self, components:list[np.ndarray]) -> GraphType:
    for component in components:
      if fastosteoid.has_cycle(component):
        return GraphType.CYCLIC

    return GraphType.TREE

  def graph_type(self) -> GraphType:
    return self._graph_type(
      fastosteoid.compute_components(self.edges, self.header.Nv)
    )

  def num_components(self) -> int:
    sentinel = np.iinfo(np.uint32).max

    if self.header.num_components < sentinel:
      return self.header.num_components

    forest = fastosteoid.compute_components(self.edges, self.header.Nv)

    if len(forest) >= sentinel:
      raise ValueError("Number of components exceeds maximum representable value.")

    self.header.num_components = len(forest)
    return len(forest)

  def change_space(self, idx:int):
    if idx < 0 or idx > len(self.spaces.spaces) or idx > 255:
      raise ValueError(f"No transform exists for this space index. {idx} (max: {len(self.spaces.spaces)})")

    if idx == self.header.space:
      return

    verts = self.vertices
    src_space = self.header.space
    identity = np.eye(4, dtype=verts.dtype)
    total_transform = identity

    if src_space != 0:
      # transform from src space to root space
      src_to_root = np.linalg.inv(self.spaces.spaces[src_space - 1].transform)
      total_transform = total_transform @ src_to_root

    if idx != 0:
      # transform from root space to target space
      root_to_dst = self.spaces.spaces[idx - 1].transform
      total_transform = total_transform @ root_to_dst

    if not np.allclose(total_transform, identity):
      ones = np.ones((verts.shape[0], 1), dtype=verts.dtype)
      verts = np.hstack([verts, ones])
      del ones
      verts = verts @ total_transform.T
      verts[:, :3] /= verts[:, 3:4]
      self.vertices = verts[:, :3]

    self.header.space = idx

  def voxel_space(self):
    self.change_space(0) # 0 is always voxel space

  def physical_space(self):
    if len(self.spaces.spaces) == 0:
      return # transform implied to be identity

    self.change_space_by_type(SpaceType.PHYSICAL)

  @property
  def unit(self):
    return self.get_space_unit(self.header.space)

  def get_space_unit(self, idx:int) -> SIUnit:
    if idx == 0:
      return self.header.length_unit
    return self.spaces.spaces[idx-1].unit

  def get_space_type(self, idx:int) -> SpaceType:
    if idx == 0:
      return self.header.space_type
    return self.spaces.spaces[idx-1].space

  def current_space_type(self) -> SpaceType:
    return self.get_space_type(self.header.space)

  def get_space_by_type(self, typ:SpaceType):
    # self.spaces : OstdTransformSection
    # self.spaces.spaces: list[OstdTransform]
    # space.space: OstdTransform.space: SpaceType

    if self.header.space_type == typ:
      return 0

    for i, space in enumerate(self.spaces.spaces):
      if space.space == typ:
        return i + 1

    raise ValueError(f"{typ} not found in space list.")

  def change_space_by_type(self, typ:SpaceType):
    i = self.get_space_by_type(typ)
    self.change_space(i)

  @classmethod
  def _decode_vertices(kls, header:OstdHeader, binary:bytes, offset:int) -> np.ndarray:
    if header.vertex_compression == CompressionType.NONE: 
      vertices = np.frombuffer(
        binary,
        offset=offset,
        count=(header.Nv * header.num_axes),
        dtype=header.vertex_dtype,
      ).reshape((header.Nv, header.num_axes), order="C")

      check_buf = vertices.view(np.uint8).reshape((vertices.nbytes,))
      off = offset + header.vertex_bytes - 4
      stored_crc32c = int.from_bytes(binary[off:off+4], 'little')
      check_crc32c(check_buf,  stored_crc32c)
      return vertices
    else:
      buffer = np.frombuffer(
        binary,
        offset=offset,
        count=header.vertex_bytes,
        dtype=np.uint8,
      )

      if header.vertex_compression == CompressionType.DRACO:
        import DracoPy
        return DracoPy.decode(buffer)
      elif header.vertex_compression == CompressionType.GZIP:
        import deflate
        vbuf = deflate.gzip_decompress(buffer)
        return np.frombuffer(
          vbuf,
          dtype=header.vertex_dtype,
        ).reshape((header.Nv, header.num_axes), order="C")
      else:
        raise ValueError(f"Compression type not supported: {header.vertex_compression}")

  @classmethod
  def from_bytes(kls, binary:bytes, offset:int = 0) -> "OstdSkeletonPart":
    header = OstdHeader.from_bytes(binary, offset=offset)

    off = offset + OstdHeader.HEADER_BYTES

    if header.has_transform:
      spaces = OstdTransformSection.from_bytes(binary, offset=off)
      off += spaces.nbytes
    else:
      spaces = OstdTransformSection([])
      off += 0

    if header.edge_compression != CompressionType.NONE:
      raise ValueError(f"Compression not yet supported.")

    vertices = kls._decode_vertices(header, binary, off)
    off += header.vertex_bytes

    edges = kls._decode_edge_representation(header, binary, off)
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
        )

        if attribute.num_components > 1:
          arr = arr.reshape((header.Nv, attribute.num_components), order="C")

        check_buf = arr.view(np.uint8).reshape((arr.nbytes,))
        off += arr.nbytes
        stored_crc32c = int.from_bytes(binary[off:off+4], 'little')
        check_crc32c(check_buf,  stored_crc32c)
        off += 4

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
    self._props = {}

  def __list__(self):
    return list(self._props.keys())

  def __contains__(self, key:str) -> bool:
    return key in self._props

  def __getattr__(self, key:str) -> np.ndarray:
    if key not in self._props:
      raise AttributeError(f"skeleton has no attribute \"{key}\"")

    val = self._props[key]
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
  def spaces(self) -> list[OstdTransform]:
    root_space = OstdTransform(
      self.parts[0].header.length_unit, 
      self.parts[0].header.space_type, 
      np.eye(4, dtype=np.float32)
    )
    return [ root_space ] + self.parts[0].spaces.spaces

  @property
  def transforms(self) -> list[npt.NDArray[np.float32]]:
    return [ np.eye(4, dtype=np.float32) ] + [ space.transform for space in self.parts[0].spaces.spaces ]

  @property
  def coordinate_frame(self):
    return self.parts[0].header.coordinate_frame

  @property
  def voxel_centered(self):
    return self.parts[0].header.voxel_centered

  def ids(self) -> list[int]:
    return [ int(part.header.id) for part in self.parts ]

  def split_by_id(self) -> dict[int,"OstdSkeleton"]:
    """
    If skeleton parts with multiple ids are contained in this file,
    split them into separate skeletons.
    """
    skels = defaultdict(list)

    for part in self.parts:
      skels[part.id].append(part)

    return { id: OstdSkeleton(parts) for id, parts in skels.items() }

  def change_space(self, idx:int):
    for part in self.parts:
      if not part.header.has_transform and idx != 0:
        raise ValueError("Skeleton part does not have a transform section in its binary.")
      part.change_space(idx)

  def change_space_by_type(self, stype:SpaceType):
    for part in self.parts:
      if not part.header.has_transform and idx != 0:
        raise ValueError("Skeleton part does not have a transform section in its binary.")
      part.change_space_by_type(stype)

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
  def attributes(self) -> list[OstdAttribute]:
    return {
      name: (unit, arr)
      for name, (unit, arr) in self.parts[0].attributes.items()
    }

  def attribute(self, key:str) -> tuple[SIUnit, np.ndarray]:
    return self.parts[0].attributes.get(key, None)

  def attribute_unit(self, key:str) -> SIUnit:
    attr = self.attribute(key)
    if attr is None:
      return None
    return attr[0]

  def attribute_value(self, key:str) -> np.ndarray:
    attr = self.attribute(key)
    if attr is None:
      return None
    return attr[1]

  @property
  def unit(self) -> SIUnit:
    if len(self.parts) == 0:
      return SIUnit()
    return self.parts[0].get_space_unit(self.parts[0].header.space)

  @property
  def human_readable_unit(self) -> str:
    return FROM_LENGTH_UNIT[self.unit.tuple()]

  @property
  def cable_length(self) -> float:
    master_unit = self.parts[0].unit
    master_space_type = self.parts[0].current_space_type()

    physical_length = 0
    for part in self.parts:
      part.change_space_by_type(master_space_type)

      part_length = part.cable_length()

      if part.unit == master_unit:
        physical_length += part_length
        continue

      conversion = 1000 ** (master_unit.si_prefix - part.unit.si_prefix)
      physical_length += part_length * conversion

    return physical_length

  @property
  def vertices(self) -> npt.NDArray[Any]:
    master_unit = self.parts[0].unit
    master_space_type = self.parts[0].current_space_type()

    verts = []
    for part in self.parts:
      part.change_space_by_type(master_space_type)
      if part.unit == master_unit:
        verts.append(part.vertices)
        continue

      # si_prefix, base_unit = part.unit.tuple()
      # factor = (master_si_value / SI_PREFIX_VALUE[si_prefix])

      # factor *= length_conversion_factor(base_unit, master_base_unit)
      factor = 1000 ** (master_unit.si_prefix - part.unit.si_prefix)
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
  def edges(self) -> npt.NDArray[np.unsignedinteger]:
    if len(self.parts) <= 1:
      return self.parts[0].edges

    Ne = np.sum([ part.header.Ne for part in self.parts ])
    combined_edge_dtype = lib.compute_dtype(Ne)

    offset = 0
    edges = []
    for part in self.parts:
      part_edges = part.edges.astype(combined_edge_dtype, copy=False)
      part_edges += offset
      edges.append(part_edges)
      offset += part.header.Nv

    return np.concatenate(edges)

  def is_multipart(self) -> bool:
    return len(self.parts) > 1

  def append(self, skel:"OstdSkeleton"):
    self.parts.extend(skel.parts)

  def save(self, filename:str):
    with open(filename, "wb") as f:
      f.write(self.to_bytes())

  @classmethod
  def load(kls, filename:str, allow_mmap:bool = False) -> "OstdSkeleton":
    from ... import util
    binary = util._load(filename, allow_mmap=allow_mmap)
    return kls.from_bytes(binary)

  @classmethod
  def create(kls, 
    vertices:npt.NDArray[np.generic], 
    edges:npt.NDArray[np.unsignedinteger],
    length_unit:str = "nm",
    space_type:SpaceType = SpaceType.UNKNOWN,
    id:Optional[int] = None,
    space:int = 0,
    spaces:list = [],
    coordinate_frame:str = "+X+Y+Z",
    voxel_centered:bool = True,
    attributes:dict[str,npt.NDArray[np.generic]] = {},
    vertex_compression:Optional[str] = None,
  ):
    Nv = vertices.shape[0]
    edge_dtype = lib.compute_dtype(Nv)

    if vertex_compression is None:
      vert_compress = CompressionType.NONE
    else:
      raise ValueError(f"Unsupported compression type: {vertex_compression}")

    header = OstdHeader(
      Nv = Nv,
      Ne = edges.shape[0],
      coordinate_frame = coordinate_frame,
      edge_data_type = TO_DATATYPE[np.dtype(edge_dtype).type],
      edge_compression = CompressionType.NONE,
      has_transform = False,
      id = id,
      length_unit = TRANSLATE_UNIT[length_unit.lower()],
      num_axes = vertices.shape[1],
      space=int(space),
      space_type=space_type,
      vertex_data_type = TO_DATATYPE[np.dtype(vertices.dtype).type],
      vertex_compression = vert_compress,
      voxel_centered = bool(voxel_centered),
    )

    spaces = OstdTransformSection([
        OstdTransform(unit=transform[0], space=transform[1], transform=transform[2])
        if isinstance(transform, tuple)
        else OstdTransform(space=SpaceType.UNKNOWN, transform=transform)
      for transform in spaces
    ])

    for attr_name, attr_value in attributes.items():
      if isinstance(attr_value, tuple):
        if isinstance(attr_value[0], str):
          attributes[attr_name] = (TRANSLATE_UNIT[attr_value[0]], attr_value[1])

    return OstdSkeleton([  
      OstdSkeletonPart(
        header=header, 
        vertices=vertices, 
        edges=edges.astype(edge_dtype, copy=False),
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
      return np.concatenate([ part.attributes[name][1] for part in skel.parts ])

    for name in parts[0].attributes:
      skel.a._props[name] = partial(getattribute, skel, name)

    return skel

  def drop_attribute(self, name:str):
    for part in self.parts:
      del part.attributes[name]

