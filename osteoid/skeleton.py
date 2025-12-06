from typing import Optional, Literal

from collections import defaultdict
import copy
import datetime
from enum import IntEnum
from io import BytesIO
import re
import os
import struct
import sys

import numpy as np
import numpy.typing as npt

import networkx as nx

import fastremap

import fastosteoid

from .exceptions import (
  SkeletonDecodeError, SkeletonEncodeError, 
  SkeletonUnassignedEdgeError, SkeletonTransformError,
  SkeletonAttributeMixingError
)
from .lib import Bbox, moving_average

IDENTITY = np.array([
  [1, 0, 0, 0],
  [0, 1, 0, 0],
  [0, 0, 1, 0],
], dtype=np.float32)

import importlib.metadata
__VERSION__ = importlib.metadata.version("osteoid")

class SWCVertexTypes(IntEnum):
  UNDEFINED = 0
  SOMA = 1
  AXON = 2
  BASAL_DENDRITE = 3
  APICAL_DENDRITE = 4

class Skeleton:
  """
  A stick figure representation of a 3D object.

  vertices: [[x,y,z], ...] float32
  edges: [[v1,v2], ...] uint32
  radii: [r1,r2,...] float32 distance from vertex to nearest boudary
  vertex_types: [t1,t2,...] uint8 SWC vertex types
  segid: numerical ID
  transform: 3x4 scaling and translation matrix (ie homogenous coordinates) 
    that represents the transformaton from voxel to physical coordinates.
    
    Example Identity Matrix:
    [
      [1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0]
    ]

  space: 'voxel', 'physical', or user choice (but other choices 
    make .physical_space() and .voxel_space() stop working as they
    become meaningless.)
  
  extra_attributes: You can specify additional per vertex
    data attributes (the most common are radii and vertex_type) 
    that are present in reading Precomputed binary skeletons using
    the following format:
    [
        {
          "id": "radius",
          "data_type": "uint8",
          "num_components": 1,
        }
    ]

    These attributes will become object properties. i.e. skel.radius

    Note that for backwards compatibility, skel.radius is treated 
    specially and is synonymous with skel.radii.

  default_attributes: when True, ensure that radius and vertex_types
    are initialized regardless of the content of extra_attributes.
    This was to conform to the SWC/Neuroglancer Precomputed
    conception of a skeleton.
  """
  def __init__(self, 
    vertices:Optional[np.ndarray] = None,
    edges:Optional[np.ndarray] = None, 
    radii:Optional[npt.NDArray[np.float32]] = None,
    vertex_types:Optional[npt.NDArray[np.uint8]] = None, 
    segid:Optional[int] = None,
    transform:Optional[npt.NDArray[np.float32]] = None,
    space:Literal['voxel', 'physical'] = 'voxel',
    extra_attributes:Optional[list] = None,
    default_attributes:bool = True,
  ):
    self.id = segid
    self.space = space

    if vertices is None:
      self.vertices = np.array([[]], dtype=np.float32).reshape(0,3)
    elif type(vertices) is list:
      self.vertices = np.array(vertices, dtype=np.float32)
    else:
      self.vertices = vertices

    if edges is None:
      self.edges = np.array([[]], dtype=np.uint32).reshape(0,2)
    elif type(edges) is list:
      self.edges = np.array(edges, dtype=np.uint32)
    else:
      self.edges = edges

    if radii is None and default_attributes:
      self.radius = np.full(shape=self.vertices.shape[0], fill_value=-1, dtype=np.float32)
    elif type(radii) is list:
      self.radius = np.asarray(radii, dtype=np.float32)
    elif radii is not None:
      self.radius = radii

    if vertex_types is None and default_attributes:
      # 0 = undefined in SWC (http://research.mssm.edu/cnic/swc.html)
      self.vertex_types = np.zeros(shape=self.vertices.shape[0], dtype=np.uint8)
    elif vertex_types is not None:
      self.vertex_types = np.asarray(vertex_types, dtype=np.uint8)

    self.default_attributes = bool(default_attributes)

    if extra_attributes is None:
      if default_attributes:
        self.extra_attributes = self._default_attributes()
      else:
        self.extra_attributes = []
    else:
      self.extra_attributes = extra_attributes

    if transform is None:
      self.transform = np.copy(IDENTITY)
    else:
      self.transform = np.asarray(transform).reshape( (3, 4) )

  @classmethod
  def _default_attributes(self):
    return [
      {
        "id": "radius",
        "data_type": "float32",
        "num_components": 1,
      }, 
      {
        "id": "vertex_types",
        "data_type": "uint8",
        "num_components": 1,
      }
    ]

  def _materialize_extra_attributes(self):
    for attr in self.extra_attributes:
      if not hasattr(self, attr["id"]):
        arr = np.zeros(
          [ self.vertices.shape[0], int(attr["num_components"]) ],
          dtype=np.dtype(attr["data_type"])
        )
        setattr(self, attr["id"], arr)

  def add_extra_attributes(self, attrs):
    for attr in attrs:
      if not hasattr(self, attr["id"]):
        self.extra_attributes.append(attr)
    self._materialize_extra_attributes()

  def add_vertex_attribute(self, name:str, arr:np.ndarray):
    """
    Add a new vertex attribute type and its data.
    """
    if not hasattr(self, name):
      self.extra_attributes.append({
        "id": name,
        "data_type": np.dtype(arr.dtype).name,
        "num_components": (1 if arr.ndim == 1 else arr.shape[2]),
      })
      setattr(self, name, arr)

  def _check_space(self):
    if self.space not in ('physical', 'voxel'):
      raise SkeletonTransformError(
        """
        Loss of coordinate frame information. If the space is not 'physical' or 'voxel',
        the meaning of applying this transform matrix is unknown.

        space: {}
        """.format(self.space)
      )

  def physical_space(self, copy:bool = True) -> "Skeleton":
    """
    Convert skeleton vertices into a physical space 
    representation if it's not already there.

    copy: if True, apply the transform to a copy

    Returns: skeleton in physical coordinates
    """
    self._check_space()

    skel = self
    if copy:
      skel = self.clone()

    if self.space != 'physical':
      skel.apply_transform()
      skel.space = 'physical'

    return skel

  def voxel_space(self, copy:bool = True) -> "Skeleton":
    """
    Convert skeleton vertices into a voxel space 
    representation if it's not already there.

    copy: if True, apply the transform to a copy

    Returns: skeleton in voxel coordinates
    """
    self._check_space()

    skel = self
    if copy:
      skel = self.clone()

    if self.space != 'voxel':
      skel.apply_inverse_transform()
      skel.space = 'voxel'
    
    return skel

  @property
  def transform(self):
    return self._transform

  @transform.setter 
  def transform(self, val):
    self._transform = np.array(val, dtype=np.float32).reshape( (3,4) )

  def transform_vertices(self, vertices, transform):
    verts = np.append(
      vertices,
      np.ones( (vertices.shape[0], 1), dtype=vertices.dtype), 
      axis=1
    )
    verts = transform.dot(verts.T).T
    return verts[:,0:3]    

  def apply_transform(self):
    self.vertices = self.transform_vertices(self.vertices, self.transform)

  def apply_inverse_transform(self, transform=None):
    if transform is None:
      transform = self.transform

    verts = np.append(
      self.vertices, 
      np.ones( (self.vertices.shape[0], 1), dtype=self.vertices.dtype), 
      axis=1
    )
    
    transform = np.zeros( (3,4), dtype=np.float32 )
    transform[:3,:3] = np.linalg.inv(self.transform[:3,:3])
    transform[:,3] = -self.transform[:,3]

    verts = transform.dot(verts.T).T
    self.vertices = verts[:,0:3]    

  @property 
  def radii(self):
    return self.radius

  @radii.setter 
  def radii(self, val):
    self.radius = val

  @classmethod
  def from_path(kls, vertices):
    """
    Given an Nx3 array of vertices that constitute a single path, 
    generate a skeleton with appropriate edges.
    """
    if vertices.shape[0] == 0:
      return Skeleton()

    skel = Skeleton(vertices)
    edges = np.zeros(shape=(skel.vertices.shape[0] - 1, 2), dtype=np.uint32)
    edges[:,0] = np.arange(skel.vertices.shape[0] - 1)
    edges[:,1] = np.arange(1, skel.vertices.shape[0])
    skel.edges = edges
    return skel

  @classmethod
  def simple_merge(kls, skeletons):
    """
    Simple concatenation of skeletons into one object 
    without adding edges between them.
    """
    if len(skeletons) == 0:
      return Skeleton()

    if type(skeletons[0]) is np.ndarray:
      skeletons = [ skeletons ]

    skeletons = [
      sk for sk in skeletons 
      if not sk.empty()
    ]

    if len(skeletons) == 0:
      return Skeleton()

    ct = 0
    edges = []
    for skel in skeletons:
      edge = skel.edges + ct
      edges.append(edge)
      ct += skel.vertices.shape[0]

    space = skeletons[0].space
    if space == "voxel":
      skeletons = [ s.voxel_space() for s in skeletons ]
    else:
      skeletons = [ s.physical_space() for s in skeletons ]

    skel = Skeleton(
      vertices=np.concatenate([ skel.vertices for skel in skeletons ], axis=0),
      edges=np.concatenate(edges, axis=0),
      segid=skeletons[0].id,
      extra_attributes=skeletons[0].extra_attributes,
      space=space,
      transform=skeletons[0].transform,
    )

    if len(skeletons) == 0:
      return skel

    first_extra_attr = skeletons[0].extra_attributes
    for skl in skeletons[1:]:
      if first_extra_attr != skl.extra_attributes:
        raise SkeletonAttributeMixingError("""
          The skeletons were unable to be merged because
          the extended vertex attributes were not uniformly
          defined.

          Template being matched against:
          {}

          Non-matching skeleton:
          {}
        """.format(first_extra_attr, skl.extra_attributes))

    for attr in skeletons[0].extra_attributes:
      setattr(skel, attr['id'], np.concatenate([
        getattr(skl, attr['id']) for skl in skeletons
      ], axis=0))

    return skel

  def terminals(self):
    """
    Returns vertex indices that correspond to terminal 
    nodes of the skeleton defined as having only one edge.
    """
    unique_nodes, unique_counts = fastremap.unique(self.edges, return_counts=True)
    return unique_nodes[ unique_counts == 1 ]

  def branches(self):
    """
    Returns vertex indices that correspond to branch points
    in the skeleton defined as having three or more edges.
    """
    unique_nodes, unique_counts = fastremap.unique(self.edges, return_counts=True)
    return unique_nodes[ unique_counts >= 3 ]

  def merge(self, skel):
    """Combine with an additional skeleton and consolidate."""
    return Skeleton.simple_merge((self, skel)).consolidate()

  def empty(self):
    return self.vertices.size == 0 or self.edges.size == 0

  def to_networkx(self):
    G = nx.Graph()
    G.add_edges_from(self.edges)
    return G

  def to_precomputed(self):
    edges = self.edges.astype(np.uint32)
    vertices = self.vertices.astype(np.float32)
    
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

    for attr in self.extra_attributes:
      arr = getattr(self, attr['id'])
      writeattr(arr, np.dtype(attr['data_type']), attr['id'])

    return result.getvalue()

  @classmethod
  def from_precomputed(kls, skelbuf, segid=None, vertex_attributes=None):
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

  @classmethod
  def equivalent(kls, first, second):
    """
    Tests that two skeletons are the same in form not merely that
    their array contents are exactly the same. This test can be
    made more sophisticated. 
    """
    if first.empty() and second.empty():
      return True
    elif first.vertices.shape[0] != second.vertices.shape[0]:
      return False
    elif first.edges.shape[0] != second.edges.shape[0]:
      return False

    EPSILON = 1e-7

    vertex1, ct1 = fastremap.unique(first.vertices, axis=0, return_counts=True)
    vertex2, ct2 = fastremap.unique(second.vertices, axis=0, return_counts=True)
    
    vertex_match = np.all(np.abs(vertex1 - vertex2) < EPSILON)
    ct_match = np.all(ct1 == ct2)
    if not (vertex_match and ct_match):
      return False

    g1 = nx.Graph()
    g1.add_edges_from(first.edges)
    g2 = nx.Graph()
    g2.add_edges_from(second.edges)
    edges_match = nx.is_isomorphic(g1, g2)
    del g1 
    del g2

    if not edges_match:
      return False

    second_verts = {}
    for i, vert in enumerate(second.vertices):
      second_verts[tuple(vert)] = i
    
    attrs = [ attr['id'] for attr in first.extra_attributes ]
    for attr in attrs:
      buf1 = getattr(first, attr)
      buf2 = getattr(second, attr)
      if len(buf1) != len(buf2):
        return False

      for i in range(len(buf1)):
        i2 = second_verts[tuple(first.vertices[i])]
        if buf1[i] != buf2[i2]:
          return False

    return True

  def crop(self, bbox):
    """
    Crop away all vertices and edges that lie outside of the given bbox.
    The edge counts as inside.

    Returns: new Skeleton
    """
    skeleton = self.clone()
    bbox = Bbox.create(bbox)

    if skeleton.empty():
      return skeleton

    nodes_valid_mask = (skeleton.vertices <= bbox.maxpt)
    nodes_valid_mask &= (bbox.minpt <= skeleton.vertices)
    nodes_valid_mask = (
      nodes_valid_mask[:,0] & nodes_valid_mask[:,1] & nodes_valid_mask[:,2]
    )

    # Set invalid vertices to be duplicates
    # so they'll be removed during consolidation
    if not np.any(nodes_valid_mask):
      empty_skel = Skeleton(extra_attributes=self.extra_attributes)
      empty_skel._materialize_extra_attributes()
      return empty_skel

    nodes_valid_idx = np.where(nodes_valid_mask)[0]

    first_node = nodes_valid_idx[0]
    skeleton.vertices[~nodes_valid_mask] = skeleton.vertices[first_node]
  
    edges_valid_mask = np.isin(skeleton.edges, nodes_valid_idx)
    edges_valid_idx = edges_valid_mask[:,0] & edges_valid_mask[:,1] 
    skeleton.edges = skeleton.edges[edges_valid_idx,:]
    return skeleton.consolidate()

  def consolidate(self, remove_disconnected_vertices=True):
    """
    Remove duplicate vertices and edges from this skeleton without
    side effects.

    Optional:
      remove_disconnected_vertices: delete vertices that have no edges
        associated with them. This does not preserve index order.

    Returns: new consolidated Skeleton 
    """
    nodes = self.vertices
    edges = self.edges 

    if self.empty():
      return Skeleton(
        segid=self.id, 
        space=self.space, 
        extra_attributes=self.extra_attributes,
        transform=self.transform,
      )
    
    eff_vertices, uniq_idx, idx_representative = fastremap.unique(
      nodes, axis=0, return_index=True, return_inverse=True
    )

    if idx_representative.ndim > 1:
      idx_representative = np.squeeze(idx_representative)

    eff_edges = idx_representative[ edges ]
    eff_edges = np.sort(eff_edges, axis=1) # sort each edge [2,1] => [1,2]
    eff_edges = eff_edges[np.lexsort(eff_edges[:,::-1].T)] # Sort rows 
    eff_edges = fastremap.unique(eff_edges, axis=0)
    eff_edges = eff_edges[ eff_edges[:,0] != eff_edges[:,1] ] # remove trivial loops

    skel = Skeleton(
      eff_vertices, eff_edges, 
      segid=self.id, 
      space=self.space,
      extra_attributes=self.extra_attributes,
      transform=self.transform,
    )

    for attr in self.extra_attributes:
      name = attr['id']
      buf = getattr(self, name)
      eff_name = buf[ uniq_idx ]
      setattr(skel, name, eff_name)

    if remove_disconnected_vertices:
      return skel.remove_disconnected_vertices()

    return skel

  def remove_disconnected_vertices(self):
    """
    Delete vertices that have no edges associated with them. 
    This does not preserve index order.

    Returns: new Skeleton
    """
    if self.empty():
      return Skeleton(
        segid=self.id, 
        space=self.space, 
        extra_attributes=self.extra_attributes,
        transform=self.transform,
      )

    all_edges = fastremap.unique(self.edges.flat)

    if len(self.vertices) == all_edges.size:
      return self.clone()

    idx_map = {}
    for i, vert in enumerate(self.vertices):
      idx_map[tuple(vert)] = i

    connected_verts = fastremap.unique(self.vertices[ self.edges.flat ], axis=0)

    edge_map = np.zeros( (len(self.vertices),), dtype=self.edges.dtype)
    vertex_remap = np.zeros( (len(self.vertices),), dtype=np.int32) - 1
    for i, vert in enumerate(connected_verts):
      reverse_idx = idx_map[tuple(vert)]
      edge_map[reverse_idx] = i
      vertex_remap[i] = reverse_idx

    edges = np.sort(edge_map[self.edges], axis=1)
    vertex_remap = vertex_remap[ vertex_remap > -1 ]

    skel = Skeleton(
      connected_verts, edges, 
      segid=self.id, 
      space=self.space,
      extra_attributes=self.extra_attributes,
      transform=self.transform,
    )

    if len(self.extra_attributes) == 0:
      return skel

    for attr in self.extra_attributes:
      name = attr['id']
      self_buf = getattr(self, name)
      skel_buf = self_buf[vertex_remap]
      setattr(skel, name, skel_buf)
        
    return skel

  def clone(self):
    vertices = np.copy(self.vertices)
    edges = np.copy(self.edges)
    radii = np.copy(self.radii)
    vertex_types = np.copy(self.vertex_types)

    skel = Skeleton(
      vertices, edges, radii, vertex_types, 
      segid=self.id, 
      space=self.space, 
      extra_attributes=self.extra_attributes,
      transform=np.copy(self.transform)
    )
    for attr in skel.extra_attributes:
      setattr(skel, attr['id'], np.copy(getattr(self, attr['id'])))

    return skel

  def cable_length(self):
    """
    Returns cable length of connected skeleton vertices in the same
    metric that this volume uses (typically nanometers).
    """
    skel = self.physical_space(copy=False)

    v1 = skel.vertices[skel.edges[:,0]]
    v2 = skel.vertices[skel.edges[:,1]]

    delta = (v2 - v1)
    delta *= delta
    dist = np.sum(delta, axis=1)
    dist = np.sqrt(dist)

    return np.sum(dist)

  def downsample(self, factor):
    """
    Compute a downsampled version of the skeleton by striding while 
    preserving endpoints.

    factor: stride length for downsampling the saved skeleton paths.

    Returns: downsampled Skeleton
    """
    if int(factor) != factor or factor < 1:
      raise ValueError("Argument `factor` must be a positive integer greater than or equal to 1. Got: <{}>({})", type(factor), factor)

    if factor == 1:
      return self.clone()

    paths = self.interjoint_paths()

    for i, path in enumerate(paths):
      paths[i] = np.concatenate(
        (path[0::factor, :], path[-1:, :]) # preserve endpoints
      )

    ds_skel = Skeleton.simple_merge(
      [ Skeleton.from_path(path) for path in paths ]
    ).consolidate()
    ds_skel.id = self.id
    ds_skel.extra_attributes = self.extra_attributes

    for attr in self.extra_attributes:
      setattr(ds_skel, attr['id'], 
        np.zeros([ds_skel.vertices.shape[0], int(attr['num_components'])], dtype=attr['data_type'])
      )

    # TODO: I'm sure this could be sped up if need be.
    index = {}
    for i, vert in enumerate(self.vertices):
      vert = tuple(vert)
      index[vert] = i

    bufs = [ getattr(ds_skel, attr['id']) for attr in self.extra_attributes ]
    orig_bufs = [ getattr(self, attr['id']) for attr in self.extra_attributes ]

    for i, vert in enumerate(ds_skel.vertices):
      reverse_i = index[tuple(vert)]
      for buf, buf_rev in zip(bufs, orig_bufs):
        buf[i] = buf_rev[reverse_i]
    
    return ds_skel

  def _single_tree_paths(self, tree, return_indices):
    """Get all traversal paths from a single tree."""
    skel = tree.consolidate()

    tree = defaultdict(list)

    for edge in skel.edges:
      svert = edge[0]
      evert = edge[1]
      tree[svert].append(evert)
      tree[evert].append(svert)

    def dfs(path, visited):
      paths = []
      stack = [ (path, visited) ]
      
      while stack:
        path, visited = stack.pop(0)

        vertex = path[-1]
        children = tree[vertex]
        
        visited[vertex] = True

        children = [ child for child in children if not visited[child] ]

        if len(children) == 0:
          paths.append(path)

        for child in children:
          stack.append( 
            (path + [child], np.copy(visited))
          )

      return paths

    size = np.max(skel.edges) + 1
    root = skel.edges[0,0]
    paths = dfs([root], np.zeros(size, dtype=bool))

    root = np.argmax([ len(_) for _ in paths ])
    root = paths[root][-1]
  
    paths = dfs([ root ], np.zeros(size, dtype=bool))

    if return_indices:
      return [ np.flip(path) for path in paths ]

    return [ np.flip(skel.vertices[path], axis=0) for path in paths ]

  def paths(self, return_indices=False):
    """
    Assuming the skeleton is structured as a single tree, return a 
    list of all traversal paths across all components. For each component, 
    start from the first vertex, find the most distant vertex by 
    hops and set that as the root. Then use depth first traversal 
    to produce paths.

    Returns: [ [(x,y,z), (x,y,z), ...], path_2, path_3, ... ]
    """
    paths = []
    for tree in self.components():
      paths += self._single_tree_paths(tree, return_indices=return_indices)
    return paths

  def _single_tree_interjoint_paths(self, skeleton, return_indices):
    vertices = skeleton.vertices
    edges = skeleton.edges

    unique_nodes, unique_counts = fastremap.unique(edges, return_counts=True)
    terminal_nodes = unique_nodes[ unique_counts == 1 ]
    branch_nodes = set(unique_nodes[ unique_counts >= 3 ])
    
    critical_points = set(terminal_nodes)
    critical_points.update(branch_nodes)

    tree = defaultdict(set)

    for e1, e2 in edges:
      tree[e1].add(e2)
      tree[e2].add(e1)

    # The below depth first search would be
    # more elegantly implemented as recursion,
    # but it quickly blows the stack, mandating
    # an iterative implementation.

    paths = []

    stack = [ terminal_nodes[0] ]
    criticals = [ terminal_nodes[0] ]
    # Saving the path stack is memory intensive
    # There might be a way to do it more linearly
    # via a DFS rather than BFS strategy.
    path_stack = [ [] ] 
    
    visited = defaultdict(bool)

    while stack:
      node = stack.pop()
      root = criticals.pop() # "root" is used v. loosely here
      path = path_stack.pop()

      path.append(node)
      visited[node] = True

      if node != root and node in critical_points:
        paths.append(path)
        path = [ node ]
        root = node

      for child in tree[node]:
        if not visited[child]:
          stack.append(child)
          criticals.append(root)
          path_stack.append(list(path))

    if return_indices:
      return paths

    return [ vertices[path] for path in paths ]

  def interjoint_paths(self, return_indices=False):
    """
    Returns paths between the adjacent critical points
    in the skeleton, where a critical point is the set of
    terminal and branch points.
    """
    paths = []
    for tree in self.components():
      subpaths = self._single_tree_interjoint_paths(
        tree, return_indices=return_indices
      )
      paths.extend(subpaths)

    return paths

  def _compute_components(self, skel):
    if skel.edges.size == 0:
      return [ skel ]

    index = defaultdict(set)
    visited = defaultdict(bool)
    for e1, e2 in skel.edges:
      index[e1].add(e2)
      index[e2].add(e1)

    def extract_component(start):
      edge_list = []
      stack = [ start ]
      parents = [ -1 ]

      while stack:
        node = stack.pop()
        parent = parents.pop()

        if node < parent:
          edge_list.append( (node, parent) )
        else:
          edge_list.append( (parent, node) )

        if visited[node]:
          continue

        visited[node] = True
        
        for child in index[node]:
          stack.append(child)
          parents.append(node)

      el = np.array(edge_list[1:], dtype=np.uint32)
      return fastremap.unique(el, axis=0)

    forest = []
    for edge in fastremap.unique(skel.edges.flat):
      if visited[edge]:
        continue

      forest.append(
        extract_component(edge)
      )

    return forest

  def average_smoothing(
    self, n:int, 
    check_boundary:bool = True,
    shrink_radii:bool = False,
  ):
    """
    Uses a moving window averaging filter to smooth
    each of the interjoint paths in the skeleton holding the
    ends fixed.

    Note: smoothing is applied in the forwards and then backwards
    direction to eliminate group delay induced by the filter.

    n: size of moving window in vertices, 1 being an exact
      copy of the oriignal object.

    check_boundary: if the skeleton has a "radii" attribute,
      check to make sure that the skeleton after smoothing is
      not outside that radius and raise an error if it does. 
      This ensures that the skeleton does not poke out of 
      the original object. 
    shrink_radii: Shrink radii by the amount the vertex is displaced
      from the original vertex. This will retain the ability to 
      reliably determine if the vertex is inside the object.
    """
    paths = self.interjoint_paths()

    if n <= 1:
      return self.clone()

    index = {}
    for i, vert in enumerate(self.vertices):
      vert = tuple(vert)
      index[vert] = i

    sub_skels = []

    for i, path in enumerate(paths):
      smooth_path = moving_average(path, n)
      smooth_path = moving_average(smooth_path[::-1], n)[::-1] # eliminate filter induced group delay
      smooth_path[0] = path[0]
      smooth_path[-1] = path[-1]

      sub_skel = Skeleton.from_path(smooth_path)
      sub_skel.add_extra_attributes(self.extra_attributes)

      bufs = [ getattr(sub_skel, attr['id']) for attr in self.extra_attributes ]
      orig_bufs = [ getattr(self, attr['id']) for attr in self.extra_attributes ]

      check_boundary = check_boundary and hasattr(self, "radii")

      for i, (vert, smooth_vert) in enumerate(zip(path, smooth_path)):
        reverse_i = index[tuple(vert)]

        dist = np.linalg.norm(vert - smooth_vert)

        if check_boundary:
          if dist > self.radii[reverse_i]:
            raise ValueError(
              f"Smoothing operation may have pushed one or more verticies "
              f"outside of the original object boundary.\n"
              f"Smoothed: {smooth_vert}\n"
              f"Original: {vert}\n"
              f"Allowed Radius: {self.radii[reverse_i]:.3f}\n"
              f"Total Radius: {dist}"
            )

        for buf, buf_rev in zip(bufs, orig_bufs):
          buf[i] = buf_rev[reverse_i]

        if shrink_radii:
          sub_skel.radii[i] -= dist

      sub_skels.append(sub_skel)

    smooth_skel = Skeleton.simple_merge(sub_skels).consolidate()
    smooth_skel.id = self.id
    return smooth_skel

  def components(self):
    """
    Extract connected components from graph. 
    Useful for ensuring that you're working with a single tree.

    Returns: [ Skeleton, Skeleton, ... ]
    """
    skel = self
    forest = fastosteoid.compute_components(skel.edges, skel.vertices.shape[0])
    
    if len(forest) == 0:
      return []
    elif len(forest) == 1:
      return [ skel.clone() ]

    skeletons = []
    for edge_list in forest:
      vert_idx = fastremap.unique(edge_list)

      vert_list = skel.vertices[vert_idx]

      if max(vert_idx[-1], vert_idx.size) < int(1e8):
        remap = np.zeros([ vert_idx[-1] + 1 ], dtype=np.uint32)
        remap[vert_idx] = np.arange(vert_idx.size)
        edge_list = remap[edge_list]
      else:
        remap = { vid: i for i, vid in enumerate(vert_idx) }
        edge_list = fastremap.remap(edge_list, remap, in_place=True)

      component_skel = Skeleton(
        vert_list, edge_list,# radii, vtypes, 
        segid=skel.id,
        extra_attributes=self.extra_attributes,
      )

      for attr in self.extra_attributes:
        vals = getattr(skel, attr['id'])[vert_idx]
        setattr(component_skel, attr['id'], vals)

      skeletons.append(component_skel)

    return skeletons

  @classmethod
  def from_navis(self, navis_skel):
    """
    Convert navis skeletons to osteoid.Skeleton. This
    should be more efficient than the SWC interchange method.
    """
    vertex_types = None
    if len(navis_skel.nodes.type):
      # 'root', 'slab', 'branch', 'end', first letter of each
      vertex_types = [ ord(t[0]) for t in navis_skel.nodes.type ] 
      mapping = np.zeros(ord('t') + 1, dtype=np.uint8)
      mapping[ord('s')] = 7 # custom
      mapping[ord('r')] = 8 # custom
      mapping[ord('b')] = 5 # fork point
      mapping[ord('t')] = 6 # end point
      vertex_types = mapping[vertex_types]

    return Skeleton(
      vertices=navis_skel.vertices, 
      edges=navis_skel.edges - 1,
      radii=navis_skel.nodes.radius.to_numpy(),
      vertex_types=vertex_types,
    )

  @classmethod
  def from_swc(self, swcstr):
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

  def to_swc(self, contributors:str = "", soma_threshold:float = np.inf) -> str:
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
    sx, sy, sz = np.diag(self.transform)[:3]

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
      if soma_threshold < np.inf and hasattr(skel, 'radii'):
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

    skels = self.components()
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

  def viewer(
    self, units='nm', 
    draw_edges=True, draw_vertices=True,
    color_by='radius', library='auto'
  ):
    """
    View the skeleton with a radius heatmap. 

    Requires either the matplotlib library or the
    microviewer+vtk libraries which are not installed
    by default.

    units: label axes with these units
    draw_edges: draw lines between vertices (more useful when skeleton is sparse)
    draw_vertices: draw each vertex colored by its radius.
    color_by: 
      'radius': color each vertex according to its radius attribute
        aliases: 'r', 'radius', 'radii'
      'component': color connected components seperately
        aliases: 'c', 'component', 'components'
      'cross_section': color each vertex according to its cross sectional area
        aliases: 'x'
      anything else: draw everything black
    library: 'auto', 'matplotlib', or 'microviewer'
      'auto' will try to detect which libraries are installed
        and use microviewer and then matplotlib in that order.
    """
    from .viewer import view_matplotlib
    if library == "auto":
      try:
        import microviewer
        import vtk
        microviewer.objects([self], skeleton_color_by=color_by)
      except ImportError:
        import matplotlib
        view_matplotlib(self, units, draw_edges, draw_vertices, color_by)
    elif library == "matplotlib":
      view_matplotlib(self, units, draw_edges, draw_vertices, color_by)
    elif library in "microviewer":
      import microviewer
      microviewer.objects([self])
    else:
      raise ValueError(f"{library} is not currently supported.")

  def __eq__(self, other):
    if self.id != other.id:
      return False
    elif self.vertices.shape[0] != other.vertices.shape[0]:
      return False
    elif self.edges.shape[0] != other.edges.shape[0]:
      return False
    elif self.extra_attributes != other.extra_attributes:
      return False

    attrs = [ attr['id'] for attr in self.extra_attributes ]
    for attr in attrs:
      buf1 = getattr(self, attr)
      buf2 = getattr(other, attr)
      if np.all(buf1 != buf2):
        return False

    return (np.all(self.vertices == other.vertices)
      and np.all(self.edges == other.edges) \
      and np.all(self.radii == other.radii) \
      and np.all(self.vertex_types == other.vertex_types))

  def chunk(
    self, 
    chunk_size:tuple[float,float,float],
    origin:Optional[np.ndarray] = None,
  ) -> dict[tuple[int,int,int], "Skeleton"]:
    """
    Cut a skeleton into a grid.
    """
    vertices = self.vertices.astype(np.float32, copy=False)
    
    if origin is None:
      origin = [None, None, None]

    chunks = fastosteoid.chunk_skeleton(
      vertices, self.edges,
      chunk_size[0], chunk_size[1], chunk_size[2],
      origin[0], origin[1], origin[2],
    )
    del vertices

    skel_chunks = {}
    for grid, (verts, edges) in chunks.items():
      sk = Skeleton(verts, edges, default_attributes=False).consolidate()
      sk.import_attributes(self)
      skel_chunks[grid] = sk

    return skel_chunks

  def import_attributes(self, other:"Skeleton"):
    """
    Copy the attributes of another skeleton
    for matching vertices.
    
    Vertices are considered matching if they have the same coordinates.
    
    Assumes importing skeleton has no extra attributes initially,
    only vertices and edges.
    
    Parameters:
    -----------
    other : Skeleton
        The skeleton to copy attributes from
    """
    if self.vertices.shape[0] == 0 or other.vertices.shape[0] == 0:
      return

    other_dict = {
      tuple(v): i
      for i, v in enumerate(other.vertices)
    }
    
    other_indices = np.array([
      other_dict[tuple(v)]
      for i, v in enumerate(self.vertices)
      if tuple(v) in other_dict
    ], dtype=np.uint64)

    if len(other_indices) == 0:
      return

    # Copy attribute data for matching vertices
    for attr in other.extra_attributes:
      attr_name = attr['id']
      other_buf = getattr(other, attr_name)
      
      buf = other_buf[other_indices]
      self.add_vertex_attribute(attr_name, buf)
      setattr(self, attr_name, buf)

  def __str__(self):
    template = "{}=({}, {})"
    attr_strings = []
    for attr in self.extra_attributes:
      attr = attr['id']
      buf = getattr(self, attr)
      attr_strings.append(
        template.format(attr, buf.shape[0], buf.dtype)
      )

    return "Skeleton(segid={}, vertices=(shape={}, {}), edges=(shape={}, {}), {}, space='{}' transform={})".format(
      self.id,
      self.vertices.shape[0], self.vertices.dtype,
      self.edges.shape[0], self.edges.dtype,
      ', '.join(attr_strings),
      self.space, self.transform.tolist()
    )

  def __repr__(self):
    return str(self)

  def __sizeof__(self):
    npy_attrs = [ getattr(self, attr['id']) for attr in self.extra_attributes ]
    npy_attrs += [ self.vertices, self.edges, self.transform ]

    npy_bytes = sum([ x.nbytes for x in npy_attrs ])

    attrs = [
      self.id, self.space, self.extra_attributes
    ]
    attr_bytes = sum(( 
      sys.getsizeof(x) for x in attrs
    ))

    return npy_bytes + attr_bytes
