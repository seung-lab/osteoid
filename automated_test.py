import pytest

import copy
import math
import numpy as np

import osteoid
from osteoid import Skeleton, Bbox
from osteoid.exceptions import SkeletonDecodeError, SkeletonAttributeMixingError

def test_consolidate():
  skel = Skeleton(
    vertices=np.array([
      (0, 0, 0),
      (1, 0, 0),
      (2, 0, 0),
      (0, 0, 0),
      (2, 1, 0),
      (2, 2, 0),
      (2, 2, 1),
      (2, 2, 2),
    ], dtype=np.float32),

    edges=np.array([
      [0, 1],
      [1, 2],
      [2, 3],
      [3, 4],
      [4, 5],
      [5, 6],
      [6, 7],
    ], dtype=np.uint32),

    radii=np.array([
      0, 1, 2, 3, 4, 5, 6, 7
    ], dtype=np.float32),

    vertex_types=np.array([
      0, 1, 2, 3, 4, 5, 6, 7
    ], dtype=np.uint8),
  )

  correct_skel = Skeleton(
    vertices=np.array([
      (0, 0, 0),
      (1, 0, 0),
      (2, 0, 0),
      (2, 1, 0),
      (2, 2, 0),
      (2, 2, 1),
      (2, 2, 2),
    ], dtype=np.float32),

    edges=np.array([
      [0, 1],
      [0, 2],
      [0, 3],
      [1, 2],
      [3, 4],
      [4, 5],
      [5, 6],
    ], dtype=np.uint32),

    radii=np.array([
      0, 1, 2, 4, 5, 6, 7
    ], dtype=np.float32),

    vertex_types=np.array([
      0, 1, 2, 4, 5, 6, 7
    ], dtype=np.uint8),
  )

  consolidated = skel.consolidate()

  assert np.all(consolidated.vertices == correct_skel.vertices)
  assert np.all(consolidated.edges == correct_skel.edges)
  assert np.all(consolidated.radii == correct_skel.radii)
  assert np.all(consolidated.vertex_types == correct_skel.vertex_types)

def test_remove_disconnected_vertices():
  skel = Skeleton(
    [ 
      (0,0,0), (1,0,0), (2,0,0),
      (0,1,0), (0,2,0), (0,3,0),
      (-1, -1, -1)
    ], 
    edges=[ 
      (0,1), (1,2), 
      (3,4), (4,5), (3,5)
    ],
    segid=666,
  )

  res = skel.remove_disconnected_vertices()
  assert res.vertices.shape[0] == 6
  assert res.edges.shape[0] == 5 
  assert res.radii.shape[0] == 6
  assert res.vertex_types.shape[0] == 6
  assert res.id == 666


def test_equivalent():
  assert Skeleton.equivalent(Skeleton(), Skeleton())

  identity = Skeleton([ (0,0,0), (1,0,0) ], [(0,1)] )
  assert Skeleton.equivalent(identity, identity)

  diffvertex = Skeleton([ (0,0,0), (0,1,0) ], [(0,1)])
  assert not Skeleton.equivalent(identity, diffvertex)

  single1 = Skeleton([ (0,0,0), (1,0,0) ], edges=[ (1,0) ])
  single2 = Skeleton([ (0,0,0), (1,0,0) ], edges=[ (0,1) ])
  assert Skeleton.equivalent(single1, single2)

  double1 = Skeleton([ (0,0,0), (1,0,0) ], edges=[ (1,0) ])
  double2 = Skeleton([ (0,0,0), (1,0,0) ], edges=[ (0,1) ])
  assert Skeleton.equivalent(double1, double2)

  double1 = Skeleton([ (0,0,0), (1,0,0), (1,1,0) ], edges=[ (1,0), (1,2) ])
  double2 = Skeleton([ (0,0,0), (1,0,0), (1,1,0) ], edges=[ (2,1), (0,1) ])
  assert Skeleton.equivalent(double1, double2)

  double1 = Skeleton([ (0,0,0), (1,0,0), (1,1,0), (1,1,3) ], edges=[ (1,0), (1,2), (1,3) ])
  double2 = Skeleton([ (0,0,0), (1,0,0), (1,1,0), (1,1,3) ], edges=[ (3,1), (2,1), (0,1) ])
  assert Skeleton.equivalent(double1, double2)

def test_cable_length():
  skel = Skeleton([ 
      (0,0,0), (1,0,0), (2,0,0), (3,0,0), (4,0,0), (5,0,0)
    ], 
    edges=[ (1,0), (1,2), (2,3), (3,4), (5,4) ],
    radii=[ 1, 2, 3, 4, 5, 6 ],
    vertex_types=[1, 2, 3, 4, 5, 6]
  )

  assert skel.cable_length() == (skel.vertices.shape[0] - 1)

  skel = Skeleton([ 
      (2,0,0), (1,0,0), (0,0,0), (0,5,0), (0,6,0), (0,7,0)
    ], 
    edges=[ (1,0), (1,2), (2,3), (3,4), (5,4) ],
    radii=[ 1, 2, 3, 4, 5, 6 ],
    vertex_types=[1, 2, 3, 4, 5, 6]
  )
  assert skel.cable_length() == 9

  skel = Skeleton([ 
      (1,1,1), (0,0,0), (1,0,0)
    ], 
    edges=[ (1,0), (1,2) ],
    radii=[ 1, 2, 3],
    vertex_types=[1, 2, 3]
  )
  assert abs(skel.cable_length() - (math.sqrt(3) + 1)) < 1e-6

def test_transform():
  skelv = Skeleton([ 
      (0,0,0), (1,0,0), (1,1,0), (1,1,3), (2,1,3), (2,2,3)
    ], 
    edges=[ (1,0), (1,2), (2,3), (3,4), (5,4) ],
    radii=[ 1, 2, 3, 4, 5, 6 ],
    vertex_types=[1, 2, 3, 4, 5, 6],
    segid=1337,
    transform=np.array([
      [2, 0, 0, 0],
      [0, 2, 0, 0],
      [0, 0, 2, 0],
    ])
  )

  skelp = skelv.physical_space()
  assert np.all(skelp.vertices == skelv.vertices * 2)
  assert np.all(skelv.vertices == skelp.voxel_space().vertices)

  skelv.transform = [
    [1, 0, 0, 1],
    [0, 1, 0, 2],
    [0, 0, 1, 3],
  ]

  skelp = skelv.physical_space()
  tmpskel = skelv.clone() 
  tmpskel.vertices[:,0] += 1
  tmpskel.vertices[:,1] += 2
  tmpskel.vertices[:,2] += 3
  assert np.all(skelp.vertices == tmpskel.vertices)
  assert np.all(skelp.voxel_space().vertices == skelv.vertices)


def test_downsample():
  skel = Skeleton([ 
      (0,0,0), (1,0,0), (1,1,0), (1,1,3), (2,1,3), (2,2,3)
    ], 
    edges=[ (1,0), (1,2), (2,3), (3,4), (5,4) ],
    radii=[ 1, 2, 3, 4, 5, 6 ],
    vertex_types=[1, 2, 3, 4, 5, 6],
    segid=1337,
  )

  def should_error(x):
    try:
      skel.downsample(x)
      assert False
    except ValueError:
      pass

  should_error(-1)
  should_error(0)
  should_error(.5)
  should_error(2.00000000000001)

  dskel = skel.downsample(1)
  assert Skeleton.equivalent(dskel, skel)
  assert dskel.id == skel.id
  assert dskel.id == 1337

  dskel = skel.downsample(2)
  dskel_gt = Skeleton(
    [ (0,0,0), (1,1,0), (2,1,3), (2,2,3) ], 
    edges=[ (1,0), (1,2), (2,3) ],
    radii=[1,3,5,6], vertex_types=[1,3,5,6] 
  )
  assert Skeleton.equivalent(dskel, dskel_gt)

  dskel = skel.downsample(3)
  dskel_gt = Skeleton(
    [ (0,0,0), (1,1,3), (2,2,3) ], edges=[ (1,0), (1,2) ],
    radii=[1,4,6], vertex_types=[1,4,6],
  )
  assert Skeleton.equivalent(dskel, dskel_gt)

  skel = Skeleton([ 
      (0,0,0), (1,0,0), (1,1,0), (1,1,3), (2,1,3), (2,2,3)
    ], 
    edges=[ (1,0), (1,2), (3,4), (5,4) ],
    radii=[ 1, 2, 3, 4, 5, 6 ],
    vertex_types=[1, 2, 3, 4, 5, 6]
  )
  dskel = skel.downsample(2)
  dskel_gt = Skeleton(
    [ (0,0,0), (1,1,0), (1,1,3), (2,2,3) ], 
    edges=[ (1,0), (2,3) ],
    radii=[1,3,4,6], vertex_types=[1,3,4,6] 
  )
  assert Skeleton.equivalent(dskel, dskel_gt)


def test_downsample_joints():
  skel = Skeleton([ 
      
                        (2, 3,0), # 0
                        (2, 2,0), # 1
                        (2, 1,0), # 2
      (0,0,0), (1,0,0), (2, 0,0), (3,0,0), (4,0,0), # 3, 4, 5, 6, 7
                        (2,-1,0), # 8
                        (2,-2,0), # 9
                        (2,-3,0), # 10

    ], 
    edges=[ 
                  (0, 1),
                  (1, 2),
                  (2, 5),
        (3,4), (4,5), (5, 6), (6,7),
                  (5, 8),
                  (8, 9),
                  (9,10)
    ],
    radii=[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
    vertex_types=[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
    segid=1337,
  )

  ds_skel = skel.downsample(2)
  ds_skel_gt = Skeleton([ 

                        (2, 3,0), # 0
                        
                        (2, 1,0), # 1
      (0,0,0),          (2, 0,0),     (4,0,0), # 2, 3, 4

                        (2,-2,0), # 5                        
                        (2,-3,0), # 6

    ], 
    edges=[ 
                  (0,1),
                  (1,3),
              (2,3),  (3,4), 
                  (3,5),
                  (5,6)
    ],
    radii=[ 0, 2, 3, 5, 7, 9, 10 ],
    vertex_types=[ 0, 2, 3, 5, 7, 9, 10 ],
    segid=1337,
  )

  assert Skeleton.equivalent(ds_skel, ds_skel_gt)


def test_read_swc():

  # From http://research.mssm.edu/cnic/swc.html
  test_file = """# ORIGINAL_SOURCE NeuronStudio 0.8.80
# CREATURE
# REGION
# FIELD/LAYER
# TYPE
# CONTRIBUTOR
# REFERENCE
# RAW
# EXTRAS
# SOMA_AREA
# SHINKAGE_CORRECTION 1.0 1.0 1.0
# VERSION_NUMBER 1.0
# VERSION_DATE 2007-07-24
# SCALE 1.0 1.0 1.0
1 1 14.566132 34.873772 7.857000 0.717830 -1
2 0 16.022520 33.760513 7.047000 0.463378 1
3 5 17.542000 32.604973 6.885001 0.638007 2
4 0 19.163984 32.022469 5.913000 0.602284 3
5 0 20.448090 30.822802 4.860000 0.436025 4
6 6 21.897903 28.881084 3.402000 0.471886 5
7 0 18.461960 30.289471 8.586000 0.447463 3
8 6 19.420759 28.730757 9.558000 0.496217 7"""

  skel = Skeleton.from_swc(test_file)
  assert skel.vertices.shape[0] == 8
  assert skel.edges.shape[0] == 7

  skel_gt = Skeleton(
    vertices=[
      [14.566132, 34.873772, 7.857000],
      [16.022520, 33.760513, 7.047000],
      [17.542000, 32.604973, 6.885001],
      [19.163984, 32.022469, 5.913000],
      [20.448090, 30.822802, 4.860000],
      [21.897903, 28.881084, 3.402000],
      [18.461960, 30.289471, 8.586000],
      [19.420759, 28.730757, 9.558000]
    ],
    edges=[ (0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (7,6) ],
    radii=[ 
      0.717830, 0.463378, 0.638007, 0.602284, 
      0.436025, 0.471886, 0.447463, 0.496217
    ],
    vertex_types=[
      1, 0, 5, 0, 0, 6, 0, 6
    ],
  )

  assert Skeleton.equivalent(skel, skel_gt)

  skel = Skeleton.from_swc(skel.to_swc())
  assert Skeleton.equivalent(skel, skel_gt)

  # sorts edges
  skel = skel.consolidate()
  skel_gt = skel_gt.consolidate()
  assert np.all(skel.edges == skel_gt.edges)
  assert np.all(np.abs(skel.radii - skel_gt.radii) < 0.00001)

  Nv = skel.vertices.shape[0]
  Ne = skel.edges.shape[0]

  for _ in range(10):
    skel = Skeleton.from_swc(skel.to_swc())
    assert skel.vertices.shape[0] == Nv 
    assert skel.edges.shape[0] == Ne

def test_read_duplicate_vertex_swc():
  test_file = """
1 0 -18.458370 23.227150 -84.035016 1.000000 -1
2 0 -18.159709 22.925778 -82.984344 1.000000 1
3 0 -17.861047 22.624407 -82.984344 1.000000 2
4 0 -17.562385 22.624407 -82.984344 1.000000 3
5 0 -16.965061 22.021663 -82.984344 1.000000 4
6 0 -16.965061 21.720292 -82.984344 1.000000 5
7 0 -16.069075 21.720292 -82.984344 1.000000 6
8 0 -16.069075 21.117548 -80.883000 1.000000 7
9 0 -15.770414 20.816176 -80.883000 1.000000 8
10 0 -15.770414 20.514805 -80.883000 1.000000 9
11 0 -15.770414 20.816176 -80.883000 1.000000 10
12 0 -16.069075 21.117548 -80.883000 1.000000 11
13 0 -16.069075 21.418920 -80.883000 1.000000 12
14 0 -16.069075 20.816176 -78.781655 1.000000 13
15 0 -15.471752 20.213433 -76.680311 1.000000 14
16 0 -15.471752 19.309318 -76.680311 1.000000 15
17 0 -15.471752 19.007946 -75.629639 1.000000 16
18 0 -15.173090 18.706574 -74.578966 1.000000 17
19 0 -14.874428 18.706574 -74.578966 1.000000 18
20 0 -14.575766 18.405202 -74.578966 1.000000 19
"""

  skel = Skeleton.from_swc(test_file)
  assert skel.vertices.shape[0] == 20
  skel2 = Skeleton.from_swc(skel.to_swc())
  assert skel2.vertices.shape[0] == 20
  assert Skeleton.equivalent(skel, skel2)

def test_components():
  skel = Skeleton(
    [ 
      (0,0,0), (1,0,0), (2,0,0),
      (0,1,0), (0,2,0), (0,3,0),
    ], 
    edges=[ 
      (0,1), (1,2), 
      (3,4), (4,5), (3,5)
    ],
    segid=666,
  )

  components = skel.components()
  assert len(components) == 2
  assert components[0].vertices.shape[0] == 3
  assert components[1].vertices.shape[0] == 3
  assert components[0].edges.shape[0] == 2
  assert components[1].edges.shape[0] == 3

  skel1_gt = Skeleton([(0,0,0), (1,0,0), (2,0,0)], [(0,1), (1,2)])
  skel2_gt = Skeleton([(0,1,0), (0,2,0), (0,3,0)], [(0,1), (0,2), (1,2)])

  assert Skeleton.equivalent(components[0], skel1_gt)
  assert Skeleton.equivalent(components[1], skel2_gt)

def test_simple_merge():
  skel1 = Skeleton(
    [ (0,0,0), (1,0,0), (2,0,0),  ], 
    edges=[ (0,1), (1,2), ],
    segid=1,
  )

  skel2 = Skeleton(
    [ (0,0,1), (1,0,2), (2,0,3),  ], 
    edges=[ (0,1), (1,2), ],
    segid=1,
  )

  result = Skeleton.simple_merge([ skel1, skel2 ])

  expected = Skeleton(
    [ (0,0,0), (1,0,0), (2,0,0), (0,0,1), (1,0,2), (2,0,3), ], 
    edges=[ (0,1), (1,2), (3,4), (4,5) ],
    segid=1,
  )

  assert result == expected

  wow_attr = {
    "id": "wow",
    "data_type": "uint8",
    "components": 1,
  }

  skel1.extra_attributes = [copy.deepcopy(wow_attr)]
  skel1.wow = np.array([1,2,3], dtype=np.uint8)

  skel2.extra_attributes = [copy.deepcopy(wow_attr)]
  skel2.wow = np.array([4,5,6], dtype=np.uint8)

  result = Skeleton.simple_merge([ skel1, skel2 ])
  expected.extra_attributes = [copy.deepcopy(wow_attr)]
  expected.wow = np.array([1,2,3,4,5,6], dtype=np.uint8)

  assert result == expected

  skel2.extra_attributes[0]['data_type'] = np.uint8

  try:
    Skeleton.simple_merge([ skel1, skel2 ])
    assert False
  except SkeletonAttributeMixingError:
    pass

  skel2.extra_attributes[0]['data_type'] = 'uint8'
  skel2.extra_attributes.append({
    "id": "amaze",
    "data_type": "float32",
    "components": 2,
  })
  skel2.amaze = np.array([])

  try:
    Skeleton.simple_merge([ skel1, skel2 ])
    assert False
  except SkeletonAttributeMixingError:
    pass

def test_paths():
  skel1 = Skeleton(
    [ (0,0,0), (1,0,0), (2,0,0), (1,1,0) ], 
    edges=[ (0,1), (1,2), (1,3) ],
    segid=1,
  )

  path = skel1.paths()

  path1 = np.array([[0,0,0], [1,0,0], [1,1,0]])
  path2 = np.array([[2,0,0], [1,0,0], [1,1,0]])

  assert np.all(path[0] == path1)
  assert np.all(path[1] == path2)

def test_crop():
  skel = Skeleton(
    [ (0,0,0), (10,0,1), (20,0,2), (30,1,3) ], 
    edges=[ (0,1), (1,2), (1,3) ],
    segid=1,
  )

  bbx = Bbox([5,-1,0], [25,5,3])
  res = skel.crop(bbx)

  assert np.all(res.vertices == np.array([(10,0,1), (20,0,2)]))

  bbx = Bbox([100,-1,0], [200,5,3])
  res = skel.crop(bbx)
  assert len(res.vertices) == 0

  bbx = Bbox([-100,-100,-100], [200,500,400])
  res = skel.crop(bbx)
  assert np.all(res.vertices == skel.vertices)

class TestSkeletonChunking:
  def test_single_chunk_small_skeleton(self):
    """Test that a small skeleton fits in a single chunk."""
    vertices = np.array([
      [0.0, 0.0, 0.0],
      [1.0, 0.0, 0.0],
      [1.0, 1.0, 0.0],
    ], dtype=np.float32)
    edges = np.array([
      [0, 1],
      [1, 2],
    ], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(10.0, 10.0, 10.0))
    
    assert len(chunks) == 1
    assert (0, 0, 0) in chunks
    assert len(chunks[(0, 0, 0)].vertices) > 0
  
  def test_multiple_chunks_2x2_grid(self):
    """Test skeleton split across a 2x2x1 grid."""
    vertices = np.array([
      [0.0, 0.0, 0.0],
      [1.5, 0.0, 0.0],
      [0.0, 1.5, 0.0],
      [1.5, 1.5, 0.0],
    ], dtype=np.float32)
    edges = np.array([
      [0, 1],
      [0, 2],
      [1, 3],
      [2, 3],
    ], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(1.0, 1.0, 1.0))
    
    # Should have chunks in multiple grid cells
    assert len(chunks) > 1
    
    # All chunks should be valid skeletons
    for chunk_skel in chunks.values():
      assert len(chunk_skel.vertices) > 0
      assert len(chunk_skel.edges) > 0
  
  def test_edge_crossing_boundary(self):
    """Test that edges crossing chunk boundaries are split."""
    vertices = np.array([
      [0.5, 0.5, 0.5],
      [1.5, 0.5, 0.5],  # Crosses boundary at x=1.0
    ], dtype=np.float32)
    edges = np.array([[0, 1]], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(
      chunk_size=(1.0, 1.0, 1.0),
      origin=(0,0,0),
    )
    
    # Should have 2 chunks
    assert len(chunks) == 2
    assert (0, 0, 0) in chunks
    assert (1, 0, 0) in chunks

    sk0 = chunks[(0, 0, 0)]
    sk1 = chunks[(1, 0, 0)]
    assert list(sk0.edges[0]) == [0,1]
    assert list(sk1.edges[0]) == [0,1]

    assert list(sk0.vertices[1]) == [1,0.5,0.5]
    assert list(sk1.vertices[0]) == [1,0.5,0.5]
    
  def test_custom_origin(self):
    """Test chunking with custom origin."""
    vertices = np.array([
      [5.0, 5.0, 5.0],
      [6.0, 5.0, 5.0],
    ], dtype=np.float32)
    edges = np.array([[0, 1]], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    
    # # With default origin (should be at min vertex)
    chunks_default = skeleton.chunk(chunk_size=(1.0, 1.0, 1.0))
    
    # With custom origin at (0, 0, 0)
    chunks_custom = skeleton.chunk(
      chunk_size=(1.0, 1.0, 1.0),
      origin=np.array([0.0, 0.0, 0.0])
    )
    
    # Keys should be different due to different grid alignment
    assert set(chunks_default.keys()) != set(chunks_custom.keys())
    
    # With custom origin, should be in grid cell (5, 5, 5)
    assert (5, 4, 4) in chunks_custom
  
  def test_empty_skeleton(self):
    """Test chunking an empty skeleton."""
    vertices = np.array([], dtype=np.float32).reshape(0, 3)
    edges = np.array([], dtype=np.uint32).reshape(0, 2)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(1.0, 1.0, 1.0))
    
    # Should return empty dict or single empty chunk
    assert len(chunks) == 0 or all(
      len(c.vertices) == 0 for c in chunks.values()
    )
  
  def test_vertex_on_boundary(self):
    """Test vertices exactly on chunk boundaries."""
    vertices = np.array([
      [0.0, 0.0, 0.0],
      [1.0, 0.0, 0.0],  # Exactly on boundary
      [2.0, 0.0, 0.0],
    ], dtype=np.float32)
    edges = np.array([
      [0, 1],
      [1, 2],
    ], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(1.0, 1.0, 1.0))
    
    # Should handle boundary vertices consistently
    assert len(chunks) >= 2
    
    # Total edges across all chunks should match or exceed original
    # (may be more due to splitting)
    total_edges = sum(len(c.edges) for c in chunks.values())
    assert total_edges >= len(edges)
  
  def test_3d_grid(self):
    """Test chunking across all three dimensions."""
    vertices = np.array([
      [0.5, 0.5, 0.5],
      [1.5, 1.5, 1.5],
    ], dtype=np.float32)
    edges = np.array([[0, 1]], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(
      chunk_size=(1.0, 1.0, 1.0),
      origin=(0,0,0),
    )
    
    # Diagonal line should cross through multiple chunks
    assert len(chunks) >= 2
    
    # Check that chunks are in different grid positions
    keys = list(chunks.keys())
    assert keys[0] != keys[-1]
  
  def test_non_uniform_chunk_size(self):
    """Test with different chunk sizes in each dimension."""
    vertices = np.array([
      [0.0, 0.0, 0.0],
      [2.0, 4.0, 8.0],
    ], dtype=np.float32)
    edges = np.array([[0, 1]], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(1.0, 2.0, 4.0))
    
    assert len(chunks) > 0
    
    # Verify grid keys make sense with non-uniform sizes
    for (gx, gy, gz) in chunks.keys():
      assert gx >= 0 and gx <= 2  # 2.0 / 1.0 = 2 chunks in x
      assert gy >= 0 and gy <= 2  # 4.0 / 2.0 = 2 chunks in y
      assert gz >= 0 and gz <= 2  # 8.0 / 4.0 = 2 chunks in z
  
  def test_large_chunk_size(self):
    """Test with chunk size larger than skeleton bounds."""
    vertices = np.array([
      [0.0, 0.0, 0.0],
      [1.0, 1.0, 1.0],
    ], dtype=np.float32)
    edges = np.array([[0, 1]], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(100.0, 100.0, 100.0))
    
    # Everything should fit in one chunk
    assert len(chunks) == 1
  
  def test_vertex_preservation(self):
    """Test that all original vertices are preserved across chunks."""
    vertices = np.array([
      [0.0, 0.0, 0.0],
      [0.5, 0.5, 0.5],
      [1.5, 0.5, 0.5],
      [2.0, 1.0, 1.0],
    ], dtype=np.float32)
    edges = np.array([
      [0, 1],
      [1, 2],
      [2, 3],
    ], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(1.0, 1.0, 1.0))
    
    # Collect all unique vertices from chunks
    all_chunk_vertices = []
    for chunk_skel in chunks.values():
      all_chunk_vertices.append(chunk_skel.vertices)
    
    # Should have at least as many vertices as original
    # (may have more due to splitting at boundaries)
    total_verts = sum(len(v) for v in all_chunk_vertices)
    assert total_verts >= len(vertices)
  
  def test_connectivity_within_chunks(self):
    """Test that edges within chunks reference valid vertices."""
    vertices = np.array([
      [0.0, 0.0, 0.0],
      [0.5, 0.0, 0.0],
      [1.0, 0.0, 0.0],
      [1.5, 0.0, 0.0],
    ], dtype=np.float32)
    edges = np.array([
      [0, 1],
      [1, 2],
      [2, 3],
    ], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(1.0, 1.0, 1.0))
    
    # Verify all edges reference valid vertex indices
    for chunk_skel in chunks.values():
      num_verts = len(chunk_skel.vertices)
      for edge in chunk_skel.edges:
        assert edge[0] < num_verts
        assert edge[1] < num_verts
        assert edge[0] >= 0
        assert edge[1] >= 0
  
  def test_different_edge_dtypes(self):
    """Test with different edge array data types."""
    vertices = np.array([
      [0.0, 0.0, 0.0],
      [1.0, 1.0, 1.0],
    ], dtype=np.float32)
    
    for dtype in [np.uint8, np.uint16, np.uint32, np.uint64]:
      edges = np.array([[0, 1]], dtype=dtype)
      skeleton = Skeleton(vertices, edges)
      chunks = skeleton.chunk(chunk_size=(1.0, 1.0, 1.0))
      
      assert len(chunks) > 0, f"Failed for dtype {dtype}"
  
  def test_negative_coordinates(self):
    """Test chunking with negative coordinates."""
    vertices = np.array([
      [-1.0, -1.0, -1.0],
      [1.0, 1.0, 1.0],
    ], dtype=np.float32)
    edges = np.array([[0, 1]], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(1.0, 1.0, 1.0))
    
    assert len(chunks) > 0
    
    # With default origin at min vertex, grid indices should be >= 0
    for key in chunks.keys():
      assert all(coord >= 0 for coord in key)
  
  def test_star_topology(self):
    """Test skeleton with star topology (one vertex connected to many)."""
    center = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
    spokes = np.array([
      [2.0, 0.0, 0.0],
      [0.0, 2.0, 0.0],
      [0.0, 0.0, 2.0],
      [-2.0, 0.0, 0.0],
      [0.0, -2.0, 0.0],
      [0.0, 0.0, -2.0],
    ], dtype=np.float32)
    
    vertices = np.vstack([center, spokes])
    edges = np.array([
      [0, i] for i in range(1, 7)
    ], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=(1.0, 1.0, 1.0))
    
    assert len(chunks) > 1
  
  @pytest.mark.parametrize("chunk_size", [
    (0.5, 0.5, 0.5),
    (1.0, 1.0, 1.0),
    (2.0, 2.0, 2.0),
    (10.0, 10.0, 10.0),
  ])
  def test_various_chunk_sizes(self, chunk_size):
    """Test with various chunk sizes."""
    vertices = np.array([
      [0.0, 0.0, 0.0],
      [5.0, 5.0, 5.0],
    ], dtype=np.float32)
    edges = np.array([[0, 1]], dtype=np.uint32)
    
    skeleton = Skeleton(vertices, edges)
    chunks = skeleton.chunk(chunk_size=chunk_size)
    
    assert len(chunks) > 0
    assert all(isinstance(k, tuple) and len(k) == 3 for k in chunks.keys())
    assert all(isinstance(v, Skeleton) for v in chunks.values())

