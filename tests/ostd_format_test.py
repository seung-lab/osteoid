import pytest
import numpy as np

# Assuming these enums and types are available from your module
from osteoid.formats.ostd import (
    OstdHeader,
    OstdSkeleton,
    OstdTransformSection,
    OstdTransform,
    DataType,
    CompressionType,
    EdgeRepresentationType,
    GraphType,
    SpaceType,
)

@pytest.fixture
def sample_header():
    return OstdHeader(
        Nv=123,
        Ne=456,
        append_mode=True,
        attribute_header_bytes=64,
        cable_length=123.456,
        coordinate_frame_orientation='-X+Y+Z',
        crc16=60279,
        edge_data_type=DataType.U16,
        edge_compression=CompressionType.ZSTD,
        edge_representation=EdgeRepresentationType.PAIR,
        edge_bytes=1024,
        format_version=OstdHeader.FORMAT_VERSION,
        graph_type=GraphType.TREE,
        has_transform=False,
        id=42,
        length_unit="um",
        num_axes=3,
        num_components=12345,
        space=SpaceType.WORLD,
        spatial_index_bytes=2048,
        total_bytes=OstdHeader.HEADER_BYTES,
        vertex_compression=CompressionType.GZIP,
        vertex_data_type=DataType.F32,
        vertex_bytes=512,
        voxel_centered=False,
    )

def test_header_roundtrip(sample_header):
    binary = sample_header.to_bytes()
    assert isinstance(binary, (bytes, bytearray))
    assert len(binary) == OstdHeader.HEADER_BYTES

    reconstructed = OstdHeader.from_bytes(binary)

    # Check that all fields match exactly
    for key, value in sample_header.__dict__.items():
        other = getattr(reconstructed, key)
        if isinstance(value, float) and np.isnan(value):
            assert np.isnan(other)
        elif isinstance(value, float):
        	assert np.isclose(value, other)
        else:
            assert other == value, f"Mismatch in field '{key}': {value} != {other}"


@pytest.fixture
def sample_transform():
    mat = np.arange(16, dtype=np.float32).reshape(4, 4)
    return OstdTransform(space=SpaceType.WORLD, transform=mat)

def test_transform_roundtrip(sample_transform):
    data = sample_transform.to_bytes()
    assert isinstance(data, (bytes, bytearray))
    assert len(data) == OstdTransform.NUM_BYTES

    recon = OstdTransform.from_bytes(data)
    assert recon.space == sample_transform.space
    np.testing.assert_array_equal(recon.transform, sample_transform.transform)

def test_section_roundtrip(sample_transform):
    section = OstdTransformSection(spaces=[sample_transform, sample_transform])
    data = section.to_bytes()
    assert isinstance(data, (bytes, bytearray))
    assert data[0] == 2  # number of spaces
    assert len(data) == 1 + 2 * OstdTransform.NUM_BYTES + 2  # header + payloads + crc16

    recon = OstdTransformSection.from_bytes(data)
    assert len(recon.spaces) == 2
    for a, b in zip(recon.spaces, section.spaces):
        assert a.space == b.space
        np.testing.assert_array_equal(a.transform, b.transform)

def test_crc_mismatch_detection(sample_transform):
    section = OstdTransformSection(spaces=[sample_transform])
    data = bytearray(section.to_bytes())
    data[-1] ^= 0xFF  # corrupt CRC
    with pytest.raises(ValueError, match="Transform header corruption detected"):
        OstdTransformSection.from_bytes(bytes(data))

@pytest.fixture
def sample_skeleton():
    vertices = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
    ], dtype=np.float32)

    edges = np.array([
        [0, 1],
        [1, 2],
    ], dtype=np.uint64)

    return OstdSkeleton.create(
        vertices, edges, 
        id=42,
        coordinate_frame_orientation="+X+Y+Z",
        voxel_centered=True
    )

def test_create_properties(sample_skeleton):
    skel = sample_skeleton
    assert skel.vertices.shape == (3, 3)
    assert skel.vertices.dtype == np.float32
    assert skel.edges.shape == (2, 2)
    assert skel.edges.dtype == np.uint64
    assert skel.id == 42
    assert skel.coordinate_frame_orientation == "+X+Y+Z"
    assert skel.voxel_centered is True

def test_serialization_roundtrip(sample_skeleton):
    skel = sample_skeleton
    data = skel.to_bytes()
    restored = OstdSkeleton.from_bytes(data)

    np.testing.assert_array_equal(restored.vertices, skel.vertices)
    np.testing.assert_array_equal(restored.edges, skel.edges)
    assert restored.id == skel.id
    assert restored.coordinate_frame_orientation == skel.coordinate_frame_orientation
    assert restored.voxel_centered == skel.voxel_centered

def test_serialization_determinism(sample_skeleton):
    data1 = sample_skeleton.to_bytes()
    data2 = sample_skeleton.to_bytes()
    assert data1 == data2

def test_from_bytes_invalid():
    with pytest.raises(Exception):
        OstdSkeleton.from_bytes(b"invalid_data")