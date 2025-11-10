import pytest
import numpy as np

# Assuming these enums and types are available from your module
from osteoid.formats.ostd import (
    OstdHeader,
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
        physical_path_length=123.456,
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
