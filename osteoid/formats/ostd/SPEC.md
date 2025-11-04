ostd ("osteoid") File Format
================

Skeletons (medial paths) may be represented as a possibly cyclic undirected graph with per vertex annotations and possibly edge annotations.

The `ostd` file format fills a gap in existing skeleton file formats by offering a self-contained, high performance, safe, binary format that supports optional vertex attributes for the serialization of skeletal structures. Incorporating metadata is not a design goal, as this core file is intended to be wrapped in a container file format for that purpose.

Typically skeleton file formats represent the geometry in either text, JSON, or XML which uses excess space and requires a comparatively slow parser. Examples include SWC, CSV/TSV, NML, OBJ. In the case of SWC, it is not easily extensible with additional attributes. Other file formats are designed to handle multiple kinds of objects. SWC only supports trees, when some skeletons may include loops.

Precomputed has most of the features one would desire, but is inflexible on the vertex and edge data types and also requires a separate info file to interpret the binary, so is not stand-alone. Furthermore, there is no indication of which physical scale to use (nanometers? micrometers?) which is a weakness of most of the other formats too. Precomputed also only supports an edge list representation, which is space inefficient. Precomputed smartly includes a 3x4 transform matrix for affine transforms to map from, e.g. voxel to physical space, but a fully forward-compatible design would use a 4x4 matrix capable of transforms using homogeneous coordinates, perspective transforms, and is more broadly compatible with graphics pipelines, especially since the matrix remains square.

Another design issue in Precomputed is that it specifies edges lists must be uint32 le, which on its face is a sensible tradeoff between space and maximum representable size, but many years later, we are finally encountering skeletons that are > 2^32 vertices (at least at certain stages of processing).

# The Design

ostd takes ideas from Precomputed, SWC, and other formats to compactly represent skeletons in a stand-alone, efficiently parsable file format.

- A 3D (XYZ) binary skeleton format allowing any data type for vertices
- Support skeletons larger than 2^32 vertices
- A header for each serialized object
- Includes a format version number to enable smooth version upgrades
- Has 128 bits for an object ID, tagging each object with a UUID4 by default but accepts uint64 ids (important for connectomics)
- Incorporates a 4x4 transform matrix and tracks which state (voxel or physical) the vertices are in
- Tracks which physical unit the vertices are in (options for common metric and imperial units)
- Blocks individually guarded against file corruption by CRCs to enable extraction of remaining good data if one block is damaged
- Supports representing edges as an edge list, parent pointers, or a path graph, saving space while retaining generality
- Advisory fields to tell you the number of connected components and the graph structure (currently GRAPH or TREE, if there are multiple components this would mean a forest of trees)
- Allows up to 2^16 - 1 additional vertex attributes (e.g. radius, cross sectional area, vertex type, etc)
- (Single-Part) Attributes header is located at the end of the file to enable efficient appending of more vertex attributes on POSIX systems
- Efficiently support both vertex and edge attributes
- Support optional octree spatial index
- Concatenate multiple ostd files together to append vertices and edges together (inhibits adding more vertex attributes)

## File Structure

An ostd skeleton file can be composed of multiple parts that have an identical structure in order to allow appending vertices and edges to an existing file.
The key difference with an appended section is that its edge list is numbered such that it includes vertices in the preceeding sections.

### Individual Part Structure

The attributes header is located at the end so that additional attributes can be easily appended to a single-part file.

| Section                  | Required | Description                                                                 |
|--------------------------|----------|-----------------------------------------------------------------------------|
| Header                   | Y        | Basic information about the file                                            |
| Transform                |          | 4x4 float32 le C order matrix describing voxel->physical space.             |
| Spatial Index            |          | Octree describing locations of vertices.                                    |
| Vertices                 | Y        | Serialized XY pairs or XYZ triples.                                         |
| Edges                    | Y        | Edge representation.                                                        |
| Vertex & Edge Attributes |          | Serialized vertex and attributes presented in order of the following table. |
| Attributes header        |          | Table of vertex and edge attributes                                         |

### Combined Structure

The parser should check if the stream is longer than the indicated content length. If so, it should check the subsequent section for the presence of another ostd file, and so forth. The vertices, edges, and attributes of each part then be combined after parsing. If an attribute is not present for a part, it should be filled with null values.

| Section | Description                  |
|---------|------------------------------|
| Part 1  | Initial ostd file.           |
| Part 2  | Appended vertices and edges. |
| ...     |                              |
| Part N  | Appended vertices and edges. |

## Header

All values are little endian except where noted. Total bytes: 137

| Field                  | Bytes | Datatype    | Value                       | Description                                                                                                       |
|------------------------|-------|-------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------|
| magic                  | 4     | string      | osdt                        | File magic number.                                                                        |
| format_version         | 1     | u8          | 0                           | Version of this file format.                                                              |
| total_bytes            | 8     | u64         | -                           | Total byte size of this part.                                                             |
| id                     | 16    | uuid4 / u64 | -                           | A wide enough field to accommodate both u64s and uuid4s                                   |
| object_flags           | 2     | bitfield    | VVVVeeeeCCCCcccc            | Vertex dtype, edge datatype, compression, see note below for definitions.                 |
| flags                  | 4     | bitfield    | EEAGGGPPPPPSSRiat*          | edge representation, graph type, physical unit, space. See note below for definitions.   |
| num_vertices (Nv)      | 8     | u64         | -                           | Number of vertices                                                                        |
| num_edges (Ne)         | 8     | u64         | -                           | Number of edges                                                                           |
| vertex_bytes           | 8     | u64         | -                           | Content length of vertices (needed for compression).                                      |
| edge_bytes             | 8     | u64         | -                           | Number of bytes encoding edges (needed for compression)                                   |
| attribute_header_bytes | 4     | u32         | -                           | Content length in bytes of the attribute header.                                                         |
| num_components         | 4     | u32         | N or (2^32-1 if unknown)    | Number of connected components in the skeleton graph. max value of uint32 is a sentinel for unknown.              |
| transform              | 64    | 4x4 f32s    | [ f32, f32, f32, f32, ... ] | Homogenous transform matrix from voxel to physical coordinates. Written in row major (C) order.                   |
| crc16                  | 2     | uint16      | -                           | crc16 using 0xFFFF init and implicit polynomial 0xd175 of header bytes excluding magic number.                    |


Flags definition:

V: Vertex data type (see Data Types)
e: Edge data type (see Data Types)
E: Edge representation
A: Append mode
	(0) This file section is completely self contained.
	(1) edges are numbered such that they include preceeding section parts.
C: Compression algorithm for vertices
c: Compression algorithm for edges
G: Graph structure (advisory)
P: Physical dimension units
S: Space (voxel or physical)
i: Spatial index present
a: axes (0) XY (1) XYZ
t: transform present
R*: RESERVED from this point forward

### Data Types

| Data Type              | Value |
|------------------------|-------|
| float8                 | 0     |
| float16                | 1     |
| float32                | 2     |
| float64                | 3     |
| uint8                  | 4     |
| uint16                 | 5     |
| uint32                 | 6     |
| uint64                 | 7     |
| int8                   | 8     |
| int16                  | 9     |
| int32                  | 10    |
| int64                  | 11    |
| boolean (1 byte)       | 12    |
| packed boolean (1 bit) | 13    |


### Edge Representation

| Type                  | Value | Description                                         | Properties                                  |
|-----------------------|-------|-----------------------------------------------------|---------------------------------------------|
| PAIR                  | 0     | edge list [(e1,e2), ... ]                           | General, no parsing. Most space.            |
| PARENT                | 1     | parent pointers [0,0,1] means (root<-node<-leaf)    | Trees only. Half space. Structure encoded.  |
| LINKED_PATHS          | 2     | set of paths linked at branch points                | General, efficient if branches are sparse.  |

Note: Only PAIR is currently supported.

### Compression Algorithm

| Algorithm | Value |
|-----------|-------|
| None      | 0     |
| gzip      | 1     |
| bzip2     | 2     |
| zstd      | 3     |
| draco     | 4     |

This indicates that each buffer field is individually compressed. e.g. vertices, edges, each attribute all have this algorithm applied to them. This way you can extract each field individually without decompressing the whole file.

Note: Only None is currently supported.

### Physical Dimension Values

| Data Type             | Value |
|-----------------------|-------|
| VOXEL (dimensionless) | 0     |
| Angstrom              | 1     |
| Femtometer            | 2     |
| Picometer             | 3     |
| Nanometer             | 4     |
| Micrometer (micron)   | 5     |
| Millimeter            | 6     |
| Centimeter            | 7     |
| Meter                 | 8     |
| Kilometer             | 9     |
| Megameter             | 10    |
| Lightyear             | 11    |
| Parsec                | 12    |
| Mil (1/1000 inch)     | 13    |
| Inch                  | 14    |
| Foot                  | 15    |
| Yard                  | 16    |
| Statute Mile          | 17    |
| Nautical Mile         | 18    |


### Graph Type

| Graph Structure       | Value |
|-----------------------|-------|
| GRAPH                 | 0     |
| TREE                  | 1     |

This value is advisory and does not control the edge representation. 
This is because a tree can be represented as an edge list. See edge representation.

### Space Type

| Space                 | Value |
|-----------------------|-------|
| Voxel                 | 0     |
| Physical              | 1     |

## Spatial Index

This optional section provides an optimized way to fetch vertices when the edge representation is LINKED_PATHS. During path generation,
the chunk size will be used to break long paths at chunk boundaries, creating more numerous paths connected by new edge pairs. This will
make it simple to fetch paths in a given region of space.

| Name       | Type        | Description                                            |
|------------|-------------|--------------------------------------------------------|
| num_bytes  | uint64      | Number of bytes in this section                        |
| minpt      | 3 x float32 | Minimum point of the bounding box.                     |
| maxpt      | 3 x float32 | Maximum point of the bounding box.                     |
| chunk_size | 3 x float32 | Chunked division of this space.                        |
| Octree     |             | Division of space terminating with a list of path ids. |


## Attribute Header

Attributes can be applied to either vertices or to edges. Vertex attributes are naturally dense since they scale linearly with the number of vertices, while edge attributes are naturally sparse because the number of possible combinations is the square of the number of vertices, but skeletons are paths with branches (and rarely loops).

| Field          | Bytes     | Datatype | Value            | Description                     |
|----------------|-----------|----------|------------------|---------------------------------|
| name           | variable? | string   | e.g. "radius"    | Name of the attribute           |
| field          | 2         | bitfield | TDDDDCCCCRRRRRRR | Packed information              |
| num_components | 1         | uint8    | -                | Number of components.           |
| content_length | 8         | uint64   | -                | Length of compressed bitstream. |

Bitfield (lsb on left):
T: (0) vertex attribute (1) edge attribute
D: Data Type
C: Compression Type
R: Reserved

## Vertex Attribute

Vertex Attributes are stored as a serialized array in C order little endian, optionally wrapped in compression.
Example vertex attributes: Vertex types (uint8), radius (float32), cross_sectional_area (float32), cross_sectional_area_contacts (uint8)

## Edge Attributes

N | (e1,e2) ... | attr12 ...

N is uint64
edge widths are based on header

## Vertex Representation

Can be XY or XYZ depending on header. Datatype set depending on header. Encoded as serialized array in C order (i.e. XYZ,XYZ,XYZ).

When the edge representation is LINKED_PATHS, the vertices will be sorted based on their connected neighbors.

## Edge Representation

### Pairs

The edges will be written explicitly as a serialized array of pairs in C order (i.e. (e1,e2),(e1,e3),...) in little endian with the data type controlled by the header.

### Parent

The edges will be written as a serialized array of parent pointers aligned with vertices. The data type is controlled by the header. This only applies to trees.

### Linked Paths

The skeleton is analyzed and deconstructed into disjoint paths that are connected at branch points. The vertices will be sorted so that each disjoint path is contiguously represented so that connected vertices are adjacent in the serialization.

The first section of the edges is then written as:

num_paths | len_1, len_2, ..., len_n

Where the lengths are the number of vertices in each path. Each path is assigned an ID numbering from 0.

The next section is the links between paths:

num_pairs | pair_1, ..., pair_n

Where the pairs are: e1,e2 with the data type controlled by the header, though typically it will be the smallest data type that encodes vertices. The edges refer to the vertices, not to the path ID.













