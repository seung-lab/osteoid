ostd ("osteoid") File Format
================


The ostd file format exists to fill a gap in existing skeleton file formats by offering a self-contained, high performance, safe, binary format that supports optional vertex attributes for the serialization of core skeletal structures. Incorporating metadata is not a design goal, as this core file is intended to be wrapped in a container file format for that purpose.

Typically skeleton files are either text, JSON, or XML based which uses excess space and requires a comparatively slow parser. Examples include SWC, CSV/TSV, NML, OBJ. In the case of SWC, it is not easily extensible with additional attributes. Other file formats are designed to handle multiple kinds of objects. SWC also only supports trees, when some skeletons may include loops.

Precomputed is pretty close, but is inflexible on the vertex and edge data types and also requires a separate info file to interpret the binary, so is not stand-alone. Furthermore, there is no indication of which physical scale to use (nanometers? micrometers?) which is a weakness of most of the other formats too. Precomputed also only supports an edge list representation, which is space inefficient when the graph is a tree. Precomputed smartly includes a 3x4 transform matrix for affine transforms to map from, e.g. voxel to physical space, but a fully forward-compatible design would use a 4x4 matrix capable of transforms using homogeneous coordinates, perspective transforms, and is more broadly compatible with graphics pipelines, especially since the matrix remains square.

# The Design

ostd takes ideas from Precomputed, SWC, and other formats to compactly represent skeletons in a stand-alone, efficiently parsable file format.

- A 3D (XYZ) skeleton format allowing any data type for vertices
- A header for each serialized object
- Has 128 bits for an object ID, tagging each object with a UUID4 by default but accepts uint64 ids (important for connectomics)
- Incorporates a 4x4 transform matrix and tracks which state (voxel or physical) the vertices are in
- Tracks which physical unit the vertices are in (options for common metric and imperial units)
- Blocks individually guarded against file corruption by CRCs to enable extraction of remaining good data if one block is damaged
- Supports representing edges as an edge list or parent pointers in the case that the graph is a tree
- Advisory fields to tell you the number of connected components and the graph structure (currently GRAPH or TREE, if there are multiple components this would mean a forest of trees)
- Allows up to 2^16 - 1 additional vertex attributes (e.g. radius, cross sectional area, vertex type, etc)
- Attributes header is located at the end of the file to enable efficient appending on POSIX systems
- Includes a format version number to enable smooth version upgrades

## Header

All values are little endian except where noted.

| Field                  | Bytes | Datatype    | Value                       | Description                                                                                                       |
|------------------------|-------|-------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------|
| magic                  | 4     | string      | osdt                        | File magic number.                                                                                                |
| format_version         | 1     | u8          | 0                           | Version of this file format.                                                                                      |
| id                     | 16    | uuid4 / u64 | -                           | A wide enough field to accommodate both u64s and uuid4s                                                           |
| flags                  | 4     | bitfield    | VVVVECCCCGGGPPPPPSSR*       | Vertex dtype, edge representation, compression, graph type, physical unit, space. See note below for definitions. |
| num_vertices (Nv)      | 8     | u64         | -                           | Number of vertices                                                                                                |
| num_edges (Ne)         | 8     | u64         | -                           | Number of edges                                                                                                   |
| num_attributes (Nattr) | 2     | u16         | -                           | Number of vertex attributes                                                                                       |
| attr_name_width        | 1     | u8          | -                           | Fixed width size of attribute id strings.                                                                         |
| num_components         | 4     | u32         | N or (2^32-1 if unknown)    | Number of connected components in the skeleton graph. max value of uint32 is a sentinel for unknown.              |
| transform              | 64    | 4x4 f32s    | [ f32, f32, f32, f32, ... ] | Homogenous transform matrix from voxel to physical coordinates. Written in row major (C) order.                   |
| crc16                  | 2     | uint16      | -                           | crc16 using 0xFFFF init and implicit polynomial 0xd175 of header bytes excluding magic number.                    |


Flags definition:

V: Vertex data type (see Data Types)
E: Edge representation
C: Compression algorithm
G: Graph structure (advisory)
P: Physical dimension units
S: Space (voxel or physical)
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

| Type                  | Value |
|-----------------------|-------|
| PAIR                  | 0     |
| PARENT                | 1     |

PAIR means an edge list like [(e1, e2), .... ] while PARENT means parent pointers like [ 0, 0, 1 ]
for a graph like (root<-node<-leaf).

### Compression Algorithm

| Algorithm | Value |
|-----------|-------|
| None      | 0     |
| gzip      | 1     |
| bzip2     | 2     |
| zstd      | 3     |

This indicates that each buffer field is individually compressed. e.g. vertices, edges, each attribute all have this algorithm applied to them. This way you can extract each field individually without decompressing the whole file.

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


