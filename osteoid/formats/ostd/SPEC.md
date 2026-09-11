ostd ("osteoid") File Format
================

Skeletons (medial paths) may be represented as a possibly cyclic undirected graph with per vertex annotations and possibly edge annotations.

The `ostd` file format fills a gap in existing skeleton file formats by offering a self-contained, high performance, small, safe, binary format that supports optional vertex attributes for the serialization of arbitrary skeletal structures. Incorporating metadata is not a design goal, as this core file is intended to be wrapped in a container file format for that purpose. It also is not intended to replace large scale chunked representations of large datasets, though it could serve as a base serialization target for those chunks.

Typically skeleton file formats represent the geometry in either text, JSON, or XML which uses excess space and requires a comparatively slow parser. Examples include SWC, CSV/TSV, NML, OBJ. In the case of SWC, it is not easily extensible with additional attributes that can be read universally (usually user-defined columns are added). Other file formats are designed to handle multiple kinds of objects. SWC only supports trees, when some skeletons may include loops (it is possible to represent loops as a forest of related trees though). Other formats, e.g. TRK, only support paths. Precomputed only supports a fully general (and space hungry) edge list.

Precomputed has most of the features one would desire, but is inflexible on the vertex and edge data types and also requires a separate info file to interpret the binary, so is not stand-alone. Furthermore, there is no indication of which physical scale to use (nanometers? micrometers?) which is a weakness of the other formats too. Precomputed only supports an edge list representation, which is space inefficient. It smartly includes a 3x4 transform matrix for affine transforms to map from, e.g. voxel to physical space, but a fully forward-compatible design would use a 4x4 matrix capable of transforms using homogeneous coordinates, perspective transforms, and is more broadly compatible with graphics pipelines, especially since the matrix remains square.

The data type inflexibility in Precomputed causes issues because edges lists must be uint32 le, which on its face is a sensible tradeoff between space and maximum representable size, but many years later, we are finally encountering skeletons that are > 2^32 vertices (at least at certain stages of processing).

# The Design Requirements

`ostd` takes ideas from Precomputed, SWC, Trk, and other formats to compactly represent skeletons in a stand-alone, efficiently parsable file format.

- A 3D (XYZ) binary skeleton format allowing any data type for vertices
- Support skeletons larger than 2^32 vertices
- A header for each serialized object
- Includes a format version number to enable smooth version upgrades
- Has 64 bits for an object ID, important for connectomics.
- Incorporates up to 255 4x4 transform matrices for space-like dimensions and tracks which state (e.g. voxel, physical) the vertices are in
- Tracks the main physical unit of the vertices.
- Tracks which orientation the coordinate frame is in.
- Blocks individually guarded against file corruption by CRCs to enable extraction of remaining good data if one block is damaged
- Represents edges as a path graph saving space while retaining generality
- Advisory fields to tell you the number of connected components, path length, and the graph structure
- (Single-Part) Attributes header is located at the end of the file to enable efficient appending of more vertex attributes on POSIX systems
- Efficiently support both vertex and edge attributes and tracks physical units.
- Support optional spatial index (in the future)
- Concatenate multiple ostd files together to append vertices and edges together (inhibits adding more vertex attributes)

## The Design Overview

A skeleton is a graph of N-dimensional vertices consisting of a set of zero or more connected components with a numerical ID. Most skeletons are anticipated to be 3d spatial coordinates, though 4-tuple (space + time) or 5D (space + time + channel) is also possible. We support up to 8 dimensions to enable use cases that were not anticipated.

An `ostd` file may contain one or more separate skeletons with potentially different IDs. Files containing one serialized skeleton are called "single-part" files, while a file containing multiple skeletons are called "multi-part" files. Parts of a skeleton that have the same ID, should have their decoded skeletons concatenated together. This enables storing multiple skeleton and/or appending to a skeleton on disk. Most files will likely be single-part.

As image derived skeletons often consist of many adjacent points, we efficiently represent the skeleton as a 3 part mathematical object that has full graph generality. First, we decompose an existing graph into a set of non-intersecting polylines with unique vertices and maintain an edge list linking these polylines. Each polyline is written into the vertex buffer one-after-another in traversal order, meaning that within each polyline, the vertex ordering implies the edges. We then write down the length of each polyline as an integer in a list of length P, where P is the number of polylines. Lastly, we write down the explicit edge list such that it indexes into the vertex buffer. In testing on a real 3D dataset, we found this reduced the size of the edge list to about 3% of the size of the file. For more information see the analysis below.

The header of the skeleton includes basic information about how to parse it, like buffer sizes, number of vertices, compression types, and also records information for ease of high speed reading like path length and number of connected components. It also provides an advisory field for whether the skeleton is cyclic, acyclic, or unknown. This enables the automatic application of appropriate algorithms without scanning the entire skeleton first.

The header indicates what space the skeleton is in from a generic but somewhat informative list (e.g. voxel space, physical space) and indicates the presence of any transform matricies that allow you to perform an affine projection to an arbitrary linear coordinate system. The units of length are noted, as well as the coordinate frame, whether vertices are voxel centered or corner centered. This allows the easy geometric reconstruction of the skeleton from files of unknown provenance.

The vertex attributes are listed at the end of the file with a table appended at the end that describes them well. For single-part files, this allows for efficient appending of vertex attributes. Each vertex attribute is annotated with information about its dimensions using SI fundemental units. For simplicity, only SI units are supported.

## Definitions

| Term | Definition |
| :--- | :--- |
| **ostd** | The name of the file format. May refer to a valid file serialized in the `ostd` format. |
| **File** | A series of bytes that could be saved on disk representing a valid `ostd` format. This file may contain one or more `ostd` segments appended end-to-end, representing one or more skeletons. |
| **Vertex** | A coordinate in space and/or time. Most frequently an `X, Y, Z` triple of `float32` numbers, though the data type and number of dimensions are configurable. |
| **Edge** | An undirected linkage between two vertices, usually represented as a pair of integers referring to the index of each vertex in the edge. |
| **Path** | A set of connected vertices that has a maximum degree of 2. Sometimes called a *polyline* or *streamline*. |
| **Graph** | A set of vertices and connected edges. |
| **Connected Component** | A set of vertices that are all mutually reachable by following their undirected edges. |
| **Vertex Attribute** | A number associated with each vertex. For example, the radius to the nearest membrane or measured signal intensity. |
| **Skeleton** | A graph of points in space and/or time representing a stick-figure representation of a biological object or traced path. It may have multiple connected components. |
| **ID** | A 64-bit integer that uniquely identifies a skeleton within the context of a particular dataset. |
| **Single-Part File** | A file that contains only a single `ostd` segment. |
| **Multi-Part File** | A file that contains multiple `ostd` segments. Segments that have the same `ID` should be considered additive to the same skeleton. |
| **Collection** | A multi-part file containing multiple skeleton `IDs`. |

## File Structure

An ostd skeleton file can be composed of multiple parts that have an identical structure in order to allow appending vertices and edges to an existing file.

Each part should have its edges numbered such that they reference the vertices within its section. If there are vertices that appear in both sections, simply duplicate the vertex in both parts. For single-part ostd files, num_components and is precise. for multi-part files, which may contain duplicated vertices, the sum of num_components across parts is an upper bound. cable length sums remain accurate.

### Individual Part Structure

The attributes header is located at the end so that additional attributes can be easily appended to a single-part file.

| Section                  | Required | Description                                                                 |
|--------------------------|----------|-----------------------------------------------------------------------------|
| Header                   | Y        | Basic information about the file                                            |
| Transform                |          | List of 4x4 float32 le C order matrix describing voxel to different spaces, such as physical or atlas space.             |
| Vertices                 | Y        | Serialized XY pairs or XYZ triples.                                         |
| Edges                    | Y        | Edge representation.                                                        |
| Vertex & Edge Attributes |          | Serialized vertex and attributes presented in order of the following table. |
| Attributes Section       |          | Table of vertex and edge attributes                                         |

### Combined Structure

The parser should check if the stream is longer than the indicated content length. If so, it should check the subsequent section for the presence of another ostd header immediately following, and so forth. The vertices, edges, and attributes of each part then be combined after parsing. If an attribute is not present for a part, it should be filled with null values.

| Section | Description                  |
|---------|------------------------------|
| Part 1  | Initial ostd file.           |
| Part 2  | Appended vertices and edges. |
| ...     |                              |
| Part N  | Appended vertices and edges. |

## Header

All values throughout this specification are little endian except where noted. Total bytes: 80

| Field                  | Bytes | Datatype    | Value                       | Description                                                                                                       |
|------------------------|-------|-------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------|
| magic                  | 4     | string      | ostd                        | File magic number.                                                                        |
| format_version         | 1     | u8          | 0                           | Version of this file format.                                                              |
| total_bytes            | 8     | u64         | -                           | Total byte size of this part.                                                             |
| id                     | 8     | u64         | -                           |                                    |
| flags                  | 8     | bitfield    | -                           | See note below for definitions.    |
| coordinate_frame       | 4     | bitfield    | -                           | See note below for definitions.    |
| current_space          | 1     | u8          | 0                           | The current transform space the vertices are in. By default 0. Every +1 means selecting the next transform from the transform list. See *Transform* | 
| num_vertices (Nv)      | 8     | u64         | -                           | Number of vertices                                                                        |
| num_edges (Ne)         | 8     | u64         | -                           | Number of edges (explicit + implicit)                                                                           |
units specified in flags.  |
| vertex_bytes           | 8     | u64         | -                           | Content length of compressed vertex stream.                                      |
| edge_bytes             | 8     | u64         | -                           | Total byte length of the edge section, comprising the polyline-offset array, the explicit-edge pair list, and the trailing CRC-32C. |                                     |
| attribute_header_bytes | 4     | u32         | -                           | Content length in bytes of the attribute header.                                                         |
| num_components         | 4     | u32         | N or (2^32-1 if unknown)    | Number of connected components in the skeleton graph. max value of uint32 is a sentinel for unknown.              |
| cable_length           | 4     | f32         | -                           | Physical path length of this object in SI prefixed meters (See flags for SI prefix). This quantity should always be set, but if it is not set, it should be NaN.                           | 
| crc16                  | 2     | uint16      | -                           | crc16, see below. The checksum is computed over bytes 4 to 77. |

Note: parsers should reject format versions above the version they were designed for.

### CRC16 Code

This CRC16 polynomial was chosen from the CRC Zoo. The header is 80 bytes (640 bits) which exceeds
the maximum capacity of a crc8 at hamming distance 2. This crc supports detection of at least 4 flipped
bits. 

using 0xFFFF init
implicit polynomial 0xd175
refin true
refout true
xorout = 0x0000 of header bytes excluding magic number.

```python
# selected from https://users.ece.cmu.edu/~koopman/crc/index.html
# based on a header that is larger than what a crc8 can handle
def crc16(data:bytes) -> int:
  # use implicit +1 representation for right shift, LSB first
  # use explicit +1 representation for left shit, MSB first
  polynomial = 0xd175 # implicit
  crc = 0xFFFF # detects zeroed data better than 0x0000
  for i in range(len(data)):
    crc ^= data[i]
    for k in range(8):
      if crc & 1:
        crc = (crc >> 1) ^ polynomial
      else:
        crc = crc >> 1

  return int(crc & 0xFFFF)

 assert crc16("123456789".encode("utf8")) == 0x97DE
```

### Flag Definitions

The least significant bit is on the left.

`VVVVeeeeCCCCccccddddDDGGGssssstR*`

| Flag   | Meaning                            | Notes                                         |
| ------ | ---------------------------------- | --------------------------------------------- |
| **V**  | Vertex data type                   | See *Data Types*                              |
| **e**  | Edge data type                     | See *Data Types*                              |
| **C**  | Compression algorithm for vertices | See *Compression Type*                        |
| **c**  | Compression algorithm for edges    | See *Compression Type*                        |
| **G**  | Graph structure (advisory)         | See *Graph Type*.                             |
| **d**  | SI Prefix                          | signed 10^(X*3) where X is the value          |
| **D**  | Scaling                            | 0: linear 1: log10 2: log2 3: ln              |
| **s**  | Default Space Type                 | Can specify what the default space (0) means. |
| **t**  | Transforms present                 | bool                                          |
| **R**  | RESERVED                           | From this point forward                       |

The SI prefix is a restricted subset compared to the attribute version for reasons of space (and you aren't going to need a scale factor of more than 10^(3\*8) meters). The exponents in the attributes section are 8-bit for ease of parsing.

## Dimension Flag Definitions

We attempt to make the geometric interpretation of voxel positions unambiguous and compactly represented. There are a number of ambiguities in presenting a list of coordinates. Firstly, how are the axes ordered? Which are space-like and which are time-like? Are the positions referring to the centroid of the grid or to the corner closest to the origin? Much of the time, this is trivial because everyone knows the convention to follow. However, sometimes you are handed a file and need to figure out what it means. See the section Common Coordinate Frames for how there are many axis orientations in common use.

We take the cartesian convention XYZT (three space-like dimensions followed by one time-like dimension) as a guiding example. When interpreting the dimension field, first consider each axis numbered from zero in the standard cartesian manner. In the case of XYZT, it would be `[0, 1, 2, 3]`. We support up to 8 axes, which hopefully is enough for most biological work which usually uses no more than 5 axes and typically 3 or 4. We store the number of dimensions, then the number of space-like dimensions. Dimensions greater than the number of space-like dimensions are considered time-like. For XYZT, we would note there are 4 dimensions and 3 are space-like, meaning 0,1,2 conventional axes are space-like and 3 is time-like. Dimensions that are neither space-like nor time-like should be stored as vertex attributes.

Next, we need to determine the coordinate frame. There are several common frames in use, but all configurations of frames are possible and it is best to be explicit about which one is in use. Starting from our convention of XYZT arranged as a standard right-handed coordinate system with X pointing to the right, Y to the top of the page, and Z out of the page, we can negate the direction of each axis (e.g. negate Y and Z to get the standard image analysis frame of Y increases down the page and Z points into the page). 

We can also establish if the axes are permuted in a non-standard way using Lehmer codes. The Axis Permutation Type tells you how the axes are permuted with respect to the standard convention. XYZT would be encoded as 0. For example, say you had ZYX instead of XYZ, this would represented as index 5 for 3 dimensions. For ZYXT it would be 14. See Axis Permutation Type for how to encode and decode the Lehmer code.

The coordinate frame information is encoded as a uint32 le.

| Property        | Bit Position | Meaning                      | Notes                    |
|-----------------|--------------|------------------------------|--------------------------|
| Num Axes        | 0-2          | Vertex dimensions - 1        | Counts from 0 to make space since 0 is useless. |
| Space-Like      | 3-5          | Number of space-like dimensions. Following dimensions are time-like.                       | This lets you set the number of spatial dimensions. e.g. XYT, XYZT, etc. |
| Voxel Centered  | 6            | If 1, vertices are located at voxel centers. If 0, top left corner. |  |
| Reserved        | 7            | Reserved for future use. Pads out to a byte.                        |   |
| Signs           | 8-(8+num_axes-1) | Whether each axis is reflected or not. 1 is positive direction, 0 is negative direction. | Unused axes should be set to positive.    |
| Reserved        | up to bit 15 | Reserved for future use.     |                          |
| Lehmer Code     | 16-31        | See Axis Permutation Type, 000000 means +X+Y+Z standard frame. |   |      

## Transform

The default space (0) is set in the header. Transforms listed below should be written such that they are a mapping from space 0 to the selected space. Transforms can then be dynamically composed to create efficient arbitrary mappings.

### Section Structure

| Field                  | Bytes | Datatype    | Value                       | Description                                                                                                       |
|------------------------|-------|-------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------|
| num_spaces             | 1     | uint8       | -                           | Number of transformations available. matrices.                   |
| Transforms             | see below | see below | -                           | A list of 4x4 transform matrices where the space number - 1 is the index (since there is a default space defined in the header).  |
| crc32c                 | 4     | uint32      | -                          |   |

### Transform Structure

| Field                  | Bytes | Datatype    | Value                       | Description                                                                                                       |
|------------------------|-------|-------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------|
| space                  | 1     | uint8       | -                           | The kind of space the transform represents. See *Space Type* |
| units                  | 8     | tuple       | See physical units.         | The physical unit this transform maps to. |
| transform              | 64    | 4x4 f32s    | [ f32, f32, f32, f32, ... ] | Homogenous transform matrix from voxel to physical coordinates. Written in row major (C) order little endian.                   |

## Attribute Section

The attributes section is located at the end of an ostd part in order for append operations to add more attributes to single part files easily.
It can be omitted if there are no attributes.

It consists of a header followed by a listing of attribute descriptions.

num_attributes | uint8
name_width     | uint8
attribute_listing
crc16

### Attribute Listing

Attributes can be applied to either vertices or to edges. Vertex attributes are naturally dense since they scale linearly with the number of vertices, while edge attributes are naturally sparse because the number of possible combinations is the square of the number of vertices, but skeletons are paths with branches (and rarely loops). Trailing right space on a name should be truncated after decoding to utf8.

| Field          | Bytes     | Datatype | Value                     | Description                     |
|----------------|-----------|----------|---------------------------|---------------------------------|
| name           | fixed<=255| string   | e.g. "radius"             | utf8, Name of the attribute     |
| flags          | 2         | bitfield | DDDDCCCCTRRRRRRR          | Packed information              |
| unit           | 8         | tuple    | See below.                | Physical unit.                  |
| num_components | 1         | uint8    | -                         | Number of components.           |
| content_length | 8         | uint64   | -                         | Length of compressed bitstream. |


#### Flags Definition

lsb on left

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

## Geometry Representation

The skeleton graph is analyzed and deconstructed into disjoint paths that are connected via an explicit edge at branch points. The vertices of each disjoint path are written in path traversal order so that connected vertices are adjacent in the serialization and imply an edge between each other.

The graph is embedded in three buffers that are written one after another.

| Buffer           | Size                                      | Note              |
|------------------|-------------------------------------------|-------------------|
| Vertices         | verted_datatype * num_axes * num_vertices | size in header    |
| Polyline Lengths | 8 + fit(Nv) * num_polylines               |                   |
| Explicit Edges   | 2 * edge_datatype * num_edges             | size in header    |

### Vertices

Vertices are written as a sequence of serialized arrays in C order (i.e. X,Y,Z,X,Y,Z,X,Y,Z,...) where the data type is specified by the header (e.g. float32). Each polyline is written in traversal order from one of its terminal points. Polylines of length 1 are permitted.

This means that the order of the vertices matters for decoding the skeleton. However, if reconstructing the edges is not important, the vertices can be read simply as a point cloud. Vertices are guaranteed to be unique.

This buffer may be compressed using e.g. gzip, zstandard, or Draco. In the case of Draco, since order matters, you must use the preserve_order flag in the encoder.

The vertex section is followed by a crc32c that is computed from the encoded stream.

### Edges

The edge binary size covers both the polyline and explicit edges. A crc32c covers both of them as well.

#### Polyline Lengths

This section describes how to extract the polylines (paths) from the vertex buffer. It consists of an array of integers indicating the length of each polyline in vertices in order of appearence in the vertex buffer.
This array is preceeded by a uint64 le number of paths. The data type used for the path length should be the smallest integer that can contain the total number of vertices in the ostd section.

1. num_paths (u64le)
2. path lengths (fit(Nv))

e.g. 

`num_paths, len_1, len_2, ..., len_n`

#### Explicit Edges

The next section is the links between paths:

pair_1, ..., pair_n

Where the pairs are positive integers: e1,e2 that indicate which vertices should be linked. The data type is controlled by the header, though typically it will be the smallest data type that can contain the number of vertices. The edges refer to the vertices numbered from 0. Duplicate edges and self-loops are disallowed.

`num_pairs = (len(edge_binary) - paths_section - 4) / edge_dtype_size / 2`

The 4 bytes refers to the size of the crc32c checksum.

The following python pseudocode algorithms can be implemented to decode or encode linked paths.

#### Decoding Linked Paths

```python
# edge_dtype is acquired from the header
# edge_binary is calculated from header information
def decode_linked_paths(num_vertices:int, edge_dtype, edge_binary:bytes) -> np.ndarray:
	crc32c = int.from_bytes(edge_binary[-4:], 'little')
	check_crc32c(edge_binary[:-4], crc32c)

	num_paths = int.from_bytes(edge_binary[:8], 'little')
	path_dtype = smallest_dtype(num_vertices, [ np.uint8, np.uint16, np.uint32, np.uint64 ]) # smallest_dtype is a pseudocode function
	path_lengths = np.frombuffer(path_lengths_binary, offset=8, count=num_path_lengths, dtype=path_dtype)

	edges = []
	edge_i = 0;
	for i in range(num_paths):
		if path_lengths[i] == 0:
			continue
		for j in range(path_lengths[i]):
			edges.append((edge_i, edge_i+1))
			edge_i += 1

	implicit_pairs = np.array(edges, dtype=edge_dtype)

	offset = 8 + num_path_lengths * np.dtype(path_dtype).itemsize
	num_pairs = (len(edge_binary) - offset - 4) // np.dtype(edge_dtype).itemsize // 2
	explicit_pairs = np.frombuffer(edge_binary, offset=offset, count=(num_pairs*2), dtype=edge_dtype)
	explicit_pairs = explicit_pairs.reshape((num_pairs, 2), 'C')

	return np.concatenate([ implicit_pairs, explicit_pairs ])
```

#### Extracting Linked Paths

```
DEF EXTRACT(edges)

	1. Given an edge list, create an adjacency index i.e. [(5,6), (6,7)] -> { 5: [6], 6: [5,7], 7: [6] }
	2. Create a visited lookup table that can contain all edges
	3. For each edge, check if it is already visited, if it is not, call extract_branches(edge, adjacency_index)
	4. return a tuple containing the paths, explicit edges, and optionally whether a cycle was detected and 
	  the number of connected components (to get those metrics for almost free)

The idea for extract branches is to traverse each edge and greedily trace a path,
recording branches that are passed as an explicit edge. Once the path terminates,
revisit those branches and trace them in the same fashion without backtracking.

DEF EXTRACT_BRANCHES(edge, adjacency_index)

	1. Create a stack for paths, parents
	2. Create an explicit edge_list, current path, and list of already found paths
	3. insert the starting position into the current path and mark it as visited
	4. insert all its neighbors into the stack and also insert all but one of them into the explicit pairs list
	5. perform a depth first search from the node that is not in the explicit pairs list
	6. each time a new node is encountered, pick the first valid node as a successor and push the rest into
		the stack and explicit pairs list. Valid means the node is not visited and does not match the parent.

Note 1: you can integrate cycle detection to avoid calculating it later
Note 2: other algorithms can be designed that would give different path lengths 
		but the same number of explicit pairs
```

## Common Coordinate Frames

To interpret the meaning of the coordinate system and axis sign,
+X+Y+Z can be interpreted differently depending on whether the 
point of view is from the object or the observer (sometimes called the 
camera perspective).

1. The orientation provided by your right-hand thumb, middle finger, and palm when your 
hand is placed parallel to the ground pointing away from your body, palm facing up with 
the thumb flexed at a 90 degree angle. This is Right, Anterior, Superior from your point 
of view in the anatomical reference frame.

2. Patient's Left, Posterior, Superior (LPS) in the anatomical reference frame from the 
point of view of a clinician looking at a front facing patient.

These are both right handed coordinate systems.

| Axes   | Anatomical Reference Axes                 | Handedness | Notes                                                                                           |
|--------|-------------------------------------------|------------|-------------------------------------------------------------------------------------------------|
| +X+Y+Z | Patient's Left, Posterior, Superior (LPS) | Right      | Standard cartesian coordinate frame.                                                            |
| +X-Y-Z | Patient's Left, Anterior, Inferior (LAI)  | Right      | Frequently used in computer image rasterization. Images are drawn left to right, top to bottom. Used in Neuroglancer. |
| -X-Y+Z | Patient's Right, Anterior, Superior (RAS) | Right      | Frequently used in Neurology.                                                                   |
| +X-Y+Z | Patient's Left, Anterior, Superior (LAS)  | Left       | Frequently used in Radiology.                                                                   |

*Thank you to Graham Wideman (http://www.grahamwideman.com/gw/brain/orientation/orientterms.htm) for the helpful information on common coordinate systems.*

## Enums

The following tables specify the meaning of various header values.

### Data Type

| Data Type              | Value |
|------------------------|-------|
| float16                | 0     |
| float32                | 1     |
| float64                | 2     |
| uint8                  | 3     |
| uint16                 | 4     |
| uint32                 | 5     |
| uint64                 | 6     |
| int8                   | 7     |
| int16                  | 8     |
| int32                  | 9    |
| int64                  | 10    |
| boolean (1 byte)       | 11    |
| packed boolean (1 bit) | 12    |

### Compression Algorithm Type

| Algorithm | Value |
|-----------|-------|
| None      | 0     |
| gzip      | 1     |
| bzip2     | 2     |
| zstd      | 3     |
| draco     | 4     |
| xz        | 5     |

Note: For Draco, preserve order must be used to preserve edge relationships to vertices. This unfortunately reduces the amount of compression possible.


### Unit Definitions

Arbitrary SI dimensions can be encoded as a tuple of fundemental units raised to a signed exponent. While
the original vision of this datastream is to only incorporate metadata that supports the interpretation of geometry, it seems strange to privilage length and time dimensions when all sorts of things might be measured along a skeleton, such as current or luminance.

For example, Joules can be expressed as W = ma x d or kg * m^2/s^2, Watts as kg * m^2/s^3 and Coulombs can be expressed as A * s. Speed is meters/sec, area is meters^2 etc, luminousity can be measured in watts or photons per a second.

Since this is designed for biological use cases, the candela which is based in human perception of light is less useful, so we reserve those bits for future use (e.g. one can imagine using them for signaling the use of US customary units).

Therefore, for attributes, we encode the dimensions as a uint64 that represents the following structure:

| Field           | Data Type        | Description                           |
|-----------------|------------------|-------------------------------------- |
| SI Prefix       | int8             | signed 10^(X*3) where X is the value  |
| amperes         | int8             | A^x                                   |
| kelvin          | int8             | K^x                                   |
| kilograms       | int8             | kg^x                                  |
| meters          | int8             | m^x                                   |
| mols            | int8             | mol^x                                 |
| seconds         | int8             | s^x                                   |
| Scaling         | uint8            | 0: linear, 1: log10, 2: log2, 3: ln   |

This allows you to specify an arbitrary SI derived unit as follows. The SI prefix is applied to the linear scaled figure and then the logarithm is applied if applicable.

All signed values are written in 2's complement. The field components are stored little endian as is the field as a whole.

For example, let's demonstrate km/s^2 as an acceleration value.

| Field           | Value            | Description                           |
|-----------------|------------------|-------------------------------------- |
| SI Prefix       | 1                | signed 10^(X*3) where X is the value  |
| amperes         | 0                | A^x                                   |
| kelvin          | 0                | K^x                                   |
| kilograms       | 0                | kg^x                                  |
| meters          | 1                | m^x                                   |
| mols            | 0                | mol^x                                 |
| seconds         | -2               | s^x                                   |
| Scaling         | 0                | 00: linear, 01: log10, 10: log2, 11: ln |

This means SI Prefix kilo (10^(1 * 3) = 1000), seconds -2 = 1/s^2, meters 1 = m. 0 scaling means the value should be read as linear.

kilo(m * 1/s^2) = km/s^2, an acceleration value

### Graph Type

3 bit

| Graph Structure       | Value | Description                |
|-----------------------|-------|----------------------------|
| graph                 | 0     | Could be any kind of graph.|
| tree                  | 1     | Acyclic graph.             |
| cyclic                | 2     | Contains one or more loops.|

This value is advisory and does not control the edge representation.
An extra bit is reserved in case additional useful categories are identified.

### Space Type

5 bits

| Space                        | Value | Description                                                       |
|------------------------------|-------|-------------------------------------------------------------------|
| Unknown                      | 0     | No semantic interpretation of the coordinate system was recorded. |
| Voxel                        | 1     | Vertices are specified in voxels.                                 |
| Physical                     | 2     | Vertices are specified in physical units like nanometers.         |
| Scanner                      | 3     | Relative to the scanner origin.                                   |
| Atlas                        | 4     | Relative to a reference atlas.                                    |
| World                        | 5     | With reference to a specific external coordinate system not described by another category. |
| Anchor                       | 6     | Vertices are relative to an anchor/landmark like a soma.          |
| Reserved                     | 7-23  | These ids are reserved for future use.
| User Defined                 | 24-31 | Users may define their own application specific spaces.           |

### Axis Permutation Type

This should work with arbitrary numbers of axes. It is the index (starting from zero) of the combination generated according to the Lehmer code, the algorithm is reproduced below. Each axis is numbered from 0 to N-1. For example, to write XYZ, you would input [0,1,2], XZY would be [0,2,1]. XYZW would be [0,1,2,3] and so on. To find the index of a given permutation, call `rank_permutation([0,2,1])`. To find the permutation corresponding to an index `k`, e.g. for three axes, call `enumerate_permutations(3,k)` and locate the index in the array containing the desired permutation.

This procedure works well with the current design limit of 8 axes which fits inside a `uint16`, but can work up to 20 axes and fit inside a `uint64`.

```python
import math

def unrank_permutation(num_axes:int, k:int) -> list[int]:
	"""Convert an index into an axis permutation for a given number of axes."""
    axes = list(range(num_axes))
    perm = []

    for i in range(num_axes, 0, -1):
        f = math.factorial(i - 1)
        idx = k // f
        k %= f
        perm.append(axes.pop(idx))

    return perm

def rank_permutation(perm:list[int]) -> int:
	"""Convert an axis permutation to an index."""
    axes = list(range(len(perm)))
    n = len(axes)
    rank = 0

    for i in range(n):
        idx = axes.index(perm[i])
        rank += idx * math.factorial(n - i - 1)
        axes.pop(idx)

    return rank
```

### Path Graph Efficiency Analysis

To give a quick analysis, in a 3D dataset using float32 vertices and uint64 edges, a vertex is 4x3 (12) bytes. An edge is 8x2 (16) bytes. In the naive approach, a polyline with Nv vertices (holding Nv > 2), would have 12xNv bytes in its vertex buffer, and 16x(Nv-1) bytes in its edge buffer. By contrast the new approach would have 12 x Nv bytes and 9 additional bytes to indicate the size of the offsets (size of buffer, length). 

Assume a polyline of Nv = 1000 vertices with no branches, vertex datatype Dv = 4 bytes, and edge datatype De = 8 bytes.

```
Naive Approach = 3 Nv Dv + 2 De (Nv - 1)
Path Graph = 3 Nv Dv + 9

Path Graph / Naive = 12009 bytes / 27984 bytes = 43%
```

In this simple example, we have produced a binary 43% the size of the original. With numerous polylines and branches, this advantage shrinks slightly. If every polyine consists of a single voxel, it becomes obvious the naive approach becomes better due to the size of the polyline offset buffer. We can calculate the crossover point at which this representation becomes more expensive. Let P again be the number of polylines and note P <= Nv. We can take our single polyline example, and arbitrarily break it into up to Nv polylines, which adds both an offset entry and an explicit edge per a polyline.

```
Path Graph = 3NvDv + (8+2*8*P) + 2 De P

Let path graph = naive and solve for P.

P = (2 De (Nv - 1) - 8) / (16 + 2 De)

Plugging in our numbers to find this crossover point.

P = 499.25
```

That's just under Nv/2. This analysis is very conservative, because if you take it too literally, you might think, well then if my data look like a binary tree, then I should go elsewhere. However, if you look at a diagram of a full binary tree, there are many nodes connected in chains! So you can elide explicitly representing edges to depth (d-1) (left) + (d-2) (right) at the first level, and so on as you progress down the tree. I think you need a pretty high depth before this becomes reasonable, but there are some possible savings here. At d=4, there are 14 edges in the graph, and 7 that can be represented implicitly (though it is offset by the polyline buffer). (Note: This is a empirically testable property.)