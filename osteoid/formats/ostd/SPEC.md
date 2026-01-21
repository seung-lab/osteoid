ostd ("osteoid") File Format
================

Skeletons (medial paths) may be represented as a possibly cyclic undirected graph with per vertex annotations and possibly edge annotations.

The `ostd` file format fills a gap in existing skeleton file formats by offering a self-contained, high performance, small, safe, binary format that supports optional vertex attributes for the serialization of arbitrary skeletal structures. Incorporating metadata is not a design goal, as this core file is intended to be wrapped in a container file format for that purpose.

Typically skeleton file formats represent the geometry in either text, JSON, or XML which uses excess space and requires a comparatively slow parser. Examples include SWC, CSV/TSV, NML, OBJ. In the case of SWC, it is not easily extensible with additional attributes. Other file formats are designed to handle multiple kinds of objects. SWC only supports trees, when some skeletons may include loops. Other formats, e.g. TRK, only support paths. Precomputed only supports a fully general (and space hungry) edge list.

Precomputed has most of the features one would desire, but is inflexible on the vertex and edge data types and also requires a separate info file to interpret the binary, so is not stand-alone. Furthermore, there is no indication of which physical scale to use (nanometers? micrometers?) which is a weakness of most of the other formats too. Precomputed also only supports an edge list representation, which is space inefficient. Precomputed smartly includes a 3x4 transform matrix for affine transforms to map from, e.g. voxel to physical space, but a fully forward-compatible design would use a 4x4 matrix capable of transforms using homogeneous coordinates, perspective transforms, and is more broadly compatible with graphics pipelines, especially since the matrix remains square.

Another design issue in Precomputed is that it specifies edges lists must be uint32 le, which on its face is a sensible tradeoff between space and maximum representable size, but many years later, we are finally encountering skeletons that are > 2^32 vertices (at least at certain stages of processing).

# The Design

ostd takes ideas from Precomputed, SWC, Trk, and other formats to compactly represent skeletons in a stand-alone, efficiently parsable file format.

- A 3D (XYZ) binary skeleton format allowing any data type for vertices
- Support skeletons larger than 2^32 vertices
- A header for each serialized object
- Includes a format version number to enable smooth version upgrades
- Has 128 bits for an object ID, tagging each object with a UUID4 by default but accepts uint64 ids (important for connectomics)
- Incorporates up to 255 4x4 transform matrices and tracks which state (e.g. voxel, physical) the vertices are in
- Tracks the main physical unit of the vertices.
- Tracks which orientation the coordinate frame is in.
- Blocks individually guarded against file corruption by CRCs to enable extraction of remaining good data if one block is damaged
- Supports representing edges as an edge list, parent pointers, or a path graph, saving space while retaining generality
- Advisory fields to tell you the number of connected components, path length, and the graph structure
- (Single-Part) Attributes header is located at the end of the file to enable efficient appending of more vertex attributes on POSIX systems
- Efficiently support both vertex and edge attributes and tracks physical units.
- Support optional spatial index
- Concatenate multiple ostd files together to append vertices and edges together (inhibits adding more vertex attributes)

## File Structure

An ostd skeleton file can be composed of multiple parts that have an identical structure in order to allow appending vertices and edges to an existing file.
The key difference with an appended section is that its edge list is numbered such that it includes vertices in the preceeding sections.

### Individual Part Structure

The attributes header is located at the end so that additional attributes can be easily appended to a single-part file.

| Section                  | Required | Description                                                                 |
|--------------------------|----------|-----------------------------------------------------------------------------|
| Header                   | Y        | Basic information about the file                                            |
| Transform                |          | List of 4x4 float32 le C order matrix describing voxel to different spaces, such as physical or atlas space.             |
| Spatial Index            |          | Octree describing locations of vertices.                                    |
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
| magic                  | 4     | string      | osdt                        | File magic number.                                                                        |
| format_version         | 1     | u8          | 0                           | Version of this file format.                                                              |
| total_bytes            | 8     | u64         | -                           | Total byte size of this part.                                                             |
| id                     | 8     | u64         | -                           |                                    |
| flags                  | 8     | bitfield    | -                           | See note below for definitions.    |
| dimension_flags        | 4     | bitfield    | -                           | See note below for definitions.    |
| num_vertices (Nv)      | 8     | u64         | -                           | Number of vertices                                                                        |
| num_edges (Ne)         | 8     | u64         | -                           | Number of edges                                                                           |
units specified in flags.  |
| vertex_bytes           | 8     | u64         | -                           | Content length of vertices (needed for compression).                                      |
| edge_bytes             | 8     | u64         | -                           | Number of bytes encoding edges (needed for compression)                                   |
| attribute_header_bytes | 4     | u32         | -                           | Content length in bytes of the attribute header.                                                         |
| num_components         | 4     | u32         | N or (2^32-1 if unknown)    | Number of connected components in the skeleton graph. max value of uint32 is a sentinel for unknown.              |
| cable_length           | 4     | f32         | -                           | Physical path length of this object in SI prefixed meters (See flags for SI prefix). This quantity should always be set, but if it is not set, it should be NaN.                           | 
| current_space          | 1     | u8          | 0                           | The current transform space the vertices are in. By default 0. Every +1 means selecting the next transform from the transform list. See *Transform* | 
| crc16                  | 2     | uint16      | -                           | crc16 using 0xFFFF init and implicit polynomial 0xd175 of header bytes excluding magic number.                    |

### Flag Definitions

The least significant bit is on the left.

`VVVVeeeeCCCCccccppppGGGEEssssstiR*`

| Flag   | Meaning                            | Notes                                         |
| ------ | ---------------------------------- | --------------------------------------------- |
| **V**  | Vertex data type                   | See *Data Types*                              |
| **e**  | Edge data type                     | See *Data Types*                              |
| **C**  | Compression algorithm for vertices | See *Compression Type*                        |
| **c**  | Compression algorithm for edges    | See *Compression Type*                        |
| **G**  | Graph structure (advisory)         | See *Graph Type*.                             |
| **p**  | SI Prefix                          | signed 10^(X*3) where X is the value          |
| **E**  | Edge representation                | See *Edge Representation*                     |
| **s**  | Default Space Type                 | Can specify what the default space (0) means. |
| **t**  | Transforms present                 | bool.                                         |
| **i**  | Spatial index present.             |                                               | 
| **R**  | RESERVED                           | From this point forward                       |

## Dimension Flag Definitions

We attempt to make the geometric interpretation of voxel positions unambiguous and compactly represented. There are a number of ambiguities in presenting a list of coordinates. Firstly, how are the axes ordered? Which are space-like and which are time-like? Are the positions referring to the centroid of the grid or to the corner closest to the origin? Much of the time, this is trivial because everyone knows the convention to follow. However, sometimes you are handed a file and need to figure out what it means. See the section Common Coordinate Frames for how there are many axis orientations in common use.

We take the cartesian convention XYZT (three space-like dimensions followed by one time-like dimension) as a guiding example. When interpreting the dimension field, first consider each axis numbered from zero in the standard cartesian manner. In the case of XYZT, it would be `[0, 1, 2, 3]`. We support up to 8 axes, which hopefully is enough for most biological work which usually uses no more than 5 axes and typically 3 or 4. We store the number of dimensions, then the number of space-like dimensions. Dimensions greater than the number of space-like dimensions are considered time-like. For XYZT, we would note there are 4 dimensions and 3 are space-like, meaning 0,1,2 conventional axes are space-like and 3 is time-like. Dimensions that are neither space-like nor time-like should be stored as vertex attributes.

Next, we need to determine the coordinate frame. There are several common frames in use, but all configurations of frames are possible and it is best to be explicit about which one is in use. Starting from our convention of XYZT arranged as a standard right-handed coordinate system with X pointing to the right, Y to the top of the page, and Z out of the page, we can negate the direction of each axis (e.g. negate Y and Z to get the standard image analysis frame of Y increases down the page and Z points into the page). 

We can also establish if the axes are permuted in a non-standard way using Lehmer codes. The Axis Permutation Type tells you how the axes are permuted with respect to the standard convention. XYZT would be encoded as 0. For example, say you had ZYX instead of XYZ, this would represented as index 5 for 3 dimensions. For ZYXT it would be 14. See Axis Permutation Type for how to encode and decode the Lehmer code.

The least significant bit is on the left.

`aaalllssssssssoooooooooooooooocR`

| Flag   | Meaning                            | Notes                                                                                |
| ------ | ---------------------------------- | ------------------------------------------------------------------------------------ |
| **a**  | Number of Axes                     | Number of axes                       |
| **l**  | Number of space-like axes.         | First l axes are space-like, the following are time-like.                       |
| **s**  | Bitfield. Sign of each axis direction compared to convention.     |  sign of X,Y,Z axes in that order (0: positive, 1: negative). Unused axes should be set to positive.                   |
| **o**  | Coordinate Frame Orientation       | Has own structure: `8s 16a`<br>s: <br>a: axis permutation<br>See Axis Permutation Type, 000000 means +X+Y+Z standard frame. By default, the signs should be positive. |
| **c** | Voxel centered or top left corner. | Describes whether voxel coordinates are interpreted as centered or in the corner closest to the origin.                                                              |
| **R*** | RESERVED                           | From this point forward               

## Transform

The default space (0) is set in the header. Transforms listed below should be written such that they are a mapping from space 0 to the selected space. Transforms can then be dynamically composed to create efficient arbitrary mappings.

| Field                  | Bytes | Datatype    | Value                       | Description                                                                                                       |
|------------------------|-------|-------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------|
| num_spaces             | 1     | uint8       | -                           | Number of transformations available. matrices.                  |
| units                  | 4     | tuple       | See physical units.         | The physical unit this transform maps to. |
| space                  | 1     | uint8       | -                           | The kind of space the transform represents. See *Space Type* |
| transform              | 64    | 4x4 f32s    | [ f32, f32, f32, f32, ... ] | Homogenous transform matrix from voxel to physical coordinates. Written in row major (C) order little endian.                   |
| crc16                  | 2     | uint16      | -                          | 16-bit CRC using 0xFF init and 0xd175 implicit polynomial                  |


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
| index      | list[(32,u32)] |  |
| paths      |             | Array of path ids. |
| crc32c     | uint32      | Checksum for spatial index.                            |


Since generally speaking, a skeleton may be elongated and curved, a dense representation of the spatial index will frequently be wasteful as many grid spaces will be empty. The grid will be numbered in fortran order (`x + cx * (y + cy * z)`) where x,y,z are grid points and cx,cy,cz are chunk sizes with 0 being the gridpoint that touches the minpt.

The format of the index is an array of `gridpoint,offset` pairs, both being uint32s sorted in eytzinger order. Offset is defined as the byte offset into the spatial index section. This section is followed by the paths section which contains lists of branch ids.


## Attribute Section

The attributes section is located at the end of an ostd part in order for append operations to add more attributes to single part files easily.
It can be omitted if there are no attributes.

It consists of a header followed by a listing of attribute descriptions.

num_attributes | uint8
name_width     | uint8
attribute_listing
crc16

### Attribute Listing

Attributes can be applied to either vertices or to edges. Vertex attributes are naturally dense since they scale linearly with the number of vertices, while edge attributes are naturally sparse because the number of possible combinations is the square of the number of vertices, but skeletons are paths with branches (and rarely loops).

| Field          | Bytes     | Datatype | Value                     | Description                     |
|----------------|-----------|----------|---------------------------|---------------------------------|
| name           | fixed<=255| string   | e.g. "radius"             | Name of the attribute           |
| flags          | 2         | bitfield | DDDDCCCCTRRRRRRR          | Packed information              |
| unit           | 4         | tuple    | See below.                | Physical unit.                  |
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

## Vertex Representation

Number of axes determined from header. Datatype set depending on header. Encoded as serialized array in C order (i.e. XYZ,XYZ,XYZ).

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

Where the lengths are the number of vertices in each path. Each path is assigned an ID numbering from 0. num_paths is a u64 little endian. Each len is the smallest data type that will hold Nv, the number of vertices in the skeleton part.

The next section is the links between paths:

pair_1, ..., pair_n

Where the pairs are: e1,e2 with the data type controlled by the header, though typically it will be the smallest data type that encodes vertices. The edges refer to the vertices, not to the path ID. 

`num_pairs = (len(edge_binary) - paths_section - 4) / edge_dtype_size / 2`

The 4 bytes refers to the size of the crc32c checksum.

The following python pseudocode algorithms can be implemented to decode or encode linked paths.

#### Decoding Linked Paths

```python
# edge_dtype is acquired from the header
# edge_binary is calculated from header information
def decode_linked_paths(edge_dtype, edge_binary:bytes) -> np.ndarray:
	crc32c = int.from_bytes(edge_binary[-4:], 'little')
	check_crc32c(edge_binary[:-4], crc32c)

	num_paths = int.from_bytes(edge_binary[:8], 'little')
	path_dtype = smallest_dtype(num_paths, [ np.uint8, np.uint16, np.uint32, np.uint64 ]) # smallest_dtype is a pseudocode function
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

### Edge Representation Type

| Type                  | Value | Description                                         | Properties                                  |
|-----------------------|-------|-----------------------------------------------------|---------------------------------------------|
| PAIR                  | 0     | edge list [(e1,e2), ... ]                           | General, no parsing. Most space.            |
| PARENT                | 1     | parent pointers [0,0,1] means (root<-node<-leaf)    | Trees only. Half space. Structure encoded.  |
| LINKED_PATHS          | 2     | set of paths linked at branch points                | General, efficient if branches are sparse.  |

### Compression Algorithm Type

| Algorithm | Value |
|-----------|-------|
| None      | 0     |
| gzip      | 1     |
| bzip2     | 2     |
| zstd      | 3     |
| draco     | 4     |
| xz        | 5     |

Note: Only None is currently supported.

Note: For Draco, preserve order must be used to preserve edge relationships to vertices. This unfortunately reduces the amount of compression possible.


### Unit Definitions

Arbitrary SI dimensions can be encoded as a tuple of fundemental units raised to a signed exponent. While
the original vision of this datastream is to only incorporate metadata that supports the interpretation of geometry, it seems strange to privilage length and time dimensions when all sorts of things might be measured along a skeleton, such as current or luminance.

For example, Joules can be expressed as W = ma x d or kg * m^2/s^2, Watts as kg * m^2/s^3 and Amperes can be expressed as C/s. Speed is meters/sec, area is meters^2 etc, luminousity can be measured in watts or photons per a second.

Since this is designed for biological use cases, the candela which is based in human perception of light is less useful, so we reserve those bits for future use (e.g. one can imagine using them for signaling the use of US customary units).

Therefore, for attributes, we encode the dimensions as a uint32 that represents the following structure:

| Field           | Bits             | Description                         |
|-----------------|------------------|-------------------------------------|
| coulombs        | 4                | c^x (-8 to 7)                       |
| grams           | 4                | g^x (-8 to 7)                       |
| kelvin          | 4                | k^x (-8 to 7)                       |
| meters          | 4                | m^x (-8 to 7)                       |
| mols            | 4                | mol^x (-8 to 7)                     |
| seconds         | 4                | s^x (-8 to 7)                       |
| SI Prefix       | 4                | signed 10^(X*3) where X is the value|
| Reserved        | 4                | Future use                          |

All signed values are written in 2's complement. The field components are stored little endian as is the field as a whole.

### Graph Type

2 bit

| Graph Structure       | Value | Description                |
|-----------------------|-------|----------------------------|
| graph                 | 0     | Could be any kind of graph.|
| tree                  | 1     | Acyclic graph.             |
| cyclic                | 2     | Contains one or more loops.|

This value is advisory and does not control the edge representation. 
This is because a tree can be represented as an edge list. See edge representation.

### Space Type

5 bits

| Space                        | Value |
|------------------------------|-------|
| Generic                      | 0     |
| Voxel                        | 1     |
| Physical                     | 2     |
| Scanner                      | 3     |
| Atlas                        | 4     |
| Aligned                      | 5     |
| World                        | 6     |
| Soma                         | 7     |
| Base                         | 8     |
| Joint                        | 9     |
| Tool                         | 10    |
| Model                        | 11    |
| Camera                       | 12    |
| Reserved                     | 13-23 |
| User Defined                 | 24-31 |

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