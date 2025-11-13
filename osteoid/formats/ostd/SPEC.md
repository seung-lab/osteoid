ostd ("osteoid") File Format
================

Skeletons (medial paths) may be represented as a possibly cyclic undirected graph with per vertex annotations and possibly edge annotations.

The `ostd` file format fills a gap in existing skeleton file formats by offering a self-contained, high performance, safe, binary format that supports optional vertex attributes for the serialization of skeletal structures. Incorporating metadata is not a design goal, as this core file is intended to be wrapped in a container file format for that purpose.

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

The parser should check if the stream is longer than the indicated content length. If so, it should check the subsequent section for the presence of another ostd file, and so forth. The vertices, edges, and attributes of each part then be combined after parsing. If an attribute is not present for a part, it should be filled with null values.

| Section | Description                  |
|---------|------------------------------|
| Part 1  | Initial ostd file.           |
| Part 2  | Appended vertices and edges. |
| ...     |                              |
| Part N  | Appended vertices and edges. |

## Header

All values are little endian except where noted. Total bytes: 88

| Field                  | Bytes | Datatype    | Value                       | Description                                                                                                       |
|------------------------|-------|-------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------|
| magic                  | 4     | string      | osdt                        | File magic number.                                                                        |
| format_version         | 1     | u8          | 0                           | Version of this file format.                                                              |
| total_bytes            | 8     | u64         | -                           | Total byte size of this part.                                                             |
| id                     | 16    | uuid4 / u64 | -                           | A wide enough field to accommodate both u64s and uuid4s                                   |
| flags                  | 8     | bitfield    | -                            | See note below for definitions.   |
| num_vertices (Nv)      | 8     | u64         | -                           | Number of vertices                                                                        |
| num_edges (Ne)         | 8     | u64         | -                           | Number of edges                                                                           |
units specified in flags.  |
| vertex_bytes           | 8     | u64         | -                           | Content length of vertices (needed for compression).                                      |
| edge_bytes             | 8     | u64         | -                           | Number of bytes encoding edges (needed for compression)                                   |
| spatial_index_bytes    | 4     | u32         | -                           | Content length in bytes of the spatial index.                                                         |
| attribute_header_bytes | 4     | u32         | -                           | Content length in bytes of the attribute header.                                                         |
| num_components         | 4     | u32         | N or (2^32-1 if unknown)    | Number of connected components in the skeleton graph. max value of uint32 is a sentinel for unknown.              |
| cable_length           | 4     | f32         | -                           | Physical path length of this object in the 
| current_space          | 1     | u8          | 0                           | The current transform space the vertices are in. By default 0. Every +1 means selecting the next transform from the transform list. See *Transform* | 
| crc16                  | 2     | uint16      | -                           | crc16 using 0xFFFF init and implicit polynomial 0xd175 of header bytes excluding magic number.                    |

### Flag Definitions

LSB on the left.

`VVVVeeeeCCCCccccGGGpppppPPPPOOOOOaaatoAEER*`

| Flag   | Meaning                            | Notes                                                                                |
| ------ | ---------------------------------- | ------------------------------------------------------------------------------------ |
| **V**  | Vertex data type                   | See *Data Types*                                                                     |
| **e**  | Edge data type                     | See *Data Types*                                                                     |
| **E**  | Edge representation                | See *Edge Representation*                                                            |
| **A**  | Append mode                        | (0) Section is self-contained<br>(1) Edges continue numbering from previous sections |
| **C**  | Compression algorithm for vertices | See *Compression Type*                                                               |
| **c**  | Compression algorithm for edges    | See *Compression Type*                                                               |
| **G**  | Graph structure (advisory)         | See *Graph Type*                                                                     |
| **p**  | SI Prefix                          | See *SIPrefixType*
| **P**  | Physical length units              | See *Length Type*                                                        |
| **a**  | Number of Axes                     | Number of axes                       |
| **t**  | Transforms present                 | bool                                                                                 |
| **O**  | Coordinate Frame Orientation       | Has own structure: `sssaaa`<br>s: sign of X,Y,Z axes in that order (0: positive, 1: negative)<br>a: axis permutation<br>See Axis Permutation Type, 000000 means +X+Y+Z standard frame |
| **o*** | Voxel centered or top left corner. | Describes whether voxel coordinates are interpreted as centered or in the top left corner.                                                              |
| **R*** | RESERVED                           | From this point forward                                                              |

## Transform

If this optional section is not present, the transform can assumed to be the 4x4 identity matrix. The reason this section is optional is to
significantly reduce the header overhead for small skeletons. 


| Field                  | Bytes | Datatype    | Value                       | Description                                                                                                       |
|------------------------|-------|-------------|-----------------------------|-------------------------------------------------------------------------------------------------------------------|
| num_spaces             | 1     | uint8       | -                           | Number of transformations available. matrices.                  |
| space                  | 1     | uint8       | -                           | The kind of space the transform represents. See *Space Type* |
| transform              | 64    | 4x4 f32s    | [ f32, f32, f32, f32, ... ] | Homogenous transform matrix from voxel to physical coordinates. Written in row major (C) order.                   |
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
| unit           | 2         | bitfield | QQQSSSSSUUUUURRR          | Physical unit.                  |
| num_components | 1         | uint8    | -                         | Number of components.           |
| content_length | 8         | uint64   | -                         | Length of compressed bitstream. |


#### Flags Definition

lsb on left

T: (0) vertex attribute (1) edge attribute  
D: Data Type  
C: Compression Type
R: Reserved

#### Unit Definition

lsb on left

Q: Quantity type (see *Physical Unit Type*)
S: SI Prefix (see *SI Prefix Type*)
U: Unit of that type

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

Note: Only None is currently supported.

### Physical Unit Type

| Unit Type         | Value |
|-------------------|-------|
| AreaType          | 0     |    
| DimensionlessType | 1     |             
| ElectricalType    | 2     |          
| EnergyType        | 3     |      
| LengthType        | 4     |      
| LuminosityType    | 5     |          
| MassType          | 6     |    
| SubstanceAmount   | 7     |           
| TemperatureType   | 8     |           
| TimeType          | 9     |    
| VolumeType        | 10    |      

### SI Prefix

5 bits

| Prefix | Value |
| ------ | ----- |
| none   | 0     |
| quecto | 1     |
| ronto  | 2     |
| yocto  | 3     |
| zepto  | 4     |
| atto   | 5     |
| femto  | 6     |
| pico   | 7     |
| nano   | 8     |
| micro  | 9     |
| milli  | 10    |
| centi  | 11    |
| deci   | 12    |
| deka   | 13    |
| hecto  | 14    |
| kilo   | 15    |
| mega   | 16    |
| giga   | 17    |
| tera   | 18    |
| peta   | 19    |
| exa    | 20    |
| zetta  | 21    |
| yotta  | 22    |
| ronna  | 23    |
| quetta | 24    |

### Length Type

4 bits

| Unit                   | Value |
|------------------------|-------|
| unknown (dimensionless)| 0     |
| voxel (dimensionless)  | 1     |
| meter                  | 2     |
| angstrom               | 3     | 
| astronomical unit      | 4     |
| lightyear              | 5     |
| parsec                 | 6     |
| mil (1/1000 inch)      | 7     |
| inch                   | 8     |
| foot                   | 9     |
| yard                   | 10    |
| statute mile           | 11    |
| nautical mile          | 12    |

### Area Type

4 bits

Inherits from length type but considers each value squared.
Additional units such as `acre` could in theory be added.

### Volume Type

4 bits

Inherits from length type but considers each value cubed.
The following values are added in addition.

| Unit                  | Value |
|-----------------------|-------|
| liter                 | 13    |

### Temperature Type

3 bits

| Unit        | Value |
|-------------|-------|
| unknown     | 0     |
| celsius     | 1     |
| fahrenheit  | 2     |
| rankine     | 3     |
| kelvin      | 4     |

### Time Type

3 bits

| Unit        | Value |
|-------------|-------|
| unknown     | 0     |
| second      | 1     |
| minute      | 2     |
| hour        | 3     |
| day         | 4     |
| month       | 5     |
| year        | 6     |
| hertz       | 7     |

### Luminosity Type

3 bits

| Unit                  | Value |
| --------------------- | ----- |
| unknown               | 0     |
| candela               | 1     |
| lumen                 | 2     |
| lux                   | 3     |
| photon                | 4     |
| photons per second    | 5     |

### Electrical Type

3 bits

| Unit            | Value |
| --------------- | ----- |
| unknown         | 0     |
| volt            | 1     |
| ampere          | 2     |
| ohm             | 3     |
| siemen          | 4     |
| farad           | 5     |
| henry           | 6     |
| coulomb         | 7     |

### Mass Type

2 bits

| Unit                  | Value |
|-----------------------|-------|
| unknown               | 0     |
| gram                  | 1     |
| dalton                | 2     |

### Substance Amount

1 bit

| Unit                  | Value |
|-----------------------|-------|
| unknown               | 0     |
| mole                  | 1     |

### Energy Type

2 bits

| Unit                  | Value |
|-----------------------|-------|
| unknown               | 0     |
| joule                 | 1     |
| watt                  | 2     |

### Graph Type

1 bit

| Graph Structure       | Value |
|-----------------------|-------|
| graph                 | 0     |
| tree                  | 1     |

This value is advisory and does not control the edge representation. 
This is because a tree can be represented as an edge list. See edge representation.

### Space Type

| Space                        | Value |
|------------------------------|-------|
| Generic                      | 0     |
| Physical                     | 1     |
| Scanner                      | 2     |
| Atlas                        | 3     |
| Aligned                      | 4     |
| World                        | 5     |
| Soma                         | 6     |
| Base                         | 7     |
| Joint                        | 8     |
| Tool                         | 9     |
| Model                        | 10    |
| Camera                       | 11    |

### Axis Permutation Type

| System | Value |
|--------|-------|
| XYZ    | 0     |
| XZY    | 1     |
| YXZ    | 2     |
| YZX    | 3     |
| ZXY    | 4     |
| ZYX    | 5     |
| XY     | 6     |
| YX     | 7     |