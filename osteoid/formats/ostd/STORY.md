
# List of Skeleton Problems

- Selecting skeletons by bounding box
- Displaying numerous skeletons, predominantly streamlines from DTI, in Neuroglancer
- Assembling multi-GB skeletons
- Selecting parts of large skeletons by bounding box
- Working with and sharing individual skeletons with attributes other than radius on disk
- Working with many skeletons on-disk

## User Story

Discussion with Forrest on 2025-11-26:

"They see a bunch of skeletons going together in a bundle, they want to see them, and select them, and then zoom out in 3d and see where they go. I of course think they should have a voxel segmentation, but the reality is that their current pipelines dont make them. So they want skeletons to do the job that voxel segmentation does for us i.e. give a sense of the dense segmentation result with spatial indicesfrom an analysis point of view, being able to do a spatial query... get me all the IDs that are in this box, or all the skeletons that are covered by this mask which again is a thing that is easy to do with a voxel segmentation, but if you don't have it, you want your skeleton store to have that capability."

"They want to see a mesoscopic kind of representation of all the streamlines as a way to orient where they are in the brain and what coarse patterns there are, and then zoom in on specific spots."

## Selecting Skeletons by Bounding Box

`f(3d bounding box) -> list of segment IDs`

There is no skeleton spatial index in Neuroglancer. People are using segmentation as a proxy.

There is a JSON spatial index in CloudVolume/Igneous used for shard construction that differs from the binary Neuroglancer annotation index. This CloudVolume index is ingestible to MySQL, Postgres, and Sqlite3, though it is not organized for bounding box queries currently in that form.

In most cases, the JSON index will be suitable, but there is no facility to generate it from raw skeletons without segmentation currently. The size of the bounding box and chunking determine if a database implementation is needed to reduce requests and data size.

**Alternatives:** sqlite or other database with R-tree, octrees, using remnant .frag files, binary format based on annotation spatial index

**Most Likely Approach:** Coordinate with Jeremy to decide on spatial format and create a method for synthesizing it from raw skeletons. For efficiency, it should be ingestable to a local sqlite database with R-tree enabled that CloudVolume can consult. The CloudVolume spatial index db table `file_lookup` could be augmented with bounding boxes to accelerate query speed, though this will cause a growth in the disk space requirement. This format should then be introduced into Neuroglancer with a control for querying by bounding box.

## Displaying Numerous Skeletons in Neuroglancer

Diffusion Tensor Imaging (DTI) is an fMRI technique to identify streamlines of information transmission within brain scans. These streamlines are single paths with no branches. They are numerous. In one file `.trx` I studied there were 87,707,950 vertices contained in 787,772 streamlines written as float16s with a total file size of 437M.

When severely downsampled, this data was able to be represented as 24 MB or less (lzip: 9.2 MB), which suggests one awkward method for large scale low-res display of downloading all skeletons at once and displaying only those that are needed. Most likely though, this level of downsampling resulted in most streamlines having only a start and end point.

Connectomics data generally has larger skeletons with intricate branching patterns, making a trx streamline format unsuitable.

There are three issues with displaying numerous skeletons. (1) How to identify which skeletons are needed? This is addressed above. (2) How to request potentially hundreds or thousands of skeletons in a reasonable time frame? (3) How to reduce the data requirement to something reasonable.

(1) is close to being addressed, at least theoretically. (3) has been addressed previously by multi-resolution meshes, and it is thought that a similar idea could work here, though no one has worked this out in detail. (2) is relatively novel, though versions of this have been solved in annotations via the "relationship index", which pre-calculates chunks of data that are related. This pre-calculation strategy is not tenable when any bounding box can be queried, but the reach of any one skeleton may cross an entire dataset in an unknown direction.

The number of skeletons needed matters. Between 1-100 skeletons, fetching by ID will be sufficient and the multi-resolution concept will be fine.

A crude idea that could work: download all skeletons hugely downsampled, and display those that are needed at low res. Store all skeletons cut up by grids, and download the grid squre that is needed as a chunk when zooming occurs. This would give a roughly 1990s/early-2000s style interaction. One can posit an extension to an octree, but this will not work because the intermediate levels will contain too much data.... or will it? Each chunk will be downsampled and its extent truncated, limiting the amount of data.

This implies two types of ways to get skeletons: by ID and by chunk. How would Neuroglancer know when to use which?

When a list of IDs is selected, fetch by ID. When a spatial query occurs or greater than N skeletons are selected, fetch by chunks.

Ideally, a format could be found that strikes a balance between its ability to represent arbitrary skeletons and data usage. Precomputed uses an explicit edge list (2x u32 = 8 bytes per edge, 3 x f32 = 12 bytes per vertex). TRX allows 3 x f16 = 6 bytes per vertex, 0 edges, making it 30% the size, but much less powerful. Experiments show that a "path graph" can represent a generic skeleton graph using ~3% the number of edges, approaching streamline data efficiecny while retaining representational power. This efficiency relies on most of the skeleton structure being representable as streamlines. Extremely dense branching will restore competitiveness to edge lists.

An alternative approach that puts potentially less strain on the client, could be to store a spatial index and a set of skeletons on a server with an SSD and let it look them up by bounding box and send them to the client. This requires more infrastructure than simply generating a dataset for cloud storage though.

It is currently unknown how many individual skeletons can be displayed by Neuroglancer.

## Assembling Multi-GB Skeletons

Multi-GB skeletons result from taking large 3D light microscopy datasets and binarizing some subset of cells. The question then is how to represent these monsters such that they are usable?

It is essential that these large skeletons can be generated in fragments and then pieced together. Edge lists should be represented using uint64 as a path graph to ensure correctness and reduce data size. When uint64 is used, edges become larger than vertices (2 x 8 = 16 bytes, vs. 12 bytes), and so their near elimination is essential.

One approach to constructing a large skeleton is to append the fragments into a single file on disk. This file can then be scanned by a smaller machine and or consolidated into a single section format by a large machine.

In order to query the skeleton for 

TBD

## Working with and Sharing Individual Skeletons with Attributes other than Radius

| Format                | Description                                                                        |
-------------------------------------------------------------------------------------------------------------| 
| SWC                   | Text based, trees only. Radius supported.                                          |
| ESWC                  | SWC + unofficial attributes indicated positionally. Used by vaa3d                  |
| Pickle                | Represents arbitrary python objects, but executes code and tied to specific class. |
| Precomputed           | Separate header file used with a directory of objects.                             |
| Parquet               | Columnar data table used with PANDAS, used by NAVis                                |
| NWB                   | Neurodata Without Borders. Open standard that represents graphs with edge lists. HDF5 based. Unable to find anything about skeletons? ChatGPT hallucination?  |

A skeleton should efficiently represent arbitrary geometry, allow vertex and edge attributes, and be transmissible as a single file that can be opened anywhere, encode a spatial index (for large skeletons), and record basic metadata about its geometry so that large or numerous skeletons can be analyzed rapidly.

Some of these features could be added to SWC/ESWC or Precomputed, others are more problematic.

### Possibilities and Problems with Existing Formats

The existing formats have 

#### SWC

SWC could potentially be extended with additional header comments, making a more flexible "ESWC". For example, adding measurements of cable length, number of components (though usually an SWC is supposed to be a single tree I think), and the meaning of additional vertex attributes could be added. In the extreme case, a spatial index could even be added as a comment.

This would have the virtue of being backwards compatible with many readers and human readable. The downside is that it remains an inefficient text format that can only represent trees. For large skeletons, this becomes infeasible.

It may still be desirable to extend SWC to preserve these features while being backwards comptaible as an additional format to generate.


#### Precomputed

Precomputed shares a common header between objects, so its fundemental design is more akin to Zarr (which uses directories) than the idea of sharing a single file between collaborators. One path could be to extend Precomputed with uint64 edges and support for additional vertex data types. CloudVolume could be augmented to emit zipped Precomputed directories for sharing single files and we could finally make it easy to load mesh or skeleton directories in isolation.

A more radical reinterpretation of Precomputed would allow specifying a path graph as the edge representation.

A spatial index per a skeleton could be saved separately into a SQLite database or tacked onto the end of the binary.


#### Parquet 

Parquet provides single files which = single tables, but a skeleton requires two tables for vertices and edges. Metadata could be a third table.

#### Pickle

The mutability and code execution capability for Pickle means that it is not suitable for sharing files between strangers and limits applicability between programming languages.

### Possibilities in New Formats

In choosing a new format, there is much more freedom to make different choices. Here are key elements:

- Path graph edge representation
- Integrated spatial index
- uint64 and float16 capable binary format
- Appendability for large skeletons
- Pre-computed statistics like cable_length
- Recording physical dimensions, axis orientation, and voxel alignment (center or corner)
- Adding edge attributes
- Optional compression of individual fields

Simplicity can mean sharing a single file, but it can also mean sharing a format that is easy to read regardless of platform.

In this project, we define a draft spect for a new binary format `ostd` that can do all this very compactly and allows sharing individual skeleton files. However, conceptually, `ostd` is a set of tables.

- vertices
- vertex attributes
- vertex attributes metadata
- transforms
- edges
	- edge list 
	- path length list
	- edge attributes (for pairs this can just be a column)
- metadata k/v store

This could be alternatively implemented in Zarr, HDF5, a delta lake, or as a sqlite database, which would enable transparent access by existing tools (though interpreting the path graph requires some custom code).

Zarr in particular would get you "transforms" for free.


























Prior art for point clouds: 

COPC - cloud optimized point cloud: https://github.com/copcio/
COGeoTIFF - Cloud Optimized GeoTIFF
Entwine Point Clouds
Sharded - Zarr, Precomputed







