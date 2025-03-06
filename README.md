# osteoid

Skeleton object used for representing neurons, adjacent cells, and organelles. 

Osteoid is the Skeleton code from [CloudVolume](https://github.com/seung-lab/cloud-volume) refactored into its own library.


## Installation

```bash
pip install osteoid
```

## Examples

```python
import osteoid

skel = osteoid.load("skeleton.swc")
osteoid.save("skeleton.swc", skel)

from osteoid import Skeleton, Bbox

skel = Skeleton(vertices, edges, radii=radii)

# you can specify a transform to e.g.
# convert the skeleton into a physical space
# with, in this example, 16x16x40 nm^3 resolution

matrix = np.array([
  [16, 0, 0, 0],
  [0, 16, 0, 0],
  [0, 0, 40, 0],
], dtype=np.float32)

skel = Skeleton(vertices, edges, radii=radii, transform=matrix)
skel = skel.physical_space() # applies transform to vertices
skel = skel.voxel_space() # removes transform from vertices

# skeleton functions

l = skel.cable_length() # physical length of the cable
comps = skel.components() # connected components
paths = skel.paths() # convert tree into a list of linear paths
skel = skel.downsample(factor) # factor must be a pos integer
skel = skel.average_smoothing(7) # smooths over a window of e.g. 7 vertices

skel2 = skel.crop([ minx, miny, minz, maxx, maxy, maxz ])
skel2 = skel.crop(Bbox([minx, miny, minz], [maxx, maxy, maxz]))

G = skel.to_networkx() # converts edges to an nx.Graph()

binary = skel.to_precomputed() # Neuroglancer compatible format
swc = skel.to_swc() # Neuroglancer compatible format
skel = Skeleton.from_swc(swc)
skel = Skeleton.from_precomputed(binary, segid=1, vertex_attributes=[
      {
        'id': 'radius',
        'num_components': 1,
        'data_type': 'float32',
      },
])

# Cross-format transform from Navis. More efficient than 
# using SWCs as an interchange medium. Navis is more fully featured,
# and written with a lot of love by Philipp Schlegel and others.
# Consider using it and citing them: https://github.com/navis-org/navis
skel = Skeleton.from_navis(navis_skel)

# remove duplicate vertices and optionally disconnected vertices
skel = skel.consolidate() 
skel2 = skel.clone()

# visualize, requires either matplotlib or microviewer+vtk installed
skel.viewer() # select library automatically

skel.viewer(color_by='radius', library='matplotlib')
skel.viewer(color_by='cross_section', library='matplotlib')
skel.viewer(color_by='component', library='matplotlib')

# gpu accelerated, fewer features, colors by component
skel.viewer(library='vtk')
```

