Ostd Format Usage
=================

This document shows how to use the ostd format implementation. See SPEC.md for controlling information on how to create or parse a file.

# Python

Most people will want to create a `osteoid.Skeleton` obect by reading or writing an ostd file just like they would with an SWC.

```python
import osteoid
skeleton = osteoid.load("example.ostd")
skeleton.save("example.ostd")
```

However, there is information in OSTD headers that could be parsed without reading in the entire file.

```python
from osteoid.formats.ostd import OstdSkeleton

skeleton = OstdSkeleton.load("example.ostd", allow_mmap=True)
skeleton.cable_length
skeleton.num_components
skeleton.num_vertices
skeleton.num_edges
```

# CLI

You can read and convert `.ostd` files with the ostd command line tool.

```bash
ostd info example.ostd # prints header and some attribute info
```

Convert between formats:

```bash
ostd convert example.swc example.ostd
ostd convert example.ostd example.swc
```

You can even view a skeleton:

```bash
ostd view example.swc
ostd view example.osdt
```



