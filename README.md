Relax with `laxpy`. Let it handle your spatial queries for you!

laxpy is a way to read and write `.lax` files generated from `lasindex` in Python. I have ambitions to merge this into
the `laspy` codebase, but for now this standalone package can handle spatial queries in conjuction with `laspy` and 
presence of `.lax` files that match the name of the `.las` file you are interested in querying.

### Why Index?

Indexing LiDAR data is only relevant if the entire file is not needed. This most frequently happens when LiDAR files need
to be clipped using some geometry smaller than the original file extent. Indexing allows loading smaller chunks
of the file into memory, rather than the whole thing, at the cost of some initial overhead.

Whether or not indexing ultimately increases computational efficiency depends on the clipping geometry, and the structure
and size of the LiDAR file.

## Release Status

Current Release: 0.3.0

Status: `laxpy` is in a stable condition.

## Installation

```{bash}
# via pip, from GitHub
pip install git+https://github.com/ForsiteConsultants/laxpy

# with .laz support
pip install "laxpy[laz] @ git+https://github.com/ForsiteConsultants/laxpy"

# from a local clone
git clone https://github.com/ForsiteConsultants/laxpy.git
cd laxpy
pip install .
```

The `laxpy` packages on PyPI and conda-forge are the upstream 0.2.3 release (laspy 1.x) and do not include these changes.

## Example

### Clipping a Polygon

`laxpy` reads only the points of the index cells that intersect a `shapely` polygon, then clips them exactly to the
polygon. The result is a regular `laspy.LasData`; the source file is never modified.

```{python}
from laxpy import IndexedLAS

ix = IndexedLAS('my_las.las')  # requires my_las.lax in the same directory; .laz works too
las = ix.query(my_shapely_polygon)

print(las.x, las.points)
las.write('clipped.las')

# All points of the intersecting index cells, without the exact clip
candidates = ix.query(my_shapely_polygon, exact=False)
```

The `.lax` file is parsed once when `IndexedLAS` is created, so querying many polygons against one file is cheap.
A query with no matching points returns a `LasData` with zero points.

Requires Python >= 3.11 and laspy >= 2.0. For laspy 1.x, use laxpy 0.2.3.
