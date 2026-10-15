## 0.4.0
- Queries read candidate intervals in chunks of `query(..., chunk_size=1_000_000)` points. Exact queries filter each
  chunk, so peak memory scales with the result rather than with all candidate points.
- Embedded LAX indexes (LASzip special EVLR written by `lasindex -append`) are supported.
  `IndexedLAS(path, index='auto' | 'embedded' | 'sidecar')`; `'auto'` prefers the embedded index and warns if a
  sidecar also exists. A malformed embedded index raises instead of falling back to the sidecar.
- LAX indexes are validated: block signatures, truncation, and point intervals against the file's point count
  (catches indexes that are stale because the file now has fewer points). Short reads raise `ValueError`.
- `LAXParser` accepts the LAX payload as bytes as well as a path; new `file.read_embedded_lax(path)`.
- Fixed: `query(..., exact=False)` returned a header (point count, bounds, return counts) describing the source file
  instead of the returned points.

## 0.3.0
- Ported to laspy 2.x; laspy 1.x is no longer supported (use 0.2.3 with laspy 1.x).
- Breaking: `IndexedLAS` no longer subclasses `laspy.file.File`. `map_polygon`/`undo_map` are replaced by
  `IndexedLAS.query(polygon, exact=True)`, which returns a new `laspy.LasData` and leaves the source file untouched.
- Queries clip exactly to the polygon by default (`shapely.contains_xy`); `exact=False` returns whole index cells.
- Empty queries return a zero-point `LasData` instead of warning and leaving all points mapped.
- `.laz` input supported via the `laz` extra.
- `.lax` cells and cell polygons are computed once per `IndexedLAS`.
- Removed `laxpy.clip` and the `numba` dependency.
- `init_lax` now raises `FileNotFoundError` when `lasindex` is missing.
- Packaging moved to `pyproject.toml`; requires Python >= 3.11, laspy >= 2.0, shapely >= 2.0.
- Tests run on GitHub Actions; Travis config and `publish.sh` removed.

## 0.2.3
- Updated deprecated aliases `np.int` to builtin `int`s

## 0.2.2

## 0.2.1

## 0.2.0

## 0.1.9
- Hotfix for empty cell problem, fixes #6

## 0.1.8
- Fixing small bug in `query` that errored on an empty query.

## 0.1.6

- Adding optional conda environment
- Adding Windows testing on Travis
- Reduced the structure of `file.init_lax` check for `lasindex` installation.

## 0.1.2a
- Added warning for empty mappings.

## 0.1.2

- Fixed issue with missing dependencies on pip install

## 0.1

Pre-release version.

### `IndexedLas`

- Removed scale argument from `.query_polygon`, it is now required to scale  (because to query a polygon requires points to be georeferneced, i.e. scaled).


