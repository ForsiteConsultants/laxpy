# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`laxpy` reads `.lax` spatial index files (produced by LAStools `lasindex`) and uses them to clip `.las`/`.laz` point clouds to a shapely polygon without loading the whole file. Targets laspy 2.x, shapely 2.x, Python >= 3.11.

## Commands

```bash
pip install .                                   # install (add [laz] for .laz support)
python -m unittest laxpytest.test_laxpy         # run all tests
python -m unittest laxpytest.test_laxpy.LAXTreeTestCase.test_get_cell_bbox   # single test
python -m build                                 # build sdist/wheel (what the publish workflow uses)
```

Test fixtures are in `laxpytest/data/` (`test.las`/`test.lax`, plus `bad_lax.lax` as a regression case for cell-count parsing). Embedded-index LAZ files are synthesized at test time (`EmbeddedLAXTestCase.write_laz`), which needs `lazrs` (`pip install ".[laz]"`). The test package is not shipped in the wheel.

## Architecture

The pipeline is `.lax` bytes -> cell intervals -> quadtree geometry -> interval reads:

- `laxpy/file.py` `LAXParser`: reads the whole `.lax` as 4-byte words (hard-coded offsets: bbox floats at words 9-12, cell count at word 15, cell records from word 16 on). Accepts a path or the payload as bytes and checks the `LASX`/`LASS`/`LASQ`/`LASV` block signatures. `cells` (cached) maps cell_index -> flat `[start, end, start, end, ...]` point intervals (inclusive). `read_embedded_lax` returns the payload of an index embedded by `lasindex -append`: it parses the LAS header and VLRs from raw bytes, takes the special-EVLR count/offset from the LASzip VLR (`laszip encoded`/22204, payload offsets 16/24, `-1` = unused) and walks that chain (60-byte headers) for `LAStools`/30. These special EVLRs are not part of the standard EVLR chain, so laspy does not expose them. `init_lax` shells out to `lasindex`.
- `laxpy/tree.py` `LAXTree`: reconstructs the quadtree implied by cell indices. Levels hold `4**i` cells; `trace_back` walks a cell up to root recording child position (1-4: SW, SE, NW, NE), and `get_cell_bbox` replays that path down from the file bbox. `cell_polygons` (cached) turns each cell into a shapely box. `plot`/`to_gdf` are debug helpers (matplotlib/geopandas are not install deps).
- `laxpy/query.py` `IndexedLAS(path, index='auto'|'embedded'|'sidecar')`: loads the embedded index or the sibling `<name>.lax` once in `__init__` (auto prefers embedded and warns if both exist; a malformed embedded index raises rather than falling back) and validates intervals against the file's point count. `query(polygon, exact=True)` merges the intervals of intersecting cells, opens the file read-only with `laspy.open` and does `seek`/`read_points` per interval, and returns a new `laspy.LasData`. With `exact=True` points are filtered by `shapely.contains_xy` (strict interior; boundary points excluded); `exact=False` returns whole cells. No match -> zero-point `LasData` with the source header.

## Constraints

- Version lives only in `laxpy/__init__.py` (`pyproject.toml` reads it dynamically); on release also update `README.md` and `CHANGELOG.md`.
- CI: `.github/workflows/tests.yml` runs the tests on Python 3.11 and 3.14; `python-publish.yml` uploads to PyPI on GitHub release.
