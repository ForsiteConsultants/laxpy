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

Test fixtures are in `laxpytest/data/` (`test.las`/`test.lax`, plus `bad_lax.lax` as a regression case for cell-count parsing). The test package is not shipped in the wheel.

## Architecture

The pipeline is `.lax` bytes -> cell intervals -> quadtree geometry -> interval reads:

- `laxpy/file.py` `LAXParser`: reads the whole `.lax` as 4-byte words (hard-coded offsets: bbox floats at words 9-12, cell count at word 15, cell records from word 16 on). `cells` (cached) maps cell_index -> flat `[start, end, start, end, ...]` point intervals (inclusive). `init_lax` shells out to `lasindex`.
- `laxpy/tree.py` `LAXTree`: reconstructs the quadtree implied by cell indices. Levels hold `4**i` cells; `trace_back` walks a cell up to root recording child position (1-4: SW, SE, NW, NE), and `get_cell_bbox` replays that path down from the file bbox. `cell_polygons` (cached) turns each cell into a shapely box. `plot`/`to_gdf` are debug helpers (matplotlib/geopandas are not install deps).
- `laxpy/query.py` `IndexedLAS`: parses the sibling `<name>.lax` once in `__init__`. `query(polygon, exact=True)` merges the intervals of intersecting cells, opens the file read-only with `laspy.open` and does `seek`/`read_points` per interval, and returns a new `laspy.LasData`. With `exact=True` points are filtered by `shapely.contains_xy` (strict interior; boundary points excluded); `exact=False` returns whole cells. No match -> zero-point `LasData` with the source header.

## Constraints

- Version lives only in `laxpy/__init__.py` (`pyproject.toml` reads it dynamically); on release also update `README.md` and `CHANGELOG.md`.
- CI: `.github/workflows/tests.yml` runs the tests on Python 3.11 and 3.14; `python-publish.yml` uploads to PyPI on GitHub release.
