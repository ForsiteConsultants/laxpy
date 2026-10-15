import copy
import os.path
import laspy
import numpy as np
import shapely
import warnings
from pathlib import Path
from laxpy import file, tree


class IndexedLAS:
    """
    The main interface for querying indexed `.las`/`.laz` files. On initialization the LAX index is loaded and
    validated once; subsequent queries read only the points of the intersecting cells.

    The index is either embedded in a LAZ file (`lasindex -append`) or a sidecar `.lax` file with the same name in
    the same directory. Use `lasindex` (or `laxpy.file.init_lax`) to create one.

    :param path: The path of a `.las` or `.laz` file to index.
    :param index: Which index to use: 'embedded', 'sidecar', or 'auto' (embedded if present, else sidecar; warns
        if both exist). A malformed embedded index raises `ValueError` and never falls back to the sidecar.
    """
    def __init__(self, path, index='auto'):
        if index not in ('auto', 'embedded', 'sidecar'):
            raise ValueError("index must be 'auto', 'embedded' or 'sidecar', not {!r}".format(index))
        self.path = path
        self.lax_path = os.path.splitext(path)[0] + '.lax'
        has_sidecar = os.path.isfile(self.lax_path)

        embedded = None if index == 'sidecar' else file.read_embedded_lax(path)
        match index, embedded, has_sidecar:
            case 'auto' | 'embedded', bytes(), _:
                self.index_source, payload = 'embedded', embedded
            case 'auto' | 'sidecar', None, True:
                self.index_source, payload = 'sidecar', Path(self.lax_path).read_bytes()
            case 'embedded', None, _:
                raise FileNotFoundError('No embedded LAX index found in {}.'.format(path))
            case 'sidecar', _, False:
                raise FileNotFoundError('A matching .lax file was not found: {}'.format(self.lax_path))
            case 'auto', None, False:
                raise FileNotFoundError('No embedded LAX index or matching .lax file found for {}.'.format(path))

        try:
            self.parser = file.LAXParser(payload)
            self.parser.cells
        except ValueError as e:
            raise ValueError('{} index of {}: {}'.format(self.index_source.capitalize(), path, e)) from e
        self.tree = tree.LAXTree(self.parser)

        with laspy.open(path) as reader:
            point_count = reader.header.point_count
        bounds = np.array([b for cell in self.parser.cells.values() for b in cell]).reshape(-1, 2)
        if (bounds[:, 0] > bounds[:, 1]).any() or (bounds >= point_count).any():
            raise ValueError('{} index of {} references points outside 0..{}; it may be stale.'
                             .format(self.index_source.capitalize(), path, point_count - 1))

        if self.index_source == 'embedded' and index == 'auto' and has_sidecar:
            warnings.warn('{} has both an embedded index and {}; using the embedded index.'.format(path, self.lax_path))

    def _query_cell(self, cell_index):
        """
        Returns the point indices of a given cell index. This is generally used internally.
        """
        return self.parser.create_point_indices(cell_index)

    def _intervals(self, q_polygon):
        """
        :return: Sorted, merged, inclusive (start, end) point intervals of all cells intersecting `q_polygon`.
        """
        intervals = []
        for cell_index, polygon in self.tree.cell_polygons.items():
            if q_polygon.intersects(polygon):
                bounds = self.parser.cells[cell_index]
                intervals.extend(zip(bounds[::2], bounds[1::2]))

        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        return merged

    def query(self, q_polygon, exact=True, chunk_size=1_000_000):
        """
        Reads the points within a polygon. The source file is opened read-only and is not modified.

        :param q_polygon: A shapely (Multi)Polygon to query.
        :param exact: If True, return only points strictly inside `q_polygon`. If False, return all points of the
            index cells that intersect `q_polygon` (a superset).
        :param chunk_size: Maximum number of points read at once. With `exact=True` each chunk is filtered before
            being kept, so this bounds peak memory; with `exact=False` all candidates are kept regardless.
        :return: A `laspy.LasData` holding the selected points. It has zero points if none were found.
        """
        if chunk_size < 1:
            raise ValueError('chunk_size must be at least 1, not {}'.format(chunk_size))
        with laspy.open(self.path) as reader:
            header = copy.deepcopy(reader.header)
            arrays = [np.empty(0, dtype=header.point_format.dtype())]
            for start, end in self._intervals(q_polygon):
                reader.seek(start)
                while start <= end:
                    n = min(end - start + 1, chunk_size)
                    points = reader.read_points(n)
                    if len(points) != n:
                        raise ValueError('{}: expected {} points at {}, read {}'.format(self.path, n, start, len(points)))
                    if exact:
                        points = points[shapely.contains_xy(q_polygon, points.x, points.y)]
                    arrays.append(points.array)
                    start += n

        las = laspy.LasData(header, laspy.PackedPointRecord(np.concatenate(arrays), header.point_format))
        las.update_header()
        return las
