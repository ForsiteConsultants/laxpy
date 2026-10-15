import copy
import os.path
import laspy
import numpy as np
import shapely
from laxpy import file, tree


class IndexedLAS:
    """
    The main interface for querying indexed `.las`/`.laz` files. On initialization the matching `.lax` file (same
    name, same directory) is parsed once; subsequent queries read only the points of the intersecting cells. If no
    such file exists please use `lasindex` (or `laxpy.file.init_lax`) to create it.

    :param path: The path of a `.las` or `.laz` file to index.
    """
    def __init__(self, path):
        self.path = path
        self.lax_path = os.path.splitext(path)[0] + '.lax'

        if not os.path.isfile(self.lax_path):
            raise FileNotFoundError('A matching .lax file was not found.')

        self.parser = file.LAXParser(self.lax_path)
        self.tree = tree.LAXTree(self.parser)

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

    def query(self, q_polygon, exact=True):
        """
        Reads the points within a polygon. The source file is opened read-only and is not modified.

        :param q_polygon: A shapely (Multi)Polygon to query.
        :param exact: If True, return only points strictly inside `q_polygon`. If False, return all points of the
            index cells that intersect `q_polygon` (a superset).
        :return: A `laspy.LasData` holding the selected points. It has zero points if none were found.
        """
        with laspy.open(self.path) as reader:
            header = copy.deepcopy(reader.header)
            arrays = [np.empty(0, dtype=header.point_format.dtype())]
            for start, end in self._intervals(q_polygon):
                reader.seek(start)
                arrays.append(reader.read_points(end - start + 1).array)

        las = laspy.LasData(header, laspy.PackedPointRecord(np.concatenate(arrays), header.point_format))
        if exact:
            las.points = las.points[shapely.contains_xy(q_polygon, las.x, las.y)]
        return las
