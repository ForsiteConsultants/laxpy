from laxpy import *
import unittest
import os
from pathlib import Path
import numpy as np
from shapely.geometry import box

data_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
test_las = os.path.join(data_dir, 'test.las')
test_lax = os.path.join(data_dir, 'test.lax')
bad_lax = os.path.join(data_dir, 'bad_lax.lax')
poly1 = box(405000, 3276400, 405020, 3276420)
poly2 = box(405021, 3276421, 405041, 3276441)

class LAXParserTestCase(unittest.TestCase):
    def setUp(self):
        self.test_parser = file.LAXParser(test_lax)

    def test_cells(self):
        self.assertEqual(len(self.test_parser.cells), 40)

    def test_create_indices(self):
        self.test_parser.create_point_indices(14)

    def test_cell_length_equal_to_parsed_length(self):
        bad_parser = file.LAXParser(bad_lax)
        bad_tree = tree.LAXTree(bad_parser)
        self.assertEqual(bad_parser.number_cells, len(bad_parser.cells.keys()))

class LAXTreeTestCase(unittest.TestCase):
    def setUp(self):
        self.test_tree = tree.LAXTree(file.LAXParser(test_lax))

    def test_tree_level_sizes(self):
        self.assertListEqual(self.test_tree.tree_level_sizes, [1, 4, 16, 64, 256, 1024])

    def test_level_edges(self):
        self.assertDictEqual(self.test_tree.level_edges, {1: (1, 4), 2: (5, 20), 3: (21, 84), 4: (85, 340), 5: (341, 1364), 6: (1365, 5460)})

    def test_get_cell_level_edges(self):
        self.assertEqual(self.test_tree.get_cell_level_edges(14), (2, (5, 20)))

    def test_get_parent_cell(self):
        self.assertEqual(self.test_tree.get_parent_cell(14), 3)

    def test_trace_back(self):
        self.assertListEqual(self.test_tree.trace_back(14), [(3,2), (0,3)])

    def test_get_cell_bbox(self):
        self.assertEqual(self.test_tree.get_cell_bbox(14), (405020.0, 405100.0,
                                                       3276400.0, 3276480.0))

    def test_cell_polygons(self):
        self.test_tree.cell_polygons

class IndexedLASTestCase(unittest.TestCase):
    def setUp(self):
        self.test_ix_las = IndexedLAS(test_las)

    def test_query_cell(self):
        self.assertEqual(len(self.test_ix_las._query_cell(14)), 67814)

    def test_query_cells(self):
        # exact=False returns every point of the intersecting index cells
        self.assertEqual(len(self.test_ix_las.query(poly1, exact=False).points), 140455)
        self.assertEqual(len(self.test_ix_las.query(poly2, exact=False).points), 67814)

    def test_query_cells_header(self):
        # Header must describe the returned subset before any write
        las = self.test_ix_las.query(poly1, exact=False)
        self.assertLess(len(las.points), 217222)
        self.assertEqual(las.header.point_count, len(las.points))
        np.testing.assert_allclose(las.header.mins, [las.x.min(), las.y.min(), las.z.min()])
        np.testing.assert_allclose(las.header.maxs, [las.x.max(), las.y.max(), las.z.max()])
        by_return = np.bincount(las.return_number, minlength=16)[1:16]
        np.testing.assert_array_equal(las.header.number_of_points_by_return, by_return)

    def test_query_cells_empty(self):
        las = self.test_ix_las.query(box(0, 0, 1, 1), exact=False)
        self.assertEqual(len(las.points), 0)
        self.assertEqual(las.header.point_count, 0)

    def test_query_exact(self):
        import laspy, shapely
        full = laspy.read(test_las)
        for poly, count in ((poly1, 1258), (poly2, 1138)):
            las = self.test_ix_las.query(poly)
            self.assertEqual(len(las.points), count)
            brute = full.points.array[shapely.contains_xy(poly, full.x, full.y)]
            self.assertTrue(np.array_equal(np.sort(las.points.array), np.sort(brute)))

    def test_query_empty(self):
        las = self.test_ix_las.query(box(0, 0, 1, 1))
        self.assertEqual(len(las.points), 0)
        self.assertEqual(las.header.point_format.id, 1)

    def test_query_leaves_source_untouched(self):
        import hashlib
        digest = lambda: hashlib.sha256(Path(test_las).read_bytes()).hexdigest()
        before = digest()
        self.test_ix_las.query(poly1)
        self.test_ix_las.query(poly2, exact=False)
        self.assertEqual(digest(), before)
