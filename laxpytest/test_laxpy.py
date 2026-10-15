from laxpy import *
import unittest
import os
import shutil
import struct
import tempfile
import warnings
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

    def test_query_chunked(self):
        expected = [self.test_ix_las.query(poly1, exact=exact) for exact in (True, False)]
        chunked = [self.test_ix_las.query(poly1, exact=exact, chunk_size=1000) for exact in (True, False)]
        for a, b in zip(chunked, expected):
            self.assertTrue(np.array_equal(a.points.array, b.points.array))
            self.assertEqual(a.header.point_count, len(b.points))
            np.testing.assert_allclose(a.header.mins, b.header.mins)
            np.testing.assert_allclose(a.header.maxs, b.header.maxs)

    def test_query_chunk_size_invalid(self):
        for chunk_size in (0, -1):
            self.assertRaises(ValueError, self.test_ix_las.query, poly1, chunk_size=chunk_size)

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


def evlr(user_id, record_id, payload):
    """A LASzip special EVLR record (60-byte header + payload)."""
    return struct.pack('<H16sHQ32s', 0, user_id, record_id, len(payload), b'') + payload


class EmbeddedLAXTestCase(unittest.TestCase):
    """Synthetic `lasindex -append` files: test.las as LAZ with special EVLRs appended."""

    @classmethod
    def setUpClass(cls):
        import io, laspy
        buf = io.BytesIO()
        laspy.read(test_las).write(buf, do_compress=True)
        cls.base_laz = buf.getvalue()
        cls.lax = Path(test_lax).read_bytes()

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.laz = os.path.join(self.tmp.name, 'tile.laz')

    def tearDown(self):
        self.tmp.cleanup()

    def write_laz(self, records=None, count=None, offset=None, truncate=0, sidecar=False):
        """Writes tile.laz with `records` appended and the LASzip VLR pointing at them."""
        records = [evlr(b'LAStools', 30, self.lax)] if records is None else records
        data = bytearray(self.base_laz)
        vlr_payload = data.index(b'laszip encoded') - 2 + 54
        struct.pack_into('<qq', data, vlr_payload + 16,
                         len(records) if count is None else count, len(data) if offset is None else offset)
        data += b''.join(records)
        Path(self.laz).write_bytes(bytes(data[:len(data) - truncate]))
        if sidecar:
            shutil.copy(test_lax, os.path.join(self.tmp.name, 'tile.lax'))

    def test_embedded_matches_sidecar(self):
        self.write_laz()
        self.assertEqual(file.read_embedded_lax(self.laz), self.lax)
        embedded, sidecar = IndexedLAS(self.laz), IndexedLAS(test_las)
        self.assertEqual(embedded.index_source, 'embedded')
        self.assertEqual(embedded.parser.cells, sidecar.parser.cells)
        for exact in (True, False):
            a, b = embedded.query(poly1, exact=exact), sidecar.query(poly1, exact=exact)
            self.assertTrue(np.array_equal(np.sort(a.points.array), np.sort(b.points.array)))

    def test_index_modes(self):
        self.write_laz(sidecar=True)
        with self.assertWarns(UserWarning):
            self.assertEqual(IndexedLAS(self.laz).index_source, 'embedded')
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            self.assertEqual(IndexedLAS(self.laz, index='embedded').index_source, 'embedded')
            self.assertEqual(IndexedLAS(self.laz, index='sidecar').index_source, 'sidecar')
        with self.assertRaises(ValueError):
            IndexedLAS(self.laz, index='other')

    def test_missing_index(self):
        self.write_laz(records=[])
        for index in ('auto', 'embedded', 'sidecar'):
            with self.assertRaises(FileNotFoundError):
                IndexedLAS(self.laz, index=index)
        self.assertRaises(FileNotFoundError, IndexedLAS, test_las, index='embedded')
        self.write_laz(sidecar=False)
        self.assertRaises(FileNotFoundError, IndexedLAS, self.laz, index='sidecar')

    def test_unrelated_records(self):
        self.write_laz(records=[evlr(b'other', 1, b'x' * 10), evlr(b'LAStools', 30, self.lax), evlr(b'other', 2, b'')])
        self.assertEqual(file.read_embedded_lax(self.laz), self.lax)
        self.write_laz(records=[evlr(b'other', 1, b'x' * 10)])
        self.assertIsNone(file.read_embedded_lax(self.laz))

    def test_unused_pointers(self):
        self.write_laz(records=[], count=-1, offset=-1)
        self.assertIsNone(file.read_embedded_lax(self.laz))
        self.assertIsNone(file.read_embedded_lax(test_las))

    def test_malformed_embedded_never_falls_back(self):
        bad_signature = b'XXXX' + self.lax[4:]
        cases = {
            'bad offset': dict(offset=len(self.base_laz) * 2),
            'bad count': dict(count=-5),
            'count past chain': dict(count=2),
            'truncated payload': dict(truncate=10),
            'duplicate': dict(records=[evlr(b'LAStools', 30, self.lax)] * 2),
            'bad signature': dict(records=[evlr(b'LAStools', 30, bad_signature)]),
            'truncated cells': dict(records=[evlr(b'LAStools', 30, self.lax[:200])]),
        }
        for name, kwargs in cases.items():
            with self.subTest(name):
                self.write_laz(sidecar=True, **kwargs)
                self.assertRaises(ValueError, IndexedLAS, self.laz)

    def test_out_of_range_interval(self):
        lax = bytearray(self.lax)
        struct.pack_into('<I', lax, 20 * 4, 10 ** 9)  # end of the first interval
        self.write_laz(records=[evlr(b'LAStools', 30, bytes(lax))])
        self.assertRaises(ValueError, IndexedLAS, self.laz)

    def test_source_unchanged_and_index_not_propagated(self):
        import laspy
        self.write_laz()
        before = Path(self.laz).read_bytes()
        ix = IndexedLAS(self.laz)
        out = os.path.join(self.tmp.name, 'out.laz')
        ix.query(poly1).write(out)
        ix.query(poly2, exact=False)
        self.assertEqual(Path(self.laz).read_bytes(), before)
        self.assertIsNone(file.read_embedded_lax(out))
        self.assertEqual(len(laspy.read(out).points), 1258)
