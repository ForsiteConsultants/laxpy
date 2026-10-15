import struct
import numpy as np
import os
from pathlib import Path
import subprocess
from functools import cached_property

# Byte offset -> block signature of the LAX layout this parser understands
SIGNATURES = {0: b'LASX', 8: b'LASS', 16: b'LASQ', 52: b'LASV'}
LASZIP_VLR = (b'laszip encoded', 22204)
LAX_EVLR = (b'LAStools', 30)


class LAXParser:
    """
    Parses a LAX index and generates point index intervals.

    :param source: The path of a `.lax` file, or the LAX payload as bytes (see `read_embedded_lax`).
    """
    def __init__(self, source):
        data = bytes(source) if isinstance(source, (bytes, bytearray, memoryview)) else Path(source).read_bytes()
        if len(data) < 64:
            raise ValueError('LAX index is truncated ({} bytes)'.format(len(data)))
        for offset, signature in SIGNATURES.items():
            if data[offset:offset + 4] != signature:
                raise ValueError('Unsupported LAX index: expected {!r} at byte {}'.format(signature, offset))

        self.bbox = list(struct.unpack_from('<4f', data, 36))
        self.number_cells = struct.unpack_from('<I', data, 60)[0]
        self.parsed_bytes = np.frombuffer(data, '<u4', count=len(data) // 4).tolist()

    def header(self):
        raise NotImplementedError('To be implemented.')
        pass

    @cached_property
    def cells(self):
        """
        :return: A dictionary of cell_index: [cell_intervals] key, value pairs.
        """

        # Index of first interval
        start_pos = 19

        cell_dict = {}

        n = self.parsed_bytes[17] # The number of intervals of the first cell
        for i in range(self.number_cells):
            if start_pos + n * 2 > len(self.parsed_bytes):
                raise ValueError('LAX index is truncated in cell {} of {}'.format(i, self.number_cells))
            cell_intervals = self.parsed_bytes[start_pos:start_pos + n * 2]
            cell_id = self.parsed_bytes[start_pos-3]
            cell_dict[cell_id] = cell_intervals

            start_pos = (start_pos + n*2 + 3)

            if (i != self.number_cells - 1):
                n = self.parsed_bytes[start_pos - 2]

        return cell_dict

    def create_point_indices(self, cell_index):
        """
        The values in the self.cells dictionary are the raw intervals. This function converts those raw intervals
        into a set of point indices that can be used to index the points of the `.las` file.

        :param cell_index: The cell index to generate point indices for.
        :return: An `ndarray` of indices of points associated with the cell specified by `cell_index`
        """

        indices = []
        i = 0
        bounds = self.cells[cell_index]
        while i < len(bounds):
            indices.append(np.arange(bounds[i], bounds[i + 1] + 1))
            i += 2

        return np.concatenate(indices)

def init_lax(path):
    """
    Creates a .lax file for a specified .las file. A simple wrapper for lasindex CLI tool.
    :return:
    """
    try:
        subprocess.run(["lasindex", "-i", path])
    except FileNotFoundError:
        raise FileNotFoundError("lasindex was not found on the system")


def _record_key(user_id, record_id):
    return user_id.split(b'\0')[0], record_id


def read_embedded_lax(path):
    """
    Reads a LAX index embedded in a LAZ file as a LASzip special EVLR (as written by `lasindex -append`).

    :param path: The path of a `.las` or `.laz` file.
    :return: The LAX payload as bytes, or None if the file has no embedded index.
    :raises ValueError: If the LASzip metadata or the special EVLR chain is malformed.
    """
    size = os.path.getsize(path)
    with open(path, 'rb') as stream:
        header = stream.read(104)
        if len(header) < 104 or header[:4] != b'LASF':
            raise ValueError('{}: not a LAS/LAZ file'.format(path))
        header_size, _, number_vlrs = struct.unpack_from('<HII', header, 94)

        stream.seek(header_size)
        for _ in range(number_vlrs):
            vlr = stream.read(54)
            if len(vlr) < 54:
                raise ValueError('{}: truncated VLR header'.format(path))
            _, user_id, record_id, length, _ = struct.unpack('<H16sHH32s', vlr)
            if _record_key(user_id, record_id) == LASZIP_VLR:
                if length < 32:
                    raise ValueError('{}: LASzip VLR too short ({} bytes)'.format(path, length))
                count, offset = struct.unpack_from('<qq', stream.read(32), 16)
                break
            stream.seek(length, os.SEEK_CUR)
        else:
            return None  # not LASzip compressed

        if count in (-1, 0):
            return None
        if count < 0 or offset < header_size or offset + 60 > size:
            raise ValueError('{}: invalid special EVLR count {} or offset {}'.format(path, count, offset))

        stream.seek(offset)
        payloads = []
        for _ in range(count):
            evlr = stream.read(60)
            if len(evlr) < 60:
                raise ValueError('{}: truncated special EVLR header'.format(path))
            _, user_id, record_id, length, _ = struct.unpack('<H16sHQ32s', evlr)
            if stream.tell() + length > size:
                raise ValueError('{}: special EVLR payload extends past the end of the file'.format(path))
            if _record_key(user_id, record_id) == LAX_EVLR:
                payloads.append(stream.read(length))
            else:
                stream.seek(length, os.SEEK_CUR)

    if len(payloads) > 1:
        raise ValueError('{}: {} embedded LAX indexes found'.format(path, len(payloads)))
    return payloads[0] if payloads else None
