#!/usr/bin/env python3
"""Reader for the binary Ac-field snapshots written by write_Ac_field_bin
(src/ssbe/multiscale_ssbe.f90, enabled with dbg_sbe_fast(2) = 1).

File layout (Fortran unformatted sequential, one record per line):
    ix1, ix2                                    (int32 x 2)
    iy1, iy2                                    (int32 x 2)
    iz1, iz2                                    (int32 x 2)
    vec_Ac(1:3, ix1:ix2, iy1:iy2, iz1:iz2)      (float64)
    vec_e(1:3, ...)                             (float64)
    vec_h(1:3, ...)                             (float64)
    vec_j_em(1:3, ...)                          (float64)
    edensity_emfield(ix1:ix2, iy1:iy2, iz1:iz2) (float64)
    edensity_absorb(ix1:ix2, iy1:iy2, iz1:iz2)  (float64)

The returned arrays keep the Fortran index order:
    data["vec_Ac"][i, ix, iy, iz]   has shape (3, nx, ny, nz)
    data["edensity_absorb"][ix, iy, iz] has shape (nx, ny, nz)
where ix = 0 corresponds to the Fortran grid index ix1 (likewise for y, z).

Usage:
    python3 read_ac_bin.py SYSNAME_sbe_RT_Ac/SYSNAME_Ac_000100.bin
"""
import sys

import numpy as np


def _read_record(fh, dtype, count):
    """Read one Fortran sequential record (4-byte length markers)."""
    head = np.fromfile(fh, dtype=np.int32, count=1)
    if head.size != 1:
        raise EOFError("unexpected end of file")
    nbyte = int(head[0])
    expected = np.dtype(dtype).itemsize * count
    if nbyte != expected:
        raise ValueError(
            f"record length mismatch: file has {nbyte} bytes, expected {expected}"
        )
    body = np.fromfile(fh, dtype=dtype, count=count)
    tail = np.fromfile(fh, dtype=np.int32, count=1)
    if body.size != count or tail.size != 1 or int(tail[0]) != nbyte:
        raise ValueError("corrupted record (trailing length marker mismatch)")
    return body


def read_ac_bin(filename):
    data = {}
    with open(filename, "rb") as fh:
        ix1, ix2 = (int(v) for v in _read_record(fh, np.int32, 2))
        iy1, iy2 = (int(v) for v in _read_record(fh, np.int32, 2))
        iz1, iz2 = (int(v) for v in _read_record(fh, np.int32, 2))
        nx = ix2 - ix1 + 1
        ny = iy2 - iy1 + 1
        nz = iz2 - iz1 + 1
        data["ix1"], data["ix2"] = ix1, ix2
        data["iy1"], data["iy2"] = iy1, iy2
        data["iz1"], data["iz2"] = iz1, iz2

        # Fortran column-major (3, nx, ny, nz) == C row-major (nz, ny, nx, 3);
        # .T restores the Fortran index order.
        for key in ("vec_Ac", "vec_e", "vec_h", "vec_j_em"):
            tmp3 = _read_record(fh, np.float64, 3 * nx * ny * nz)
            data[key] = tmp3.reshape([nz, ny, nx, 3]).T
        for key in ("edensity_emfield", "edensity_absorb"):
            tmp1 = _read_record(fh, np.float64, nx * ny * nz)
            data[key] = tmp1.reshape([nz, ny, nx]).T
    return data


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    for fname in sys.argv[1:]:
        d = read_ac_bin(fname)
        print(f"# {fname}")
        print(f"#   ix: {d['ix1']}..{d['ix2']}  iy: {d['iy1']}..{d['iy2']}  iz: {d['iz1']}..{d['iz2']}")
        for k in ("vec_Ac", "vec_e", "vec_h", "vec_j_em", "edensity_emfield", "edensity_absorb"):
            print(f"#   {k:17s} shape={d[k].shape}  min={d[k].min(): .6e}  max={d[k].max(): .6e}")
