"""Write metadata JSON for a table produced by exact/bfs_fast.c."""

import hashlib
import json
import os
import sys

from bfs import factorials

TABLES = os.path.join(os.path.dirname(__file__), "..", "data", "tables")


def main():
    n = int(sys.argv[1])
    note = sys.argv[2] if len(sys.argv) > 2 else ""
    path = os.path.join(TABLES, f"dist_n{n}.bin")
    with open(path, "rb") as f:
        data = f.read()
    assert len(data) == factorials(n)[n]
    assert 0xFF not in data, "incomplete table"
    meta = {
        "n": n,
        "generators": "L (left cyclic shift), R (right cyclic shift), "
                      "X (swap first two); PROBLEM.md conventions",
        "indexing": "Lehmer rank (rank_perm in exact/bfs.py)",
        "dtype": "uint8, 0xFF = unvisited",
        "core_version": "oracle-1.0",
        "builder": "exact/bfs_fast.c (byte-identical to exact/bfs.py "
                   "for 4<=n<=10)",
        "complete": True,
        "max_distance": max(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        "note": note,
    }
    with open(os.path.join(TABLES, f"dist_n{n}.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
