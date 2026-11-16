"""Check vertex manifoldness by requiring each STL vertex link to be one cycle."""
import argparse
import collections
import json
from pathlib import Path
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stl', type=Path)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    data = a.stl.read_bytes()
    dtype = np.dtype([('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('attr', '<u2')])
    tri = np.frombuffer(data, dtype=dtype, offset=84)['v']
    vertices, inv = np.unique(tri.reshape(-1, 3), axis=0, return_inverse=True)
    faces = inv.reshape(-1, 3)
    links = [collections.defaultdict(list) for _ in vertices]
    for x, y, z in faces:
        for v, b, c in ((x, y, z), (y, z, x), (z, x, y)):
            links[v][b].append(c)
            links[v][c].append(b)
    bad = []
    for i, link in enumerate(links):
        if not link or any(len(n) != 2 for n in link.values()):
            bad.append(i)
            continue
        seen, todo = set(), [next(iter(link))]
        while todo:
            v = todo.pop()
            if v not in seen:
                seen.add(v)
                todo.extend(link[v])
        if len(seen) != len(link):
            bad.append(i)
    result = {'checked_vertices': len(vertices), 'non_manifold_vertex_links': len(bad),
              'bad_vertex_indices': bad,
              'criterion': 'Every vertex link is one connected cycle with degree two at every link vertex'}
    a.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'bad_vertex_indices'}))
    return bool(bad)


if __name__ == '__main__':
    raise SystemExit(main())
