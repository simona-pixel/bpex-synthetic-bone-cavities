"""Independent continuous triangle-surface clearance audit for two STL meshes.

Run inside Blender 5.x.  The script does not modify the historical generation
pipeline.  It computes the global minimum Euclidean distance between two
triangulated surfaces by:

1. obtaining a conservative upper bound from bidirectional vertex-to-surface
   queries;
2. using epsilon-expanded Blender BVHs to enumerate every triangle pair whose
   axis-aligned boxes can contain a distance no larger than that upper bound;
3. evaluating the exact distance for each candidate triangle pair as the
   minimum of all vertex-to-triangle and edge-to-edge distances.

The output is a JSON audit record suitable for an independent post-generation
clearance check.  Numerical results remain specific to the input triangulation
and Blender floating-point arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import bpy
import numpy as np
from mathutils.bvhtree import BVHTree


TINY = 1.0e-24


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def mul(a, value):
    return (a[0] * value, a[1] * value, a[2] * value)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm_sq(a):
    return dot(a, a)


def point_triangle_closest(point, first, second, third):
    """Return the closest point on a triangle (Ericson region tests)."""

    ab = sub(second, first)
    ac = sub(third, first)
    ap = sub(point, first)
    d1 = dot(ab, ap)
    d2 = dot(ac, ap)
    if d1 <= 0.0 and d2 <= 0.0:
        return first

    bp = sub(point, second)
    d3 = dot(ab, bp)
    d4 = dot(ac, bp)
    if d3 >= 0.0 and d4 <= d3:
        return second

    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        fraction = d1 / (d1 - d3)
        return add(first, mul(ab, fraction))

    cp = sub(point, third)
    d5 = dot(ab, cp)
    d6 = dot(ac, cp)
    if d6 >= 0.0 and d5 <= d6:
        return third

    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        fraction = d2 / (d2 - d6)
        return add(first, mul(ac, fraction))

    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        edge = sub(third, second)
        fraction = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return add(second, mul(edge, fraction))

    denominator = va + vb + vc
    if abs(denominator) <= TINY:
        candidates = (first, second, third)
        return min(candidates, key=lambda item: norm_sq(sub(point, item)))
    inverse = 1.0 / denominator
    v = vb * inverse
    w = vc * inverse
    return add(first, add(mul(ab, v), mul(ac, w)))


def segment_segment_closest(first_a, first_b, second_a, second_b):
    """Return closest points on two closed 3-D segments."""

    first_direction = sub(first_b, first_a)
    second_direction = sub(second_b, second_a)
    offset = sub(first_a, second_a)
    aa = dot(first_direction, first_direction)
    ee = dot(second_direction, second_direction)
    ff = dot(second_direction, offset)

    if aa <= TINY and ee <= TINY:
        return first_a, second_a
    if aa <= TINY:
        first_fraction = 0.0
        second_fraction = min(1.0, max(0.0, ff / ee))
    else:
        cc = dot(first_direction, offset)
        if ee <= TINY:
            second_fraction = 0.0
            first_fraction = min(1.0, max(0.0, -cc / aa))
        else:
            bb = dot(first_direction, second_direction)
            denominator = aa * ee - bb * bb
            if abs(denominator) > TINY:
                first_fraction = min(1.0, max(0.0, (bb * ff - cc * ee) / denominator))
            else:
                first_fraction = 0.0
            second_fraction = (bb * first_fraction + ff) / ee
            if second_fraction < 0.0:
                second_fraction = 0.0
                first_fraction = min(1.0, max(0.0, -cc / aa))
            elif second_fraction > 1.0:
                second_fraction = 1.0
                first_fraction = min(1.0, max(0.0, (bb - cc) / aa))

    first_point = add(first_a, mul(first_direction, first_fraction))
    second_point = add(second_a, mul(second_direction, second_fraction))
    return first_point, second_point


def triangle_triangle_closest(outer_triangle, inner_triangle):
    best_sq = math.inf
    best_outer = None
    best_inner = None
    best_feature = None

    for index, outer_point in enumerate(outer_triangle):
        inner_point = point_triangle_closest(outer_point, *inner_triangle)
        distance_sq = norm_sq(sub(outer_point, inner_point))
        if distance_sq < best_sq:
            best_sq = distance_sq
            best_outer = outer_point
            best_inner = inner_point
            best_feature = f"outer_vertex_{index}-inner_triangle"

    for index, inner_point in enumerate(inner_triangle):
        outer_point = point_triangle_closest(inner_point, *outer_triangle)
        distance_sq = norm_sq(sub(outer_point, inner_point))
        if distance_sq < best_sq:
            best_sq = distance_sq
            best_outer = outer_point
            best_inner = inner_point
            best_feature = f"inner_vertex_{index}-outer_triangle"

    outer_edges = ((0, 1), (1, 2), (2, 0))
    inner_edges = ((0, 1), (1, 2), (2, 0))
    for outer_index, outer_edge in enumerate(outer_edges):
        for inner_index, inner_edge in enumerate(inner_edges):
            outer_point, inner_point = segment_segment_closest(
                outer_triangle[outer_edge[0]],
                outer_triangle[outer_edge[1]],
                inner_triangle[inner_edge[0]],
                inner_triangle[inner_edge[1]],
            )
            distance_sq = norm_sq(sub(outer_point, inner_point))
            if distance_sq < best_sq:
                best_sq = distance_sq
                best_outer = outer_point
                best_inner = inner_point
                best_feature = f"outer_edge_{outer_index}-inner_edge_{inner_index}"

    return best_sq, best_outer, best_inner, best_feature


def reset_scene():
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)


def read_stl(path: Path):
    before = set(bpy.data.objects)
    bpy.ops.wm.stl_import(filepath=str(path))
    imported = [item for item in bpy.data.objects if item not in before and item.type == "MESH"]
    if len(imported) != 1:
        raise RuntimeError(f"Expected one imported mesh for {path}, observed {len(imported)}")
    obj = imported[0]
    matrix = obj.matrix_world
    vertices = [tuple(matrix @ vertex.co) for vertex in obj.data.vertices]
    polygons = []
    for polygon in obj.data.polygons:
        indices = tuple(polygon.vertices)
        if len(indices) != 3:
            raise RuntimeError(f"Non-triangle polygon in STL {path}: {len(indices)} vertices")
        polygons.append(indices)
    bpy.data.objects.remove(obj, do_unlink=True)
    return vertices, polygons


def triangles(vertices, polygons):
    return [
        (vertices[first], vertices[second], vertices[third])
        for first, second, third in polygons
    ]


class AABBNode:
    __slots__ = ("minimum", "maximum", "indices", "left", "right")

    def __init__(self, minimum, maximum, indices=None, left=None, right=None):
        self.minimum = minimum
        self.maximum = maximum
        self.indices = indices
        self.left = left
        self.right = right


def build_aabb_tree(minimums, maximums, centers, indices, leaf_size=24):
    node_minimum = minimums[indices].min(axis=0)
    node_maximum = maximums[indices].max(axis=0)
    if len(indices) <= leaf_size:
        return AABBNode(node_minimum, node_maximum, indices=indices)
    axis = int(np.argmax(node_maximum - node_minimum))
    midpoint = len(indices) // 2
    order = np.argpartition(centers[indices, axis], midpoint)
    ordered = indices[order]
    left = build_aabb_tree(minimums, maximums, centers, ordered[:midpoint], leaf_size)
    right = build_aabb_tree(minimums, maximums, centers, ordered[midpoint:], leaf_size)
    return AABBNode(node_minimum, node_maximum, left=left, right=right)


def aabb_distance_sq(first_minimum, first_maximum, second_minimum, second_maximum):
    gap = np.maximum(0.0, np.maximum(second_minimum - first_maximum, first_minimum - second_maximum))
    return float(np.dot(gap, gap))


def query_aabb_tree(node, query_minimum, query_maximum, radius_sq, minimums, maximums):
    if aabb_distance_sq(query_minimum, query_maximum, node.minimum, node.maximum) > radius_sq:
        return ()
    if node.indices is not None:
        candidate_minimums = minimums[node.indices]
        candidate_maximums = maximums[node.indices]
        gaps = np.maximum(
            0.0,
            np.maximum(candidate_minimums - query_maximum, query_minimum - candidate_maximums),
        )
        distances_sq = np.einsum("ij,ij->i", gaps, gaps)
        return node.indices[distances_sq <= radius_sq + 1.0e-18]
    left = query_aabb_tree(
        node.left, query_minimum, query_maximum, radius_sq, minimums, maximums
    )
    right = query_aabb_tree(
        node.right, query_minimum, query_maximum, radius_sq, minimums, maximums
    )
    if len(left) == 0:
        return right
    if len(right) == 0:
        return left
    return np.concatenate((left, right))


def vertex_surface_upper_bound(source_vertices, target_tree):
    minimum = math.inf
    source_point = None
    target_point = None
    for point in source_vertices:
        hit = target_tree.find_nearest(point)
        if hit[0] is None:
            raise RuntimeError("BVH nearest query unexpectedly returned no hit")
        distance = float(hit[3])
        if distance < minimum:
            minimum = distance
            source_point = point
            target_point = tuple(hit[0])
    return minimum, source_point, target_point


def vector_list(value):
    return [float(component) for component in value]


def audit(outer_path: Path, inner_path: Path, case_id: str, required_clearance: float):
    started = time.perf_counter()
    reset_scene()
    outer_vertices, outer_polygons = read_stl(outer_path)
    inner_vertices, inner_polygons = read_stl(inner_path)
    outer_tree = BVHTree.FromPolygons(outer_vertices, outer_polygons, all_triangles=True)
    inner_tree = BVHTree.FromPolygons(inner_vertices, inner_polygons, all_triangles=True)

    inner_to_outer = vertex_surface_upper_bound(inner_vertices, outer_tree)
    outer_to_inner = vertex_surface_upper_bound(outer_vertices, inner_tree)
    upper_bound = min(inner_to_outer[0], outer_to_inner[0])
    if not math.isfinite(upper_bound):
        raise RuntimeError("Could not establish a finite vertex-to-surface upper bound")

    outer_triangles = triangles(outer_vertices, outer_polygons)
    inner_triangles = triangles(inner_vertices, inner_polygons)
    outer_array = np.asarray(outer_triangles, dtype=np.float64)
    inner_array = np.asarray(inner_triangles, dtype=np.float64)
    outer_minimums = outer_array.min(axis=1)
    outer_maximums = outer_array.max(axis=1)
    inner_minimums = inner_array.min(axis=1)
    inner_maximums = inner_array.max(axis=1)
    inner_centers = 0.5 * (inner_minimums + inner_maximums)
    inner_tree_aabb = build_aabb_tree(
        inner_minimums,
        inner_maximums,
        inner_centers,
        np.arange(len(inner_triangles), dtype=np.int64),
    )

    best_sq = upper_bound * upper_bound
    best = None
    candidate_pair_count = 0
    for outer_index, outer_triangle in enumerate(outer_triangles):
        candidate_indices = query_aabb_tree(
            inner_tree_aabb,
            outer_minimums[outer_index],
            outer_maximums[outer_index],
            best_sq + 1.0e-18,
            inner_minimums,
            inner_maximums,
        )
        candidate_pair_count += len(candidate_indices)
        for inner_index_value in candidate_indices:
            inner_index = int(inner_index_value)
            distance_sq, outer_point, inner_point, feature = triangle_triangle_closest(
                outer_triangle, inner_triangles[inner_index]
            )
            if distance_sq <= best_sq:
                best_sq = distance_sq
                best = (outer_index, inner_index, outer_point, inner_point, feature)

    if best is None:
        # The upper-bound pair can sit exactly on a pruning boundary.  Retain
        # its known closest points so the result remains conservative.
        if inner_to_outer[0] <= outer_to_inner[0]:
            best = (-1, -1, inner_to_outer[2], inner_to_outer[1], "vertex-surface upper bound")
        else:
            best = (-1, -1, outer_to_inner[1], outer_to_inner[2], "vertex-surface upper bound")

    continuous_minimum = math.sqrt(max(0.0, best_sq))
    if continuous_minimum > upper_bound + 1.0e-5:
        raise RuntimeError(
            "Continuous minimum exceeded its vertex-derived upper bound: "
            f"{continuous_minimum:.9f} > {upper_bound:.9f}"
        )

    outer_index, inner_index, outer_point, inner_point, feature = best
    return {
        "schema_version": 1,
        "case_id": case_id,
        "outer_path": str(outer_path.resolve()),
        "inner_path": str(inner_path.resolve()),
        "outer_sha256": sha256(outer_path),
        "inner_sha256": sha256(inner_path),
        "outer_vertex_count": len(outer_vertices),
        "inner_vertex_count": len(inner_vertices),
        "outer_triangle_count": len(outer_polygons),
        "inner_triangle_count": len(inner_polygons),
        "bidirectional_vertex_surface_upper_bound_mm": upper_bound,
        "continuous_triangle_surface_minimum_mm": continuous_minimum,
        "required_clearance_mm": required_clearance,
        "continuous_clearance_pass": continuous_minimum + 1.0e-9 >= required_clearance,
        "candidate_triangle_pair_count": candidate_pair_count,
        "closest_outer_triangle_index": outer_index,
        "closest_inner_triangle_index": inner_index,
        "closest_outer_point_mm": vector_list(outer_point),
        "closest_inner_point_mm": vector_list(inner_point),
        "closest_feature_pair": feature,
        "method": (
            "bidirectional vertex-surface upper bound; hierarchical AABB candidate "
            "enumeration; exact vertex-triangle and edge-edge triangle-pair distance"
        ),
        "blender_version": bpy.app.version_string,
        "python_version": platform.python_version(),
        "elapsed_seconds": time.perf_counter() - started,
    }


def parse_args():
    values = sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else []
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer", required=True, type=Path)
    parser.add_argument("--inner", required=True, type=Path)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--required-clearance", type=float, default=0.15)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args(values)


def main():
    args = parse_args()
    result = audit(args.outer, args.inner, args.case_id, args.required_clearance)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
