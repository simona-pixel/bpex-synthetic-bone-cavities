"""Configurable Blender pipeline for the medullary-canal benchmark (v3).

This module consolidates the outer/inner algorithm that was duplicated in the
three frozen July V4 scripts.  It deliberately excludes the independent physis
workflow.  Version 3 replaces nearest-triangle-normal containment signs with a
three-direction ray-parity classification while retaining unsigned clearance.
Run it with Blender, not with the system Python interpreter.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import struct
import sys
import time
import traceback

import bmesh
import bpy
import mathutils
from mathutils.bvhtree import BVHTree


METHODS = ("regional", "uniform_normal", "solidify")


def arguments():
    argv = sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else []
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--bone", required=True)
    parser.add_argument("--method", required=True, choices=METHODS)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--thickness", type=float)
    parser.add_argument(
        "--regional-postprocess",
        choices=("historical", "decimate_only", "smooth_only", "none"),
        default="historical",
    )
    parser.add_argument(
        "--regional-kernel",
        choices=("legacy_deform", "validated_loft"),
        default="legacy_deform",
    )
    parser.add_argument("--allow-invalid-source", action="store_true")
    parser.add_argument("--allow-invalid-output", action="store_true")
    parser.add_argument("--containment-guard", action="store_true")
    return parser.parse_args(argv)


def clamp(value, lower, upper):
    return max(lower, min(upper, value))


def smoothstep01(value):
    value = clamp(value, 0.0, 1.0)
    return value * value * (3.0 - 2.0 * value)


def smootherstep01(value):
    value = clamp(value, 0.0, 1.0)
    return value * value * value * (value * (value * 6.0 - 15.0) + 10.0)


def lerp(first, second, fraction):
    return first + (second - first) * fraction


def select_only(obj):
    bpy.ops.object.select_all(action="DESELECT")
    obj.hide_set(False)
    obj.hide_viewport = False
    obj.hide_render = False
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj


def clear_scene():
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)


def duplicate_object(obj, name):
    duplicate = obj.copy()
    duplicate.data = obj.data.copy()
    duplicate.name = name
    duplicate.data.name = name + "_mesh"
    bpy.context.collection.objects.link(duplicate)
    return duplicate


def require_file(path):
    if not os.path.isfile(path):
        raise RuntimeError(f"Required file not found: {path}")


def import_stl(path, name):
    require_file(path)
    before = set(bpy.data.objects)
    try:
        result = bpy.ops.import_mesh.stl(filepath=path)
    except Exception:
        result = bpy.ops.wm.stl_import(filepath=path)
    if "FINISHED" not in result:
        raise RuntimeError(f"STL import failed: {result}")
    created = [obj for obj in set(bpy.data.objects) - before if obj.type == "MESH"]
    if len(created) != 1:
        raise RuntimeError(f"Expected one imported mesh, received {len(created)}: {path}")
    obj = created[0]
    obj.name = name
    obj.data.name = name + "_mesh"
    select_only(obj)
    bpy.ops.object.transform_apply(location=False, rotation=True, scale=True)
    return obj


def export_stl(obj, path):
    path = os.path.abspath(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    partial = path + ".partial.stl"
    if os.path.exists(partial):
        os.remove(partial)
    select_only(obj)
    try:
        result = bpy.ops.export_mesh.stl(filepath=partial, use_selection=True)
    except Exception:
        result = bpy.ops.wm.stl_export(
            filepath=partial,
            export_selected_objects=True,
            ascii_format=False,
        )
    if "FINISHED" not in result or not os.path.isfile(partial) or os.path.getsize(partial) <= 0:
        raise RuntimeError(f"STL export failed: {path}")
    os.replace(partial, path)


def recalc_normals_outside(obj):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bmesh.ops.recalc_face_normals(bm, faces=list(bm.faces))
    bm.normal_update()
    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()


def merge_by_distance(obj, distance):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bmesh.ops.remove_doubles(bm, verts=list(bm.verts), dist=distance)
    bm.normal_update()
    bm.to_mesh(obj.data)
    bm.free()
    obj.data.update()


def winding_conflicts(bm):
    conflicts = 0
    for edge in bm.edges:
        loops = list(edge.link_loops)
        if len(loops) == 2:
            starts = [loop.vert is edge.verts[0] for loop in loops]
            if starts[0] == starts[1]:
                conflicts += 1
    return conflicts


def connected_components(bm):
    unseen = set(bm.verts)
    count = 0
    while unseen:
        count += 1
        stack = [unseen.pop()]
        while stack:
            vertex = stack.pop()
            for edge in vertex.link_edges:
                other = edge.other_vert(vertex)
                if other in unseen:
                    unseen.remove(other)
                    stack.append(other)
    return count


def mesh_diagnostics(obj):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bm.normal_update()
    boundary = sum(1 for edge in bm.edges if len(edge.link_faces) == 1)
    non_manifold = sum(1 for edge in bm.edges if len(edge.link_faces) != 2)
    loose_vertices = sum(1 for vertex in bm.verts if not vertex.link_edges)
    loose_edges = sum(1 for edge in bm.edges if not edge.link_faces)
    degenerate_edges = sum(1 for edge in bm.edges if edge.calc_length() <= 1.0e-12)
    degenerate_faces = sum(1 for face in bm.faces if face.calc_area() <= 1.0e-12)
    components = connected_components(bm) if bm.verts else 0
    conflicts = winding_conflicts(bm)
    closed = bool(bm.faces) and non_manifold == 0 and loose_vertices == 0 and loose_edges == 0
    signed_volume = float(bm.calc_volume(signed=True)) if closed else math.nan
    result = {
        "vertices": len(bm.verts),
        "edges": len(bm.edges),
        "faces": len(bm.faces),
        "components": components,
        "boundary_edges": boundary,
        "non_manifold_edges": non_manifold,
        "loose_vertices": loose_vertices,
        "loose_edges": loose_edges,
        "winding_conflicts": conflicts,
        "degenerate_edges": degenerate_edges,
        "degenerate_faces": degenerate_faces,
        "closed": closed,
        "signed_volume_mm3": signed_volume,
        "volume_mm3": abs(signed_volume) if closed else math.nan,
    }
    bm.free()
    return result


def diagnostic_failures(diagnostics):
    failures = []
    for key in (
        "boundary_edges",
        "non_manifold_edges",
        "loose_vertices",
        "loose_edges",
        "winding_conflicts",
        "degenerate_edges",
        "degenerate_faces",
    ):
        if diagnostics[key] != 0:
            failures.append(f"{key}={diagnostics[key]}")
    if diagnostics["components"] != 1:
        failures.append(f"components={diagnostics['components']}")
    if not diagnostics["closed"]:
        failures.append("closed=false")
    if not math.isfinite(diagnostics["volume_mm3"]) or diagnostics["volume_mm3"] <= 1.0e-9:
        failures.append(f"volume_mm3={diagnostics['volume_mm3']}")
    return failures


def mesh_signature(obj):
    digest = hashlib.sha256()
    mesh = obj.data
    digest.update(struct.pack("<III", len(mesh.vertices), len(mesh.edges), len(mesh.polygons)))
    digest.update(struct.pack("<" + "d" * 16, *(float(value) for row in obj.matrix_world for value in row)))
    for vertex in mesh.vertices:
        digest.update(struct.pack("<ddd", *map(float, vertex.co)))
    for edge in mesh.edges:
        digest.update(struct.pack("<II", *(int(index) for index in edge.vertices)))
    for polygon in mesh.polygons:
        indices = tuple(int(index) for index in polygon.vertices)
        digest.update(struct.pack("<I", len(indices)))
        if indices:
            digest.update(struct.pack("<" + "I" * len(indices), *indices))
    return digest.hexdigest()


def safe_inner_distance(outer_distance, target_distance):
    if outer_distance <= 0.10:
        return 0.05
    return clamp(target_distance, 0.05, outer_distance - 0.05)


def canal_ratio_profile(u, params):
    distance = abs(u - 0.5) / 0.5
    if distance < params["half_isthme"]:
        ratio = params["canal_ratio_isthme"]
    elif distance < params["half_diaphyse"]:
        fraction = (distance - params["half_isthme"]) / (
            params["half_diaphyse"] - params["half_isthme"]
        )
        ratio = lerp(
            params["canal_ratio_isthme"],
            params["canal_ratio_diaphyse"],
            smoothstep01(fraction),
        )
    else:
        fraction = (distance - params["half_diaphyse"]) / (1.0 - params["half_diaphyse"])
        ratio = lerp(
            params["canal_ratio_diaphyse"],
            params["canal_ratio_metaphyse"],
            smoothstep01(fraction),
        )
    return clamp(ratio, 0.05, 0.98)


def presence_from_one_end(distance_mm, end_zone_mm, flare_mm):
    if distance_mm <= end_zone_mm:
        return 0.0
    if distance_mm < end_zone_mm + flare_mm:
        return smootherstep01((distance_mm - end_zone_mm) / flare_mm)
    return 1.0


def shell_floor_from_end(distance_mm, shell_epi_mm, shell_meta_mm, plateau_mm, epi_fade_mm, zero_fade_mm):
    if distance_mm <= plateau_mm:
        return shell_epi_mm
    first_end = plateau_mm + epi_fade_mm
    if distance_mm <= first_end:
        fraction = (distance_mm - plateau_mm) / max(epi_fade_mm, 1.0e-9)
        return lerp(shell_epi_mm, shell_meta_mm, smootherstep01(fraction))
    second_end = first_end + zero_fade_mm
    if distance_mm <= second_end:
        fraction = (distance_mm - first_end) / max(zero_fade_mm, 1.0e-9)
        return lerp(shell_meta_mm, 0.0, smootherstep01(fraction))
    return 0.0


def build_spine(obj, params):
    vertices = obj.data.vertices
    z_values = [vertex.co.z for vertex in vertices]
    if not z_values:
        raise RuntimeError("Source STL contains no vertices")
    z_min, z_max = min(z_values), max(z_values)
    length = z_max - z_min
    if length <= 1.0e-6:
        raise RuntimeError("Longitudinal Z extent is nearly zero")
    section_count = int(params["nb_sections"])
    z_step = length / (section_count - 1)
    z_window = length * params["z_window_ratio"]
    raw = []
    for section in range(section_count):
        z_target = z_min + section * z_step
        sample = [vertex.co for vertex in vertices if abs(vertex.co.z - z_target) <= z_window]
        if sample:
            center_x = statistics.median(co.x for co in sample)
            center_y = statistics.median(co.y for co in sample)
        elif raw:
            center_x, center_y = raw[-1].x, raw[-1].y
        else:
            center_x = center_y = 0.0
        raw.append(mathutils.Vector((center_x, center_y, z_target)))
    half_window = int(params["spine_smooth_window"])
    spine = []
    for index in range(section_count):
        chunk = raw[max(0, index - half_window) : min(section_count, index + half_window + 1)]
        spine.append(
            mathutils.Vector(
                (
                    sum(point.x for point in chunk) / len(chunk),
                    sum(point.y for point in chunk) / len(chunk),
                    raw[index].z,
                )
            )
        )
    return spine, z_min, z_max, length, z_step


def _ray_intersection_count(tree, origin, direction, maximum_distance, epsilon=1.0e-3):
    count = 0
    travelled = 0.0
    current = origin.copy()
    while travelled < maximum_distance:
        hit = tree.ray_cast(current, direction, maximum_distance - travelled)
        if not hit or hit[0] is None:
            break
        distance = float(hit[3])
        step = distance + epsilon if distance > epsilon else 10.0 * epsilon
        current = current + direction * step
        travelled += step
        count += 1
        if count > 1000:
            raise RuntimeError("Unexpected ray-intersection count")
    return count


def point_inside_votes(tree, point, maximum_distance):
    directions = (
        mathutils.Vector((1.0, 0.173, 0.071)).normalized(),
        mathutils.Vector((-0.119, 1.0, 0.227)).normalized(),
        mathutils.Vector((0.193, -0.137, 1.0)).normalized(),
    )
    return sum(
        _ray_intersection_count(tree, point, direction, maximum_distance) % 2 == 1
        for direction in directions
    )


def point_inside(tree, point, maximum_distance):
    return point_inside_votes(tree, point, maximum_distance) >= 2


def build_interior_spine(
    obj,
    params,
    section_count,
    grid_size,
    z_margin_mm,
    z_margin_fraction=0.0,
):
    """Build a longitudinal line whose sampled centers are inside the outer volume."""
    vertices = [obj.matrix_world @ vertex.co for vertex in obj.data.vertices]
    z_min = min(point.z for point in vertices)
    z_max = max(point.z for point in vertices)
    length = z_max - z_min
    if length <= 1.0e-6:
        raise RuntimeError("Longitudinal Z extent is nearly zero")
    margin = max(float(z_margin_mm), float(z_margin_fraction) * length)
    if margin <= 0.0 or 2.0 * margin >= length:
        raise RuntimeError("Invalid validated-loft Z margin")
    z_start, z_end = z_min + margin, z_max - margin
    z_step = (z_end - z_start) / (section_count - 1)
    local_window = max(length * params["z_window_ratio"], 1.5 * z_step)
    outer_bm, tree = bvh_for(obj)
    minimum = [min(point[index] for point in vertices) for index in range(3)]
    maximum = [max(point[index] for point in vertices) for index in range(3)]
    maximum_distance = math.sqrt(sum((maximum[index] - minimum[index]) ** 2 for index in range(3))) * 4.0
    raw = []
    for section in range(section_count):
        z = z_start + section * z_step
        local = [point for point in vertices if abs(point.z - z) <= local_window]
        if not local:
            outer_bm.free()
            raise RuntimeError(f"No outer vertices near validated spine section {section}")
        xs = [point.x for point in local]
        ys = [point.y for point in local]
        x_min, x_max = min(xs), max(xs)
        y_min, y_max = min(ys), max(ys)
        candidates = []
        if raw:
            candidates.append(mathutils.Vector((raw[-1].x, raw[-1].y, z)))
        candidates.extend(
            [
                mathutils.Vector((statistics.median(xs), statistics.median(ys), z)),
                mathutils.Vector((statistics.fmean(xs), statistics.fmean(ys), z)),
            ]
        )
        for x_index in range(grid_size):
            x = x_min + (x_index + 0.5) / grid_size * (x_max - x_min)
            for y_index in range(grid_size):
                y = y_min + (y_index + 0.5) / grid_size * (y_max - y_min)
                candidates.append(mathutils.Vector((x, y, z)))
        best = None
        best_score = -math.inf
        previous = raw[-1] if raw else None
        for candidate in candidates:
            if not point_inside(tree, candidate, maximum_distance):
                continue
            hit = tree.find_nearest(candidate)
            if not hit or hit[0] is None:
                continue
            depth = float(hit[3])
            continuity_penalty = 0.02 * float((candidate - previous).length) if previous else 0.0
            score = depth - continuity_penalty
            if score > best_score:
                best = candidate
                best_score = score
        if best is None:
            outer_bm.free()
            raise RuntimeError(f"Could not find an interior center at validated spine section {section}")
        raw.append(best)
    half_window = max(1, int(params["spine_smooth_window"]))
    spine = []
    for index, raw_point in enumerate(raw):
        chunk = raw[max(0, index - half_window) : min(len(raw), index + half_window + 1)]
        smooth = mathutils.Vector(
            (
                sum(point.x for point in chunk) / len(chunk),
                sum(point.y for point in chunk) / len(chunk),
                raw_point.z,
            )
        )
        chosen = raw_point
        for fraction in (1.0, 0.75, 0.50, 0.25):
            candidate = raw_point.lerp(smooth, fraction)
            candidate.z = raw_point.z
            if point_inside(tree, candidate, maximum_distance):
                chosen = candidate
                break
        spine.append(chosen)
    outer_bm.free()
    return spine, z_min, z_max, length, z_start, z_end, z_step


def required_shell_mm(z, z_min, z_max, params):
    proximal = shell_floor_from_end(
        z_max - z,
        params["prox_epi_shell_min_mm"],
        params["prox_meta_shell_min_mm"],
        params["prox_epi_plateau_mm"],
        params["prox_epi_to_meta_fade_mm"],
        params["prox_meta_to_zero_fade_mm"],
    )
    distal = shell_floor_from_end(
        z - z_min,
        params["dist_epi_shell_min_mm"],
        params["dist_meta_shell_min_mm"],
        params["dist_epi_plateau_mm"],
        params["dist_epi_to_meta_fade_mm"],
        params["dist_meta_to_zero_fade_mm"],
    )
    return max(proximal, distal)


def expected_inner_radius(outer_radius, z, z_min, z_max, params):
    length = z_max - z_min
    u = (z - z_min) / length
    distal_presence = presence_from_one_end(
        z - z_min, params["end_zone_dist_mm"], params["flare_dist_mm"]
    )
    proximal_presence = presence_from_one_end(
        z_max - z, params["end_zone_prox_mm"], params["flare_prox_mm"]
    )
    presence = min(distal_presence, proximal_presence)
    terminal = max(0.05, outer_radius - params["terminal_backing_offset_mm"])
    canal = max(params["rmin_absolute_mm"], outer_radius * canal_ratio_profile(u, params))
    canal = safe_inner_distance(outer_radius, canal)
    target = lerp(terminal, canal, presence)
    shell_floor = required_shell_mm(z, z_min, z_max, params)
    if shell_floor > 0.0:
        target = min(target, safe_inner_distance(outer_radius, outer_radius - shell_floor))
    return safe_inner_distance(outer_radius, target)


def generate_regional(source, bone_cfg, postprocess="historical"):
    params = bone_cfg["parameters"]
    outer = duplicate_object(source, Path(bone_cfg["outer_filename"]).stem)
    inner = duplicate_object(source, Path(bone_cfg["inner_filename"]).stem)
    recalc_normals_outside(outer)
    recalc_normals_outside(inner)
    spine, z_min, z_max, _, z_step = build_spine(inner, params)
    bm = bmesh.new()
    bm.from_mesh(inner.data)
    for vertex in bm.verts:
        index = int((vertex.co.z - z_min) / z_step)
        index = max(0, min(len(spine) - 1, index))
        center = spine[index]
        radial = mathutils.Vector((vertex.co.x - center.x, vertex.co.y - center.y, 0.0))
        outer_radius = radial.length
        if outer_radius < 1.0e-9:
            continue
        target_radius = expected_inner_radius(outer_radius, vertex.co.z, z_min, z_max, params)
        vertex.co.x = center.x + radial.x / outer_radius * target_radius
        vertex.co.y = center.y + radial.y / outer_radius * target_radius
    bm.to_mesh(inner.data)
    bm.free()
    inner.data.update()
    merge_by_distance(inner, params["merge_distance_mm"])
    recalc_normals_outside(inner)
    select_only(inner)
    if params["use_decimate"] and postprocess in {"historical", "decimate_only"}:
        modifier = inner.modifiers.new(name="Decimate", type="DECIMATE")
        modifier.ratio = params["decimate_ratio"]
        bpy.ops.object.modifier_apply(modifier=modifier.name)
    if postprocess in {"historical", "smooth_only"}:
        modifier = inner.modifiers.new(name="Smooth", type="SMOOTH")
        modifier.iterations = int(params["smooth_iterations"])
        bpy.ops.object.modifier_apply(modifier=modifier.name)
    merge_by_distance(inner, params["merge_distance_mm"])
    recalc_normals_outside(inner)
    return outer, inner


def generate_regional_loft(source, bone_cfg, benchmark_cfg):
    params = bone_cfg["parameters"]
    outer = duplicate_object(source, Path(bone_cfg["outer_filename"]).stem)
    recalc_normals_outside(outer)
    section_count = int(benchmark_cfg["validated_loft_sections"])
    angle_count = int(benchmark_cfg["validated_loft_angles"])
    grid_size = int(benchmark_cfg["validated_spine_grid_size"])
    z_margin = float(benchmark_cfg["validated_loft_z_margin_mm"])
    z_margin_fraction = float(benchmark_cfg["validated_loft_z_margin_fraction"])
    safety = float(benchmark_cfg["validated_loft_radial_safety_mm"])
    spine, z_min, z_max, _, _, _, _ = build_interior_spine(
        outer,
        params,
        section_count,
        grid_size,
        z_margin,
        z_margin_fraction,
    )
    outer_bm, outer_tree = bvh_for(outer)
    outer_bounds = [outer.matrix_world @ vertex.co for vertex in outer.data.vertices]
    maximum_distance = math.sqrt(
        sum(
            (max(point[index] for point in outer_bounds) - min(point[index] for point in outer_bounds)) ** 2
            for index in range(3)
        )
    ) * 2.0
    vertices = []
    for center in spine:
        for angle_index in range(angle_count):
            angle = 2.0 * math.pi * angle_index / angle_count
            direction = mathutils.Vector((math.cos(angle), math.sin(angle), 0.0))
            hit = outer_tree.ray_cast(center, direction, maximum_distance)
            if not hit or hit[0] is None:
                outer_bm.free()
                raise RuntimeError(
                    f"Outer ray miss at z={center.z:.6f}, angle={math.degrees(angle):.1f} degrees"
                )
            outer_radius = float(hit[3])
            target = expected_inner_radius(outer_radius, center.z, z_min, z_max, params)
            target = max(0.05, target - safety)
            vertices.append(
                (
                    center.x + direction.x * target,
                    center.y + direction.y * target,
                    center.z,
                )
            )
    distal_center_index = len(vertices)
    vertices.append(tuple(spine[0]))
    proximal_center_index = len(vertices)
    vertices.append(tuple(spine[-1]))
    faces = []
    for section in range(section_count - 1):
        lower = section * angle_count
        upper = (section + 1) * angle_count
        for angle_index in range(angle_count):
            next_angle = (angle_index + 1) % angle_count
            faces.append(
                (
                    lower + angle_index,
                    lower + next_angle,
                    upper + next_angle,
                    upper + angle_index,
                )
            )
    for angle_index in range(angle_count):
        next_angle = (angle_index + 1) % angle_count
        faces.append((distal_center_index, next_angle, angle_index))
        base = (section_count - 1) * angle_count
        faces.append((proximal_center_index, base + angle_index, base + next_angle))
    mesh = bpy.data.meshes.new(Path(bone_cfg["inner_filename"]).stem + "_mesh")
    mesh.from_pydata(vertices, [], faces)
    mesh.update()
    inner = bpy.data.objects.new(Path(bone_cfg["inner_filename"]).stem, mesh)
    bpy.context.collection.objects.link(inner)
    triangulated = bmesh.new()
    triangulated.from_mesh(inner.data)
    bmesh.ops.triangulate(triangulated, faces=list(triangulated.faces))
    triangulated.to_mesh(inner.data)
    triangulated.free()
    inner.data.update()
    recalc_normals_outside(inner)
    outer_bm.free()
    return outer, inner, spine


def generate_uniform_normal(source, thickness, bone):
    inner = duplicate_object(source, f"{bone}_uniform_normal_inner")
    recalc_normals_outside(inner)
    bm = bmesh.new()
    bm.from_mesh(inner.data)
    bm.normal_update()
    for vertex in bm.verts:
        vertex.co -= vertex.normal.normalized() * thickness
    bm.to_mesh(inner.data)
    bm.free()
    inner.data.update()
    recalc_normals_outside(inner)
    return inner


def bvh_for(obj):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    bm.transform(obj.matrix_world)
    bm.normal_update()
    return bm, BVHTree.FromBMesh(bm)


def polygon_tree(obj, epsilon=0.0):
    vertices = [obj.matrix_world @ vertex.co for vertex in obj.data.vertices]
    polygons = [tuple(int(index) for index in polygon.vertices) for polygon in obj.data.polygons]
    if any(len(polygon) != 3 for polygon in polygons):
        raise RuntimeError(f"Expected triangulated STL mesh: {obj.name}")
    return vertices, polygons, BVHTree.FromPolygons(
        vertices, polygons, all_triangles=True, epsilon=epsilon
    )


def self_intersection_count(obj):
    _, polygons, tree = polygon_tree(obj)
    pairs = set()
    for first, second in tree.overlap(tree):
        if first == second:
            continue
        if first > second:
            first, second = second, first
        if set(polygons[first]).intersection(polygons[second]):
            continue
        pairs.add((first, second))
    return len(pairs)


def relation_diagnostics(outer, inner, tolerance=1.0e-4, required_clearance=0.0):
    outer_vertices, _, outer_tree = polygon_tree(outer)
    _, _, inner_tree = polygon_tree(inner)
    overlap_pairs = outer_tree.overlap(inner_tree)
    outside = []
    clearance_violations = []
    minimum_signed = math.inf
    maximum_signed = -math.inf
    minimum_unsigned = math.inf
    maximum_unsigned = -math.inf
    parity_vote_histogram = {0: 0, 1: 0, 2: 0, 3: 0}
    minimum = [min(point[index] for point in outer_vertices) for index in range(3)]
    maximum = [max(point[index] for point in outer_vertices) for index in range(3)]
    maximum_distance = math.sqrt(
        sum((maximum[index] - minimum[index]) ** 2 for index in range(3))
    ) * 4.0
    for vertex in inner.data.vertices:
        point = inner.matrix_world @ vertex.co
        hit = outer_tree.find_nearest(point)
        if not hit or hit[0] is None:
            continue
        location, normal, _, _ = hit
        signed = float((point - location).dot(normal))
        unsigned = float((point - location).length)
        minimum_signed = min(minimum_signed, signed)
        maximum_signed = max(maximum_signed, signed)
        minimum_unsigned = min(minimum_unsigned, unsigned)
        maximum_unsigned = max(maximum_unsigned, unsigned)
        inside_votes = point_inside_votes(outer_tree, point, maximum_distance)
        parity_vote_histogram[inside_votes] += 1
        is_inside = inside_votes >= 2
        if not is_inside:
            outside.append(unsigned)
        if not is_inside or unsigned < required_clearance - tolerance:
            clearance_violations.append(unsigned)
    self_pairs = self_intersection_count(inner)
    return {
        "containment_classifier": "majority vote of three deterministic ray-parity tests",
        "outer_inner_intersection_pairs": len(overlap_pairs),
        "outside_vertices_count": len(outside),
        "outside_max_mm": max(outside) if outside else 0.0,
        "required_clearance_mm": required_clearance,
        "clearance_violation_count": len(clearance_violations),
        "nearest_unsigned_min_mm": minimum_unsigned,
        "nearest_unsigned_max_mm": maximum_unsigned,
        "nearest_signed_min_mm": minimum_signed,
        "nearest_signed_max_mm": maximum_signed,
        "nearest_signed_role": "diagnostic only; not used for containment or clearance decisions",
        "parity_vote_histogram": parity_vote_histogram,
        "parity_disagreement_vertices_count": (
            parity_vote_histogram[1] + parity_vote_histogram[2]
        ),
        "self_intersection_pairs": self_pairs,
        "geometrically_valid": len(overlap_pairs) == 0 and not outside and self_pairs == 0,
        "clearance_protected_valid": (
            len(overlap_pairs) == 0
            and not outside
            and not clearance_violations
            and self_pairs == 0
        ),
    }


def enforce_containment(
    outer,
    inner,
    params,
    minimum_clearance,
    minimum_scale,
    search_step,
    binary_iterations,
):
    """Find the least global radial contraction that guarantees confinement."""
    original = [vertex.co.copy() for vertex in inner.data.vertices]
    spine, z_min, _, _, z_step = build_spine(outer, params)

    def apply_scale(scale):
        for index, vertex in enumerate(inner.data.vertices):
            point = original[index]
            spine_index = int((point.z - z_min) / z_step)
            spine_index = max(0, min(len(spine) - 1, spine_index))
            center = spine[spine_index]
            vertex.co.x = center.x + (point.x - center.x) * scale
            vertex.co.y = center.y + (point.y - center.y) * scale
            vertex.co.z = point.z
        inner.data.update()
        recalc_normals_outside(inner)
        relation = relation_diagnostics(
            outer,
            inner,
            required_clearance=minimum_clearance,
        )
        return relation

    evaluations = []
    upper = 1.0
    upper_relation = apply_scale(upper)
    evaluations.append({"radial_scale": upper, "relation": upper_relation})
    if upper_relation["clearance_protected_valid"]:
        selected_scale = upper
        selected_relation = upper_relation
    else:
        lower = None
        trial = 1.0 - search_step
        while trial >= minimum_scale - 1.0e-12:
            relation = apply_scale(trial)
            evaluations.append({"radial_scale": trial, "relation": relation})
            if relation["clearance_protected_valid"]:
                lower = trial
                selected_relation = relation
                break
            upper = trial
            trial -= search_step
        if lower is None:
            raise RuntimeError(
                "Containment guard could not find a valid radial scale: "
                + json.dumps(json_safe(evaluations[-5:]), sort_keys=True)
            )
        for _ in range(binary_iterations):
            middle = (lower + upper) / 2.0
            relation = apply_scale(middle)
            evaluations.append({"radial_scale": middle, "relation": relation})
            if relation["clearance_protected_valid"]:
                lower = middle
                selected_relation = relation
            else:
                upper = middle
        selected_scale = lower
        selected_relation = apply_scale(selected_scale)
    displacements = [
        float((vertex.co - original[index]).length)
        for index, vertex in enumerate(inner.data.vertices)
    ]
    after_diagnostics = mesh_diagnostics(inner)
    report = {
        "enabled": True,
        "strategy": "minimal global radial contraction about the frozen spine",
        "minimum_clearance_mm": minimum_clearance,
        "minimum_allowed_radial_scale": minimum_scale,
        "selected_radial_scale": selected_scale,
        "radial_reduction_percent": (1.0 - selected_scale) * 100.0,
        "surface_displacement_max_mm": max(displacements) if displacements else 0.0,
        "surface_displacement_mean_mm": statistics.fmean(displacements) if displacements else 0.0,
        "evaluations": evaluations,
        "after_topology": after_diagnostics,
        "final_relation": selected_relation,
    }
    if diagnostic_failures(after_diagnostics) or not selected_relation["clearance_protected_valid"]:
        raise RuntimeError(
            "Containment guard did not produce a valid inner surface: "
            + json.dumps(json_safe(report), sort_keys=True)
        )
    return report


def enforce_loft_containment(
    outer,
    inner,
    spine,
    angle_count,
    minimum_clearance,
    minimum_scale,
    search_step,
    binary_iterations,
):
    """Contract validated-loft rings about their proven interior centers.

    The loft generator creates ``len(spine)`` ordered rings followed by two cap
    centers.  Contracting those rings about the same frozen centers preserves
    the requested angular profile while providing a deterministic, fail-closed
    containment guard for interpolation between adjacent sections.
    """
    ring_vertex_count = len(spine) * int(angle_count)
    if len(inner.data.vertices) != ring_vertex_count + 2:
        raise RuntimeError(
            "Unexpected validated-loft vertex layout: "
            f"{len(inner.data.vertices)} vertices for {len(spine)} rings x {angle_count} angles"
        )
    original = [vertex.co.copy() for vertex in inner.data.vertices]

    def apply_scale(scale):
        for section_index, center in enumerate(spine):
            start = section_index * angle_count
            for angle_index in range(angle_count):
                vertex_index = start + angle_index
                point = original[vertex_index]
                vertex = inner.data.vertices[vertex_index]
                vertex.co.x = center.x + (point.x - center.x) * scale
                vertex.co.y = center.y + (point.y - center.y) * scale
                vertex.co.z = point.z
        for vertex_index in range(ring_vertex_count, len(inner.data.vertices)):
            inner.data.vertices[vertex_index].co = original[vertex_index]
        inner.data.update()
        recalc_normals_outside(inner)
        return relation_diagnostics(
            outer,
            inner,
            required_clearance=minimum_clearance,
        )

    evaluations = []
    upper = 1.0
    upper_relation = apply_scale(upper)
    evaluations.append({"radial_scale": upper, "relation": upper_relation})
    if upper_relation["clearance_protected_valid"]:
        selected_scale = upper
        selected_relation = upper_relation
    else:
        lower = None
        trial = 1.0 - search_step
        while trial >= minimum_scale - 1.0e-12:
            relation = apply_scale(trial)
            evaluations.append({"radial_scale": trial, "relation": relation})
            if relation["clearance_protected_valid"]:
                lower = trial
                selected_relation = relation
                break
            upper = trial
            trial -= search_step
        if lower is None:
            raise RuntimeError(
                "Validated-loft containment guard could not find a valid radial scale: "
                + json.dumps(json_safe(evaluations[-5:]), sort_keys=True)
            )
        for _ in range(binary_iterations):
            middle = (lower + upper) / 2.0
            relation = apply_scale(middle)
            evaluations.append({"radial_scale": middle, "relation": relation})
            if relation["clearance_protected_valid"]:
                lower = middle
                selected_relation = relation
            else:
                upper = middle
        selected_scale = lower
        selected_relation = apply_scale(selected_scale)

    displacements = [
        float((vertex.co - original[index]).length)
        for index, vertex in enumerate(inner.data.vertices)
    ]
    after_diagnostics = mesh_diagnostics(inner)
    report = {
        "enabled": True,
        "strategy": "minimal global radial contraction about validated-loft interior centers",
        "minimum_clearance_mm": minimum_clearance,
        "minimum_allowed_radial_scale": minimum_scale,
        "selected_radial_scale": selected_scale,
        "radial_reduction_percent": (1.0 - selected_scale) * 100.0,
        "surface_displacement_max_mm": max(displacements) if displacements else 0.0,
        "surface_displacement_mean_mm": statistics.fmean(displacements) if displacements else 0.0,
        "evaluations": evaluations,
        "after_topology": after_diagnostics,
        "final_relation": selected_relation,
    }
    if diagnostic_failures(after_diagnostics) or not selected_relation["clearance_protected_valid"]:
        raise RuntimeError(
            "Validated-loft containment guard did not produce a valid inner surface: "
            + json.dumps(json_safe(report), sort_keys=True)
        )
    return report


def mean_nearest_distance(obj, tree):
    values = []
    stride = max(1, len(obj.data.vertices) // 5000)
    for index in range(0, len(obj.data.vertices), stride):
        hit = tree.find_nearest(obj.matrix_world @ obj.data.vertices[index].co)
        if hit and hit[0] is not None:
            values.append(float(hit[3]))
    return statistics.fmean(values) if values else math.inf


def generate_solidify(source, thickness, bone, method_cfg):
    source_bm, source_tree = bvh_for(source)
    shell = duplicate_object(source, f"{bone}_solidify_shell")
    select_only(shell)
    modifier = shell.modifiers.new(name="Solidify", type="SOLIDIFY")
    modifier.thickness = thickness
    modifier.solidify_mode = str(method_cfg["solidify_mode"])
    modifier.offset = float(method_cfg["offset"])
    modifier.use_rim = bool(method_cfg["use_rim"])
    modifier.use_even_offset = bool(method_cfg["use_even_offset"])
    modifier.use_quality_normals = bool(method_cfg["use_quality_normals"])
    modifier.thickness_clamp = float(method_cfg["thickness_clamp"])
    modifier.use_thickness_angle_clamp = bool(method_cfg["use_thickness_angle_clamp"])
    bpy.ops.object.modifier_apply(modifier=modifier.name)
    select_only(shell)
    bpy.ops.object.mode_set(mode="EDIT")
    bpy.ops.mesh.select_all(action="SELECT")
    bpy.ops.mesh.separate(type="LOOSE")
    bpy.ops.object.mode_set(mode="OBJECT")
    components = [obj for obj in bpy.context.selected_objects if obj.type == "MESH"]
    if len(components) != 2:
        source_bm.free()
        raise RuntimeError(f"Solidify yielded {len(components)} loose component(s), expected exactly two")
    ranked = sorted(
        ((mean_nearest_distance(component, source_tree), component) for component in components),
        key=lambda item: item[0],
    )
    inner_components = [component for distance, component in ranked if distance > 0.5 * thickness]
    outer_components = [component for distance, component in ranked if distance <= 0.5 * thickness]
    if len(inner_components) != 1 or len(outer_components) != 1:
        source_bm.free()
        raise RuntimeError(
            "Could not classify Solidify components: "
            + ", ".join(f"{distance:.6f}" for distance, _ in ranked)
        )
    inner = inner_components[0]
    inner.name = f"{bone}_solidify_inner"
    for component in outer_components:
        bpy.data.objects.remove(component, do_unlink=True)
    source_bm.free()
    recalc_normals_outside(inner)
    return inner


def json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def run(args):
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    if args.bone not in config["bones"]:
        raise RuntimeError(f"Unknown bone in configuration: {args.bone}")
    bone_cfg = config["bones"][args.bone]
    if args.method != "regional" and (args.thickness is None or args.thickness <= 0.0):
        raise RuntimeError("A positive --thickness is required for constant-thickness comparators")
    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)
    clear_scene()
    total_started = time.perf_counter()
    import_started = time.perf_counter()
    source = import_stl(args.source, f"{args.bone}_source_imported")
    import_seconds = time.perf_counter() - import_started
    source_raw = mesh_diagnostics(source)
    historical_raw_reproduction = bool(
        args.method == "regional"
        and args.regional_kernel == "legacy_deform"
        and args.allow_invalid_source
    )
    if not historical_raw_reproduction:
        recalc_normals_outside(source)
    source_prepared = mesh_diagnostics(source)
    source_failures = diagnostic_failures(source_prepared)
    if source_failures and not args.allow_invalid_source:
        raise RuntimeError("Fail-closed source validation: " + ", ".join(source_failures))

    generation_started = time.perf_counter()
    outer = None
    validated_spine = None
    if args.method == "regional":
        if args.regional_kernel == "validated_loft":
            outer, inner, validated_spine = generate_regional_loft(
                source, bone_cfg, config["benchmark"]
            )
        else:
            outer, inner = generate_regional(source, bone_cfg, postprocess=args.regional_postprocess)
    elif args.method == "uniform_normal":
        outer = source
        inner = generate_uniform_normal(source, args.thickness, args.bone)
    else:
        outer = source
        inner = generate_solidify(source, args.thickness, args.bone, config["methods"]["solidify"])
    generation_seconds = time.perf_counter() - generation_started

    outer_diagnostics = mesh_diagnostics(outer)
    inner_diagnostics = mesh_diagnostics(inner)
    containment_guard = {"enabled": False}
    if args.method == "regional" and args.containment_guard:
        benchmark_cfg = config["benchmark"]
        if validated_spine is not None:
            containment_guard = enforce_loft_containment(
                outer,
                inner,
                validated_spine,
                int(benchmark_cfg["validated_loft_angles"]),
                float(benchmark_cfg["containment_guard_minimum_clearance_mm"]),
                float(benchmark_cfg["containment_guard_minimum_radial_scale"]),
                float(benchmark_cfg["containment_guard_search_step"]),
                int(benchmark_cfg["containment_guard_binary_iterations"]),
            )
        else:
            containment_guard = enforce_containment(
                outer,
                inner,
                bone_cfg["parameters"],
                float(benchmark_cfg["containment_guard_minimum_clearance_mm"]),
                float(benchmark_cfg["containment_guard_minimum_radial_scale"]),
                float(benchmark_cfg["containment_guard_search_step"]),
                int(benchmark_cfg["containment_guard_binary_iterations"]),
            )
        inner_diagnostics = mesh_diagnostics(inner)
    benchmark_cfg = config["benchmark"]
    relation = relation_diagnostics(
        outer,
        inner,
        tolerance=float(benchmark_cfg["geometry_tolerance_mm"]),
        required_clearance=float(benchmark_cfg["containment_guard_minimum_clearance_mm"]),
    )
    topology_failures = diagnostic_failures(inner_diagnostics)
    geometry_failures = []
    if relation["self_intersection_pairs"]:
        geometry_failures.append(f"self_intersection_pairs={relation['self_intersection_pairs']}")
    if relation["outer_inner_intersection_pairs"]:
        geometry_failures.append(
            f"outer_inner_intersection_pairs={relation['outer_inner_intersection_pairs']}"
        )
    if relation["outside_vertices_count"]:
        geometry_failures.append(f"outside_vertices_count={relation['outside_vertices_count']}")
    clearance_failures = []
    if relation["clearance_violation_count"]:
        clearance_failures.append(
            f"clearance_violation_count={relation['clearance_violation_count']}"
        )
    output_failures = topology_failures + geometry_failures + clearance_failures
    if output_failures and not args.allow_invalid_output:
        raise RuntimeError("Fail-closed output validation: " + ", ".join(output_failures))

    export_started = time.perf_counter()
    if args.method == "regional":
        outer_path = os.path.join(output_dir, bone_cfg["outer_filename"])
        inner_path = os.path.join(output_dir, bone_cfg["inner_filename"])
        export_stl(outer, outer_path)
    else:
        outer_path = None
        inner_path = os.path.join(output_dir, f"{args.bone}_{args.method}_inner.stl")
    export_stl(inner, inner_path)
    export_seconds = time.perf_counter() - export_started

    metadata = {
        "schema_version": 2,
        "bone": args.bone,
        "method": args.method,
        "analysis_scope": bone_cfg["analysis_scope"],
        "blender_version": bpy.app.version_string,
        "source_stl": os.path.abspath(args.source),
        "source_raw": source_raw,
        "source_prepared": source_prepared,
        "source_normals_recalculated_before_generation": not historical_raw_reproduction,
        "source_topology_valid": not source_failures,
        "source_validation_failures": source_failures,
        "allow_invalid_source": bool(args.allow_invalid_source),
        "outer": outer_diagnostics,
        "inner": inner_diagnostics,
        "output_topology_valid": not topology_failures,
        "output_geometry_valid": not topology_failures and not geometry_failures,
        "output_clearance_valid": not topology_failures and not geometry_failures and not clearance_failures,
        "output_all_valid": not output_failures,
        "topology_validation_failures": topology_failures,
        "geometry_validation_failures": geometry_failures,
        "clearance_validation_failures": clearance_failures,
        "output_validation_failures": output_failures,
        "allow_invalid_output": bool(args.allow_invalid_output),
        "containment_guard": containment_guard,
        "relation_diagnostics": relation,
        "thickness_mm": args.thickness,
        "regional_postprocess": args.regional_postprocess if args.method == "regional" else None,
        "regional_kernel": args.regional_kernel if args.method == "regional" else None,
        "outer_path": outer_path,
        "inner_path": inner_path,
        "outer_mesh_signature": mesh_signature(outer),
        "inner_mesh_signature": mesh_signature(inner),
        "timing_seconds": {
            "import": import_seconds,
            "generation": generation_seconds,
            "export": export_seconds,
            "total_internal": time.perf_counter() - total_started,
        },
    }
    metadata_path = os.path.join(output_dir, f"{args.method}_run_metadata.json")
    Path(metadata_path).write_text(json.dumps(json_safe(metadata), indent=2, sort_keys=True), encoding="utf-8")
    print("FORARM3D_RUN=" + json.dumps(json_safe(metadata), sort_keys=True))
    return metadata


def main():
    args = arguments()
    try:
        run(args)
    except Exception:
        traceback.print_exc()
        sys.exit(2)


if __name__ == "__main__":
    main()
