#!/usr/bin/env python3
"""Prepare one external long-bone envelope without geometric repair.

Allowed operations are deliberately narrow: load, concatenate source objects,
merge coincident vertices, select the largest face-connected component, orient
normals, center, and rigidly align the principal longitudinal axis to Z.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import trimesh


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--units", default="mm", choices=("mm",))
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def geometry_sha(mesh: trimesh.Trimesh) -> str:
    digest = hashlib.sha256()
    digest.update(np.asarray(mesh.vertices, dtype="<f8").tobytes(order="C"))
    digest.update(np.asarray(mesh.faces, dtype="<i8").tobytes(order="C"))
    return digest.hexdigest()


def diagnostics(mesh: trimesh.Trimesh) -> dict:
    extents = np.asarray(mesh.extents, dtype=float)
    degenerate = int(np.count_nonzero(mesh.area_faces <= 1.0e-12))
    return {
        "vertices": int(len(mesh.vertices)),
        "faces": int(len(mesh.faces)),
        "body_count": int(mesh.body_count),
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "is_volume": bool(mesh.is_volume),
        "euler_number": int(mesh.euler_number),
        "degenerate_faces": degenerate,
        "volume_mm3": float(mesh.volume),
        "surface_area_mm2": float(mesh.area),
        "bbox_extents_mm": extents.tolist(),
        "bbox_volume_mm3": float(np.prod(extents)),
        "geometry_sha256": geometry_sha(mesh),
    }


def load_mesh(path: Path) -> trimesh.Trimesh:
    loaded = trimesh.load(path, process=False, maintain_order=True)
    if isinstance(loaded, trimesh.Scene):
        geometries = [
            geometry.copy()
            for _, geometry in sorted(loaded.geometry.items(), key=lambda item: item[0])
            if isinstance(geometry, trimesh.Trimesh) and len(geometry.faces)
        ]
        if not geometries:
            raise RuntimeError("No triangular mesh geometry found in source scene")
        return trimesh.util.concatenate(geometries)
    if not isinstance(loaded, trimesh.Trimesh) or not len(loaded.faces):
        raise RuntimeError("Source does not contain a non-empty triangular mesh")
    return loaded.copy()


def face_components(mesh: trimesh.Trimesh) -> list[np.ndarray]:
    adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
    nodes = np.arange(len(mesh.faces), dtype=np.int64)
    return [
        np.asarray(component, dtype=np.int64)
        for component in trimesh.graph.connected_components(
            adjacency, nodes=nodes, min_len=1, engine=None
        )
    ]


def select_primary_component(mesh: trimesh.Trimesh) -> tuple[trimesh.Trimesh, list[dict]]:
    records = []
    for index, component_faces in enumerate(face_components(mesh)):
        component = mesh.submesh([component_faces], append=True, repair=False)
        record = diagnostics(component)
        record["component_index"] = index
        records.append((record, component))
    if not records:
        raise RuntimeError("No face-connected component found")
    records.sort(
        key=lambda item: (
            -item[0]["faces"],
            -abs(item[0]["volume_mm3"]),
            -item[0]["bbox_volume_mm3"],
            item[0]["geometry_sha256"],
        )
    )
    for rank, (record, _) in enumerate(records, start=1):
        record["selection_rank"] = rank
    return records[0][1].copy(), [record for record, _ in records]


def canonical_sign(vector: np.ndarray) -> np.ndarray:
    index = int(np.argmax(np.abs(vector)))
    return vector if vector[index] >= 0.0 else -vector


def pca_align(mesh: trimesh.Trimesh) -> tuple[trimesh.Trimesh, dict]:
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    centroid = vertices.mean(axis=0)
    centered = vertices - centroid
    covariance = centered.T @ centered / float(len(centered))
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    z_axis = canonical_sign(eigenvectors[:, 0])
    x_axis = eigenvectors[:, 1] - z_axis * float(np.dot(z_axis, eigenvectors[:, 1]))
    x_axis /= np.linalg.norm(x_axis)
    x_axis = canonical_sign(x_axis)
    y_axis = np.cross(z_axis, x_axis)
    y_axis /= np.linalg.norm(y_axis)
    x_axis = np.cross(y_axis, z_axis)
    x_axis /= np.linalg.norm(x_axis)
    basis = np.column_stack((x_axis, y_axis, z_axis))
    if np.linalg.det(basis) <= 0.0:
        raise RuntimeError("PCA basis is not right-handed")
    aligned = mesh.copy()
    aligned.vertices = centered @ basis
    transform = {
        "centroid_original_mm": centroid.tolist(),
        "covariance_eigenvalues_mm2": eigenvalues.tolist(),
        "basis_columns_original_coordinates": basis.tolist(),
        "determinant": float(np.linalg.det(basis)),
        "formula": "prepared_row = (original_row - centroid_row) @ basis",
    }
    return aligned, transform


def main() -> None:
    args = arguments()
    source = Path(args.source).resolve()
    output = Path(args.output).resolve()
    report_path = Path(args.report).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    raw = load_mesh(source)
    raw_loaded = diagnostics(raw)
    before_merge_vertices = len(raw.vertices)
    raw.merge_vertices(merge_tex=True, merge_norm=True)
    after_merge = diagnostics(raw)
    primary, components = select_primary_component(raw)
    selected_before_normals = diagnostics(primary)
    primary.fix_normals(multibody=True)
    selected_after_normals = diagnostics(primary)
    source_failures = []
    if not selected_after_normals["watertight"]:
        source_failures.append("source_not_watertight")
    if not selected_after_normals["winding_consistent"]:
        source_failures.append("source_winding_inconsistent")
    if selected_after_normals["degenerate_faces"]:
        source_failures.append(
            f"source_degenerate_faces={selected_after_normals['degenerate_faces']}"
        )
    if source_failures:
        raise RuntimeError("Fail-closed external source: " + ", ".join(source_failures))
    aligned, transformation = pca_align(primary)
    aligned.fix_normals(multibody=True)
    aligned_diagnostics = diagnostics(aligned)
    output.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(aligned.export(file_type="stl"))
    exported = trimesh.load_mesh(output, process=True)
    exported_diagnostics = diagnostics(exported)
    if not exported_diagnostics["watertight"] or exported_diagnostics["degenerate_faces"]:
        raise RuntimeError("Prepared STL failed post-export validation")
    length_mm = float(exported.bounds[1, 2] - exported.bounds[0, 2])
    report = {
        "schema_version": 1,
        "case_id": args.case_id,
        "units": args.units,
        "source": str(source),
        "source_sha256": sha256(source),
        "output": str(output),
        "output_sha256": sha256(output),
        "allowed_operations": [
            "source object concatenation",
            "coincident vertex merge without smoothing",
            "largest face-connected component selection",
            "normal orientation",
            "rigid PCA centering and alignment",
            "binary STL serialization",
        ],
        "forbidden_operations_applied": [],
        "raw_loaded": raw_loaded,
        "coincident_vertex_merge": {
            "vertices_before": int(before_merge_vertices),
            "vertices_after": int(len(raw.vertices)),
        },
        "after_vertex_merge": after_merge,
        "components_ranked": components,
        "selected_before_normals": selected_before_normals,
        "selected_after_normals": selected_after_normals,
        "pca_transformation": transformation,
        "aligned_before_export": aligned_diagnostics,
        "prepared_after_export": exported_diagnostics,
        "longitudinal_length_mm": length_mm,
    }
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"case_id": args.case_id, "length_mm": length_mm, "output": str(output)}))


if __name__ == "__main__":
    main()

