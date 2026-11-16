#!/usr/bin/env python3
"""Materialize one v2-compatible case config from the frozen normalized profile."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def clamp(value: float, bounds: list[float]) -> float:
    return max(float(bounds[0]), min(float(bounds[1]), value))


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--prepared-report", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> None:
    args = arguments()
    report = json.loads(Path(args.prepared_report).read_text(encoding="utf-8"))
    frozen = json.loads(Path(args.profile).read_text(encoding="utf-8"))
    profile = frozen["profile"]
    length = float(report["longitudinal_length_mm"])
    terminal = clamp(
        length * profile["terminal_backing_fraction_of_length"],
        profile["terminal_backing_bounds_mm"],
    )
    minimum_inner = clamp(
        length * profile["minimum_inner_radius_fraction_of_length"],
        profile["minimum_inner_radius_bounds_mm"],
    )
    epi_shell = clamp(
        length * profile["epiphyseal_shell_fraction_of_length"],
        profile["epiphyseal_shell_bounds_mm"],
    )
    meta_shell = clamp(
        length * profile["metaphyseal_shell_fraction_of_length"],
        profile["metaphyseal_shell_bounds_mm"],
    )
    end_zone = length * profile["end_zone_fraction_of_length"]
    flare = length * profile["flare_fraction_of_length"]
    plateau = length * profile["epiphyseal_plateau_fraction_of_length"]
    epi_fade = length * profile["epiphyseal_to_metaphyseal_fade_fraction_of_length"]
    zero_fade = length * profile["metaphyseal_to_zero_fade_fraction_of_length"]
    parameters = {
        "nb_sections": int(profile["nb_sections"]),
        "z_window_ratio": float(profile["z_window_ratio"]),
        "spine_smooth_window": int(profile["spine_smooth_window"]),
        "canal_ratio_isthme": float(profile["canal_ratio_isthme"]),
        "canal_ratio_diaphyse": float(profile["canal_ratio_diaphyse"]),
        "canal_ratio_metaphyse": float(profile["canal_ratio_metaphyse"]),
        "half_isthme": float(profile["half_isthme"]),
        "half_diaphyse": float(profile["half_diaphyse"]),
        "end_zone_dist_mm": end_zone,
        "end_zone_prox_mm": end_zone,
        "flare_dist_mm": flare,
        "flare_prox_mm": flare,
        "terminal_backing_offset_mm": terminal,
        "rmin_absolute_mm": minimum_inner,
        "prox_epi_shell_min_mm": epi_shell,
        "prox_meta_shell_min_mm": meta_shell,
        "prox_epi_plateau_mm": plateau,
        "prox_epi_to_meta_fade_mm": epi_fade,
        "prox_meta_to_zero_fade_mm": zero_fade,
        "dist_epi_shell_min_mm": epi_shell,
        "dist_meta_shell_min_mm": meta_shell,
        "dist_epi_plateau_mm": plateau,
        "dist_epi_to_meta_fade_mm": epi_fade,
        "dist_meta_to_zero_fade_mm": zero_fade,
        "use_decimate": bool(profile["use_decimate"]),
        "decimate_ratio": float(profile["decimate_ratio"]),
        "smooth_iterations": int(profile["smooth_iterations"]),
        "merge_distance_mm": float(profile["merge_distance_mm"]),
    }
    loft = frozen["validated_loft"]
    guard = frozen["containment_guard"]
    validation = frozen["validation"]
    config = {
        "schema_version": 3,
        "external_profile_id": frozen["profile_id"],
        "case_materialization": {
            "case_id": args.case_id,
            "prepared_report": str(Path(args.prepared_report).resolve()),
            "longitudinal_length_mm": length,
        },
        "benchmark": {
            "geometry_tolerance_mm": float(validation["geometry_tolerance_mm"]),
            "self_intersection_epsilon_mm": float(validation["self_intersection_epsilon_mm"]),
            "containment_guard_minimum_clearance_mm": float(guard["minimum_clearance_mm"]),
            "containment_guard_minimum_radial_scale": float(guard["minimum_radial_scale"]),
            "containment_guard_search_step": float(guard["search_step"]),
            "containment_guard_binary_iterations": int(guard["binary_iterations"]),
            "validated_loft_sections": int(loft["sections"]),
            "validated_loft_angles": int(loft["angles_per_section"]),
            "validated_spine_grid_size": int(loft["spine_grid_size"]),
            "validated_loft_z_margin_mm": float(loft["z_margin_mm"]),
            "validated_loft_z_margin_fraction": float(loft["z_margin_fraction"]),
            "validated_loft_radial_safety_mm": float(loft["radial_safety_mm"]),
            "random_seed": int(validation["random_seed"]),
        },
        "methods": {"solidify": {}},
        "bones": {
            args.case_id: {
                "label": args.case_id,
                "analysis_scope": "external_validation_portability",
                "outer_filename": f"{args.case_id}_outer.stl",
                "inner_filename": f"{args.case_id}_inner_negative.stl",
                "parameters": parameters,
            }
        },
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(config, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"case_id": args.case_id, "length_mm": length, "output": str(output)}))


if __name__ == "__main__":
    main()

