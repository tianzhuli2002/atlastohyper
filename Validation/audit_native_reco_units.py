#!/usr/bin/env python3
"""Audit native dilepton constituent units against objects and MET in ROOT."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


METHODS = ("DA", "EM", "NW", "VyPER")


def native_p4(tree, name: str, stop: int) -> np.ndarray:
    values = tree[name].array(entry_stop=stop, library="np")["fCoordinates"]
    pt = np.asarray(values["fCoordinates.fPt"], dtype=np.float64)
    eta = np.asarray(values["fCoordinates.fEta"], dtype=np.float64)
    phi = np.asarray(values["fCoordinates.fPhi"], dtype=np.float64)
    mass = np.asarray(values["fCoordinates.fM"], dtype=np.float64)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(np.maximum(px * px + py * py + pz * pz + mass * mass, 0.0))
    return np.stack((px, py, pz, energy), axis=-1)


def max_component_difference(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.max(np.abs(left - right)))


def invariant_mass(p4: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(p4[:, 3] ** 2 - np.sum(p4[:, :3] ** 2, axis=1), 0.0))


def scale_summary(values: list[float]) -> dict[str, float | int]:
    array = np.asarray(values, dtype=np.float64)
    if not len(array):
        return {"matched": 0}
    return {
        "matched": int(len(array)),
        "median_native_over_object_pt": float(np.median(array)),
        "p16_native_over_object_pt": float(np.quantile(array, 0.16)),
        "p84_native_over_object_pt": float(np.quantile(array, 0.84)),
        "fraction_scale_1": float(np.mean(np.isclose(array, 1.0, rtol=1e-5, atol=1e-5))),
        "fraction_scale_0p001": float(np.mean(np.isclose(array, 0.001, rtol=1e-5, atol=1e-7))),
    }


def object_match_scales(tree, component: np.ndarray, stop: int, prefix: str) -> dict:
    collection_names = ("jet",) if prefix == "jet" else ("el", "mu")
    object_arrays = {}
    for collection in collection_names:
        object_arrays[collection] = tuple(
            ak.to_list(tree[f"{collection}_{field}"].array(entry_stop=stop, library="ak"))
            for field in ("pt_NOSYS", "eta", "phi")
        )

    scales: list[float] = []
    matched_drs: list[float] = []
    for row, p4 in enumerate(component):
        pt = float(np.hypot(p4[0], p4[1]))
        eta = float(np.arcsinh(p4[2] / pt)) if pt > 0 else np.nan
        phi = float(np.arctan2(p4[1], p4[0]))
        best = None
        for collection, (pts, etas, phis) in object_arrays.items():
            pts, etas, phis = (np.asarray(field[row], dtype=np.float64) for field in (pts, etas, phis))
            if not len(pts):
                continue
            dphi = (phis - phi + np.pi) % (2.0 * np.pi) - np.pi
            dr = np.hypot(etas - eta, dphi)
            index = int(np.argmin(dr))
            candidate = (float(dr[index]), float(pts[index]))
            if best is None or candidate[0] < best[0]:
                best = candidate
        if best is not None and best[0] < 1e-5 and best[1] > 0:
            matched_drs.append(best[0])
            scales.append(pt / best[1])

    result = scale_summary(scales)
    result["max_matched_deltaR"] = max(matched_drs, default=None)
    return result


def met_balance(tree, nu: np.ndarray, nubar: np.ndarray, met: np.ndarray, met_phi: np.ndarray) -> dict:
    met_px = met * np.cos(met_phi)
    met_py = met * np.sin(met_phi)
    result = {}
    for scale in (1.0, 1000.0):
        residual = np.hypot(scale * (nu[:, 0] + nubar[:, 0]) - met_px,
                            scale * (nu[:, 1] + nubar[:, 1]) - met_py) / 1000.0
        result[f"component_scale_{int(scale)}"] = {
            "median_MET_residual_GeV": float(np.median(residual)),
            "p16_MET_residual_GeV": float(np.quantile(residual, 0.16)),
            "p84_MET_residual_GeV": float(np.quantile(residual, 0.84)),
        }
    return result


def audit(path: Path, max_events: int) -> dict:
    with uproot.open(path) as source:
        if "reco" not in source:
            raise ValueError(f"{path}: required tree 'reco' is absent")
        tree = source["reco"]
        entries = min(int(tree.num_entries), max_events)
        if entries == 0:
            raise ValueError(f"{path}: tree 'reco' is empty")
        required = ["met_met_NOSYS", "met_phi_NOSYS"]
        for method in METHODS:
            required.extend(
                f"TtbarDiLepton_{method}_{suffix}_NOSYS"
                for suffix in (
                    "Top_b_p4", "Top_lep_p4", "Top_nu_p4", "Top",
                    "Topbar_bbar_p4", "Topbar_lepbar_p4", "Topbar_nubar_p4", "Tbar",
                )
            )
        required.extend(f"{collection}_{field}" for collection in ("jet", "el", "mu")
                        for field in ("pt_NOSYS", "eta", "phi"))
        missing = [name for name in required if name not in tree]
        if missing:
            raise ValueError(f"{path}: missing unit-audit branches: {missing}")

        met = tree["met_met_NOSYS"].array(entry_stop=entries, library="np").astype(np.float64)
        met_phi = tree["met_phi_NOSYS"].array(entry_stop=entries, library="np").astype(np.float64)
        results = {}
        for method in METHODS:
            stem = f"TtbarDiLepton_{method}_"
            top_b = native_p4(tree, stem + "Top_b_p4_NOSYS", entries)
            top_lep = native_p4(tree, stem + "Top_lep_p4_NOSYS", entries)
            top_nu = native_p4(tree, stem + "Top_nu_p4_NOSYS", entries)
            tbar_b = native_p4(tree, stem + "Topbar_bbar_p4_NOSYS", entries)
            tbar_lep = native_p4(tree, stem + "Topbar_lepbar_p4_NOSYS", entries)
            tbar_nu = native_p4(tree, stem + "Topbar_nubar_p4_NOSYS", entries)
            stored_top = native_p4(tree, stem + "Top_NOSYS", entries)
            stored_tbar = native_p4(tree, stem + "Tbar_NOSYS", entries)
            raw_top_sum = top_b + top_lep + top_nu
            raw_tbar_sum = tbar_b + tbar_lep + tbar_nu
            corrected_top_sum = top_b + 1000.0 * (top_lep + top_nu)
            corrected_tbar_sum = tbar_b + 1000.0 * (tbar_lep + tbar_nu)
            results[method] = {
                "top_lepton_pt_unit_match": object_match_scales(tree, top_lep, entries, "lepton"),
                "antitop_lepton_pt_unit_match": object_match_scales(tree, tbar_lep, entries, "lepton"),
                "top_b_pt_unit_match": object_match_scales(tree, top_b, entries, "jet"),
                "antitop_b_pt_unit_match": object_match_scales(tree, tbar_b, entries, "jet"),
                "neutrino_MET_balance": met_balance(tree, top_nu, tbar_nu, met, met_phi),
                "stored_top_closure_using_native_numeric_scales": {
                    "max_abs_component_difference": max_component_difference(raw_top_sum, stored_top),
                    "median_four_vector_norm_difference": float(np.median(np.linalg.norm(raw_top_sum - stored_top, axis=1))),
                },
                "stored_antitop_closure_using_native_numeric_scales": {
                    "max_abs_component_difference": max_component_difference(raw_tbar_sum, stored_tbar),
                    "median_four_vector_norm_difference": float(np.median(np.linalg.norm(raw_tbar_sum - stored_tbar, axis=1))),
                },
                "top_sum_if_lepton_and_neutrino_scaled_by_1000_vs_stored_top": {
                    "median_four_vector_norm_difference": float(np.median(np.linalg.norm(corrected_top_sum - stored_top, axis=1))),
                    "median_mass_GeV": float(np.median(invariant_mass(corrected_top_sum)) / 1000.0),
                },
                "antitop_sum_if_lepton_and_neutrino_scaled_by_1000_vs_stored_antitop": {
                    "median_four_vector_norm_difference": float(np.median(np.linalg.norm(corrected_tbar_sum - stored_tbar, axis=1))),
                    "median_mass_GeV": float(np.median(invariant_mass(corrected_tbar_sum)) / 1000.0),
                },
                "stored_top_mass_GeV": {
                    "median": float(np.median(invariant_mass(stored_top)) / 1000.0),
                },
                "stored_antitop_mass_GeV": {
                    "median": float(np.median(invariant_mass(stored_tbar)) / 1000.0),
                },
            }

        return {
            "input": str(path),
            "tree": "reco",
            "input_entries": int(tree.num_entries),
            "entries_audited": entries,
            "assumed_object_and_MET_units": "MeV",
            "neutrino_MET_diagnostic": "vector sum of the two reconstructed neutrinos compared with MET; a diagnostic, not a strict equality requirement",
            "scaled_sum_diagnostic": "The factor 1000 comparison tests whether lepton/neutrino constituent p4 values use GeV while b jets use MeV; it does not rewrite native outputs.",
            "methods": results,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="representative raw or rich FastFrames ROOT file")
    parser.add_argument("--max-events", type=int, default=500)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = audit(args.root, args.max_events)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
