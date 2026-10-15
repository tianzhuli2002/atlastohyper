#!/usr/bin/env python3
"""Compare raw HyPER+NuFlows top vectors with constituent sums in pilot events."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import uproot


def p4_from_pt_eta_phi_m(tree, prefix: str, rows: np.ndarray) -> np.ndarray:
    values = []
    for suffix in ("fPt", "fEta", "fPhi", "fM"):
        key = f"{prefix}/fCoordinates/fCoordinates.{suffix}"
        values.append(tree[key].array(library="np")[rows].astype(np.float64))
    pt, eta, phi, mass = values
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(np.maximum(px * px + py * py + pz * pz + mass * mass, 0.0))
    return np.stack((px, py, pz, energy), axis=1)


def scalar_p4(tree, prefix: str) -> np.ndarray:
    pt = tree[f"{prefix}_pt_HYPER_NUFLOWS_NOSYS"].array(library="np").astype(np.float64)
    eta = tree[f"{prefix}_eta_HYPER_NUFLOWS_NOSYS"].array(library="np").astype(np.float64)
    phi = tree[f"{prefix}_phi_HYPER_NUFLOWS_NOSYS"].array(library="np").astype(np.float64)
    mass = tree[f"{prefix}_m_HYPER_NUFLOWS_NOSYS"].array(library="np").astype(np.float64)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(np.maximum(px * px + py * py + pz * pz + mass * mass, 0.0))
    return np.stack((px, py, pz, energy), axis=1)


def residual(left: np.ndarray, right: np.ndarray) -> dict:
    differences = np.abs(left - right)
    return {
        "events": int(len(left)),
        "median_abs_by_cartesian_component": np.median(differences, axis=0).tolist(),
        "max_abs_by_cartesian_component": np.max(differences, axis=0).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", required=True, type=Path)
    parser.add_argument("--fastframes", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    with uproot.open(args.raw) as raw_file, uproot.open(args.fastframes) as output_file:
        raw = raw_file["reco"]
        output = output_file["reco"]
        raw_ids = raw["eventNumber"].array(library="np").astype(np.uint64)
        output_ids = output["eventNumber"].array(library="np").astype(np.uint64)
        raw_index = {int(event): index for index, event in enumerate(raw_ids)}
        if len(raw_index) != len(raw_ids) or len(np.unique(output_ids)) != len(output_ids):
            raise AssertionError("eventNumber must be unique in the raw and FastFrames files")
        missing = [int(event) for event in output_ids if int(event) not in raw_index]
        if missing:
            raise AssertionError(f"FastFrames events absent from raw file: {missing[:5]}")
        raw_rows = np.asarray([raw_index[int(event)] for event in output_ids], dtype=np.int64)
        result = {
            "raw": str(args.raw),
            "fastframes": str(args.fastframes),
            "events_joined_by_eventNumber": int(len(output_ids)),
            "units_observation": "HyPER constituent p4 branches are on a MeV scale; NuFlows p4 branches are on a GeV scale",
            "custom_native_note": "The custom path re-decodes the HyPER assignment and charge-orders leptons. Its difference from the native combined branch remains unclassified; do not treat the constituent residual as a validated custom-versus-native comparison.",
            "methods": {},
        }

        for side, raw_parts, native_prefix, custom_prefix in (
            (
                "top",
                (
                    "TtbarDiLepton_HyPER_Top_b_p4_NOSYS",
                    "TtbarDiLepton_HyPER_Top_lep_p4_NOSYS",
                    "TtbarDiLepton_NuFlows_Top_nu_p4_NOSYS",
                ),
                "TtbarDiLepton_HyPER_plus_NuFlows_Top_NOSYS",
                "top",
            ),
            (
                "antitop",
                (
                    "TtbarDiLepton_HyPER_Topbar_bbar_p4_NOSYS",
                    "TtbarDiLepton_HyPER_Topbar_lepbar_p4_NOSYS",
                    "TtbarDiLepton_NuFlows_Topbar_nubar_p4_NOSYS",
                ),
                "TtbarDiLepton_HyPER_plus_NuFlows_Tbar_NOSYS",
                "antitop",
            ),
        ):
            b, lepton, neutrino = [p4_from_pt_eta_phi_m(raw, name, raw_rows) for name in raw_parts]
            native = p4_from_pt_eta_phi_m(raw, native_prefix, raw_rows)
            custom = scalar_p4(output, custom_prefix)
            result["methods"][side] = {
                "native_vs_b_plus_lepton_plus_unscaled_GeV_nu": residual(native, b + lepton + neutrino),
                "native_vs_b_plus_lepton_plus_GeV_to_MeV_nu": residual(
                    native, b + lepton + 1000.0 * neutrino
                ),
                "custom_vs_native": residual(custom, native),
            }
    text = json.dumps(result, indent=2)
    print(text)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text + "\n")


if __name__ == "__main__":
    main()
