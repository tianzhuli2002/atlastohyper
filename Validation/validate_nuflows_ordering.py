#!/usr/bin/env python3
"""Validate NuFlows output-slot order against its native parent-side vectors."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


def p4_branch(tree, name: str, stop: int) -> np.ndarray:
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


def validate(path: Path, max_events: int, tolerance_gev: float) -> dict:
    with uproot.open(path) as source:
        if "reco" not in source:
            raise ValueError(f"{path}: required tree 'reco' is absent")
        tree = source["reco"]
        entries = min(int(tree.num_entries), max_events)
        required = [
            "TtbarDiLepton_NuFlows_nuflows_nu_out_NOSYS",
            "TtbarDiLepton_NuFlows_Top_nu_p4_NOSYS",
            "TtbarDiLepton_NuFlows_Topbar_nubar_p4_NOSYS",
        ] + [
            f"{prefix}_{field}"
            for prefix in ("el", "mu")
            for field in ("pt_NOSYS", "charge", "select_tight_NOSYS")
        ]
        missing = [name for name in required if name not in tree]
        if missing:
            raise ValueError(f"{path}: missing NuFlows ordering branches: {missing}")

        outputs = ak.to_list(
            tree["TtbarDiLepton_NuFlows_nuflows_nu_out_NOSYS"].array(entry_stop=entries, library="ak")
        )
        top_nu = p4_branch(tree, "TtbarDiLepton_NuFlows_Top_nu_p4_NOSYS", entries)
        antitop_nubar = p4_branch(tree, "TtbarDiLepton_NuFlows_Topbar_nubar_p4_NOSYS", entries)
        collections = {}
        for prefix in ("el", "mu"):
            collections[prefix] = {
                field: ak.to_list(tree[f"{prefix}_{field}"].array(entry_stop=entries, library="ak"))
                for field in ("pt_NOSYS", "charge", "select_tight_NOSYS")
            }

        checked = 0
        malformed_output = 0
        wrong_lepton_count = 0
        invalid_native_vectors = 0
        exact_charge_order_matches = 0
        slot0_parent_counts = {"top": 0, "antitop": 0}
        channel_counts: dict[str, dict[str, int]] = {}
        max_expected_component_difference = 0.0
        max_best_component_difference = 0.0

        for row in range(entries):
            values = np.asarray(outputs[row], dtype=np.float64)
            if values.shape != (6,) or not np.all(np.isfinite(values)):
                malformed_output += 1
                continue

            # The native NuFlows input contract uses the selected electron
            # collection followed by the selected muon collection, preserving
            # each collection's existing event order.
            ordered_lepton_charges: list[int] = []
            channel = []
            for prefix in ("el", "mu"):
                pts = collections[prefix]["pt_NOSYS"][row]
                charges = collections[prefix]["charge"][row]
                selected = collections[prefix]["select_tight_NOSYS"][row]
                if not (len(pts) == len(charges) == len(selected)):
                    raise AssertionError(f"{path}: {prefix} object fields differ in length at row {row}")
                for charge, keep in zip(charges, selected):
                    if bool(keep):
                        ordered_lepton_charges.append(int(np.sign(charge)))
                        channel.append("e" if prefix == "el" else "m")
            if len(ordered_lepton_charges) != 2 or 0 in ordered_lepton_charges:
                wrong_lepton_count += 1
                continue
            if ordered_lepton_charges[0] == ordered_lepton_charges[1]:
                wrong_lepton_count += 1
                continue

            vectors = []
            for offset in (0, 3):
                px, py, pz = values[offset : offset + 3]
                energy = float(np.sqrt(px * px + py * py + pz * pz))
                vectors.append(np.asarray((px, py, pz, energy), dtype=np.float64))
            if not np.all(np.isfinite(top_nu[row])) or not np.all(np.isfinite(antitop_nubar[row])):
                invalid_native_vectors += 1
                continue
            if np.hypot(top_nu[row, 0], top_nu[row, 1]) <= 0 or np.hypot(antitop_nubar[row, 0], antitop_nubar[row, 1]) <= 0:
                invalid_native_vectors += 1
                continue

            identity = [
                np.max(np.abs(vectors[0] - top_nu[row])),
                np.max(np.abs(vectors[1] - antitop_nubar[row])),
            ]
            swapped = [
                np.max(np.abs(vectors[0] - antitop_nubar[row])),
                np.max(np.abs(vectors[1] - top_nu[row])),
            ]
            identity_cost = float(sum(identity))
            swapped_cost = float(sum(swapped))
            best_is_identity = identity_cost <= swapped_cost
            expected_slot0_top = ordered_lepton_charges[0] > 0
            expected_is_identity = expected_slot0_top
            expected_diffs = identity if expected_is_identity else swapped
            best_diffs = identity if best_is_identity else swapped
            max_expected_component_difference = max(max_expected_component_difference, *map(float, expected_diffs))
            max_best_component_difference = max(max_best_component_difference, *map(float, best_diffs))
            checked += 1

            slot0_top = best_is_identity
            slot0_parent_counts["top" if slot0_top else "antitop"] += 1
            channel_key = "".join(channel)
            counts = channel_counts.setdefault(channel_key, {"events": 0, "charge_order_matches": 0})
            counts["events"] += 1
            matches = (best_is_identity == expected_is_identity) and max(expected_diffs) <= tolerance_gev
            if matches:
                exact_charge_order_matches += 1
                counts["charge_order_matches"] += 1

        if checked == 0:
            raise AssertionError(f"{path}: no well-formed two-lepton NuFlows events were checked")
        if max_best_component_difference > tolerance_gev:
            raise AssertionError(
                f"{path}: NuFlows output slots do not match the native top/antitop neutrino branches; "
                f"maximum residual {max_best_component_difference:g} GeV exceeds {tolerance_gev:g} GeV"
            )
        if exact_charge_order_matches != checked:
            raise AssertionError(
                f"{path}: only {exact_charge_order_matches}/{checked} NuFlows slot pairs match the selected-lepton charge order"
            )
        if max_expected_component_difference > tolerance_gev:
            raise AssertionError(
                f"{path}: charge-mapped NuFlows components exceed tolerance ({max_expected_component_difference:g} GeV)"
            )

        return {
            "input": str(path),
            "tree": "reco",
            "input_entries": int(tree.num_entries),
            "entries_read": entries,
            "events_with_well_formed_two_lepton_order": checked,
            "malformed_six_component_outputs": malformed_output,
            "events_without_exactly_two_opposite_charge_selected_leptons": wrong_lepton_count,
            "events_with_invalid_native_nu_vectors": invalid_native_vectors,
            "output_units": "nuflows_nu_out components and native NuFlows constituent p4 branches are GeV",
            "input_slot_order": "selected electron collection in event order, followed by selected muon collection in event order",
            "parent_mapping": "positive charged lepton -> nu from top; negative charged lepton -> nubar from antitop",
            "slot0_native_parent_counts": slot0_parent_counts,
            "charge_order_matches": exact_charge_order_matches,
            "per_channel": channel_counts,
            "max_expected_cartesian_component_difference_GeV": max_expected_component_difference,
            "max_best_permutation_component_difference_GeV": max_best_component_difference,
            "comparison_tolerance_GeV": tolerance_gev,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="representative raw or rich FastFrames ROOT input")
    parser.add_argument("--max-events", type=int, default=500)
    parser.add_argument("--tolerance-gev", type=float, default=0.01)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = validate(args.root, args.max_events, args.tolerance_gev)
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
