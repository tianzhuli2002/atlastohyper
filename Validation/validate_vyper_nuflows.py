#!/usr/bin/env python3
"""Validate the fixed VyPER-assignment + native NuFlows reconstruction pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


CONSTITUENT_BRANCHES = {
    "b": "TtbarDiLepton_VyPER_Top_b_p4_NOSYS",
    "lep": "TtbarDiLepton_VyPER_Top_lep_p4_NOSYS",
    "bbar": "TtbarDiLepton_VyPER_Topbar_bbar_p4_NOSYS",
    "antilep": "TtbarDiLepton_VyPER_Topbar_lepbar_p4_NOSYS",
}


def read_native_p4(tree, name: str, entries: int) -> np.ndarray:
    values = tree[name].array(entry_stop=entries, library="np")["fCoordinates"]
    pt = np.asarray(values["fCoordinates.fPt"], dtype=np.float64)
    eta = np.asarray(values["fCoordinates.fEta"], dtype=np.float64)
    phi = np.asarray(values["fCoordinates.fPhi"], dtype=np.float64)
    mass = np.asarray(values["fCoordinates.fM"], dtype=np.float64)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(np.maximum(px * px + py * py + pz * pz + mass * mass, 0.0))
    return np.stack((px, py, pz, energy), axis=-1)


def invariant_mass(vectors: list[np.ndarray]) -> float:
    total = np.sum(np.stack(vectors), axis=0)
    mass2 = total[3] * total[3] - float(np.dot(total[:3], total[:3]))
    return float(np.sqrt(max(mass2, 0.0)))


def validate(path: Path, max_events: int, tolerance_mev: float) -> dict:
    with uproot.open(path) as source:
        if "reco" not in source:
            raise ValueError(f"{path}: required tree 'reco' is absent")
        tree = source["reco"]
        entries = min(int(tree.num_entries), max_events)
        required = [
            "eventNumber",
            "TtbarDiLepton_NuFlows_nuflows_nu_out_NOSYS",
            "el_charge",
            "el_select_tight_NOSYS",
            "mu_charge",
            "mu_select_tight_NOSYS",
            "VYPER_NUFLOWS_assignment_valid_NOSYS",
            "VYPER_NUFLOWS_neutrino_valid_NOSYS",
            "VYPER_NUFLOWS_reco_valid_NOSYS",
            "VYPER_NUFLOWS_failure_code_NOSYS",
            "VYPER_NUFLOWS_top_components_NOSYS",
            "mttbar_VYPER_NUFLOWS_NOSYS",
            "nu1_pt_VYPER_NUFLOWS_NOSYS",
            "nu1_pz_VYPER_NUFLOWS_NOSYS",
            "nu2_pt_VYPER_NUFLOWS_NOSYS",
            "nu2_pz_VYPER_NUFLOWS_NOSYS",
        ] + list(CONSTITUENT_BRANCHES.values())
        missing = [name for name in required if name not in tree]
        if missing:
            raise ValueError(f"{path}: missing required VyPER+NuFlows branches: {missing}")

        event_numbers = np.asarray(tree["eventNumber"].array(entry_stop=entries, library="np"))
        if len(np.unique(event_numbers)) != entries:
            raise AssertionError(f"{path}: duplicate eventNumber keys in first {entries} rows")

        nuflows = ak.to_list(
            tree["TtbarDiLepton_NuFlows_nuflows_nu_out_NOSYS"].array(entry_stop=entries, library="ak")
        )
        charges = {
            prefix: ak.to_list(tree[f"{prefix}_charge"].array(entry_stop=entries, library="ak"))
            for prefix in ("el", "mu")
        }
        selected = {
            prefix: ak.to_list(tree[f"{prefix}_select_tight_NOSYS"].array(entry_stop=entries, library="ak"))
            for prefix in ("el", "mu")
        }
        constituents = {
            key: read_native_p4(tree, name, entries)
            for key, name in CONSTITUENT_BRANCHES.items()
        }
        assignment_valid = np.asarray(
            tree["VYPER_NUFLOWS_assignment_valid_NOSYS"].array(entry_stop=entries, library="np"), dtype=int
        )
        neutrino_valid = np.asarray(
            tree["VYPER_NUFLOWS_neutrino_valid_NOSYS"].array(entry_stop=entries, library="np"), dtype=int
        )
        reco_valid = np.asarray(
            tree["VYPER_NUFLOWS_reco_valid_NOSYS"].array(entry_stop=entries, library="np"), dtype=int
        )
        failure_code = np.asarray(
            tree["VYPER_NUFLOWS_failure_code_NOSYS"].array(entry_stop=entries, library="np"), dtype=int
        )
        mtt = np.asarray(tree["mttbar_VYPER_NUFLOWS_NOSYS"].array(entry_stop=entries, library="np"), dtype=float)
        components = ak.to_list(
            tree["VYPER_NUFLOWS_top_components_NOSYS"].array(entry_stop=entries, library="ak")
        )
        output_nu = {
            branch: np.asarray(tree[branch].array(entry_stop=entries, library="np"), dtype=float)
            for branch in (
                "nu1_pt_VYPER_NUFLOWS_NOSYS",
                "nu1_pz_VYPER_NUFLOWS_NOSYS",
                "nu2_pt_VYPER_NUFLOWS_NOSYS",
                "nu2_pz_VYPER_NUFLOWS_NOSYS",
            )
        }

    selected_rows = 0
    invalid_rows = 0
    maximum_mtt_difference = 0.0
    maximum_component_mtt_difference = 0.0
    maximum_neutrino_kinematic_difference = 0.0
    mtt_values_gev = []

    for row in range(entries):
        if not (assignment_valid[row] and neutrino_valid[row] and reco_valid[row]):
            invalid_rows += 1
            if failure_code[row] == 0:
                raise AssertionError(f"{path}: row {row} has invalid reconstruction but failure code 0")
            continue
        if failure_code[row] != 0:
            raise AssertionError(f"{path}: row {row} is valid but has failure code {failure_code[row]}")
        flow = np.asarray(nuflows[row], dtype=np.float64)
        if flow.shape != (6,) or not np.all(np.isfinite(flow)):
            raise AssertionError(f"{path}: row {row} has malformed or non-finite NuFlows six-vector")

        selected_charges: list[float] = []
        for prefix in ("el", "mu"):
            row_charges = charges[prefix][row]
            row_selected = selected[prefix][row]
            if len(row_charges) != len(row_selected):
                raise AssertionError(f"{path}: row {row} has unequal {prefix} charge/selection lengths")
            selected_charges.extend(
                float(charge) for charge, keep in zip(row_charges, row_selected) if bool(keep)
            )
        if len(selected_charges) != 2 or selected_charges[0] * selected_charges[1] >= 0:
            raise AssertionError(f"{path}: valid row {row} does not have two opposite-charge selected leptons")

        slot_vectors = []
        for offset in (0, 3):
            px, py, pz = flow[offset : offset + 3] * 1000.0  # NuFlows input is GeV; output is MeV.
            energy = float(np.sqrt(px * px + py * py + pz * pz))
            slot_vectors.append(np.asarray((px, py, pz, energy), dtype=np.float64))
        top_nu = slot_vectors[0 if selected_charges[0] > 0 else 1]
        antitop_nubar = slot_vectors[1 if selected_charges[0] > 0 else 0]
        independent_mtt_mev = invariant_mass(
            [
                constituents["b"][row],
                constituents["lep"][row],
                top_nu,
                constituents["bbar"][row],
                constituents["antilep"][row],
                antitop_nubar,
            ]
        )
        difference = abs(independent_mtt_mev - float(mtt[row]))
        maximum_mtt_difference = max(maximum_mtt_difference, difference)

        component_row = np.asarray(components[row], dtype=np.float64)
        if component_row.shape != (12,) or not np.all(np.isfinite(component_row)):
            raise AssertionError(f"{path}: row {row} has malformed or non-finite top component vector")
        component_difference = abs(float(component_row[11]) - independent_mtt_mev)
        maximum_component_mtt_difference = max(maximum_component_mtt_difference, component_difference)

        expected_kinematics = (
            float(np.hypot(top_nu[0], top_nu[1])),
            float(top_nu[2]),
            float(np.hypot(antitop_nubar[0], antitop_nubar[1])),
            float(antitop_nubar[2]),
        )
        observed_kinematics = tuple(
            float(output_nu[name][row])
            for name in (
                "nu1_pt_VYPER_NUFLOWS_NOSYS",
                "nu1_pz_VYPER_NUFLOWS_NOSYS",
                "nu2_pt_VYPER_NUFLOWS_NOSYS",
                "nu2_pz_VYPER_NUFLOWS_NOSYS",
            )
        )
        kinematic_difference = max(
            abs(observed - expected)
            for observed, expected in zip(observed_kinematics, expected_kinematics)
        )
        maximum_neutrino_kinematic_difference = max(maximum_neutrino_kinematic_difference, kinematic_difference)
        mtt_values_gev.append(float(mtt[row]) / 1000.0)
        selected_rows += 1

    if selected_rows == 0:
        raise AssertionError(f"{path}: no valid VyPER+NuFlows rows were checked")
    if maximum_mtt_difference > tolerance_mev:
        raise AssertionError(
            f"{path}: independent mtt differs from VYPER_NUFLOWS output by {maximum_mtt_difference:g} MeV "
            f"(tolerance {tolerance_mev:g} MeV)"
        )
    if maximum_component_mtt_difference > tolerance_mev:
        raise AssertionError(
            f"{path}: top-components mtt differs from independent sum by "
            f"{maximum_component_mtt_difference:g} MeV (tolerance {tolerance_mev:g} MeV)"
        )
    if maximum_neutrino_kinematic_difference > tolerance_mev:
        raise AssertionError(
            f"{path}: charge-reordered NuFlows kinematics differ by "
            f"{maximum_neutrino_kinematic_difference:g} MeV"
        )

    return {
        "input": str(path),
        "tree": "reco",
        "entries_read": entries,
        "valid_rows_checked": selected_rows,
        "invalid_rows": invalid_rows,
        "assignment_valid_fraction": float(np.mean(assignment_valid)),
        "neutrino_valid_fraction": float(np.mean(neutrino_valid)),
        "reco_valid_fraction": float(np.mean(reco_valid)),
        "failure_code_counts": {str(code): int(np.count_nonzero(failure_code == code)) for code in np.unique(failure_code)},
        "flow_input_units": "GeV",
        "hybrid_kinematic_output_units": "MeV for pT, pz and mtt; radians for angular quantities",
        "neutrino_ordering": "positive charged lepton -> top nu; negative charged lepton -> antitop nubar",
        "maximum_independent_mtt_difference_MeV": maximum_mtt_difference,
        "maximum_top_component_mtt_difference_MeV": maximum_component_mtt_difference,
        "maximum_charge_ordered_neutrino_kinematic_difference_MeV": maximum_neutrino_kinematic_difference,
        "comparison_tolerance_MeV": tolerance_mev,
        "mtt_median_GeV": float(np.median(mtt_values_gev)),
        "mtt_min_GeV": float(np.min(mtt_values_gev)),
        "mtt_max_GeV": float(np.max(mtt_values_gev)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="rich FastFrames pilot ROOT output")
    parser.add_argument("--max-events", type=int, default=1000)
    parser.add_argument("--tolerance-mev", type=float, default=1.0)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = validate(args.root, args.max_events, args.tolerance_mev)
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
