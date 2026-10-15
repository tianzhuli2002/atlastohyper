#!/usr/bin/env python3
"""Validate the VyPER assignment + parent-truth-neutrino ROOT branches."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


def branch(tree, name: str):
    return tree[name].array(library="np")


def p4_from_pt_eta_phi_m(pt, eta, phi, mass):
    pt = np.asarray(pt, dtype=np.float64)
    eta = np.asarray(eta, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    mass = np.asarray(mass, dtype=np.float64)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    energy = np.sqrt(np.maximum(px * px + py * py + pz * pz + mass * mass, 0.0))
    return np.stack((px, py, pz, energy), axis=-1)


def p4_branch(tree, prefix: str):
    return p4_from_pt_eta_phi_m(
        branch(tree, f"{prefix}_fCoordinates_fPt"),
        branch(tree, f"{prefix}_fCoordinates_fEta"),
        branch(tree, f"{prefix}_fCoordinates_fPhi"),
        branch(tree, f"{prefix}_fCoordinates_fM"),
    )


def mass(p4):
    m2 = p4[..., 3] ** 2 - np.sum(p4[..., :3] ** 2, axis=-1)
    return np.sqrt(np.maximum(m2, 0.0))


def max_abs_diff(left, right):
    return float(np.max(np.abs(np.asarray(left) - np.asarray(right))))


def matched_jet_indices(candidate_p4, jet_pt, jet_eta, jet_phi):
    matched = []
    for row, (pts, etas, phis) in enumerate(zip(jet_pt, jet_eta, jet_phi)):
        pt = float(np.hypot(candidate_p4[row, 0], candidate_p4[row, 1]))
        eta = float(np.arcsinh(candidate_p4[row, 2] / pt)) if pt else np.nan
        phi = float(np.arctan2(candidate_p4[row, 1], candidate_p4[row, 0]))
        pts = np.asarray(pts, dtype=np.float64)
        etas = np.asarray(etas, dtype=np.float64)
        phis = np.asarray(phis, dtype=np.float64)
        dphi = (phis - phi + np.pi) % (2.0 * np.pi) - np.pi
        close = (
            np.isclose(pts, pt, rtol=1e-6, atol=0.05)
            & np.isclose(etas, eta, rtol=0.0, atol=1e-6)
            & np.isclose(dphi, 0.0, rtol=0.0, atol=1e-6)
        )
        indices = np.flatnonzero(close)
        matched.append(int(indices[0]) if len(indices) == 1 else -1)
    return matched


def validate(root_path: Path):
    tree = uproot.open(root_path)["reco"]
    entries = tree.num_entries
    event_number = branch(tree, "eventNumber")
    if len(np.unique(event_number)) != entries:
        raise AssertionError("eventNumber is not unique in the pilot output")

    required = [
        "has_native_VYPER_NOSYS",
        "VYPER_NATIVE_assignment_valid_NOSYS",
        "VYPER_NATIVE_neutrino_valid_NOSYS",
        "VYPER_NATIVE_reco_valid_NOSYS",
        "VYPER_NATIVE_failure_code_NOSYS",
        "VYPER_TRUTHNU_assignment_valid_NOSYS",
        "VYPER_TRUTHNU_neutrino_valid_NOSYS",
        "VYPER_TRUTHNU_reco_valid_NOSYS",
        "VYPER_TRUTHNU_failure_code_NOSYS",
        "mttbar_VYPER_TRUTHNU_NOSYS",
    ]
    missing = [name for name in required if name not in tree]
    if missing:
        raise AssertionError(f"Missing required VyPER truth-neutrino branches: {missing}")

    native_b = p4_branch(tree, "TtbarDiLepton_VyPER_Top_b_p4_NOSYS")
    native_lep = p4_branch(tree, "TtbarDiLepton_VyPER_Top_lep_p4_NOSYS")
    native_nu = p4_branch(tree, "TtbarDiLepton_VyPER_Top_nu_p4_NOSYS")
    native_bbar = p4_branch(tree, "TtbarDiLepton_VyPER_Topbar_bbar_p4_NOSYS")
    native_lepbar = p4_branch(tree, "TtbarDiLepton_VyPER_Topbar_lepbar_p4_NOSYS")
    native_nubar = p4_branch(tree, "TtbarDiLepton_VyPER_Topbar_nubar_p4_NOSYS")
    native_top = p4_branch(tree, "TtbarDiLepton_VyPER_Top_NOSYS")
    native_tbar = p4_branch(tree, "TtbarDiLepton_VyPER_Tbar_NOSYS")
    jet_pt = ak.to_list(tree["jet_pt_NOSYS"].array(library="ak"))
    jet_eta = ak.to_list(tree["jet_eta"].array(library="ak"))
    jet_phi = ak.to_list(tree["jet_phi"].array(library="ak"))
    native_b_jet_indices = matched_jet_indices(native_b, jet_pt, jet_eta, jet_phi)
    native_bbar_jet_indices = matched_jet_indices(native_bbar, jet_pt, jet_eta, jet_phi)

    truth_nu = p4_from_pt_eta_phi_m(
        branch(tree, "truth_nu1_pt_NOSYS"),
        branch(tree, "truth_nu1_eta_NOSYS"),
        branch(tree, "truth_nu1_phi_NOSYS"),
        branch(tree, "truth_nu1_m_NOSYS"),
    )
    truth_nubar = p4_from_pt_eta_phi_m(
        branch(tree, "truth_nu2_pt_NOSYS"),
        branch(tree, "truth_nu2_eta_NOSYS"),
        branch(tree, "truth_nu2_phi_NOSYS"),
        branch(tree, "truth_nu2_m_NOSYS"),
    )

    top_native_composed = native_b + native_lep + native_nu
    tbar_native_composed = native_bbar + native_lepbar + native_nubar
    top_truthnu_composed = native_b + native_lep + truth_nu
    tbar_truthnu_composed = native_bbar + native_lepbar + truth_nubar
    ttbar_truthnu_composed = top_truthnu_composed + tbar_truthnu_composed
    independent_mtt = mass(ttbar_truthnu_composed)
    stored_mtt = branch(tree, "mttbar_VYPER_TRUTHNU_NOSYS").astype(np.float64)
    stored_top_pt = branch(tree, "top_pt_VYPER_TRUTHNU_NOSYS")
    stored_tbar_pt = branch(tree, "antitop_pt_VYPER_TRUTHNU_NOSYS")
    stored_top_m = branch(tree, "top_m_VYPER_TRUTHNU_NOSYS")
    stored_tbar_m = branch(tree, "antitop_m_VYPER_TRUTHNU_NOSYS")

    flag_counts = {
        name: int(np.count_nonzero(branch(tree, name) == 1))
        for name in (
            "has_native_VYPER_NOSYS",
            "VYPER_NATIVE_assignment_valid_NOSYS",
            "VYPER_NATIVE_neutrino_valid_NOSYS",
            "VYPER_NATIVE_reco_valid_NOSYS",
            "VYPER_TRUTHNU_assignment_valid_NOSYS",
            "VYPER_TRUTHNU_neutrino_valid_NOSYS",
            "VYPER_TRUTHNU_reco_valid_NOSYS",
        )
    }
    if any(count != entries for count in flag_counts.values()):
        raise AssertionError(f"Pilot validity flags are not all valid: {flag_counts}")
    failure_codes = np.unique(branch(tree, "VYPER_TRUTHNU_failure_code_NOSYS"))
    if not np.array_equal(failure_codes, np.array([0])):
        raise AssertionError(f"Unexpected VyPER truth-neutrino failure codes: {failure_codes}")

    # Native component and lepton-charge conventions are checked directly
    # against their summed native top vectors and selected object branches.
    native_top_sum_diff = max_abs_diff(top_native_composed, native_top)
    native_tbar_sum_diff = max_abs_diff(tbar_native_composed, native_tbar)
    if native_top_sum_diff > 1.0 or native_tbar_sum_diff > 1.0:
        raise AssertionError(
            "Native VyPER top vectors do not close to their constituent vectors "
            f"within 1 MeV: top={native_top_sum_diff}, antitop={native_tbar_sum_diff}"
        )

    positive_lepton_match_rows = []
    negative_lepton_match_rows = []
    for row in range(entries):
        matches_positive = []
        matches_negative = []
        for flavour in ("el", "mu"):
            pts = ak.to_list(tree[f"{flavour}_pt_NOSYS"].array(library="ak")[row])
            etas = ak.to_list(tree[f"{flavour}_eta"].array(library="ak")[row])
            phis = ak.to_list(tree[f"{flavour}_phi"].array(library="ak")[row])
            charges = ak.to_list(tree[f"{flavour}_charge"].array(library="ak")[row])
            for idx, charge in enumerate(charges):
                p4 = native_lep[row] if charge > 0 else native_lepbar[row]
                dphi = (phis[idx] - np.arctan2(p4[1], p4[0]) + np.pi) % (2 * np.pi) - np.pi
                pt = np.hypot(p4[0], p4[1])
                eta = np.arcsinh(p4[2] / pt) if pt else np.nan
                is_match = (
                    np.isclose(pts[idx], pt, rtol=1e-6, atol=0.05)
                    and np.isclose(etas[idx], eta, rtol=0.0, atol=1e-6)
                    and np.isclose(dphi, 0.0, rtol=0.0, atol=1e-6)
                )
                if is_match:
                    (matches_positive if charge > 0 else matches_negative).append(row)
        if len(matches_positive) == 1:
            positive_lepton_match_rows.append(row)
        if len(matches_negative) == 1:
            negative_lepton_match_rows.append(row)
    if len(positive_lepton_match_rows) != entries or len(negative_lepton_match_rows) != entries:
        raise AssertionError("VyPER Top/Topbar leptons are not positive/negative-charge ordered")
    if any(index < 0 for index in native_b_jet_indices + native_bbar_jet_indices):
        raise AssertionError("VyPER b constituent vectors do not uniquely match selected jets")
    if any(b == bbar for b, bbar in zip(native_b_jet_indices, native_bbar_jet_indices)):
        raise AssertionError("VyPER b and bbar constituents map to the same selected jet")

    max_mtt_diff = max_abs_diff(stored_mtt, independent_mtt)
    if not np.allclose(stored_mtt, independent_mtt, rtol=1e-6, atol=1.0):
        raise AssertionError(f"Stored and independent mtt differ by {max_mtt_diff} MeV")
    top_pt_diff = max_abs_diff(stored_top_pt, np.hypot(top_truthnu_composed[:, 0], top_truthnu_composed[:, 1]))
    tbar_pt_diff = max_abs_diff(stored_tbar_pt, np.hypot(tbar_truthnu_composed[:, 0], tbar_truthnu_composed[:, 1]))
    top_m_diff = max_abs_diff(stored_top_m, mass(top_truthnu_composed))
    tbar_m_diff = max_abs_diff(stored_tbar_m, mass(tbar_truthnu_composed))
    if max(top_pt_diff, tbar_pt_diff, top_m_diff, tbar_m_diff) > 1.0:
        raise AssertionError("Stored top components disagree with independent four-vector sums")

    return {
        "input": str(root_path),
        "tree": "reco",
        "events": entries,
        "unique_eventNumber": len(np.unique(event_number)),
        "units": "MeV for four-vectors, mtt, and top kinematics; betaStar is dimensionless",
        "truth_neutrino_convention": "truth_nu1 = from top; truth_nu2 = from antitop",
        "native_vyper_constituents": "Top=(b, positive lepton, nu), Tbar=(bbar, negative lepton, nubar)",
        "validity_counts": flag_counts,
        "failure_codes": failure_codes.tolist(),
        "positive_lepton_charge_matches": len(positive_lepton_match_rows),
        "negative_lepton_charge_matches": len(negative_lepton_match_rows),
        "top_b_unique_selected_jet_matches": sum(index >= 0 for index in native_b_jet_indices),
        "antitop_bbar_unique_selected_jet_matches": sum(index >= 0 for index in native_bbar_jet_indices),
        "distinct_b_and_bbar_jet_matches": sum(
            b != bbar for b, bbar in zip(native_b_jet_indices, native_bbar_jet_indices)
        ),
        "max_native_top_constituent_closure_MeV": native_top_sum_diff,
        "max_native_antitop_constituent_closure_MeV": native_tbar_sum_diff,
        "max_independent_mtt_difference_MeV": max_mtt_diff,
        "max_independent_top_pt_difference_MeV": top_pt_diff,
        "max_independent_antitop_pt_difference_MeV": tbar_pt_diff,
        "max_independent_top_mass_difference_MeV": top_m_diff,
        "max_independent_antitop_mass_difference_MeV": tbar_m_diff,
        "mtt_VYPER_TRUTHNU_median_GeV": float(np.median(stored_mtt) / 1000.0),
        "mtt_VYPER_NATIVE_constituent_sum_median_GeV": float(
            np.median(mass(native_top + native_tbar)) / 1000.0
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path, help="rich FastFrames reconstruction ROOT output")
    parser.add_argument("--report", type=Path, help="optional JSON report output")
    args = parser.parse_args()
    report = validate(args.root)
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
