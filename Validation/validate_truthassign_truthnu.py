#!/usr/bin/env python3
"""Independently close the stored truth-assignment + truth-neutrino mtt oracle."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


def p4(pt, eta, phi, energy=None):
    pt = float(pt)
    eta = float(eta)
    phi = float(phi)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    e = pt * np.cosh(eta) if energy is None else float(energy)
    return np.array([px, py, pz, e], dtype=np.float64)


def invariant_mass(vector):
    mass_squared = vector[3] ** 2 - np.dot(vector[:3], vector[:3])
    return float(np.sqrt(max(mass_squared, 0.0)))


def array(tree, name):
    return tree[name].array(library="ak")


def validate(path: Path, max_events: int | None = None):
    with uproot.open(path) as root_file:
        if "reco" not in root_file:
            raise KeyError(f"{path}: required tree 'reco' is absent")
        tree = root_file["reco"]
        required = (
            "eventNumber", "mttbar_truth_assignment_valid_NOSYS",
            "mttbar_truth_bjet_index_NOSYS", "mttbar_truth_bbarjet_index_NOSYS",
            "mttbar_TRUTH_NOSYS", "jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS",
            "el_pt_NOSYS", "el_eta", "el_phi", "el_e_NOSYS", "el_charge", "el_select_tight_NOSYS",
            "mu_pt_NOSYS", "mu_eta", "mu_phi", "mu_e_NOSYS", "mu_charge", "mu_select_tight_NOSYS",
            *(f"{flavour}_dRmatched_to_Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
              for flavour in ("el", "mu") for daughter in (1, 2) for side in ("t", "tbar")),
            *(f"n_lep_matched_to_Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
              for daughter in (1, 2) for side in ("t", "tbar")),
            "truth_nu1_pt_NOSYS", "truth_nu1_eta_NOSYS", "truth_nu1_phi_NOSYS",
            "truth_nu2_pt_NOSYS", "truth_nu2_eta_NOSYS", "truth_nu2_phi_NOSYS",
            "Ttbar_MC_ttbar_afterFSR_m",
        )
        missing = [name for name in required if name not in tree]
        if missing:
            raise KeyError(f"{path}: missing oracle validation branches: {missing}")
        stop = min(tree.num_entries, max_events) if max_events is not None else tree.num_entries
        values = {name: array(tree, name)[:stop] for name in required}

    event_numbers = ak.to_numpy(values["eventNumber"])
    stored_valid = ak.to_numpy(values["mttbar_truth_assignment_valid_NOSYS"]) == 1
    b_indices = ak.to_numpy(values["mttbar_truth_bjet_index_NOSYS"]).astype(int)
    bbar_indices = ak.to_numpy(values["mttbar_truth_bbarjet_index_NOSYS"]).astype(int)
    stored_mtt = ak.to_numpy(values["mttbar_TRUTH_NOSYS"]).astype(float)
    jet_pt = ak.to_list(values["jet_pt_NOSYS"])
    jet_eta = ak.to_list(values["jet_eta"])
    jet_phi = ak.to_list(values["jet_phi"])
    jet_e = ak.to_list(values["jet_e_NOSYS"])

    leptons_by_charge = []
    for i in range(stop):
        selected = []
        for flavour in ("el", "mu"):
            pts = ak.to_list(values[f"{flavour}_pt_NOSYS"][i])
            etas = ak.to_list(values[f"{flavour}_eta"][i])
            phis = ak.to_list(values[f"{flavour}_phi"][i])
            energies = ak.to_list(values[f"{flavour}_e_NOSYS"][i])
            charges = ak.to_list(values[f"{flavour}_charge"][i])
            tight = ak.to_list(values[f"{flavour}_select_tight_NOSYS"][i])
            for j, keep in enumerate(tight):
                if keep:
                    selected.append((float(pts[j]), float(charges[j]),
                                     p4(pts[j], etas[j], phis[j], energies[j])))
        selected.sort(key=lambda row: -row[0])
        positive = [row[2] for row in selected if row[1] > 0]
        negative = [row[2] for row in selected if row[1] < 0]
        leptons_by_charge.append((positive, negative))

    nu1_pt = ak.to_numpy(values["truth_nu1_pt_NOSYS"]).astype(float)
    nu1_eta = ak.to_numpy(values["truth_nu1_eta_NOSYS"]).astype(float)
    nu1_phi = ak.to_numpy(values["truth_nu1_phi_NOSYS"]).astype(float)
    nu2_pt = ak.to_numpy(values["truth_nu2_pt_NOSYS"]).astype(float)
    nu2_eta = ak.to_numpy(values["truth_nu2_eta_NOSYS"]).astype(float)
    nu2_phi = ak.to_numpy(values["truth_nu2_phi_NOSYS"]).astype(float)

    independently_valid = np.zeros(stop, dtype=bool)
    calculated = np.full(stop, np.nan)
    selected_lepton_count = np.zeros(stop, dtype=int)
    distinct_jets = np.zeros(stop, dtype=bool)
    truth_lepton_pair_matched = np.zeros(stop, dtype=bool)
    truth_parent_mtt = ak.to_numpy(values["Ttbar_MC_ttbar_afterFSR_m"]).astype(float)
    for i in range(stop):
        selected_ids = {}
        selected_charges = []
        for flavour in ("el", "mu"):
            charges = ak.to_list(values[f"{flavour}_charge"][i])
            tight = ak.to_list(values[f"{flavour}_select_tight_NOSYS"][i])
            for index, (charge, keep) in enumerate(zip(charges, tight)):
                if keep:
                    selected_ids[(flavour, index)] = float(charge)
                    selected_charges.append(float(charge))
        good = (
            len(selected_ids) == 2
            and sum(charge > 0.0 for charge in selected_charges) == 1
            and sum(charge < 0.0 for charge in selected_charges) == 1
        )
        matched_ids = []
        for side, expected_charge in (("t", 1), ("tbar", -1)):
            side_ids = []
            for daughter in (1, 2):
                count = int(ak.to_numpy(values[
                    f"n_lep_matched_to_Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
                ])[i])
                candidates = []
                for flavour in ("el", "mu"):
                    index = int(ak.to_numpy(values[
                        f"{flavour}_dRmatched_to_Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
                    ])[i])
                    if index >= 0:
                        candidates.append((flavour, index))
                if count == 0 and not candidates:
                    continue
                if count != 1 or len(candidates) != 1:
                    good = False
                    continue
                side_ids.extend(candidates)
            if len(side_ids) != 1:
                good = False
                continue
            flavour, index = side_ids[0]
            charge = selected_ids.get((flavour, index))
            if charge is None or np.sign(charge) != expected_charge:
                good = False
            matched_ids.append((flavour, index))
        if len(set(matched_ids)) != 2 or set(matched_ids) != set(selected_ids):
            good = False
        truth_lepton_pair_matched[i] = good

    for i in range(stop):
        positive, negative = leptons_by_charge[i]
        selected_lepton_count[i] = len(positive) + len(negative)
        b, bbar = b_indices[i], bbar_indices[i]
        in_range = 0 <= b < len(jet_pt[i]) and 0 <= bbar < len(jet_pt[i])
        distinct_jets[i] = in_range and b != bbar
        valid = (
            stored_valid[i] and distinct_jets[i] and len(positive) == 1 and len(negative) == 1
            and np.all(np.isfinite([nu1_pt[i], nu1_eta[i], nu1_phi[i], nu2_pt[i], nu2_eta[i], nu2_phi[i]]))
        )
        if not valid:
            continue
        b_p4 = p4(jet_pt[i][b], jet_eta[i][b], jet_phi[i][b], jet_e[i][b])
        bbar_p4 = p4(jet_pt[i][bbar], jet_eta[i][bbar], jet_phi[i][bbar], jet_e[i][bbar])
        nu = p4(nu1_pt[i], nu1_eta[i], nu1_phi[i])
        nubar = p4(nu2_pt[i], nu2_eta[i], nu2_phi[i])
        total = b_p4 + bbar_p4 + positive[0] + negative[0] + nu + nubar
        calculated[i] = invariant_mass(total)
        independently_valid[i] = np.isfinite(calculated[i]) and calculated[i] > 0.0

    if np.any(stored_valid & (selected_lepton_count != 2)):
        bad = np.flatnonzero(stored_valid & (selected_lepton_count != 2))[:5]
        raise AssertionError(f"Oracle-valid events do not contain two tight selected leptons: rows {bad.tolist()}")
    if np.any(stored_valid & ~distinct_jets):
        bad = np.flatnonzero(stored_valid & ~distinct_jets)[:5]
        raise AssertionError(f"Oracle-valid events do not have two distinct matched jets: rows {bad.tolist()}")
    if not np.array_equal(stored_valid, independently_valid):
        raise AssertionError(
            "Stored and independently reconstructed oracle validity differs: "
            f"stored={int(stored_valid.sum())}, independent={int(independently_valid.sum())}"
        )
    truth_neutrinos_valid = (
        np.isfinite(nu1_pt) & (nu1_pt > 0.0)
        & np.isfinite(nu1_eta) & np.isfinite(nu1_phi)
        & np.isfinite(nu2_pt) & (nu2_pt > 0.0)
        & np.isfinite(nu2_eta) & np.isfinite(nu2_phi)
    )
    partonic_truth_valid = np.isfinite(truth_parent_mtt) & (truth_parent_mtt > 0.0)
    study_truth_reconstructible = (
        independently_valid & truth_lepton_pair_matched
        & truth_neutrinos_valid & partonic_truth_valid
    )
    residual = calculated[independently_valid] - stored_mtt[independently_valid]
    if residual.size == 0:
        raise AssertionError(f"{path}: no truth-reconstructible events found in inspected entries")
    summary = {
        "input": str(path.resolve()),
        "tree": "reco",
        "entries_inspected": int(stop),
        "unique_event_numbers": int(np.unique(event_numbers).size),
        "truth_assignment_oracle_valid": int(independently_valid.sum()),
        "truth_lepton_pair_reconstructible": int(truth_lepton_pair_matched.sum()),
        "truth_reconstructible": int(study_truth_reconstructible.sum()),
        "invalid_oracle_rows": int((~independently_valid).sum()),
        "truth_b_mapping": "mttbar_truth_bjet_index_NOSYS and bbar counterpart index raw jet_TLV_NOSYS ordering",
        "leptons": "exactly two el/mu tight-selected objects, charge ordered positive then negative; each selected object uniquely matches one W daughter on its corresponding parent side, with positive charge for top and negative for antitop; input four-vectors in MeV",
        "neutrinos": "truth_nu1 from parent top and truth_nu2 from parent antitop, massless PtEtaPhi four-vectors in MeV",
        "oracle_definition": "matched reconstructed b+bbar jets + selected positive/negative leptons + parent-truth neutrinos",
        "stored_branch_unit": "MeV",
        "closure": {
            "mean_difference_MeV": float(np.mean(residual)),
            "median_abs_difference_MeV": float(np.median(np.abs(residual))),
            "max_abs_difference_MeV": float(np.max(np.abs(residual))),
            "allclose_rtol_1e-6_atol_1MeV": bool(np.allclose(
                calculated[independently_valid], stored_mtt[independently_valid], rtol=1e-6, atol=1.0
            )),
        },
    }
    if not summary["closure"]["allclose_rtol_1e-6_atol_1MeV"]:
        raise AssertionError(f"Truth-assignment oracle failed closure: {summary['closure']}")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rootfile", type=Path)
    parser.add_argument("--max-events", type=int)
    parser.add_argument("--json", type=Path, help="Optional path for a JSON report")
    args = parser.parse_args()
    result = validate(args.rootfile, args.max_events)
    payload = json.dumps(result, indent=2, allow_nan=False)
    print(payload)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload + "\n")


if __name__ == "__main__":
    main()
