#!/usr/bin/env python3
"""Audit rich ROOT files before running the expensive VyPER study suite."""

import argparse
import json
from pathlib import Path

import numpy as np
import uproot


LEPTON_MATCH = [
    f"{flavour}_dRmatched_to_Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
    for side in ("t", "tbar") for daughter in (1, 2) for flavour in ("el", "mu")
] + [
    f"n_lep_matched_to_Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
    for side in ("t", "tbar") for daughter in (1, 2)
]
REQUIRED = {
    "identity": ("mcChannelNumber", "eventNumber"),
    "truth_and_matching": (
        "Ttbar_MC_t_afterFSR_pt", "Ttbar_MC_tbar_afterFSR_pt",
        "Ttbar_MC_ttbar_afterFSR_m", "Truth_had_b_index", "Truth_had_bbar_index",
        "mttbar_truth_assignment_valid_NOSYS", "mttbar_truth_bjet_index_NOSYS",
        "mttbar_truth_bbarjet_index_NOSYS", "mttbar_TRUTH_NOSYS",
        "truth_nu1_pz_NOSYS", "truth_nu2_pz_NOSYS",
    ),
    "objects_and_candidates": (
        "jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS",
        "jet_select_baselineJvt_NOSYS", "el_select_tight_NOSYS",
        "mu_select_tight_NOSYS", "n_jets_NOSYS",
    ),
    "native_and_mechanism": (
        "VYPER_NATIVE_reco_valid_NOSYS", "TtbarDiLepton_VyPER_Top_NOSYS",
        "TtbarDiLepton_VyPER_Tbar_NOSYS",
        "TtbarDiLepton_VyPER_Top_Indices_NOSYS",
        "TtbarDiLepton_VyPER_AntiTop_Indices_NOSYS",
        "nuflows_nu1_pz_NOSYS", "nuflows_nu2_pz_NOSYS",
        "mttbar_VYPER_NUFLOWS_NOSYS", "mttbar_VYPER_TRUTHNU_NOSYS",
    ),
    "strict_lepton_match": LEPTON_MATCH,
}


def audit_file(path):
    path = Path(path)
    file = uproot.open(path)
    if "reco" not in file:
        raise KeyError(f"{path}: reco tree is absent")
    tree = file["reco"]
    available = set(tree.keys())
    groups = {group: {"present": len(set(branches) & available),
                      "required": len(branches),
                      "missing": sorted(set(branches) - available)}
              for group, branches in REQUIRED.items()}
    identities = None
    if not groups["identity"]["missing"]:
        arrays = tree.arrays(REQUIRED["identity"], library="np")
        keys = list(zip(arrays["mcChannelNumber"].astype(int),
                        arrays["eventNumber"].astype(int)))
        identities = {"n_rows": len(keys), "n_unique_channel_event": len(set(keys)),
                      "duplicate_channel_event": len(keys) - len(set(keys)),
                      "n_dsid": int(np.unique(arrays["mcChannelNumber"]).size)}
    usable = (tree.num_entries > 0 and identities is not None
              and identities["duplicate_channel_event"] == 0
              and identities["n_dsid"] == 1
              and all(not result["missing"] for group, result in groups.items()
                      if group != "strict_lepton_match"))
    return {"path": str(path.resolve()), "entries": tree.num_entries,
            "groups": groups, "identity": identities,
            "usable_core_schema": usable,
            "strict_lepton_match_present": not groups["strict_lepton_match"]["missing"],
            "sidecar_required": usable and bool(groups["strict_lepton_match"]["missing"])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rootfiles", nargs="+")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = [audit_file(path) for path in args.rootfiles]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"inputs": rows}, indent=2) + "\n")
    for row in rows:
        print(row["path"], "entries", row["entries"],
              "core", row["usable_core_schema"],
              "lepton_match", row["strict_lepton_match_present"],
              "sidecar_required", row["sidecar_required"])
    print(output)


if __name__ == "__main__":
    main()
