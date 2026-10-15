#!/usr/bin/env python3
"""Audit raw truth-mtt support in the nominal MC20a SM ttbar ntuple set.

This reads only the after-FSR truth-mtt and W-daughter PDG-ID branches. It is a
raw generated-event support audit, not the selected/reconstructible method
population used in the reconstruction study.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import uproot


DEFAULT_INPUT_DIR = Path(
    "/data/leet/ntuples/run2_ttbar_dilep/"
    "user.leet.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p7266."
    "TtNT_2L_2Lsys_b71171e_output"
)
TREE_NAME = "reco"
MTT_BRANCH = "Ttbar_MC_ttbar_afterFSR_m"
MASS_EDGES_GEV = np.asarray((0, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 10000), dtype=float)
THRESHOLDS_GEV = (1000, 1600, 2000, 3000, 4000, 5000, 6000)


def audit(input_dir: Path) -> dict:
    if not input_dir.is_dir():
        raise FileNotFoundError(f"SM ntuple directory does not exist: {input_dir}")
    files = sorted(input_dir.glob("*.root"))
    if not files:
        raise FileNotFoundError(f"No ROOT files in exact SM ntuple directory: {input_dir}")
    branch_names = [MTT_BRANCH]
    for side in ("t", "tbar"):
        branch_names.extend(
            f"Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}_pdgId"
            for daughter in (1, 2)
        )

    total_entries = 0
    finite_entries = 0
    direct_entries = 0
    max_mtt = None
    max_direct_mtt = None
    hist_all = np.zeros(len(MASS_EDGES_GEV) - 1, dtype=np.int64)
    hist_direct = np.zeros_like(hist_all)
    above_all = {str(threshold): 0 for threshold in THRESHOLDS_GEV}
    above_direct = {str(threshold): 0 for threshold in THRESHOLDS_GEV}

    for index, path in enumerate(files, start=1):
        with uproot.open(path) as root_file:
            if TREE_NAME not in root_file:
                raise KeyError(f"{path}: required tree {TREE_NAME!r} is absent")
            tree = root_file[TREE_NAME]
            missing = sorted(set(branch_names) - set(tree.keys()))
            if missing:
                raise KeyError(f"{path}: missing required support-audit branches: {missing}")
            arrays = tree.arrays(branch_names, library="np")
            mtt = np.asarray(arrays[MTT_BRANCH], dtype=float) * 0.001
            if mtt.ndim != 1 or mtt.size != tree.num_entries:
                raise ValueError(f"{path}: {MTT_BRANCH} must be one scalar per raw tree entry")
            direct = np.isfinite(mtt) & (mtt > 0)
            for side, sign in (("t", -1), ("tbar", 1)):
                charged = np.zeros(mtt.size, dtype=bool)
                for daughter in (1, 2):
                    pdg = np.asarray(
                        arrays[f"Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}_pdgId"], dtype=int
                    )
                    if pdg.shape != mtt.shape:
                        raise ValueError(f"{path}: W daughter PDG-ID branch shape differs from event mtt")
                    charged |= (np.sign(pdg) == sign) & np.isin(np.abs(pdg), (11, 13))
                direct &= charged

        finite = np.isfinite(mtt) & (mtt > 0)
        total_entries += int(mtt.size)
        finite_entries += int(np.count_nonzero(finite))
        direct_entries += int(np.count_nonzero(direct))
        if np.any(finite):
            file_max = float(np.max(mtt[finite]))
            max_mtt = file_max if max_mtt is None else max(max_mtt, file_max)
        if np.any(direct):
            file_max = float(np.max(mtt[direct]))
            max_direct_mtt = file_max if max_direct_mtt is None else max(max_direct_mtt, file_max)
        hist_all += np.histogram(mtt[finite], bins=MASS_EDGES_GEV)[0]
        hist_direct += np.histogram(mtt[direct], bins=MASS_EDGES_GEV)[0]
        for threshold in THRESHOLDS_GEV:
            above_all[str(threshold)] += int(np.count_nonzero(finite & (mtt > threshold)))
            above_direct[str(threshold)] += int(np.count_nonzero(direct & (mtt > threshold)))
        if index % 100 == 0 or index == len(files):
            print(f"Scanned {index}/{len(files)} SM ROOT files", flush=True)

    return {
        "dataset": ("MC20a nominal Powheg+Pythia8 ttbar, DSID 410472, r13167"
                    if input_dir.resolve() == DEFAULT_INPUT_DIR.resolve() else input_dir.name),
        "scope": "one exact local ntuple directory; raw generated events before analysis selection",
        "input_directory": str(input_dir.resolve()),
        "file_count": len(files),
        "tree": TREE_NAME,
        "truth_mtt_branch": MTT_BRANCH,
        "truth_mtt_units_in_file": "MeV; converted to GeV for counts",
        "direct_2l_definition": "one after-FSR e/mu charged daughter on each top side using daughter PDG IDs",
        "raw_tree_entries": total_entries,
        "finite_positive_mtt_entries": finite_entries,
        "direct_2l_entries": direct_entries,
        "max_finite_mtt_GeV": max_mtt,
        "max_direct_2l_mtt_GeV": max_direct_mtt,
        "mtt_bin_edges_GeV": MASS_EDGES_GEV.tolist(),
        "all_event_counts_by_mtt_bin": hist_all.tolist(),
        "direct_2l_counts_by_mtt_bin": hist_direct.tolist(),
        "thresholds_GeV": THRESHOLDS_GEV,
        "all_event_counts_above_threshold": above_all,
        "direct_2l_counts_above_threshold": above_direct,
        "warning": "Counts are raw support only; analysis-selected, truth-reconstructible and candidate-supported counts are lower.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output", type=Path, default=Path("Validation/high_mass_sm_support_r13167.json"))
    args = parser.parse_args()
    result = audit(args.input_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
