#!/usr/bin/env python3
"""Export conservative pilot metrics for a rich dilepton comparison ROOT file."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Plot.reco_plot.io import RecoData
from Plot.reco_plot.units import display_values


ASSIGNMENT_BRANCHES = {
    "HYPER_adapter": "HYPER_correct_NOSYS",
    "CHI2": "CHI2_correct_NOSYS",
    "DELTAR": "DELTAR_correct_NOSYS",
    "PT": "PT_correct_NOSYS",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--region", default="ttll_analysis")
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()

    with uproot.open(args.root) as root_file:
        tree = root_file["reco"]
        entries = int(tree.num_entries)
        event_numbers = tree["eventNumber"].array(library="np").astype(np.uint64)
        if len(np.unique(event_numbers)) != len(event_numbers):
            raise AssertionError("eventNumber is not unique in pilot ROOT output")
        assignment = {}
        truth_assignment_valid_events = None
        if "mttbar_truth_assignment_valid_NOSYS" in tree:
            truth_validity = np.asarray(
                tree["mttbar_truth_assignment_valid_NOSYS"].array(library="np")
            )
            truth_assignment_valid_events = int(np.count_nonzero(truth_validity == 1))
        for method, branch in ASSIGNMENT_BRANCHES.items():
            if branch not in tree:
                assignment[method] = {"available": False, "branch": branch}
                continue
            values = np.asarray(tree[branch].array(library="np"), dtype=float)
            valid = np.isin(values, (0.0, 1.0))
            denominator = int(np.count_nonzero(valid))
            correct = int(np.count_nonzero(values == 1.0))
            assignment[method] = {
                "available": True,
                "branch": branch,
                "analysis_selected_events": entries,
                "valid_assignment_events": denominator,
                "correct": correct,
                "incorrect": int(np.count_nonzero(values == 0.0)),
                "invalid_or_unmatched": int(entries - denominator),
                "efficiency_on_valid_assignment_events": correct / denominator if denominator else None,
            }

    data = RecoData(str(args.root), region=args.region)
    reco_raw = data.vyper_array("mttbar", "pred", source="VYPER_NATIVE")
    truth_raw = data.truth_array("mttbar", "partonic_after_fsr")
    if reco_raw is None or truth_raw is None:
        raise RuntimeError("pilot ROOT lacks the native VyPER or after-FSR truth input")
    reco, reco_unit, reco_display, _ = display_values(reco_raw, "mttbar")
    truth, truth_unit, truth_display, _ = display_values(truth_raw, "mttbar")
    if len(reco) != entries or len(truth) != entries:
        raise AssertionError("native VyPER, truth, and event counts differ")
    valid = np.isfinite(reco) & np.isfinite(truth) & (truth > 0.0)
    residual = reco[valid] - truth[valid]
    relative = residual / truth[valid]
    if not valid.any():
        raise AssertionError("no valid native VyPER / after-FSR truth pairs")

    mtt = {
        "analysis_selected_events": entries,
        "valid_native_vyper_and_truth_pairs": int(np.count_nonzero(valid)),
        "invalid_pairs": int(entries - np.count_nonzero(valid)),
        "reco_unit": reco_unit,
        "truth_unit": truth_unit,
        "display_unit": reco_display or truth_display,
        "reco_median_GeV": float(np.median(reco[valid])),
        "truth_median_GeV": float(np.median(truth[valid])),
        "absolute_residual_GeV": {
            "median": float(np.median(residual)),
            "p16": float(np.percentile(residual, 16)),
            "p84": float(np.percentile(residual, 84)),
            "central_68_width": float(np.percentile(residual, 84) - np.percentile(residual, 16)),
            "rms": float(np.sqrt(np.mean(residual**2))),
        },
        "relative_residual": {
            "median": float(np.median(relative)),
            "p16": float(np.percentile(relative, 16)),
            "p84": float(np.percentile(relative, 84)),
            "rms": float(np.sqrt(np.mean(relative**2))),
            "fraction_abs_gt_0p2": float(np.mean(np.abs(relative) > 0.2)),
        },
        "truth_definition": "partonic top and antitop after FSR from Ttbar_MC_* branches",
    }
    result = {
        "sample": "KK candidate DSID 574797, mc20a, fullsim; mass/benchmark not AMI-verified",
        "root": str(args.root),
        "region": args.region,
        "event_identity": "eventNumber (unique in this single-DSID pilot)",
        "population": {
            "analysis_selected_events": entries,
            "truth_assignment_valid_events": truth_assignment_valid_events,
        },
        "assignment": assignment,
        "native_vyper_mttbar": mtt,
        "interpretation_limit": "Single low-statistics pilot file; not a resonance-mass result or checkpoint zero-shot evaluation.",
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n")
    rows = []
    for method, values in assignment.items():
        rows.append({"category": "assignment", "method": method, **values})
    rows.extend(
        {"category": f"mttbar_{group}", "method": "VYPER_NATIVE", "metric": name, "value": value}
        for group, metrics in mtt.items()
        if isinstance(metrics, dict)
        for name, value in metrics.items()
        if isinstance(value, (int, float))
    )
    with args.csv.open("w", newline="") as stream:
        fieldnames = sorted({key for row in rows for key in row})
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
