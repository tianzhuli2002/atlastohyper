#!/usr/bin/env python3
"""Recompute matched-domain summaries and plots from strict event caches."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from Plot.reco_plot.bias_diagnostics import (
    ASSIGNMENT_CATEGORY_CODE,
    _matched_2d_rows,
    _matched_3d_nu_rows,
    _matched_mtt_only_rows,
)
from Plot.reco_plot.method_study import _write_machine_csv, render_method_study


def _load_cached_samples(cache_dir):
    inverse_categories = {value: key for key, value in ASSIGNMENT_CATEGORY_CODE.items()}
    samples = []
    seen_cache_names = set()
    for sidecar_path in sorted(Path(cache_dir).glob("*.json")):
        metadata = json.loads(sidecar_path.read_text())
        cache_path = Path(metadata["cache_path"])
        if not cache_path.exists():
            cache_path = sidecar_path.with_suffix(".npz")
        if cache_path.name in seen_cache_names:
            raise ValueError(f"duplicate event cache for {cache_path.name}")
        seen_cache_names.add(cache_path.name)
        with np.load(cache_path, allow_pickle=False) as cache:
            n = len(cache["eventNumber"])
            if n != int(metadata["n_events"]):
                raise ValueError(f"{cache_path}: cached rows disagree with sidecar n_events")
            identities = set(zip(cache["mcChannelNumber"].tolist(), cache["eventNumber"].tolist()))
            if len(identities) != n:
                raise ValueError(f"{cache_path}: duplicate (mcChannelNumber,eventNumber) identities")
            arrays = {
                "eventNumber": cache["eventNumber"],
                "mcChannelNumber": cache["mcChannelNumber"],
                "reco": {"PARTONIC_TRUTH": {}},
                "assignment_variables": {
                    "truth_ttbar_pt": cache["truth_variable_truth_ttbar_pt"],
                    "truth_nu_pair_energy": cache["truth_variable_truth_nu_pair_energy"],
                },
                "method_valid": {},
                "assignments": {},
                "categories": {},
            }
            for observable in ("mttbar", "top_pt"):
                arrays["reco"]["PARTONIC_TRUTH"][observable] = cache[f"reco_PARTONIC_TRUTH_{observable}"]
            for method in ("VYPER_NATIVE", "EM", "VYPER_NUFLOWS", "TRUTHASSIGN_TRUTHNU"):
                arrays["reco"][method] = {
                    observable: cache[f"reco_{method}_{observable}"]
                    for observable in ("mttbar", "top_pt")
                }
            for method in ("VYPER", "HYPER", "EM", "NW", "DA"):
                arrays["assignments"][method] = cache[f"assignment_{method}_indices"]
                codes = cache[f"assignment_{method}_category"]
                try:
                    arrays["categories"][method] = np.asarray(
                        [inverse_categories[int(code)] for code in codes]
                    )
                except KeyError as exc:
                    raise ValueError(f"{cache_path}: unknown assignment category code {exc.args[0]}") from exc
            for method in ("VYPER_NATIVE", "EM", "NW", "VYPER_NUFLOWS", "TRUTHASSIGN_TRUTHNU"):
                arrays["method_valid"][method] = cache[f"method_valid_{method}"]
        samples.append({
            "sample": {
                "label": metadata["sample"], "domain": metadata["domain"],
                "dsid": metadata["dsid"], "mass_GeV": metadata.get("nominal_MG"),
                "role": metadata.get("role", "representative"),
            },
            "_population": np.ones(n, dtype=bool),
            "_arrays": arrays,
        })
    if not samples:
        raise ValueError(f"no cache sidecars found under {cache_dir}")
    return samples


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", required=True, help="Directory containing NPZ caches and JSON sidecars")
    parser.add_argument("--summary-json", required=True, help="Existing method_study_metrics.json to update")
    parser.add_argument("--outdir", required=True, help="Study output directory for summary CSVs and figures")
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    args = parser.parse_args(argv)
    if args.bootstrap_resamples < 1:
        parser.error("--bootstrap-resamples must be at least 1")

    summary_path = Path(args.summary_json)
    outdir = Path(args.outdir)
    summary = json.loads(summary_path.read_text())
    samples = _load_cached_samples(args.cache_dir)
    cells_1d, rows_1d = _matched_mtt_only_rows(
        samples, minimum=50, bootstrap_resamples=args.bootstrap_resamples,
    )
    cells_2d, rows_2d = _matched_2d_rows(
        samples, minimum=50, bootstrap_resamples=args.bootstrap_resamples,
    )
    cells_3d, rows_3d = _matched_3d_nu_rows(
        samples, minimum=50, bootstrap_resamples=args.bootstrap_resamples,
    )
    summary.update({
        "matched_mtt_only_cells": cells_1d,
        "matched_mtt_only_response": rows_1d,
        "matched_mtt_ttpt_cells": cells_2d,
        "matched_mtt_ttpt_response": rows_2d,
        "matched_mtt_ttpt_nupair_cells": cells_3d,
        "matched_mtt_ttpt_nupair_response": rows_3d,
        "matched_domain_final_bootstrap_resamples": int(args.bootstrap_resamples),
    })
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    _write_machine_csv(summary_path.parent / "matched_mtt_only_cells.csv", cells_1d)
    _write_machine_csv(summary_path.parent / "matched_mtt_only_response.csv", rows_1d)
    _write_machine_csv(summary_path.parent / "matched_mtt_ttpt_cells.csv", cells_2d)
    _write_machine_csv(summary_path.parent / "matched_mtt_ttpt_response.csv", rows_2d)
    _write_machine_csv(summary_path.parent / "matched_mtt_ttpt_nupair_cells.csv", cells_3d)
    _write_machine_csv(summary_path.parent / "matched_mtt_ttpt_nupair_response.csv", rows_3d)

    render_method_study(summary, outdir)
    summary["figures"] = sorted(str(path.relative_to(outdir)) for path in outdir.rglob("*.png"))
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(f"Updated matched summaries from {len(samples)} strict event caches with "
          f"{args.bootstrap_resamples} weighted bootstrap resamples")
    print(f"Wrote {len(summary['figures'])} PNG figures under {outdir}")


if __name__ == "__main__":
    main()
