#!/usr/bin/env python3
"""Aggregate raw acceptance diagnostics across representative files and masses."""

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


MASS_BY_DSID = {574792: 350, 574793: 400, 574794: 500, 574795: 750,
                574796: 850, 574797: 1000, 574798: 1250, 574799: 1500,
                574800: 1750, 574801: 2000, 574802: 2500, 574803: 3000,
                574804: 3500, 574805: 4000, 574806: 4500, 574807: 5000,
                574808: 6000}


def summarize(report, outdir):
    inputs = json.loads(Path(report).read_text())["inputs"]
    input_rows = []
    input_cutflow_rows = []
    for item in inputs:
        c = item["stage_counts"]
        ntruth = c["truth_direct_2l"]
        nselected_direct = c["analysis_selected"]
        nselected_all = c["all_analysis_selected"]
        input_rows.append({
            "dsid": item["dsid"], "mass_GeV": MASS_BY_DSID.get(item["dsid"]),
            "raw_file": item["raw_file"],
            "derivation_tag": (re.search(r"(r\d+)_p", Path(item["raw_file"]).parent.name).group(1)
                               if re.search(r"(r\d+)_p", Path(item["raw_file"]).parent.name) else "unknown"),
            "n_raw": c["raw_entries"],
            "n_truth_direct_2l": ntruth,
            "n_truth_2l_reco_lepton_pair": c["reco_os_tight_lepton_pair"],
            "n_truth_2l_triggered_after_lepton_pair": c["single_lepton_trigger"],
            "n_truth_2l_at_least_two_baseline_jvt_jets": c["at_least_two_baseline_jvt_jets"],
            "n_truth_2l_at_least_one_btagged_jet": c["at_least_one_baseline_77pct_btag_jet_pt30"],
            "n_truth_2l_selected": nselected_direct,
            "n_all_selected": nselected_all,
            "n_selected_b_pair_matched": c["unique_matched_b_pair"],
            "n_selected_lepton_pair_matched": c["selected_lepton_pair_truth_matched"],
            "n_selected_combined_match": c["strict_reco_object_match"],
            "f_reco_lepton_pair_given_truth_2l": c["reco_os_tight_lepton_pair"] / ntruth if ntruth else None,
            "f_prebtag_jet_stage_given_truth_2l": c["at_least_two_baseline_jvt_jets"] / ntruth if ntruth else None,
            "f_reco_object_eligible_given_truth_2l": c["at_least_one_baseline_77pct_btag_jet_pt30"] / ntruth if ntruth else None,
            "f_btag_given_prebtag_jet_stage": (
                c["at_least_one_baseline_77pct_btag_jet_pt30"] / c["at_least_two_baseline_jvt_jets"]
                if c["at_least_two_baseline_jvt_jets"] else None
            ),
            "f_analysis_selected_given_truth_2l": nselected_direct / ntruth if ntruth else None,
            "f_selected_b_pair_matched": c["unique_matched_b_pair"] / nselected_direct if nselected_direct else None,
            "f_selected_lepton_pair_matched": c["selected_lepton_pair_truth_matched"] / nselected_direct if nselected_direct else None,
            "f_selected_combined_match": c["strict_reco_object_match"] / nselected_direct if nselected_direct else None,
            "n_tau_w_daughter_selected": item["tau_category"]["n_tau_w_daughter_selected"],
            "f_tau_w_daughter_given_all_selected": (
                item["tau_category"]["n_tau_w_daughter_selected"] / nselected_all if nselected_all else None
            ),
            "n_reconstructed_topology_matched": item["reconstructed_topology_resolution"]["n_events"],
            "median_reco_minus_truth_min_dr": item["reconstructed_topology_resolution"]["median_delta_r"],
            "sigma68_reco_minus_truth_min_dr": item["reconstructed_topology_resolution"]["sigma68_delta_r"],
        })
        input_cutflow_rows.extend({"dsid": item["dsid"], "raw_file": item["raw_file"], **stage}
                                  for stage in item["cutflow"])

    # Aggregate independent derivation-tag files by DSID using event counts,
    # never by averaging per-file efficiencies.
    rows = []
    by_dsid = {}
    for item in inputs:
        by_dsid.setdefault(item["dsid"], []).append(item)
    count_fields = (
        "n_raw", "n_truth_direct_2l", "n_truth_2l_reco_lepton_pair",
        "n_truth_2l_triggered_after_lepton_pair", "n_truth_2l_at_least_two_baseline_jvt_jets",
        "n_truth_2l_at_least_one_btagged_jet", "n_truth_2l_selected", "n_all_selected",
        "n_selected_b_pair_matched", "n_selected_lepton_pair_matched",
        "n_selected_combined_match", "n_tau_w_daughter_selected",
    )
    for dsid, group in sorted(by_dsid.items()):
        group_rows = [row for row in input_rows if row["dsid"] == dsid]
        combined = {key: sum(row[key] for row in group_rows) for key in count_fields}
        ntruth, nselected_direct = combined["n_truth_direct_2l"], combined["n_truth_2l_selected"]
        nselected_all = combined["n_all_selected"]
        aggregate = {
            "dsid": dsid, "mass_GeV": MASS_BY_DSID.get(dsid), "n_input_files": len(group),
            "raw_files": ";".join(item["raw_file"] for item in group), **combined,
            "f_reco_lepton_pair_given_truth_2l": combined["n_truth_2l_reco_lepton_pair"] / ntruth if ntruth else None,
            "f_prebtag_jet_stage_given_truth_2l": combined["n_truth_2l_at_least_two_baseline_jvt_jets"] / ntruth if ntruth else None,
            "f_reco_object_eligible_given_truth_2l": combined["n_truth_2l_at_least_one_btagged_jet"] / ntruth if ntruth else None,
            "f_btag_given_prebtag_jet_stage": (
                combined["n_truth_2l_at_least_one_btagged_jet"]
                / combined["n_truth_2l_at_least_two_baseline_jvt_jets"]
                if combined["n_truth_2l_at_least_two_baseline_jvt_jets"] else None
            ),
            "f_analysis_selected_given_truth_2l": combined["n_truth_2l_selected"] / ntruth if ntruth else None,
            "f_selected_b_pair_matched": combined["n_selected_b_pair_matched"] / nselected_direct if nselected_direct else None,
            "f_selected_lepton_pair_matched": combined["n_selected_lepton_pair_matched"] / nselected_direct if nselected_direct else None,
            "f_selected_combined_match": combined["n_selected_combined_match"] / nselected_direct if nselected_direct else None,
            "f_tau_w_daughter_given_all_selected": (
                combined["n_tau_w_daughter_selected"] / nselected_all if nselected_all else None
            ),
        }
        weight_names = sorted(set().union(*(item["raw_weight_summary"].keys() for item in group)))
        for branch in weight_names:
            values = [item["raw_weight_summary"].get(branch) for item in group]
            numerator = sum(value["sum_analysis_selected"] for value in values if value is not None)
            denominator = sum(value["sum_truth_direct_2l"] for value in values if value is not None)
            aggregate[f"weighted_efficiency_{branch}"] = numerator / denominator if denominator else None
        rows.append(aggregate)

    stage_groups = {}
    for item in inputs:
        for stage in item["cutflow"]:
            key = (item["dsid"], stage["stage"])
            target = stage_groups.setdefault(key, {**stage, "dsid": item["dsid"],
                "raw_file": None, "n_events": 0, "denominator": 0,
                "weighted_numerator": 0.0, "weighted_denominator": 0.0, "n_inputs": 0})
            target["n_events"] += stage["n_events"]
            target["denominator"] += stage["denominator"]
            target["weighted_numerator"] += stage["weighted_numerator"] or 0.0
            target["weighted_denominator"] += stage["weighted_denominator"] or 0.0
            target["n_inputs"] += 1
    cutflow_rows = []
    for row in stage_groups.values():
        row["efficiency"] = row["n_events"] / row["denominator"] if row["denominator"] else None
        row["cumulative_fraction_of_raw"] = row["n_events"] / (
            sum(item["stage_counts"]["raw_entries"] for item in by_dsid[row["dsid"]])
        )
        row["weighted_efficiency"] = (row["weighted_numerator"] / row["weighted_denominator"]
                                      if row["weighted_denominator"] else None)
        row["weighting"] = "weight_mc_NOSYS"
        cutflow_rows.append(row)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "raw_acceptance_representative.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    input_csv_path = outdir / "raw_acceptance_by_input.csv"
    with input_csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(input_rows[0]))
        writer.writeheader()
        writer.writerows(input_rows)
    cutflow_path = outdir / "raw_acceptance_cutflow.csv"
    with cutflow_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cutflow_rows[0]))
        writer.writeheader()
        writer.writerows(cutflow_rows)
    input_differential_rows = []
    for item in inputs:
        for variable, bins in item.get("differential", {}).items():
            input_differential_rows.extend({"dsid": item["dsid"],
                                            "mass_GeV": MASS_BY_DSID.get(item["dsid"]),
                                            "raw_file": item["raw_file"],
                                            "x_variable": variable, **row} for row in bins)
    differential_groups = {}
    for row in input_differential_rows:
        key = (row["dsid"], row["x_variable"], row["x_low"], row["x_high"])
        target = differential_groups.setdefault(key, {
            "dsid": row["dsid"], "mass_GeV": row["mass_GeV"], "x_variable": row["x_variable"],
            "x_low": row["x_low"], "x_high": row["x_high"], "n_truth_direct_2l": 0,
            "n_analysis_selected": 0, "n_selected_b_pair_matched": 0,
            "n_reconstructed_topology_matched": 0,
            "reco_truth_dr_shift_histogram": [0] * len(row.get("reco_truth_dr_shift_histogram") or []),
        })
        for name in ("n_truth_direct_2l", "n_analysis_selected", "n_selected_b_pair_matched",
                     "n_reconstructed_topology_matched"):
            target[name] += row.get(name, 0) or 0
        hist = row.get("reco_truth_dr_shift_histogram")
        if hist is not None:
            target["reco_truth_dr_shift_histogram"] = [a + b for a, b in zip(target["reco_truth_dr_shift_histogram"], hist)]
    differential_rows = []
    for row in differential_groups.values():
        den, selected, matched = (row["n_truth_direct_2l"], row["n_analysis_selected"],
                                  row["n_selected_b_pair_matched"])
        row.update({
            "f_selected_given_truth_2l": selected / den if den else None,
            "f_b_pair_matched_given_truth_2l": matched / den if den else None,
            "f_complete_selection_and_b_pair_given_truth_2l": matched / den if den else None,
            "f_b_pair_matched_given_selected": matched / selected if selected else None,
        })
        differential_rows.append(row)
    if differential_rows:
        differential_path = outdir / "raw_acceptance_differential.csv"
        with differential_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(differential_rows[0]))
            writer.writeheader()
            writer.writerows(differential_rows)
        input_differential_path = outdir / "raw_acceptance_differential_by_input.csv"
        with input_differential_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(input_differential_rows[0]))
            writer.writeheader()
            writer.writerows(input_differential_rows)
    kk = sorted((row for row in rows if row["mass_GeV"] is not None),
                key=lambda row: row["mass_GeV"])
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 7.0), sharex=True, constrained_layout=True)
    x = [row["mass_GeV"] for row in kk]
    for field, label, color in (
        ("f_reco_lepton_pair_given_truth_2l", "Reco OS lepton pair / direct truth 2L", "#7389a1"),
        ("f_btag_given_prebtag_jet_stage", "≥1 b-tag jet / prior lepton+trigger+jet stage", "#dfb56c"),
        ("f_analysis_selected_given_truth_2l", "Analysis selected / direct truth 2L", "#826fa6"),
    ):
        axes[0].plot(x, [row[field] for row in kk], "o-", label=label, color=color)
    for field, label, color in (
        ("f_selected_b_pair_matched", "Unique reco b pair / direct-2L selected", "#7389a1"),
        ("f_selected_lepton_pair_matched", "Matched selected leptons / direct-2L selected", "#dfb56c"),
        ("f_selected_combined_match", "Combined match / direct-2L selected", "#826fa6"),
    ):
        axes[1].plot(x, [row[field] for row in kk], "o-", label=label, color=color)
    for ax in axes:
        ax.set(ylabel="Unweighted fraction", ylim=(-0.02, 1.02))
        ax.grid(axis="y", linestyle=":", alpha=0.3)
        ax.legend(fontsize=8)
    axes[1].set_xlabel(r"Generated benchmark mass $M_G$ [GeV]")
    fig.suptitle("Raw acceptance stages — independent derivation-tag slices (diagnostic sample)")
    plot_path = outdir / "raw_acceptance_representative.png"
    fig.savefig(plot_path, bbox_inches="tight")
    plt.close(fig)
    topology_rows = [row for row in input_rows
                     if row["mass_GeV"] is not None and row["median_reco_minus_truth_min_dr"] is not None]
    if topology_rows:
        fig, ax = plt.subplots(figsize=(8.0, 4.6), constrained_layout=True)
        tag_colors = {"r13144": "#7389a1", "r13145": "#dfb56c", "r13167": "#826fa6"}
        for tag in sorted({row["derivation_tag"] for row in topology_rows}):
            points = sorted((row for row in topology_rows if row["derivation_tag"] == tag),
                            key=lambda row: row["mass_GeV"])
            ax.plot([row["mass_GeV"] for row in points],
                    [row["median_reco_minus_truth_min_dr"] for row in points],
                    "o", color=tag_colors.get(tag, "#589c91"), label=tag)
        ax.axhline(0.0, color="#555555", linestyle=":", linewidth=1)
        ax.set(xlabel=r"Generated benchmark mass $M_G$ [GeV]",
               ylabel=r"Median $\Delta R_{b\ell}^{reco}-\Delta R_{b\ell}^{truth}$")
        ax.grid(axis="y", linestyle=":", alpha=0.3)
        ax.legend(title="Derivation tag", fontsize=8)
        ax.set_title("Matched reconstructed topology response (per independent raw input)")
        topology_plot = outdir / "raw_topology_resolution_by_input.png"
        fig.savefig(topology_plot, bbox_inches="tight")
        plt.close(fig)
        print(topology_plot)
    if differential_rows:
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.8), constrained_layout=True)
        panels = (
            ("truth_min_dr_bl", "f_selected_given_truth_2l", "Truth min $\\Delta R(b,\\ell)$",
             "Selected / direct truth 2L"),
            ("truth_min_dr_bl", "f_b_pair_matched_given_truth_2l", "Truth min $\\Delta R(b,\\ell)$",
             "Selected unique b pair / direct truth 2L"),
            ("truth_min_dr_bl", "f_b_pair_matched_given_selected", "Truth min $\\Delta R(b,\\ell)$",
             "Unique reco b pair / selected"),
            ("truth_top_pt_GeV", "f_selected_given_truth_2l", "Truth $p_T(t)$ [GeV]",
             "Selected / direct truth 2L"),
        )
        for ax, (variable, field, xlabel, ylabel) in zip(axes.ravel(), panels):
            for mass, color in ((500, "#7389a1"), (1000, "#826fa6"),
                                (3000, "#dfb56c"), (5000, "#589c91")):
                bins = [row for row in differential_rows if row["mass_GeV"] == mass
                        and row["x_variable"] == variable
                         and row["n_truth_direct_2l"] >= 20
                        and (field != "f_b_pair_matched_given_selected" or row["n_analysis_selected"] >= 20)]
                ax.plot([(row["x_low"] + row["x_high"]) / 2 for row in bins],
                        [row[field] for row in bins], "o-", color=color, label=f"{mass} GeV")
            ax.set(xlabel=xlabel, ylabel=ylabel, ylim=(-0.02, 1.02))
            ax.grid(axis="y", linestyle=":", alpha=0.3)
        axes[0, 0].legend(fontsize=8, ncol=2)
        fig.suptitle("Raw acceptance vs truth kinematics — independent derivation-tag slices")
        differential_plot = outdir / "raw_acceptance_vs_truth_kinematics.png"
        fig.savefig(differential_plot, bbox_inches="tight")
        plt.close(fig)
        print(differential_path)
        print(input_differential_path)
        print(differential_plot)
    print(csv_path)
    print(input_csv_path)
    print(cutflow_path)
    print(plot_path)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report_json")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    summarize(args.report_json, args.outdir)


if __name__ == "__main__":
    main()
