#!/usr/bin/env python3
"""Count the Run-2 dilepton stages before the selected rich ntuple exists.

This reads only raw branches used by the current TtbarDileptonFrameReco region
or by the explicit after-FSR truth definition.  Counts are unweighted event
counts; one raw file is a diagnostic slice, not a campaign acceptance result.
"""

import argparse
import json
import sys
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Plot.reco_plot.study import _truth_lepton_pair_reconstructible


TRIGGERS = (
    "trigPassed_HLT_e60_lhmedium_L1EM22VHI",
    "trigPassed_HLT_e140_lhloose_L1EM22VHI",
    "trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI",
    "trigPassed_HLT_e26_lhtight_nod0_ivarloose",
    "trigPassed_HLT_e140_lhloose_nod0",
    "trigPassed_HLT_e60_lhmedium_nod0",
    "trigPassed_HLT_mu24_ivarmedium_L1MU14FCH",
    "trigPassed_HLT_mu50_L1MU14FCH",
    "trigPassed_HLT_mu50",
    "trigPassed_HLT_mu26_ivarmedium",
)


def stage_cutflow(stages, denominator_names, weights=None):
    """Return nested/conditional event counts and efficiencies with fixed denominators."""
    weights = None if weights is None else np.asarray(weights, dtype=float).reshape(-1)
    arrays = {name: np.asarray(mask, dtype=bool).reshape(-1) for name, mask in stages.items()}
    sizes = {mask.size for mask in arrays.values()}
    if len(sizes) != 1 or (weights is not None and weights.size not in sizes):
        raise ValueError("cutflow masks and weights must have the same event count")
    rows = []
    for name, numerator_mask in arrays.items():
        denominator_name = denominator_names[name]
        denominator_mask = arrays[denominator_name]
        if np.any(numerator_mask & ~denominator_mask):
            raise ValueError(f"cutflow stage {name!r} is not nested in denominator {denominator_name!r}")
        numerator = int(np.count_nonzero(numerator_mask))
        denominator = int(np.count_nonzero(denominator_mask))
        row = {"stage": name, "denominator_stage": denominator_name,
               "n_events": numerator, "denominator": denominator,
               "efficiency": numerator / denominator if denominator else None,
               "cumulative_fraction_of_raw": numerator / arrays["raw_entries"].size}
        if weights is not None:
            finite = np.isfinite(weights)
            numerator_weight = float(np.sum(weights[numerator_mask & finite]))
            denominator_weight = float(np.sum(weights[denominator_mask & finite]))
            row.update({"weighting": "weight_mc_NOSYS", "weighted_numerator": numerator_weight,
                        "weighted_denominator": denominator_weight,
                        "weighted_efficiency": numerator_weight / denominator_weight
                        if denominator_weight != 0.0 else None})
        else:
            row.update({"weighting": "unweighted", "weighted_numerator": None,
                        "weighted_denominator": None, "weighted_efficiency": None})
        rows.append(row)
    return rows


def _raw_file(path, max_events=None, selected_root=None):
    file = uproot.open(path)
    if "reco" not in file:
        raise KeyError(f"{path}: raw reco tree is absent")
    tree = file["reco"]
    required = ["mcChannelNumber", "eventNumber", "Ttbar_MC_ttbar_afterFSR_m",
                "Ttbar_MC_t_afterFSR_pt",
                "Ttbar_MC_b_afterFSR_from_t_eta", "Ttbar_MC_b_afterFSR_from_t_phi",
                "Ttbar_MC_bbar_afterFSR_from_tbar_eta", "Ttbar_MC_bbar_afterFSR_from_tbar_phi",
                "jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_select_baselineJvt_NOSYS",
                "jet_select_GN2v01_FixedCutBEff_77",
                "n_jet_matched_to_Ttbar_MC_b_afterFSR_from_t",
                "n_jet_matched_to_Ttbar_MC_bbar_afterFSR_from_tbar",
                "jet_dRmatched_to_Ttbar_MC_b_afterFSR_from_t",
                "jet_dRmatched_to_Ttbar_MC_bbar_afterFSR_from_tbar", *TRIGGERS]
    for flavour in ("el", "mu"):
        required += [f"{flavour}_pt_NOSYS", f"{flavour}_eta", f"{flavour}_phi",
                     f"{flavour}_charge", f"{flavour}_select_tight_NOSYS"]
    for side in ("t", "tbar"):
        for daughter in (1, 2):
            prefix = f"Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
            required += [f"{prefix}_pdgId", f"{prefix}_pt", f"{prefix}_eta", f"{prefix}_phi",
                         f"n_lep_matched_to_{prefix}"]
            for flavour in ("el", "mu"):
                required.append(f"{flavour}_dRmatched_to_{prefix}")
    tree_keys = set(tree.keys())
    absent = sorted(set(required) - tree_keys)
    if absent:
        raise KeyError(f"{path}: required raw acceptance branches absent: {absent}")
    optional = [name for name in ("weight_mc_NOSYS", "weight_total_NOSYS") if name in tree_keys]
    arrays = tree.arrays(required + optional, library="ak", entry_stop=max_events)
    n = len(arrays["eventNumber"])
    arr = lambda name: np.asarray(ak.to_numpy(arrays[name]))
    truth_mtt = arr("Ttbar_MC_ttbar_afterFSR_m") * 0.001
    truth_2l = np.isfinite(truth_mtt) & (truth_mtt > 0)
    charged_daughters = {}
    tau_daughters = {}
    for side, expected_sign in (("t", -1), ("tbar", 1)):
        daughters = [arr(f"Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}_pdgId") for daughter in (1, 2)]
        charged_daughters[side] = daughters
        has_charged = np.logical_or.reduce([
            (np.sign(pdg) == expected_sign) & np.isin(np.abs(pdg), (11, 13))
            for pdg in daughters
        ])
        truth_2l &= has_charged
        tau_daughters[side] = np.logical_or.reduce([np.abs(pdg) == 15 for pdg in daughters])
    truth_tau_w = tau_daughters["t"] | tau_daughters["tbar"]

    leptons = {}
    for flavour in ("el", "mu"):
        leptons[flavour] = {
            key: ak.to_list(arrays[f"{flavour}_{branch}"])
            for key, branch in (("pt", "pt_NOSYS"), ("charge", "charge"),
                                ("eta", "eta"), ("phi", "phi"),
                                ("tight", "select_tight_NOSYS"))
        }
    lepton_pair = np.zeros(n, dtype=bool)
    lepton_match = np.zeros(n, dtype=bool)
    for i in range(n):
        charges = [float(q) for flavour in ("el", "mu")
                   for pt, q, keep in zip(leptons[flavour]["pt"][i],
                                          leptons[flavour]["charge"][i],
                                          leptons[flavour]["tight"][i])
                   if keep and pt > 25000]
        lepton_pair[i] = len(charges) == 2 and np.isfinite(charges).all() and sum(charges) == 0
        if not lepton_pair[i]:
            continue
        event_leptons = {flavour: {
            "charge": leptons[flavour]["charge"][i],
            "tight": [keep and pt > 25000 for pt, keep in
                      zip(leptons[flavour]["pt"][i], leptons[flavour]["tight"][i])],
        } for flavour in ("el", "mu")}
        match_indices = {}
        match_counts = {}
        for side in ("t", "tbar"):
            for daughter in (1, 2):
                prefix = f"Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
                match_counts[(side, daughter)] = int(arrays[f"n_lep_matched_to_{prefix}"][i])
                for flavour in ("el", "mu"):
                    match_indices[(side, daughter, flavour)] = int(arrays[f"{flavour}_dRmatched_to_{prefix}"][i])
        lepton_match[i] = _truth_lepton_pair_reconstructible(event_leptons, match_indices, match_counts)

    jet_pt = arrays["jet_pt_NOSYS"]
    baseline = arrays["jet_select_baselineJvt_NOSYS"] != 0
    btag = arrays["jet_select_GN2v01_FixedCutBEff_77"] != 0
    # The configured region uses n_jets_NOSYS, whose frame definition counts
    # baseline-JVT selected jets with a zero-GeV extra threshold.
    jet_eligible = np.asarray(ak.to_numpy(ak.sum(baseline, axis=1) >= 2), dtype=bool)
    btag_eligible = np.asarray(ak.to_numpy(ak.sum(baseline & btag & (jet_pt > 30000), axis=1) >= 1), dtype=bool)
    triggered = np.logical_or.reduce([arr(name) != 0 for name in TRIGGERS])
    analysis_selected = lepton_pair & jet_eligible & btag_eligible & triggered
    b_count = arr("n_jet_matched_to_Ttbar_MC_b_afterFSR_from_t")
    bbar_count = arr("n_jet_matched_to_Ttbar_MC_bbar_afterFSR_from_tbar")
    b_index = arr("jet_dRmatched_to_Ttbar_MC_b_afterFSR_from_t")
    bbar_index = arr("jet_dRmatched_to_Ttbar_MC_bbar_afterFSR_from_tbar")
    raw_njet = np.asarray(ak.to_numpy(ak.num(jet_pt)), dtype=int)
    b_pair = ((b_count == 1) & (bbar_count == 1) & (b_index >= 0) & (bbar_index >= 0)
              & (b_index < raw_njet) & (bbar_index < raw_njet) & (b_index != bbar_index))

    stages = {
        "raw_entries": np.ones(n, dtype=bool),
        "truth_direct_2l": truth_2l,
        "reco_os_tight_lepton_pair": truth_2l & lepton_pair,
        "single_lepton_trigger": truth_2l & lepton_pair & triggered,
        "at_least_two_baseline_jvt_jets": truth_2l & lepton_pair & triggered & jet_eligible,
        "at_least_one_baseline_77pct_btag_jet_pt30": truth_2l & lepton_pair & triggered & jet_eligible & btag_eligible,
        "analysis_selected": truth_2l & analysis_selected,
        "all_analysis_selected": analysis_selected,
        "unique_matched_b_pair": truth_2l & analysis_selected & b_pair,
        "selected_lepton_pair_truth_matched": truth_2l & analysis_selected & lepton_match,
        "strict_reco_object_match": truth_2l & analysis_selected & b_pair & lepton_match,
    }
    stage_denominators = {
        "raw_entries": "raw_entries",
        "truth_direct_2l": "raw_entries",
        "reco_os_tight_lepton_pair": "truth_direct_2l",
        "single_lepton_trigger": "reco_os_tight_lepton_pair",
        "at_least_two_baseline_jvt_jets": "single_lepton_trigger",
        "at_least_one_baseline_77pct_btag_jet_pt30": "at_least_two_baseline_jvt_jets",
        "analysis_selected": "at_least_one_baseline_77pct_btag_jet_pt30",
        "all_analysis_selected": "raw_entries",
        "unique_matched_b_pair": "analysis_selected",
        "selected_lepton_pair_truth_matched": "analysis_selected",
        "strict_reco_object_match": "analysis_selected",
    }
    cutflow = stage_cutflow(stages, stage_denominators,
                            arr("weight_mc_NOSYS") if "weight_mc_NOSYS" in optional else None)

    truth_dr = np.full(n, np.inf)
    for side, bprefix, sign in (
        ("t", "Ttbar_MC_b_afterFSR_from_t", -1),
        ("tbar", "Ttbar_MC_bbar_afterFSR_from_tbar", 1),
    ):
        lepton_eta = np.full(n, np.nan)
        lepton_phi = np.full(n, np.nan)
        for daughter in (1, 2):
            prefix = f"Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
            pdg = arr(prefix + "_pdgId")
            charged = (np.sign(pdg) == sign) & np.isin(np.abs(pdg), (11, 13))
            lepton_eta[charged] = arr(prefix + "_eta")[charged]
            lepton_phi[charged] = arr(prefix + "_phi")[charged]
        dphi = (arr(bprefix + "_phi") - lepton_phi + np.pi) % (2 * np.pi) - np.pi
        dr = np.hypot(arr(bprefix + "_eta") - lepton_eta, dphi)
        truth_dr = np.minimum(truth_dr, dr)
    truth_dr[~np.isfinite(truth_dr)] = np.nan

    reco_dr = np.full(n, np.nan)
    reco_dr_shift = np.full(n, np.nan)
    jet_eta = ak.to_list(arrays["jet_eta"])
    jet_phi = ak.to_list(arrays["jet_phi"])
    for i in np.flatnonzero(analysis_selected & b_pair & lepton_match & truth_2l):
        matched_dr = []
        for side, jet_index, expected_sign in (
            ("t", int(b_index[i]), -1), ("tbar", int(bbar_index[i]), 1),
        ):
            matched_leptons = []
            for daughter in (1, 2):
                pdg = charged_daughters[side][daughter - 1][i]
                if np.sign(pdg) != expected_sign or abs(pdg) not in (11, 13):
                    continue
                prefix = f"Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
                count = int(arrays[f"n_lep_matched_to_{prefix}"][i])
                candidates = []
                for flavour in ("el", "mu"):
                    index = int(arrays[f"{flavour}_dRmatched_to_{prefix}"][i])
                    if index >= 0:
                        candidates.append((flavour, index))
                if count == 1 and len(candidates) == 1:
                    matched_leptons = candidates
            if len(matched_leptons) != 1:
                matched_dr = []
                break
            flavour, lepton_index = matched_leptons[0]
            if lepton_index >= len(leptons[flavour]["eta"][i]) or jet_index >= len(jet_eta[i]):
                matched_dr = []
                break
            dphi = (float(jet_phi[i][jet_index]) -
                    float(leptons[flavour]["phi"][i][lepton_index]) + np.pi) % (2 * np.pi) - np.pi
            matched_dr.append(float(np.hypot(
                float(jet_eta[i][jet_index]) - float(leptons[flavour]["eta"][i][lepton_index]),
                dphi,
            )))
        if len(matched_dr) == 2:
            reco_dr[i] = min(matched_dr)
            reco_dr_shift[i] = reco_dr[i] - truth_dr[i]
    matched_topology = analysis_selected & b_pair & lepton_match & truth_2l & np.isfinite(reco_dr_shift)

    differential = {}
    axes = {
        "truth_mtt_GeV": (truth_mtt, [300, 400, 500, 650, 800, 1000, 1250, 1500,
                                      2000, 2500, 3000, 4000, 5000, 6500, 10000]),
        "truth_top_pt_GeV": (arr("Ttbar_MC_t_afterFSR_pt") * 0.001,
                              [0, 100, 200, 300, 400, 600, 800, 1200, 1800, 2500, 4000]),
        "truth_min_dr_bl": (truth_dr,
                             [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]),
    }
    for name, (values, edges) in axes.items():
        rows = []
        for low, high in zip(edges[:-1], edges[1:]):
            in_bin = np.isfinite(values) & (values >= low) & (values < high)
            denominator = int(np.count_nonzero(in_bin & truth_2l))
            selected_n = int(np.count_nonzero(in_bin & truth_2l & analysis_selected))
            matched_n = int(np.count_nonzero(in_bin & truth_2l & analysis_selected & b_pair))
            topology_values = reco_dr_shift[in_bin & matched_topology]
            rows.append({"x_low": low, "x_high": high, "n_truth_direct_2l": denominator,
                         "n_analysis_selected": selected_n, "n_selected_b_pair_matched": matched_n,
                         "f_selected_given_truth_2l": selected_n / denominator if denominator else None,
                         "f_b_pair_matched_given_selected": matched_n / selected_n if selected_n else None,
                         "f_b_pair_matched_given_truth_2l": matched_n / denominator if denominator else None,
                         "f_complete_selection_and_b_pair_given_truth_2l": matched_n / denominator if denominator else None,
                         "n_reconstructed_topology_matched": int(topology_values.size),
                         "median_reco_minus_truth_min_dr": float(np.median(topology_values)) if topology_values.size else None,
                         "sigma68_reco_minus_truth_min_dr": float(np.percentile(topology_values, 84) - np.percentile(topology_values, 16)) if topology_values.size else None,
                         "reco_truth_dr_shift_histogram": np.histogram(
                             topology_values, bins=[-2, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 2, 5]
                         )[0].astype(int).tolist() if name == "truth_min_dr_bl" else None})
        differential[name] = rows
    keys = list(zip(arr("mcChannelNumber").astype(int), arr("eventNumber").astype(int)))
    if len(set(keys)) != len(keys):
        raise ValueError(f"{path}: duplicate (mcChannelNumber,eventNumber) in raw input")
    key_validation = None
    if selected_root:
        output = uproot.open(selected_root)["reco"]
        out = output.arrays(["mcChannelNumber", "eventNumber"], library="np")
        output_keys = set(zip(out["mcChannelNumber"].astype(int), out["eventNumber"].astype(int)))
        selected_keys = {key for key, keep in zip(keys, analysis_selected) if keep}
        key_validation = {"selected_raw": len(selected_keys), "selected_output": len(output_keys),
                          "missing_from_output": len(selected_keys - output_keys),
                          "extra_in_output": len(output_keys - selected_keys)}
        if selected_keys != output_keys:
            raise ValueError(f"{path}: raw selection replay does not match {selected_root}: {key_validation}")
    raw_weights = {}
    for branch in optional:
        weights = np.asarray(ak.to_numpy(arrays[branch]), dtype=float)
        finite = np.isfinite(weights)
        raw_weights[branch] = {
            "n_finite": int(np.count_nonzero(finite)),
            "sum_truth_direct_2l": float(np.sum(weights[truth_2l & finite])),
            "sum_analysis_selected": float(np.sum(weights[analysis_selected & truth_2l & finite])),
            "weighted_selected_given_truth_direct_2l": (
                float(np.sum(weights[analysis_selected & truth_2l & finite]) / np.sum(weights[truth_2l & finite]))
                if np.sum(weights[truth_2l & finite]) != 0 else None
            ),
            "minimum": float(np.min(weights[finite])) if np.any(finite) else None,
            "maximum": float(np.max(weights[finite])) if np.any(finite) else None,
            "standard_deviation": float(np.std(weights[finite])) if np.any(finite) else None,
        }
    dr_shift = reco_dr_shift[matched_topology]
    topology_resolution = {
        "definition": "minimum reconstructed matched b-lepton DeltaR minus truth minimum DeltaR; requires selected event, unique b-pair match, and selected-lepton truth match",
        "n_events": int(dr_shift.size),
        "median_delta_r": float(np.median(dr_shift)) if dr_shift.size else None,
        "sigma68_delta_r": float(np.percentile(dr_shift, 84) - np.percentile(dr_shift, 16)) if dr_shift.size else None,
        "histogram_edges": [-2, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 2, 5],
        "histogram_counts": np.histogram(dr_shift, bins=[-2, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 2, 5])[0].astype(int).tolist(),
    }
    tau_selected = int(np.count_nonzero(truth_tau_w & analysis_selected))
    report = {"raw_file": str(Path(path).resolve()), "entries_read": n,
              "dsid": int(np.unique(arr("mcChannelNumber"))[0]),
              "denominators": {"truth_stages": "direct after-FSR W->e/mu from both top parents within entries read",
                               "selected_stages": "all ttll_analysis selected events within entries read",
                               "weights": "unweighted counts plus nominal generator weight where available"},
              "stage_counts": {name: int(np.count_nonzero(mask)) for name, mask in stages.items()},
              "n_analysis_selected_all": int(np.count_nonzero(analysis_selected)),
              "cutflow": cutflow,
              "tau_category": {"definition": "either top-parent W has an after-FSR tau daughter; separate from direct-e/mu denominator",
                               "n_with_tau_w_daughter": int(np.count_nonzero(truth_tau_w)),
                               "n_tau_w_daughter_selected": tau_selected,
                               "tau_selected_fraction_of_selected": tau_selected / int(np.count_nonzero(analysis_selected))
                               if np.count_nonzero(analysis_selected) else None},
              "raw_weight_summary": raw_weights,
              "reconstructed_topology_resolution": topology_resolution,
              "differential": differential,
              "key_validation": key_validation,
              "limitations": "Truth-direct-2l excludes W->tau decays. The raw replay uses exactly the configured region's two tight pT>25 GeV OS leptons, trigger OR, n_jets_NOSYS>=2 baseline-JVT jets, and n_bjets_NOSYS>=1 selected baseline 77%-tagged jet with pT>30 GeV. Raw match stages do not include the rich-ntuple neutrino/oracle validity checks."}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_root", nargs="*", help="Raw ntuple ROOT files (reco trees)")
    parser.add_argument("--filelist", help="Text file with one explicit raw ROOT path per line")
    parser.add_argument("--max-events", type=int)
    parser.add_argument("--selected-root", help="Validate exact selected event keys against one rich output")
    parser.add_argument("--output", required=True, help="JSON report")
    args = parser.parse_args()
    paths = list(args.raw_root)
    if args.filelist:
        paths.extend(line.strip() for line in Path(args.filelist).read_text().splitlines()
                     if line.strip() and not line.lstrip().startswith("#"))
    if not paths:
        parser.error("give raw_root paths or --filelist")
    if args.selected_root and len(paths) != 1:
        parser.error("--selected-root requires exactly one raw ROOT input")
    reports = [_raw_file(path, args.max_events, args.selected_root) for path in paths]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"inputs": reports}, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
