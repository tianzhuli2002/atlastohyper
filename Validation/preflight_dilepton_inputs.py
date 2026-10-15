#!/usr/bin/env python3
"""Preflight Run-2 dilepton ROOT inputs and KK-graviton metadata.

Reads ROOT trees with uproot. It does not use MakeClass-generated readers, whose
split Lorentz-vector leaves can collide in C++ branch naming.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import awkward as ak
import numpy as np
import uproot


DEFAULT_METADATA = Path("/home/leet/FF/metadata_kkgraviton")
SPIN_OBSERVABLES = (
    "cos_phi",
    "cos_theta_helicity_m",
    "cos_theta_helicity_p",
    "cos_theta_raxis_m",
    "cos_theta_raxis_p",
    "cos_theta_transverse_m",
    "cos_theta_transverse_p",
    "cos_han",
)


# Exact object-level aliases found in the raw Run-2 ntuple. FastFrames output
# aliases (Jet_pt, Electron_pt, etc.) are not raw-input names.
OBJECT_BRANCH_TYPES = {
    "eventNumber": "scalar_int",
    "mcChannelNumber": "scalar_int",
    "runNumber": "scalar_int",
    "jet_pt_NOSYS": "vector_float",
    "jet_eta": "vector_float",
    "jet_phi": "vector_float",
    "jet_e_NOSYS": "vector_float",
    "jet_GN2v01_Continuous_quantile": "vector_int",
    "el_pt_NOSYS": "vector_float",
    "el_eta": "vector_float",
    "el_phi": "vector_float",
    "el_e_NOSYS": "vector_float",
    "el_charge": "vector_float",
    "mu_pt_NOSYS": "vector_float",
    "mu_eta": "vector_float",
    "mu_phi": "vector_float",
    "mu_e_NOSYS": "vector_float",
    "mu_charge": "vector_float",
    "met_met_NOSYS": "scalar_float",
    "met_phi_NOSYS": "scalar_float",
}


def _truth_specs() -> dict[str, str]:
    result: dict[str, str] = {}
    truth_particles = (
        "Ttbar_MC_t_afterFSR",
        "Ttbar_MC_tbar_afterFSR",
        "Ttbar_MC_t_beforeFSR",
        "Ttbar_MC_tbar_beforeFSR",
        "Ttbar_MC_b_afterFSR_from_t",
        "Ttbar_MC_b_beforeFSR_from_t",
        "Ttbar_MC_bbar_afterFSR_from_tbar",
        "Ttbar_MC_bbar_beforeFSR_from_tbar",
        "Ttbar_MC_W_afterFSR_from_t",
        "Ttbar_MC_W_afterFSR_from_tbar",
        "Ttbar_MC_W_beforeFSR_from_t",
        "Ttbar_MC_W_beforeFSR_from_tbar",
        "Ttbar_MC_Wdecay1_afterFSR_from_t",
        "Ttbar_MC_Wdecay2_afterFSR_from_t",
        "Ttbar_MC_Wdecay1_afterFSR_from_tbar",
        "Ttbar_MC_Wdecay2_afterFSR_from_tbar",
        "Ttbar_MC_Wdecay1_beforeFSR_from_t",
        "Ttbar_MC_Wdecay2_beforeFSR_from_t",
        "Ttbar_MC_Wdecay1_beforeFSR_from_tbar",
        "Ttbar_MC_Wdecay2_beforeFSR_from_tbar",
        "Ttbar_MC_ttbar_afterFSR",
        "Ttbar_MC_ttbar_beforeFSR",
    )
    for particle in truth_particles:
        for component in ("pt", "eta", "phi", "m"):
            result[f"{particle}_{component}"] = "scalar_float"
        if "ttbar" not in particle:
            result[f"{particle}_pdgId"] = "scalar_int"

    # Raw per-object nearest-match indices and match multiplicities. These are
    # raw-schema ID/matching branches; they are not reinterpreted as the
    # FastFrames Truth_had_b_index labels.
    match_targets = (
        "Ttbar_MC_Wdecay1_afterFSR_from_t",
        "Ttbar_MC_Wdecay2_afterFSR_from_t",
        "Ttbar_MC_Wdecay1_afterFSR_from_tbar",
        "Ttbar_MC_Wdecay2_afterFSR_from_tbar",
        "Ttbar_MC_b_afterFSR_from_t",
        "Ttbar_MC_bbar_afterFSR_from_tbar",
    )
    for target in match_targets:
        result[f"jet_dRmatched_to_{target}"] = "scalar_int"
    for prefix in ("el", "mu", "tau"):
        for target in match_targets[:4]:
            result[f"{prefix}_dRmatched_to_{target}"] = "scalar_int"
    for target in match_targets:
        result[f"n_jet_matched_to_{target}"] = "scalar_int"
    for target in match_targets[:4]:
        result[f"n_lep_matched_to_{target}"] = "scalar_int"
    return result


TRUTH_BRANCH_TYPES = _truth_specs()


def _p4(prefix: str, names: tuple[str, ...]) -> dict[str, str]:
    return {f"{prefix}_{name}_NOSYS": "p4" for name in names}


def _scalars(prefix: str, names: tuple[str, ...]) -> dict[str, str]:
    return {f"{prefix}_{name}_NOSYS": "scalar_float" for name in names}


def _method_specs() -> dict[str, dict[str, str]]:
    da = {
        **_p4("TtbarDiLepton_DA", ("Top", "Tbar", "Top_b_p4", "Top_lep_p4", "Top_nu_p4", "Topbar_bbar_p4", "Topbar_lepbar_p4", "Topbar_nubar_p4")),
        "TtbarDiLepton_DA_jet_indices_NOSYS": "vector_int",
        "TtbarDiLepton_DA_jet_probs_NOSYS": "vector_float",
        "TtbarDiLepton_DA_neutrino_momenta_NOSYS": "vector_float",
        **_scalars("TtbarDiLepton_DA", ("nu_pt", "nu_eta", "nu_phi", "nu_e", "nubar_pt", "nubar_eta", "nubar_phi", "nubar_e")),
    }
    em = {
        **_p4("TtbarDiLepton_EM", ("Top", "Tbar", "Top_b_p4", "Top_lep_p4", "Top_nu_p4", "Topbar_bbar_p4", "Topbar_lepbar_p4", "Topbar_nubar_p4")),
        "TtbarDiLepton_EM_b_idx_NOSYS": "scalar_int",
        "TtbarDiLepton_EM_bbar_idx_NOSYS": "scalar_int",
        **_scalars("TtbarDiLepton_EM", ("nu_pt", "nu_eta", "nu_phi", "nu_e", "nubar_pt", "nubar_eta", "nubar_phi", "nubar_e")),
    }
    nw = {
        **_p4("TtbarDiLepton_NW", ("Top", "Tbar", "Top_b_p4", "Top_lep_p4", "Top_nu_p4", "Topbar_bbar_p4", "Topbar_lepbar_p4", "Topbar_nubar_p4")),
        "TtbarDiLepton_NW_b_idx_NOSYS": "scalar_int",
        "TtbarDiLepton_NW_bbar_idx_NOSYS": "scalar_int",
        "TtbarDiLepton_NW_weight_max_NOSYS": "scalar_float",
        **_scalars("TtbarDiLepton_NW", ("nu_pt", "nu_eta", "nu_phi", "nu_e", "nubar_pt", "nubar_eta", "nubar_phi", "nubar_e")),
    }
    nuflows = {
        **_p4("TtbarDiLepton_NuFlows", ("Top_nu_p4", "Topbar_nubar_p4")),
        "TtbarDiLepton_NuFlows_nuflows_nu_out_NOSYS": "vector_float",
        "TtbarDiLepton_NuFlows_nuflows_loglik_NOSYS": "scalar_float",
        **_scalars("TtbarDiLepton_NuFlows", ("nu_pt", "nu_eta", "nu_phi", "nu_e", "nubar_pt", "nubar_eta", "nubar_phi", "nubar_e")),
    }
    hyper = {
        "TtbarDiLepton_HyPER_Top1_IDs_NOSYS": "vector_int",
        "TtbarDiLepton_HyPER_Top1_Indices_NOSYS": "vector_int",
        "TtbarDiLepton_HyPER_Top1_Score_NOSYS": "scalar_float",
        "TtbarDiLepton_HyPER_Top2_IDs_NOSYS": "vector_int",
        "TtbarDiLepton_HyPER_Top2_Indices_NOSYS": "vector_int",
        "TtbarDiLepton_HyPER_Top2_Score_NOSYS": "scalar_float",
        "TtbarDiLepton_HyPER_Classification_Score_NOSYS": "scalar_float",
        "TtbarDiLepton_HyPER_HE_Score_NOSYS": "scalar_float",
        **_p4("TtbarDiLepton_HyPER", ("Top_b_p4", "Top_lep_p4", "Topbar_bbar_p4", "Topbar_lepbar_p4")),
        **_p4("TtbarDiLepton_HyPER_plus_EM", ("Top", "Tbar")),
        **_p4("TtbarDiLepton_HyPER_plus_NW", ("Top", "Tbar")),
        **_p4("TtbarDiLepton_HyPER_plus_NuFlows", ("Top", "Tbar")),
    }
    vyper = {
        "TtbarDiLepton_VyPER_Top_IDs_NOSYS": "vector_int",
        "TtbarDiLepton_VyPER_Top_Indices_NOSYS": "vector_int",
        "TtbarDiLepton_VyPER_Top_Score_NOSYS": "scalar_float",
        "TtbarDiLepton_VyPER_AntiTop_IDs_NOSYS": "vector_int",
        "TtbarDiLepton_VyPER_AntiTop_Indices_NOSYS": "vector_int",
        "TtbarDiLepton_VyPER_AntiTop_Score_NOSYS": "scalar_float",
        "TtbarDiLepton_VyPER_Top1_IDs_NOSYS": "vector_int",
        "TtbarDiLepton_VyPER_Top1_Indices_NOSYS": "vector_int",
        "TtbarDiLepton_VyPER_Top1_Score_NOSYS": "scalar_float",
        "TtbarDiLepton_VyPER_Top1_Neutrino_NOSYS": "vector_float",
        "TtbarDiLepton_VyPER_Top2_IDs_NOSYS": "vector_int",
        "TtbarDiLepton_VyPER_Top2_Indices_NOSYS": "vector_int",
        "TtbarDiLepton_VyPER_Top2_Score_NOSYS": "scalar_float",
        "TtbarDiLepton_VyPER_Top2_Neutrino_NOSYS": "vector_float",
        **_p4("TtbarDiLepton_VyPER", ("Top", "Tbar", "Top_b_p4", "Top_lep_p4", "Top_nu_p4", "Topbar_bbar_p4", "Topbar_lepbar_p4", "Topbar_nubar_p4")),
        **_scalars("TtbarDiLepton_VyPER", ("nu_pt", "nu_eta", "nu_phi", "nu_e", "nubar_pt", "nubar_eta", "nubar_phi", "nubar_e")),
    }
    return {
        "DA": da,
        "EM": em,
        "NW": nw,
        "NuFlows": nuflows,
        "HyPER": hyper,
        "VyPER": vyper,
    }


METHOD_BRANCH_TYPES = _method_specs()
SPIN_BRANCH_TYPES = {
    f"spin_TtbarDiLepton_{method}_{observable}_NOSYS": "scalar_float"
    for method in ("DA", "EM", "NW", "VyPER", "HyPER_plus_EM", "HyPER_plus_NW", "HyPER_plus_NuFlows")
    for observable in SPIN_OBSERVABLES
}


def parse_table(path: Path, columns: int, label: str) -> tuple[list[list[str]], list[str]]:
    rows: list[list[str]] = []
    errors: list[str] = []
    if not path.exists():
        return rows, [f"{label} file does not exist: {path}"]
    for line_number, raw in enumerate(path.read_text().splitlines(), start=1):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) != columns:
            errors.append(f"{label}:{line_number}: expected {columns} fields, found {len(fields)}")
            continue
        rows.append(fields)
    return rows, errors


def inspect_metadata(metadata_dir: Path, shape_only: bool, check_paths: bool) -> tuple[dict[str, Any], list[str]]:
    issues: list[str] = []
    file_rows, errors = parse_table(metadata_dir / "filelist.txt", 4, "filelist")
    issues.extend(errors)
    sow_rows, errors = parse_table(metadata_dir / "sum_of_weights.txt", 5, "sum_of_weights")
    issues.extend(errors)

    path_counts = Counter(row[3] for row in file_rows)
    duplicate_paths = sorted(path for path, count in path_counts.items() if count > 1)
    if duplicate_paths:
        issues.append(f"{len(duplicate_paths)} duplicated file paths in filelist")

    group_counts = Counter(tuple(row[:3]) for row in file_rows)
    campaigns_by_dsid: dict[str, list[str]] = defaultdict(list)
    for dsid, campaign, _simulation in group_counts:
        campaigns_by_dsid[dsid].append(campaign)
    for dsid in campaigns_by_dsid:
        campaigns_by_dsid[dsid] = sorted(set(campaigns_by_dsid[dsid]))

    sow_map: dict[tuple[str, str, str, str], float] = {}
    sow_key_counts: Counter[tuple[str, str, str, str]] = Counter()
    nonfinite_sow: list[tuple[int, str, str]] = []
    zero_sow: list[tuple[int, str, str]] = []
    negative_sow: list[tuple[int, str, str]] = []
    invalid_sow: list[str] = []
    for line_number, row in enumerate(sow_rows, start=1):
        key = tuple(row[:4])
        sow_key_counts[key] += 1
        try:
            value = float(row[4])
        except ValueError:
            invalid_sow.append(f"sum_of_weights row {line_number}: {row[4]!r} is not numeric")
            continue
        if not math.isfinite(value):
            nonfinite_sow.append((line_number, row[0], row[1]))
        elif value == 0:
            zero_sow.append((line_number, row[0], row[1]))
        else:
            sow_map[key] = value
            if value < 0:
                negative_sow.append((line_number, row[0], row[1]))
    if invalid_sow:
        issues.extend(invalid_sow)
    duplicate_sow_keys = sorted(key for key, count in sow_key_counts.items() if count > 1)
    if duplicate_sow_keys:
        issues.append(f"{len(duplicate_sow_keys)} duplicate DSID/campaign/simulation/systematic sum-of-weights keys")
    if nonfinite_sow:
        issues.append(f"{len(nonfinite_sow)} non-finite sum-of-weights values")
    if zero_sow:
        issues.append(f"{len(zero_sow)} zero sum-of-weights values")

    file_groups = set(group_counts)
    sow_groups = {key[:3] for key, value in sow_map.items() if math.isfinite(value)}
    missing_sow_groups = sorted(file_groups - sow_groups)
    extra_sow_groups = sorted(sow_groups - file_groups)
    # A duplicate SOW key is still listed as an error even though the set view
    # above is useful for detecting missing groups.
    if missing_sow_groups and not shape_only:
        issues.append(f"{len(missing_sow_groups)} filelist DSID/campaign/simulation groups have no finite SOW entry")

    missing_paths: list[str] = []
    if check_paths:
        for row in file_rows:
            if not os.path.isfile(row[3]):
                missing_paths.append(row[3])
        if missing_paths:
            issues.append(f"{len(missing_paths)} filelist paths are missing or not regular files")

    kk_groups = sorted(group for group in file_groups if group[0].startswith(("700", "701")))
    kk_missing_sow = sorted(set(kk_groups) - sow_groups)
    kk_dsids = sorted({group[0] for group in kk_groups})
    expected_kk_campaigns = ("mc20a", "mc20d", "mc20e")
    observed_kk_pairs = {(group[0], group[1]) for group in kk_groups}
    kk_campaign_gaps = sorted(
        (dsid, campaign)
        for dsid in kk_dsids
        for campaign in expected_kk_campaigns
        if (dsid, campaign) not in observed_kk_pairs
    )
    metadata = {
        "filelist_rows": len(file_rows),
        "filelist_groups": len(file_groups),
        "filelist_dsids": len({row[0] for row in file_rows}),
        "campaign_rows": dict(sorted(Counter(row[1] for row in file_rows).items())),
        "duplicate_file_paths": duplicate_paths,
        "missing_file_paths": missing_paths,
        "missing_file_paths_count": len(missing_paths),
        "file_groups": [list(group) for group in sorted(file_groups)],
        "campaigns_by_dsid": campaigns_by_dsid,
        "sow_rows": len(sow_rows),
        "sow_unique_keys": len(sow_key_counts),
        "sow_duplicate_keys": [list(key) for key in duplicate_sow_keys],
        "sow_nonfinite_count": len(nonfinite_sow),
        "sow_zero_count": len(zero_sow),
        "sow_negative_count": len(negative_sow),
        "sow_min_finite_nonzero": min(sow_map.values()) if sow_map else None,
        "sow_max_finite_nonzero": max(sow_map.values()) if sow_map else None,
        "missing_sow_groups": [list(group) for group in missing_sow_groups],
        "extra_sow_groups": [list(group) for group in extra_sow_groups],
        "kk_style_dsid_campaign_groups": [list(group) for group in kk_groups],
        "kk_style_missing_sow_groups": [list(group) for group in kk_missing_sow],
        "kk_expected_campaigns": list(expected_kk_campaigns),
        "kk_campaign_gaps": [list(pair) for pair in kk_campaign_gaps],
        "shape_only": shape_only,
    }
    return metadata, issues


def choose_root_file(
    metadata_dir: Path,
    dsid: str,
    campaign: str,
    override: str | None,
    tree_name: str,
) -> tuple[str | None, int | None, list[str]]:
    candidates: list[str] = []
    if override:
        candidates = [override]
    else:
        rows, errors = parse_table(metadata_dir / "filelist.txt", 4, "filelist")
        if errors:
            return None, None, errors
        candidates = [row[3] for row in rows if row[0] == dsid and row[1] == campaign]
    problems: list[str] = []
    if not candidates:
        return None, None, [f"no filelist entry for DSID {dsid}, campaign {campaign}"]
    for path in candidates:
        if not os.path.isfile(path):
            problems.append(f"missing file: {path}")
            continue
        try:
            with uproot.open(path) as root_file:
                if tree_name not in root_file:
                    problems.append(f"{path}: tree {tree_name!r} absent")
                    continue
                tree = root_file[tree_name]
                entries = int(tree.num_entries)
                if entries <= 0:
                    problems.append(f"{path}: tree {tree_name!r} has zero entries")
                    continue
                return path, entries, problems
        except Exception as error:  # report candidate failures while trying next file
            problems.append(f"{path}: uproot open failed: {error}")
    return None, None, problems + [f"no nonempty readable {tree_name!r} tree for DSID {dsid}, campaign {campaign}"]


def type_matches(typename: str, expected: str) -> bool:
    name = typename.replace(" ", "")
    lower = name.lower()
    if expected == "p4":
        return "lorentzvector" in lower and "ptetaphim" in lower
    if expected.startswith("vector_"):
        if not (lower.startswith("std::vector<") or lower.startswith("vector<")):
            return False
        element = lower.split("<", 1)[1].split(">", 1)[0]
        if expected == "vector_float":
            return any(token in element for token in ("float", "double"))
        if expected == "vector_int":
            return any(token in element for token in ("int", "short", "char", "long"))
        return False
    if "vector<" in lower or "lorentzvector" in lower:
        return False
    if expected == "scalar_float":
        return any(token in lower for token in ("float", "double"))
    if expected == "scalar_int":
        return any(token in lower for token in ("int", "short", "char", "long", "bool"))
    return False


def p4_coordinate_path(branch_name: str, component: str) -> str:
    return f"{branch_name}/fCoordinates/fCoordinates.f{component}"


def _finite_count(array: Any) -> tuple[int, int]:
    if isinstance(array, ak.Array):
        flat = ak.flatten(array, axis=None)
        values = np.asarray(ak.to_numpy(flat))
    else:
        values = np.asarray(array)
    if values.dtype.kind not in "fciub":
        return int(values.size), 0
    finite = np.isfinite(values)
    return int(values.size), int(np.count_nonzero(~finite))


def inspect_schema(path: str, dsid: str, campaign: str, tree_name: str, sample_events: int) -> tuple[dict[str, Any], list[str]]:
    issues: list[str] = []
    report: dict[str, Any] = {
        "dsid": dsid,
        "campaign": campaign,
        "path": path,
        "tree": tree_name,
    }
    try:
        root_file = uproot.open(path)
    except Exception as error:
        return report, [f"{path}: unable to open ROOT file: {error}"]
    with root_file:
        if tree_name not in root_file:
            return report, [f"{path}: required tree {tree_name!r} is absent"]
        tree = root_file[tree_name]
        if not isinstance(tree, uproot.behaviors.TTree.TTree):
            return report, [f"{path}: {tree_name!r} is {type(tree).__name__}, expected a TTree"]
        entries = int(tree.num_entries)
        report["entries"] = entries
        if entries <= 0:
            issues.append(f"{path}: tree {tree_name!r} has zero entries")

        # Use only actual top-level branches for type checks. Recursive leaf
        # paths are inspected explicitly for split Lorentz-vector contents.
        branch_names = set(tree.keys(recursive=False))
        recursive_names = set(tree.keys(recursive=True))
        type_specs: dict[str, str] = {
            **OBJECT_BRANCH_TYPES,
            **TRUTH_BRANCH_TYPES,
            **{branch: kind for family in METHOD_BRANCH_TYPES.values() for branch, kind in family.items()},
            **SPIN_BRANCH_TYPES,
        }
        missing = sorted(set(type_specs) - branch_names)
        mismatched_types: list[dict[str, str]] = []
        present_specs = {name: expected for name, expected in type_specs.items() if name in branch_names}
        for name, expected in present_specs.items():
            actual = str(tree[name].typename)
            if not type_matches(actual, expected):
                mismatched_types.append({"branch": name, "expected": expected, "actual": actual})
        report["required_branch_count"] = len(type_specs)
        report["missing_branches"] = missing
        report["type_mismatches"] = mismatched_types
        if missing:
            issues.append(f"{len(missing)} required raw branches missing in {path}")
        if mismatched_types:
            issues.append(f"{len(mismatched_types)} required branches have incompatible types in {path}")

        families = {
            "objects": set(OBJECT_BRANCH_TYPES),
            "truth_and_ids": set(TRUTH_BRANCH_TYPES),
            **{f"native_{name}": set(specs) for name, specs in METHOD_BRANCH_TYPES.items()},
            "spin": set(SPIN_BRANCH_TYPES),
        }
        family_report: dict[str, dict[str, Any]] = {}
        for family, required in families.items():
            absent = sorted(required - branch_names)
            family_report[family] = {
                "required": len(required),
                "present": len(required) - len(absent),
                "complete": not absent,
                "missing": absent,
            }
        report["families"] = family_report

        # Every required branch should have one entry per tree event.
        entry_mismatch: list[dict[str, int | str]] = []
        for name in present_specs:
            branch_entries = int(tree[name].num_entries)
            if branch_entries != entries:
                entry_mismatch.append({"branch": name, "branch_entries": branch_entries, "tree_entries": entries})
        report["entry_mismatches"] = entry_mismatch
        if entry_mismatch:
            issues.append(f"{len(entry_mismatch)} required branches have event-count mismatch in {path}")

        nread = min(max(entries, 0), max(1, sample_events))
        report["sampled_entries"] = nread
        finite_report: dict[str, dict[str, int]] = {}
        finite_checks: dict[str, str] = {}
        for name, expected in present_specs.items():
            if expected == "p4":
                for component in ("pt", "eta", "phi", "m"):
                    leaf = p4_coordinate_path(name, component)
                    if leaf in recursive_names:
                        finite_checks[leaf] = "p4_leaf"
            elif expected in ("scalar_float", "vector_float"):
                finite_checks[name] = expected

        for name, expected in finite_checks.items():
            try:
                arr = tree[name].array(entry_stop=nread, library="ak" if expected in ("vector_float", "vector_int") else "np")
                total, bad = _finite_count(arr)
                finite_report[name] = {"values": total, "nonfinite": bad}
            except Exception as error:
                finite_report[name] = {"values": 0, "nonfinite": -1}
                issues.append(f"{path}: cannot read sample values from {name}: {error}")
        report["finite_sample"] = finite_report

        core_names = set(OBJECT_BRANCH_TYPES) | set(TRUTH_BRANCH_TYPES)
        core_nonfinite = [name for name in core_names if name in finite_report and finite_report[name]["nonfinite"] > 0]
        report["core_nonfinite_branches"] = sorted(core_nonfinite)
        if core_nonfinite:
            issues.append(f"{len(core_nonfinite)} core object/truth branches contain non-finite sampled values in {path}")

        # Check jagged length agreement event-by-event for each object record.
        collection_fields = {
            "jet": ("jet_pt_NOSYS", "jet_eta", "jet_phi", "jet_e_NOSYS", "jet_GN2v01_Continuous_quantile"),
            "electron": ("el_pt_NOSYS", "el_eta", "el_phi", "el_e_NOSYS", "el_charge"),
            "muon": ("mu_pt_NOSYS", "mu_eta", "mu_phi", "mu_e_NOSYS", "mu_charge"),
        }
        length_report: dict[str, dict[str, int]] = {}
        for collection, fields in collection_fields.items():
            if not all(name in branch_names for name in fields):
                continue
            arrays = tree.arrays(list(fields), entry_stop=nread, library="ak")
            lengths = [ak.num(arrays[name], axis=1) for name in fields]
            mismatch_mask = ak.Array(np.zeros(nread, dtype=np.bool_))
            for size in lengths[1:]:
                mismatch_mask = mismatch_mask | (size != lengths[0])
            mismatched_events = int(ak.sum(mismatch_mask))
            length_report[collection] = {"events_checked": nread, "events_with_length_mismatch": mismatched_events}
            if mismatched_events:
                issues.append(f"{path}: {collection} fields have differing jagged lengths in {mismatched_events}/{nread} sampled events")
        report["jagged_length_checks"] = length_report

        spin_out_of_range: dict[str, dict[str, float | int]] = {}
        spin_sentinel_counts: dict[str, int] = {}
        for name in SPIN_BRANCH_TYPES:
            if name not in branch_names:
                continue
            try:
                values = np.asarray(tree[name].array(entry_stop=nread, library="np"), dtype=np.float64)
                finite = values[np.isfinite(values)]
                if finite.size:
                    sentinels = finite[finite == -999.0]
                    if sentinels.size:
                        spin_sentinel_counts[name] = int(sentinels.size)
                    physics_values = finite[finite != -999.0]
                    if not physics_values.size:
                        continue
                    low, high = float(np.min(physics_values)), float(np.max(physics_values))
                    if low < -1.0001 or high > 1.0001:
                        spin_out_of_range[name] = {"min": low, "max": high, "out_of_range_count": int(np.count_nonzero((physics_values < -1.0001) | (physics_values > 1.0001)))}
            except Exception as error:
                issues.append(f"{path}: cannot read spin branch {name}: {error}")
        report["spin_out_of_range"] = spin_out_of_range
        report["spin_sentinel_counts"] = spin_sentinel_counts

    return report, issues


def print_report(metadata: dict[str, Any], samples: list[dict[str, Any]], issues: list[str], shape_only: bool) -> None:
    print("METADATA PREFLIGHT")
    print(f"  filelist: {metadata['filelist_rows']:,} rows, {metadata['filelist_dsids']} DSIDs, {metadata['filelist_groups']} DSID/campaign/simulation groups")
    print(f"  campaign rows: {metadata['campaign_rows']}")
    print(f"  duplicate file paths: {len(metadata['duplicate_file_paths'])}")
    print(f"  paths checked: {metadata['filelist_rows']:,}; missing/non-file paths: {metadata['missing_file_paths_count']:,}")
    print(f"  sum-of-weights: {metadata['sow_rows']} rows, {metadata['sow_unique_keys']} unique keys, duplicate keys: {len(metadata['sow_duplicate_keys'])}")
    print(f"  SOW finite nonzero range: {metadata['sow_min_finite_nonzero']} .. {metadata['sow_max_finite_nonzero']}")
    print(f"  SOW zero/negative/non-finite values: {metadata['sow_zero_count']}/{metadata['sow_negative_count']}/{metadata['sow_nonfinite_count']}")
    print(f"  missing DSID/campaign/simulation SOW groups: {len(metadata['missing_sow_groups'])}")
    print(f"  extra SOW groups without filelist: {len(metadata['extra_sow_groups'])}")
    print(f"  700/701 DSID campaign groups: {len(metadata['kk_style_dsid_campaign_groups'])}; missing SOW: {len(metadata['kk_style_missing_sow_groups'])}")
    print(f"  700/701 DSIDs with an absent mc20a/d/e filelist combination: {len(metadata['kk_campaign_gaps'])}")
    for dsid, campaign in metadata["kk_campaign_gaps"]:
        print(f"    {dsid}: no {campaign} filelist group")
    print("  observed campaign coverage for 700/701 DSIDs:")
    for dsid in sorted(k for k in metadata["campaigns_by_dsid"] if k.startswith(("700", "701"))):
        print(f"    {dsid}: {','.join(metadata['campaigns_by_dsid'][dsid])}")
    if metadata["missing_sow_groups"]:
        print("  incomplete SOW groups (DSID campaign simulation):")
        for group in metadata["missing_sow_groups"]:
            print("    " + " ".join(group))
    if metadata["missing_file_paths"]:
        print("  missing paths (first 20):")
        for path in metadata["missing_file_paths"][:20]:
            print(f"    {path}")
    if shape_only:
        print("  shape-only: incomplete SOW groups are reported but do not block shape-only reconstruction studies")

    for sample in samples:
        print(f"\nROOT SCHEMA: DSID {sample.get('dsid')} {sample.get('campaign')}")
        print(f"  file: {sample.get('path')}")
        print(f"  tree: {sample.get('tree')} entries: {sample.get('entries')}")
        print(f"  required branch presence: {sample.get('required_branch_count', 0) - len(sample.get('missing_branches', []))}/{sample.get('required_branch_count', 0)}")
        for family, detail in sample.get("families", {}).items():
            print(f"  {family}: {detail['present']}/{detail['required']} {'OK' if detail['complete'] else 'INCOMPLETE'}")
            if detail["missing"]:
                print("    missing: " + ", ".join(detail["missing"]))
        if sample.get("type_mismatches"):
            print("  type mismatches:")
            for item in sample["type_mismatches"]:
                print(f"    {item['branch']}: expected {item['expected']}, got {item['actual']}")
        if sample.get("core_nonfinite_branches"):
            print("  core non-finite sampled branches: " + ", ".join(sample["core_nonfinite_branches"]))
        if sample.get("jagged_length_checks"):
            print(f"  object collection length checks: {sample['jagged_length_checks']}")
        nonfinite = {name: detail for name, detail in sample.get("finite_sample", {}).items() if detail["nonfinite"] > 0}
        if nonfinite:
            print("  sampled non-finite values (branch: count/total):")
            for name, detail in sorted(nonfinite.items()):
                print(f"    {name}: {detail['nonfinite']}/{detail['values']}")
        if sample.get("spin_out_of_range"):
            print("  spin values outside [-1,1] (+/-1e-4):")
            for name, limits in sorted(sample["spin_out_of_range"].items()):
                print(f"    {name}: {limits['min']:.6g} .. {limits['max']:.6g} ({limits['out_of_range_count']} values)")
        if sample.get("spin_sentinel_counts"):
            total = sum(sample["spin_sentinel_counts"].values())
            print(f"  spin -999 sentinel values: {total} across {len(sample['spin_sentinel_counts'])} branches")

    if issues:
        print("\nISSUES")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print("\nNo blocking metadata/schema issues detected.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-dir", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--sm-dsid", default="410472", help="SM representative DSID (default: 410472)")
    parser.add_argument("--sm-campaign", default="mc20d", help="SM campaign (default: mc20d)")
    parser.add_argument("--sm-file", help="Optional direct SM ROOT input override")
    parser.add_argument("--kk-dsid", default="700320", help="KK representative DSID (default: 700320)")
    parser.add_argument("--kk-campaign", default="mc20a", help="KK campaign (default: mc20a)")
    parser.add_argument("--kk-file", help="Optional direct KK ROOT input override")
    parser.add_argument("--tree", default="reco", help="Required raw event tree name (default: reco)")
    parser.add_argument("--sample-events", type=int, default=2000, help="Entries read for finite/shape checks")
    parser.add_argument("--shape-only", action="store_true", help="Report incomplete SOW metadata without requiring it for shape-only studies")
    parser.add_argument("--skip-path-check", action="store_true", help="Skip checking all filelist paths on disk")
    parser.add_argument("--strict-method-finite", action="store_true", help="Treat sampled non-finite native reconstruction/spin values as blocking issues")
    parser.add_argument("--json-report", type=Path, help="Write machine-readable report to this path")
    args = parser.parse_args()

    if args.sample_events <= 0:
        parser.error("--sample-events must be positive")

    metadata, issues = inspect_metadata(args.metadata_dir, args.shape_only, not args.skip_path_check)
    samples: list[dict[str, Any]] = []
    for label, dsid, campaign, override in (
        ("SM", args.sm_dsid, args.sm_campaign, args.sm_file),
        ("KK", args.kk_dsid, args.kk_campaign, args.kk_file),
    ):
        path, _entries, selection_problems = choose_root_file(args.metadata_dir, dsid, campaign, override, args.tree)
        if selection_problems:
            for problem in selection_problems:
                print(f"{label} candidate: {problem}")
        if path is None:
            issues.append(f"could not select a readable nonempty {label} ROOT file for DSID {dsid}, campaign {campaign}")
            samples.append({"dsid": dsid, "campaign": campaign, "path": override, "tree": args.tree})
            continue
        sample_report, sample_issues = inspect_schema(path, dsid, campaign, args.tree, args.sample_events)
        sample_report["label"] = label
        samples.append(sample_report)
        issues.extend(sample_issues)
        if args.strict_method_finite:
            nonfinite = {
                name: data
                for name, data in sample_report.get("finite_sample", {}).items()
                if name not in OBJECT_BRANCH_TYPES and name not in TRUTH_BRANCH_TYPES and data["nonfinite"] > 0
            }
            if nonfinite:
                issues.append(f"{len(nonfinite)} native reconstruction/spin branches contain sampled non-finite values in {path}")

    print_report(metadata, samples, issues, args.shape_only)
    if args.json_report:
        args.json_report.parent.mkdir(parents=True, exist_ok=True)
        payload = {"metadata": metadata, "samples": samples, "issues": issues, "shape_only": args.shape_only}
        args.json_report.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"\nJSON report: {args.json_report}")

    return 1 if issues else 0


if __name__ == "__main__":
    sys.exit(main())
