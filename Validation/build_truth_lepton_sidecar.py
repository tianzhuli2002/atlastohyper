#!/usr/bin/env python3
"""Parent-ordered truth charged leptons from raw ``reco`` ntuples (evaluation only).

The truth definition follows the generator W decay record, never reconstructed
matching.  For each top side the two after-FSR W daughters
``Ttbar_MC_Wdecay{1,2}_afterFSR_from_{t,tbar}_{pt,eta,phi,m,pdgId}`` are read and

* top side:     direct l+ means pdgId in {-11, -13} with sibling nu in {+12, +14},
* antitop side: direct l- means pdgId in {+11, +13} with sibling nubar in {-12, -14}.

Taus (|pdgId| = 15) and hadronic W decays are not direct e/mu and carry
``truth_spin_support = False``.  Outputs are E,px,py,pz in GeV, ordered
``truth_lplus`` (from t) and ``truth_lminus`` (from tbar); there is no
electron-slot/muon-slot ordering.

Subcommands
  inspect   one-file schema, PDG and unit report (no output written)
  build     chunked extraction for a list of files -> one .npz + .manifest.json
  merge     concatenate per-file outputs with a global duplicate-key check
  validate  recompute every counter and invariant from a written .npz

Event identity stored per row: mcChannelNumber, runNumber (if present),
eventNumber, campaign, source_file_index, source_entry.  Nothing is joined by
position anywhere in this script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

SIDES = ("t", "tbar")
FIELDS = ("pt", "eta", "phi", "m", "pdgId")
SCHEMA = "atlastohyper_truth_lepton_sidecar_v1"
TRUTH_DEFINITION = "after_fsr_W_daughter_direct_emu_parent_charge_ordered_v1"
OUTPUT_FIELDS = (
    "truth_lplus_E", "truth_lplus_px", "truth_lplus_py", "truth_lplus_pz", "truth_lplus_pdgId",
    "truth_lminus_E", "truth_lminus_px", "truth_lminus_py", "truth_lminus_pz", "truth_lminus_pdgId",
    "truth_nu_top_pdgId", "truth_nu_tbar_pdgId",
    "truth_direct_emu_top", "truth_direct_emu_tbar", "truth_spin_support",
    "truth_tau_top", "truth_tau_tbar",
)
COUNTERS = (
    "n_raw",
    "n_direct_top_e", "n_direct_top_mu", "n_direct_tbar_e", "n_direct_tbar_mu",
    "n_truth_spin_support",
    "n_tau_top", "n_tau_tbar", "n_hadronic_top", "n_hadronic_tbar",
    "n_wrong_sign", "n_missing_charged_daughter", "n_multiple_direct_charged_daughters",
    "n_sibling_neutrino_flavour_mismatch", "n_nonfinite_kinematics",
)
# Run-number convention for ATLAS MC campaigns (randomRunNumber is not used here).
CAMPAIGN_BY_RUN = {284500: "mc20a", 300000: "mc20d", 310000: "mc20e",
                   410000: "mc23a", 450000: "mc23d", 470000: "mc23e"}
CAMPAIGN_BY_TAG = {"r13167": "mc20a", "r13144": "mc20d", "r13145": "mc20e",
                   "r14622": "mc23a", "r15224": "mc23d", "r16083": "mc23e"}


def truth_branches() -> list[str]:
    return [f"Ttbar_MC_Wdecay{d}_afterFSR_from_{side}_{field}"
            for side in SIDES for d in (1, 2) for field in FIELDS]


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parent,
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unknown"


def _sha256(path: Path, limit: int | None = None) -> str:
    digest = hashlib.sha256()
    read = 0
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
            read += len(block)
            if limit is not None and read >= limit:
                break
    return digest.hexdigest()


def campaign_from(path: str, run_numbers: np.ndarray | None) -> np.ndarray:
    """Per-row campaign label; runNumber wins, the path r-tag is the fallback."""
    count = 0 if run_numbers is None else len(run_numbers)
    if run_numbers is not None and count:
        labels = np.array([CAMPAIGN_BY_RUN.get(int(r), "") for r in run_numbers], dtype="U8")
        if (labels != "").all():
            return labels
    for tag, value in CAMPAIGN_BY_TAG.items():
        if tag in str(path):
            return np.full(count, value, dtype="U8")
    match = re.search(r"mc2[03][ade]", str(path))
    if match:
        return np.full(count, match.group(0), dtype="U8")
    return np.full(count, "unknown", dtype="U8")


def classify_side(pdg1: np.ndarray, pdg2: np.ndarray, side: str) -> dict[str, np.ndarray]:
    """Classify one W decay.  Returns which daughter (0/1, -1 if none) is the direct e/mu."""
    pdg = np.stack([np.asarray(pdg1, dtype=np.int64), np.asarray(pdg2, dtype=np.int64)], axis=1)
    apdg = np.abs(pdg)
    expected_sign = -1 if side == "t" else +1  # e+ = -11 from t (W+), e- = +11 from tbar (W-)
    is_emu = np.isin(apdg, (11, 13))
    is_tau = apdg == 15
    is_nu = np.isin(apdg, (12, 14, 16))
    n_emu = is_emu.sum(axis=1)
    single = n_emu == 1
    charged_index = np.where(single, np.argmax(is_emu, axis=1), -1)
    rows = np.arange(len(pdg))
    safe = np.clip(charged_index, 0, 1)
    lepton_pdg = np.where(single, pdg[rows, safe], 0)
    sibling_pdg = np.where(single, pdg[rows, 1 - safe], 0)
    right_sign = np.sign(lepton_pdg) == expected_sign
    sibling_ok = (np.abs(sibling_pdg) == np.abs(lepton_pdg) + 1) & (np.sign(sibling_pdg) == -np.sign(lepton_pdg))
    direct = single & right_sign & sibling_ok
    any_tau = is_tau.any(axis=1)
    leptonic_any = (is_emu | is_tau).any(axis=1)
    hadronic = ~leptonic_any & ~is_nu.any(axis=1)
    nu_pdg = np.where(is_nu[:, 0], pdg[:, 0], np.where(is_nu[:, 1], pdg[:, 1], 0))
    return {
        "direct": direct,
        "index": np.where(direct, charged_index, -1),
        "lepton_pdg": np.where(direct, lepton_pdg, 0),
        "nu_pdg": nu_pdg,
        "electron": direct & (np.abs(lepton_pdg) == 11),
        "muon": direct & (np.abs(lepton_pdg) == 13),
        "tau": any_tau,
        "hadronic": hadronic,
        "wrong_sign": single & ~right_sign,
        "sibling_mismatch": single & right_sign & ~sibling_ok,
        "multiple": n_emu > 1,
        # a leptonic-looking W (neutrino present) with no e/mu/tau daughter
        "missing_charged": is_nu.any(axis=1) & ~leptonic_any,
    }


def four_vectors(pt, eta, phi, mass, scale: float) -> np.ndarray:
    pt = np.asarray(pt, dtype=np.float64) * scale
    mass = np.asarray(mass, dtype=np.float64) * scale
    eta = np.asarray(eta, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    px, py, pz = pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)
    energy = np.sqrt(px * px + py * py + pz * pz + np.clip(mass, 0.0, None) ** 2)
    return np.stack([energy, px, py, pz], axis=-1)


def extract(arrays: dict[str, np.ndarray], scale: float) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    """Pure extraction from flat per-event arrays; ``scale`` converts input energy units to GeV."""
    count = len(np.asarray(arrays["eventNumber"]))
    out: dict[str, np.ndarray] = {}
    counters = {name: 0 for name in COUNTERS}
    counters["n_raw"] = count
    sides = {}
    for side in SIDES:
        info = classify_side(arrays[f"Ttbar_MC_Wdecay1_afterFSR_from_{side}_pdgId"],
                             arrays[f"Ttbar_MC_Wdecay2_afterFSR_from_{side}_pdgId"], side)
        p4 = np.stack([four_vectors(*(arrays[f"Ttbar_MC_Wdecay{d}_afterFSR_from_{side}_{f}"] for f in FIELDS[:4]),
                                    scale=scale) for d in (1, 2)], axis=1)  # [N,2,4]
        chosen = p4[np.arange(count), np.clip(info["index"], 0, 1)]
        finite = np.isfinite(chosen).all(axis=1)
        counters["n_nonfinite_kinematics"] += int((info["direct"] & ~finite).sum())
        info["direct"] = info["direct"] & finite
        chosen[~info["direct"]] = np.nan
        info["p4"] = chosen
        sides[side] = info
    for side, name in (("t", "lplus"), ("tbar", "lminus")):
        info = sides[side]
        for k, comp in enumerate(("E", "px", "py", "pz")):
            out[f"truth_{name}_{comp}"] = info["p4"][:, k].astype(np.float32)
        out[f"truth_{name}_pdgId"] = info["lepton_pdg"].astype(np.int32)
    out["truth_nu_top_pdgId"] = sides["t"]["nu_pdg"].astype(np.int32)
    out["truth_nu_tbar_pdgId"] = sides["tbar"]["nu_pdg"].astype(np.int32)
    out["truth_direct_emu_top"] = sides["t"]["direct"]
    out["truth_direct_emu_tbar"] = sides["tbar"]["direct"]
    out["truth_spin_support"] = sides["t"]["direct"] & sides["tbar"]["direct"]
    out["truth_tau_top"] = sides["t"]["tau"]
    out["truth_tau_tbar"] = sides["tbar"]["tau"]
    for side in SIDES:
        info = sides[side]
        counters[f"n_direct_{side}_e"] = int(info["electron"].sum())
        counters[f"n_direct_{side}_mu"] = int(info["muon"].sum())
        counters[f"n_tau_{'top' if side == 't' else 'tbar'}"] = int(info["tau"].sum())
        counters[f"n_hadronic_{'top' if side == 't' else 'tbar'}"] = int(info["hadronic"].sum())
        counters["n_wrong_sign"] += int(info["wrong_sign"].sum())
        counters["n_missing_charged_daughter"] += int(info["missing_charged"].sum())
        counters["n_multiple_direct_charged_daughters"] += int(info["multiple"].sum())
        counters["n_sibling_neutrino_flavour_mismatch"] += int(info["sibling_mismatch"].sum())
    counters["n_truth_spin_support"] = int(out["truth_spin_support"].sum())
    return out, counters


def detect_scale(pt_values: np.ndarray, unit: str) -> tuple[float, str]:
    if unit == "GeV":
        return 1.0, "declared GeV"
    if unit == "MeV":
        return 1.0e-3, "declared MeV"
    finite = np.asarray(pt_values, dtype=np.float64)
    finite = finite[np.isfinite(finite) & (finite > 0)]
    if not len(finite):
        raise ValueError("Cannot infer units: no finite positive W-daughter pT")
    median = float(np.median(finite))
    return (1.0e-3, f"auto MeV (median pT {median:.1f})") if median > 1000.0 else (1.0, f"auto GeV (median pT {median:.2f})")


def _open_tree(path: str, tree_name: str):
    import uproot
    handle = uproot.open(path)
    if tree_name not in handle:
        raise KeyError(f"{path}: tree {tree_name!r} absent; trees={[k for k, c in handle.classnames().items() if c == 'TTree']}")
    return handle[tree_name]


def inspect(path: str, tree_name: str, entries: int) -> dict:
    import uproot
    handle = uproot.open(path)
    trees = {key: cls for key, cls in handle.classnames().items() if cls == "TTree"}
    tree = handle[tree_name]
    keys = set(tree.keys())
    wanted = ["mcChannelNumber", "eventNumber", "runNumber", *truth_branches()]
    present = {name: (name in keys) for name in wanted}
    report = {"file": path, "trees": trees, "tree": tree_name, "entries": int(tree.num_entries),
              "branches_present": present,
              "missing_required": [n for n, ok in present.items() if not ok and n != "runNumber"],
              "related_truth_branches": sorted(k for k in keys if "Wdecay" in k)[:80]}
    readable = [n for n in wanted if present[n]]
    arrays = tree.arrays(readable, entry_stop=min(entries, tree.num_entries), library="np")
    for side in SIDES:
        for d in (1, 2):
            name = f"Ttbar_MC_Wdecay{d}_afterFSR_from_{side}_pdgId"
            if name in arrays:
                values, counts = np.unique(np.asarray(arrays[name]).astype(np.int64), return_counts=True)
                report[f"pdg_{side}_daughter{d}"] = {int(v): int(c) for v, c in zip(values, counts)}
            name = f"Ttbar_MC_Wdecay{d}_afterFSR_from_{side}_pt"
            if name in arrays:
                report[f"pt_{side}_daughter{d}_quantiles"] = np.nanquantile(np.asarray(arrays[name], dtype=float), [0.1, 0.5, 0.9]).tolist()
    if "runNumber" in arrays:
        report["runNumber_values"] = sorted({int(v) for v in np.asarray(arrays["runNumber"])})[:10]
    if not report["missing_required"]:
        scale, note = detect_scale(arrays["Ttbar_MC_Wdecay1_afterFSR_from_t_pt"], "auto")
        _, counters = extract({k: np.asarray(v) for k, v in arrays.items()}, scale)
        report["unit_detection"] = note
        report["counters_first_entries"] = counters
    return report


def build(paths: list[str], output: Path, tree_name: str, unit: str, step_size: str,
          campaign: str | None, hash_bytes: int | None) -> dict:
    import uproot
    output = output.expanduser().resolve()
    manifest_path = output.with_suffix(".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    columns: dict[str, list[np.ndarray]] = {}
    counters = {name: 0 for name in COUNTERS}
    files = []
    scale = note = None
    start = time.monotonic()
    for file_index, path in enumerate(paths):
        tree = _open_tree(path, tree_name)
        keys = set(tree.keys())
        required = ["mcChannelNumber", "eventNumber", *truth_branches()]
        missing = [n for n in required if n not in keys]
        if missing:
            raise KeyError(f"{path}: required branches absent: {missing}")
        read = required + (["runNumber"] if "runNumber" in keys else [])
        entry_offset = 0
        for chunk in uproot.iterate(f"{path}:{tree_name}", read, step_size=step_size, library="np"):
            chunk = {k: np.asarray(v) for k, v in chunk.items()}
            if scale is None:
                scale, note = detect_scale(chunk["Ttbar_MC_Wdecay1_afterFSR_from_t_pt"], unit)
            n = len(chunk["eventNumber"])
            values, chunk_counters = extract(chunk, scale)
            for name, value in chunk_counters.items():
                counters[name] += value
            run = chunk.get("runNumber")
            labels = np.full(n, campaign, dtype="U8") if campaign else campaign_from(path, run)
            ids = {
                "mcChannelNumber": chunk["mcChannelNumber"].astype(np.int64),
                "runNumber": (run.astype(np.int64) if run is not None else np.full(n, -1, dtype=np.int64)),
                "eventNumber": chunk["eventNumber"].astype(np.uint64).astype(np.int64),
                "campaign": labels,
                "source_file_index": np.full(n, file_index, dtype=np.int32),
                "source_entry": np.arange(entry_offset, entry_offset + n, dtype=np.int64),
            }
            entry_offset += n
            for name, value in {**ids, **values}.items():
                columns.setdefault(name, []).append(value)
        if entry_offset != tree.num_entries:
            raise ValueError(f"{path}: iterated {entry_offset} rows but tree has {tree.num_entries}")
        files.append({"index": file_index, "path": str(path), "entries": int(tree.num_entries),
                      "size_bytes": os.path.getsize(path) if os.path.exists(path) else None,
                      "sha256" + ("_first_bytes" if hash_bytes else ""):
                          (_sha256(Path(path), hash_bytes) if os.path.exists(path) else None)})
        print(json.dumps({"file": file_index + 1, "of": len(paths), "rows": entry_offset,
                          "elapsed_s": round(time.monotonic() - start, 1)}), flush=True)
    arrays = {name: np.concatenate(parts) for name, parts in columns.items()}
    duplicates = _duplicate_count(arrays)
    if duplicates:
        raise ValueError(f"{duplicates} duplicate (mcChannelNumber, runNumber/campaign, eventNumber) keys; join would be unsafe")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp-{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    os.replace(temporary, output)
    manifest = {
        "schema": SCHEMA, "truth_definition": TRUTH_DEFINITION,
        "units": "GeV; E,px,py,pz", "unit_detection": note, "input_scale_to_gev": scale,
        "ordering": "truth_lplus from top W+, truth_lminus from antitop W-",
        "tau_policy": "taus are not direct e/mu; truth_spin_support=False",
        "repository_sha": _git_sha(), "script_sha256": _sha256(Path(__file__)),
        "tree": tree_name, "files": files, "rows": int(len(arrays["eventNumber"])),
        "counters": counters, "duplicate_keys": duplicates,
        "output": str(output), "output_sha256": _sha256(output),
        "event_identity": ["mcChannelNumber", "runNumber", "campaign", "eventNumber"],
        "created_unix": int(time.time()),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def _duplicate_count(arrays: dict[str, np.ndarray]) -> int:
    _, campaign_code = np.unique(arrays["campaign"].astype("U8"), return_inverse=True)
    key = np.stack([np.asarray(arrays["mcChannelNumber"], dtype=np.int64), np.asarray(arrays["runNumber"], dtype=np.int64),
                    campaign_code.reshape(-1).astype(np.int64), np.asarray(arrays["eventNumber"], dtype=np.int64)], axis=1)
    if not len(key):
        return 0
    _, counts = np.unique(key, axis=0, return_counts=True)
    return int((counts > 1).sum())


def merge(inputs: list[Path], output: Path) -> dict:
    output = output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(output)
    parts, manifests = [], []
    for path in inputs:
        with np.load(path, allow_pickle=False) as data:
            parts.append({k: data[k] for k in data.files})
        manifests.append(json.loads(Path(path).with_suffix(".manifest.json").read_text()))
    offset = 0
    for part, manifest in zip(parts, manifests):
        part["source_file_index"] = part["source_file_index"] + offset
        offset += len(manifest["files"])
    arrays = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    duplicates = _duplicate_count(arrays)
    if duplicates:
        raise ValueError(f"Merged sidecar has {duplicates} duplicate keys")
    temporary = output.with_name(f".{output.name}.tmp-{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    os.replace(temporary, output)
    counters = {name: sum(m["counters"][name] for m in manifests) for name in COUNTERS}
    files = [dict(f, index=f["index"] + sum(len(m["files"]) for m in manifests[:i]))
             for i, m in enumerate(manifests) for f in m["files"]]
    manifest = {**{k: manifests[0][k] for k in ("schema", "truth_definition", "units", "ordering", "tau_policy", "tree")},
                "repository_sha": _git_sha(), "script_sha256": _sha256(Path(__file__)),
                "merged_from": [str(p) for p in inputs], "files": files, "rows": int(len(arrays["eventNumber"])),
                "counters": counters, "duplicate_keys": 0, "output": str(output), "output_sha256": _sha256(output),
                "event_identity": manifests[0]["event_identity"], "created_unix": int(time.time())}
    output.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def validate(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        arrays = {k: data[k] for k in data.files}
    manifest = json.loads(Path(path).with_suffix(".manifest.json").read_text())
    support = arrays["truth_spin_support"]
    lplus = np.stack([arrays[f"truth_lplus_{c}"] for c in ("E", "px", "py", "pz")], axis=1).astype(np.float64)
    lminus = np.stack([arrays[f"truth_lminus_{c}"] for c in ("E", "px", "py", "pz")], axis=1).astype(np.float64)
    mass2 = lambda v: v[:, 0] ** 2 - (v[:, 1:] ** 2).sum(axis=1)
    checks = {
        "rows": int(len(support)),
        "rows_match_manifest": int(len(support)) == manifest["rows"],
        "support_count_matches_manifest": int(support.sum()) == manifest["counters"]["n_truth_spin_support"],
        "support_fraction": float(support.mean()) if len(support) else None,
        "finite_on_support": bool(np.isfinite(lplus[support]).all() and np.isfinite(lminus[support]).all()),
        "nan_off_support_lplus": bool(np.isnan(lplus[~arrays["truth_direct_emu_top"]]).all()),
        "lplus_pdg_negative_on_support": bool((arrays["truth_lplus_pdgId"][support] < 0).all()),
        "lminus_pdg_positive_on_support": bool((arrays["truth_lminus_pdgId"][support] > 0).all()),
        "sibling_flavour_on_support": bool(
            (np.abs(arrays["truth_nu_top_pdgId"][support]) == np.abs(arrays["truth_lplus_pdgId"][support]) + 1).all()
            and (np.abs(arrays["truth_nu_tbar_pdgId"][support]) == np.abs(arrays["truth_lminus_pdgId"][support]) + 1).all()),
        "max_abs_lepton_mass_gev": float(np.sqrt(np.abs(np.concatenate([mass2(lplus[support]), mass2(lminus[support])]))).max()) if support.any() else None,
        "duplicate_keys": _duplicate_count(arrays),
        "flavour_pairs_on_support": {f"{a},{b}": int(((arrays['truth_lplus_pdgId'] == a) & (arrays['truth_lminus_pdgId'] == b) & support).sum())
                                     for a in (-11, -13) for b in (11, 13)},
        "campaigns": {str(k): int(v) for k, v in zip(*np.unique(arrays["campaign"], return_counts=True))},
        "dsids": {str(k): int(v) for k, v in zip(*np.unique(arrays["mcChannelNumber"], return_counts=True))},
        "counters": manifest["counters"],
    }
    must = ["rows_match_manifest", "support_count_matches_manifest", "finite_on_support", "nan_off_support_lplus",
            "lplus_pdg_negative_on_support", "lminus_pdg_positive_on_support", "sibling_flavour_on_support"]
    checks["passed"] = all(checks[name] for name in must) and checks["duplicate_keys"] == 0 \
        and manifest["counters"]["n_wrong_sign"] == 0 and manifest["counters"]["n_multiple_direct_charged_daughters"] == 0 \
        and (checks["max_abs_lepton_mass_gev"] is None or checks["max_abs_lepton_mass_gev"] < 0.2)
    return checks


def _read_list(values: list[str]) -> list[str]:
    paths = []
    for value in values:
        if value.endswith(".txt") and Path(value).is_file():
            paths += [line.split()[0] for line in Path(value).read_text().splitlines() if line.strip() and not line.startswith("#")]
        else:
            paths.append(value)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("inspect"); p.add_argument("root"); p.add_argument("--tree", default="reco"); p.add_argument("--entries", type=int, default=20000)
    p = sub.add_parser("build"); p.add_argument("inputs", nargs="+", help="ROOT files or .txt filelists")
    p.add_argument("--output", type=Path, required=True); p.add_argument("--tree", default="reco")
    p.add_argument("--unit", choices=("auto", "MeV", "GeV"), default="auto"); p.add_argument("--step-size", default="200 MB")
    p.add_argument("--campaign", default=None, help="override; default from runNumber, then path r-tag")
    p.add_argument("--hash-bytes", type=int, default=64 << 20, help="hash only the first N bytes of each input (0 = full file)")
    p = sub.add_parser("merge"); p.add_argument("inputs", nargs="+", type=Path); p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("validate"); p.add_argument("sidecar", type=Path); p.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.command == "inspect":
        result = inspect(args.root, args.tree, args.entries)
    elif args.command == "build":
        result = build(_read_list(args.inputs), args.output, args.tree, args.unit, args.step_size, args.campaign,
                       args.hash_bytes or None)
    elif args.command == "merge":
        result = merge(args.inputs, args.output)
    else:
        result = validate(args.sidecar)
        if args.output:
            args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    if args.command == "validate" and not result["passed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
