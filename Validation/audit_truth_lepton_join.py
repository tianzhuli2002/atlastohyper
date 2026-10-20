#!/usr/bin/env python3
"""Gate-0A identity audit: attach truth-lepton sidecar rows to an H5 population.

The canonical JER release uses ``EVENT_ID = META/eventNumber`` and records
``META/source_entry_index`` = the entry number in the single ROOT file that the
H5 converter (``H5Converters/TtbarDiLeptonTopRegressionCombined_ParentNu.py``)
read.  Neither DSID nor campaign is stored in the H5, so the join is done in two
explicit steps; nothing is ever joined by row position.

  population-keys  read the converter-input ROOT at every recorded
                   source_entry_index, require eventNumber equality for every
                   row, and emit (eventNumber, mcChannelNumber, runNumber)
                   keys in H5 row order.
  join             match those keys one-to-one against a truth-lepton sidecar
                   (from build_truth_lepton_sidecar.py).  Uses every identity
                   column both sides have; fails if any key is missing or
                   ambiguous.  Writes an EVENT_ID-keyed aligned sidecar.

The aligned output is what JER's build_complete_state_targets.py consumes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ALIGNED_SCHEMA = "truth_lepton_sidecar_aligned_v1"
TRUTH_COLUMNS = (
    "truth_lplus_E", "truth_lplus_px", "truth_lplus_py", "truth_lplus_pz", "truth_lplus_pdgId",
    "truth_lminus_E", "truth_lminus_px", "truth_lminus_py", "truth_lminus_pz", "truth_lminus_pdgId",
    "truth_direct_emu_top", "truth_direct_emu_tbar", "truth_spin_support", "truth_tau_top", "truth_tau_tbar",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_npz(path: Path, **arrays) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}.npz")
    np.savez_compressed(temporary, **arrays)
    os.replace(temporary, path)


def _integer_keys(left: list[np.ndarray], right: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Map multi-column int keys of both tables to shared dense integer codes."""
    stacked = np.concatenate([np.stack([np.asarray(c, dtype=np.int64) for c in left], axis=1),
                              np.stack([np.asarray(c, dtype=np.int64) for c in right], axis=1)])
    _, codes = np.unique(stacked, axis=0, return_inverse=True)
    codes = codes.reshape(-1)
    return codes[:len(left[0])], codes[len(left[0]):]


def population_keys(h5_path: Path, root_path: str, tree_name: str, output: Path, chunk: int) -> dict:
    import h5py
    import uproot
    if output.exists():
        raise FileExistsError(output)
    with h5py.File(h5_path, "r") as h5:
        event_number = np.asarray(h5["META/eventNumber"][:], dtype=np.int64)
        entries = np.asarray(h5["META/source_entry_index"][:], dtype=np.int64)
        attrs = {k: (v.decode() if isinstance(v, bytes) else (v.item() if hasattr(v, "item") else v)) for k, v in h5.attrs.items()}
    if len(np.unique(event_number)) != len(event_number):
        raise ValueError(f"{h5_path}: META/eventNumber is not unique; EVENT_ID cannot be eventNumber")
    if (np.diff(entries) <= 0).any():
        raise ValueError(f"{h5_path}: source_entry_index is not strictly increasing; converter provenance differs from expectation")
    tree = uproot.open(root_path)[tree_name]
    keys = set(tree.keys())
    read = ["eventNumber"] + [k for k in ("mcChannelNumber", "runNumber") if k in keys]
    if entries.max() >= tree.num_entries:
        raise ValueError(f"source_entry_index {entries.max()} beyond {root_path} entries {tree.num_entries}")
    found = {name: np.empty(len(entries), dtype=np.int64) for name in read}
    for start in range(0, len(entries), chunk):
        stop = min(len(entries), start + chunk)
        lo, hi = int(entries[start]), int(entries[stop - 1]) + 1
        block = tree.arrays(read, entry_start=lo, entry_stop=hi, library="np")
        local = entries[start:stop] - lo
        for name in read:
            found[name][start:stop] = np.asarray(block[name]).astype(np.uint64).astype(np.int64)[local] if name == "eventNumber" \
                else np.asarray(block[name]).astype(np.int64)[local]
    mismatch = int((found["eventNumber"] != event_number).sum())
    result = {
        "h5": str(h5_path), "h5_sha256": _sha256(h5_path), "h5_attrs": attrs,
        "converter_input_root": root_path, "tree": tree_name, "root_entries": int(tree.num_entries),
        "rows": int(len(entries)), "eventNumber_mismatches": mismatch,
        "identity_columns_available": read,
        "passed": mismatch == 0,
    }
    if mismatch:
        first = int(np.flatnonzero(found["eventNumber"] != event_number)[0])
        result["first_mismatch"] = {"row": first, "h5_eventNumber": int(event_number[first]),
                                    "root_eventNumber": int(found["eventNumber"][first]), "entry": int(entries[first])}
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        _atomic_npz(output, row=np.arange(len(entries), dtype=np.int64), event_id=event_number,
                    source_entry_index=entries, **{k: v for k, v in found.items() if k != "eventNumber"})
        result["output"] = str(output)
        result["output_sha256"] = _sha256(output)
    output.with_suffix(".audit.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n")
    return result


def join(keys_path: Path, sidecar_paths: list[Path], output: Path, restrict_dsid: list[int] | None) -> dict:
    if output.exists():
        raise FileExistsError(output)
    with np.load(keys_path, allow_pickle=False) as data:
        pop = {k: data[k] for k in data.files}
    side_parts, manifests = [], []
    for path in sidecar_paths:
        with np.load(path, allow_pickle=False) as data:
            side_parts.append({k: data[k] for k in data.files})
        manifests.append(json.loads(Path(path).with_suffix(".manifest.json").read_text()))
    side = {k: np.concatenate([p[k] for p in side_parts]) for k in side_parts[0]}
    if restrict_dsid:
        keep = np.isin(side["mcChannelNumber"], restrict_dsid)
        side = {k: v[keep] for k, v in side.items()}
    columns = ["eventNumber"] + [c for c in ("mcChannelNumber", "runNumber") if c in pop and (side[c] >= 0).all()]
    pop_key, side_key = _integer_keys([pop["event_id"]] + [pop[c] for c in columns[1:]],
                                      [side["eventNumber"]] + [side[c] for c in columns[1:]])
    side_counts = np.bincount(side_key, minlength=int(max(pop_key.max(initial=0), side_key.max(initial=0))) + 1)
    order = np.argsort(side_key, kind="stable")
    sorted_keys = side_key[order]
    position = np.clip(np.searchsorted(sorted_keys, pop_key), 0, max(len(sorted_keys) - 1, 0))
    found = (sorted_keys[position] == pop_key) if len(sorted_keys) else np.zeros(len(pop_key), dtype=bool)
    ambiguous = side_counts[pop_key] > 1
    report = {
        "keys": str(keys_path), "sidecars": [str(p) for p in sidecar_paths],
        "sidecar_sha256": [m.get("output_sha256") for m in manifests],
        "join_columns": columns, "restrict_dsid": restrict_dsid,
        "population_rows": int(len(pop_key)), "sidecar_rows": int(len(side_key)),
        "matched": int((found & ~ambiguous).sum()), "unmatched": int((~found).sum()),
        "ambiguous": int((found & ambiguous).sum()),
        "sidecar_dsids": {str(k): int(v) for k, v in zip(*np.unique(side["mcChannelNumber"], return_counts=True))},
        "sidecar_campaigns": {str(k): int(v) for k, v in zip(*np.unique(side["campaign"], return_counts=True))},
    }
    report["passed"] = report["unmatched"] == 0 and report["ambiguous"] == 0
    if report["passed"]:
        rows = order[position]
        aligned = {name: side[name][rows] for name in TRUTH_COLUMNS}
        matched_ids = {f"sidecar_{c}": side[c][rows] for c in ("mcChannelNumber", "runNumber", "campaign", "source_file_index", "source_entry")}
        support = aligned["truth_spin_support"]
        report["truth_spin_support_fraction"] = float(support.mean())
        report["matched_dsids"] = {str(k): int(v) for k, v in zip(*np.unique(matched_ids["sidecar_mcChannelNumber"], return_counts=True))}
        report["matched_campaigns"] = {str(k): int(v) for k, v in zip(*np.unique(matched_ids["sidecar_campaign"], return_counts=True))}
        output.parent.mkdir(parents=True, exist_ok=True)
        _atomic_npz(output, event_id=pop["event_id"], source_entry_index=pop["source_entry_index"], **aligned, **matched_ids)
        report["output"] = str(output)
        report["output_sha256"] = _sha256(output)
        report["schema"] = ALIGNED_SCHEMA
        report["event_id_hash"] = hashlib.sha256(np.sort(pop["event_id"]).astype("<i8").tobytes()).hexdigest()
    report["created_unix"] = int(time.time())
    output.with_suffix(".audit.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("population-keys")
    p.add_argument("--h5", type=Path, required=True); p.add_argument("--converter-input-root", required=True)
    p.add_argument("--tree", default="reco"); p.add_argument("--output", type=Path, required=True)
    p.add_argument("--chunk", type=int, default=200000)
    p = sub.add_parser("join")
    p.add_argument("--keys", type=Path, required=True); p.add_argument("--sidecar", type=Path, nargs="+", required=True)
    p.add_argument("--output", type=Path, required=True); p.add_argument("--restrict-dsid", type=int, nargs="*")
    args = parser.parse_args()
    if args.command == "population-keys":
        result = population_keys(args.h5, args.converter_input_root, args.tree, args.output, args.chunk)
    else:
        result = join(args.keys, args.sidecar, args.output, args.restrict_dsid)
    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    sys.exit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
