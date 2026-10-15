#!/usr/bin/env python3
"""Check event-keyed preservation of native Run-2 reconstruction branches.

This verifies that the rich FastFrames ntuple carries the selected raw native
method outputs without row reordering or serialization changes. It does not
compare the native algorithms with independent legacy implementations.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


def canonical(value):
    if isinstance(value, dict):
        return tuple((key, canonical(item)) for key, item in sorted(value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(canonical(item) for item in value)
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and np.isnan(value):
        return ("NaN",)
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", required=True, type=Path)
    parser.add_argument("--fastframes", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    with uproot.open(args.raw) as source, uproot.open(args.fastframes) as output:
        src = source["reco"]
        dst = output["reco"]
        src_keys = set(src.keys())
        dst_keys = set(dst.keys())
        source_ids = src["eventNumber"].array(library="np").astype(np.uint64)
        output_ids = dst["eventNumber"].array(library="np").astype(np.uint64)
        if len(np.unique(source_ids)) != len(source_ids):
            raise AssertionError("Raw input has duplicate eventNumber values")
        if len(np.unique(output_ids)) != len(output_ids):
            raise AssertionError("FastFrames output has duplicate eventNumber values")
        source_index = {int(event): index for index, event in enumerate(source_ids)}
        missing = [int(event) for event in output_ids if int(event) not in source_index]
        if missing:
            raise AssertionError(f"FastFrames event IDs missing from raw input: {missing[:10]}")
        source_rows = np.asarray([source_index[int(event)] for event in output_ids], dtype=int)

        candidates = sorted(
            name for name in (src_keys & dst_keys)
            if "/" not in name
            and name.endswith("_NOSYS")
            and (name.startswith("spin_TtbarDiLepton_") or name.startswith("TtbarDiLepton_"))
        )
        compared = []
        mismatches = {}
        for name in candidates:
            src_values = ak.to_list(src[name].array(library="ak"))
            dst_values = ak.to_list(dst[name].array(library="ak"))
            bad = []
            for out_row, source_row in enumerate(source_rows):
                if canonical(dst_values[out_row]) != canonical(src_values[source_row]):
                    bad.append(int(output_ids[out_row]))
                    if len(bad) == 3:
                        break
            if bad:
                mismatches[name] = bad
            compared.append(name)
        if mismatches:
            raise AssertionError(f"Native branch values differ after FastFrames copy: {mismatches}")

        report = {
            "raw_path": str(args.raw),
            "fastframes_path": str(args.fastframes),
            "raw_entries": int(src.num_entries),
            "fastframes_entries": int(dst.num_entries),
            "event_identity": "eventNumber; unique in both files",
            "matched_event_ids": int(len(output_ids)),
            "native_top_level_branches_compared": len(compared),
            "spin_branches_compared": sum(name.startswith("spin_") for name in compared),
            "method_branches_compared": sum(name.startswith("TtbarDiLepton_") for name in compared),
            "value_mismatches": 0,
            "branches": compared,
            "scope_note": "This is a raw-to-FastFrames copy/identity check, not a native-versus-legacy algorithm comparison.",
        }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
