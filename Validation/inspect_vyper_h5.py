#!/usr/bin/env python3
"""Summarise H5 schema, masks, labels, and feature ranges for pilot files."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def _array_summary(array: np.ndarray) -> dict:
    values = np.asarray(array)
    finite = np.isfinite(values) if np.issubdtype(values.dtype, np.number) else np.ones(values.shape, bool)
    finite_values = values[finite]
    return {
        "shape": list(values.shape),
        "dtype": str(values.dtype),
        "finite": int(np.count_nonzero(finite)),
        "nonfinite": int(values.size - np.count_nonzero(finite)),
        "min": float(np.min(finite_values)) if finite_values.size else None,
        "median": float(np.median(finite_values)) if finite_values.size else None,
        "max": float(np.max(finite_values)) if finite_values.size else None,
    }


def inspect(path: Path) -> dict:
    result = {"file": str(path), "groups": {}, "labels": {}, "metadata": {}}
    with h5py.File(path, "r") as h5:
        for group_name in ("INPUTS", "LABELS"):
            if group_name not in h5:
                raise KeyError(f"{path}: missing required group {group_name}")
            group = h5[group_name]
            result["groups"][group_name] = {}
            for name, dataset in group.items():
                values = dataset[:]
                entry = {"shape": list(dataset.shape), "dtype": str(dataset.dtype)}
                if dataset.dtype.names:
                    entry["dtype_fields"] = list(dataset.dtype.names)
                    entry["fields"] = {
                        field: _array_summary(values[field]) for field in dataset.dtype.names
                    }
                else:
                    entry.update(_array_summary(values))
                result["groups"][group_name][name] = entry

        if "META" in h5:
            meta = h5["META"]
            if "eventNumber" in meta:
                event_numbers = meta["eventNumber"][:].astype(np.uint64, copy=False)
                result["metadata"]["eventNumber"] = {
                    "count": int(len(event_numbers)),
                    "unique": int(len(np.unique(event_numbers))),
                    "even": int(np.count_nonzero(event_numbers % 2 == 0)),
                    "odd": int(np.count_nonzero(event_numbers % 2 == 1)),
                }
                if len(np.unique(event_numbers)) != len(event_numbers):
                    raise AssertionError(f"{path}: duplicated eventNumber values")
            for key, value in meta.attrs.items():
                result["metadata"][key] = value.decode() if isinstance(value, bytes) else str(value)

        jet_labels = h5["LABELS/JET"][:]
        result["labels"]["jet_class_counts"] = {
            str(label): int(np.count_nonzero(jet_labels == label)) for label in (0, 1, 2)
        }
        result["labels"]["events_with_b_and_bbar_labels"] = int(
            np.count_nonzero(np.any(jet_labels == 1, axis=1) & np.any(jet_labels == 2, axis=1))
        )
        for name in ("ELECTRON", "MUON"):
            values = h5[f"LABELS/{name}"][:]
            result["labels"][f"{name.lower()}_class_counts"] = {
                str(label): int(np.count_nonzero(values == label)) for label in (1, 2)
            }
        neutrinos = h5["LABELS/NEUTRINO"][:]
        result["labels"]["neutrinos"] = {
            field: _array_summary(neutrinos[field]) for field in neutrinos.dtype.names
        }
        global_values = h5["INPUTS/GLOBAL"][:]
        if global_values.dtype.names:
            result["event_object_counts"] = {
                field: {
                    "min": float(np.nanmin(global_values[field])),
                    "median": float(np.nanmedian(global_values[field])),
                    "max": float(np.nanmax(global_values[field])),
                }
                for field in ("njet", "nelectron", "nmuon")
            }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = {str(path): inspect(path) for path in args.files}
    output = json.dumps(result, indent=2, allow_nan=False)
    print(output)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(output + "\n")


if __name__ == "__main__":
    main()
