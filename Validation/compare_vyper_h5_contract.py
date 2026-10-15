#!/usr/bin/env python3
"""Compare a converter H5 output with an established VyPER H5 reference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


CONTRACT_GROUPS = ("INPUTS", "LABELS")


def _compare_dataset(current: h5py.Dataset, reference: h5py.Dataset, tolerance: float) -> dict:
    if current.shape != reference.shape:
        raise AssertionError(
            f"{current.name}: shape differs: current={current.shape}, reference={reference.shape}"
        )
    if current.dtype != reference.dtype:
        raise AssertionError(
            f"{current.name}: dtype differs: current={current.dtype}, reference={reference.dtype}"
        )
    left = current[:]
    right = reference[:]
    if left.dtype.names:
        fields = left.dtype.names
        if fields != right.dtype.names:
            raise AssertionError(f"{current.name}: structured field order differs")
        comparisons = {}
        for field in fields:
            a = left[field]
            b = right[field]
            if np.issubdtype(a.dtype, np.floating):
                if not np.array_equal(np.isnan(a), np.isnan(b)):
                    raise AssertionError(f"{current.name}.{field}: NaN masks differ")
                valid = np.isfinite(a) & np.isfinite(b)
                max_abs = float(np.max(np.abs(a[valid] - b[valid]))) if valid.any() else 0.0
                if not np.allclose(a, b, atol=tolerance, rtol=0, equal_nan=True):
                    raise AssertionError(f"{current.name}.{field}: values differ; max_abs={max_abs}")
                comparisons[field] = max_abs
            elif not np.array_equal(a, b):
                raise AssertionError(f"{current.name}.{field}: values differ")
        return {"shape": list(left.shape), "dtype_fields": list(fields), "max_abs_by_field": comparisons}

    if np.issubdtype(left.dtype, np.floating):
        if not np.array_equal(np.isnan(left), np.isnan(right)):
            raise AssertionError(f"{current.name}: NaN masks differ")
        valid = np.isfinite(left) & np.isfinite(right)
        max_abs = float(np.max(np.abs(left[valid] - right[valid]))) if valid.any() else 0.0
        if not np.allclose(left, right, atol=tolerance, rtol=0, equal_nan=True):
            raise AssertionError(f"{current.name}: values differ; max_abs={max_abs}")
        return {"shape": list(left.shape), "dtype": str(left.dtype), "max_abs": max_abs}
    if not np.array_equal(left, right):
        raise AssertionError(f"{current.name}: values differ")
    return {"shape": list(left.shape), "dtype": str(left.dtype), "max_abs": 0.0}


def compare(current_path: Path, reference_path: Path, tolerance: float) -> dict:
    with h5py.File(current_path, "r") as current, h5py.File(reference_path, "r") as reference:
        for group in CONTRACT_GROUPS:
            if group not in current or group not in reference:
                raise AssertionError(f"required H5 group {group} missing")
            if set(current[group]) != set(reference[group]):
                raise AssertionError(
                    f"{group} datasets differ: current={sorted(current[group])}, "
                    f"reference={sorted(reference[group])}"
                )
        result = {
            "current": str(current_path),
            "reference": str(reference_path),
            "reference_has_event_identity": "META" in reference,
            "current_event_count": int(current["INPUTS/JET"].shape[0]),
            "reference_event_count": int(reference["INPUTS/JET"].shape[0]),
            "datasets": {},
        }
        for group in CONTRACT_GROUPS:
            for name in sorted(current[group]):
                path = f"/{group}/{name}"
                result["datasets"][path] = _compare_dataset(
                    current[f"{group}/{name}"], reference[f"{group}/{name}"], tolerance
                )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--tolerance", type=float, default=1.0e-6)
    args = parser.parse_args()
    result = compare(args.current, args.reference, args.tolerance)
    print(json.dumps(result, indent=2))
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
