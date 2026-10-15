#!/usr/bin/env python3
import argparse

import numpy as np
import uproot

BRANCHES = [
    "truth_nu1_px_NOSYS", "truth_nu1_py_NOSYS", "truth_nu1_pz_NOSYS",
    "truth_nu2_px_NOSYS", "truth_nu2_py_NOSYS", "truth_nu2_pz_NOSYS",
    "nuflows_nu1_px_NOSYS", "nuflows_nu1_py_NOSYS", "nuflows_nu1_pz_NOSYS",
    "nuflows_nu2_px_NOSYS", "nuflows_nu2_py_NOSYS", "nuflows_nu2_pz_NOSYS",
    "dilepassigner_nu1_px_NOSYS", "dilepassigner_nu1_py_NOSYS", "dilepassigner_nu1_pz_NOSYS",
    "dilepassigner_nu2_px_NOSYS", "dilepassigner_nu2_py_NOSYS", "dilepassigner_nu2_pz_NOSYS",
]


def pick_tree(root_file):
    for name in ("reco", "NOSYS"):
        if name in root_file:
            return root_file[name]
    for _, obj in root_file.items():
        if isinstance(obj, uproot.behaviors.TTree.TTree):
            return obj
    return root_file[root_file.keys()[0]]


def median_or_nan(values):
    if values.size == 0:
        return float("nan")
    return float(np.median(values))


def main():
    parser = argparse.ArgumentParser(description="Validate neutrino Cartesian components in reco ntuples.")
    parser.add_argument("path", help="Path to the FastFrames output ROOT file")
    args = parser.parse_args()

    with uproot.open(args.path) as root_file:
        tree = pick_tree(root_file)
        keys = set(tree.keys())

        print("Tree:", tree.name)
        print("Branch presence:")
        for branch in BRANCHES:
            print(f"  {branch}: {'OK' if branch in keys else 'MISSING'}")

        print("\nCartesian consistency:")
        for prefix in ("truth", "nuflows", "dilepassigner"):
            for slot in (1, 2):
                px_b = f"{prefix}_nu{slot}_px_NOSYS"
                py_b = f"{prefix}_nu{slot}_py_NOSYS"
                pz_b = f"{prefix}_nu{slot}_pz_NOSYS"
                pt_b = f"{prefix}_nu{slot}_pt_NOSYS"

                if not all(b in keys for b in (px_b, py_b, pz_b, pt_b)):
                    print(f"  {prefix} nu{slot}: missing one or more branches")
                    continue

                px = tree[px_b].array(library="np")
                py = tree[py_b].array(library="np")
                pz = tree[pz_b].array(library="np")
                pt = tree[pt_b].array(library="np")

                calc_pt = np.sqrt(px * px + py * py)
                mask = np.isfinite(calc_pt) & np.isfinite(pt) & (pt > 0)
                ratio = calc_pt[mask] / pt[mask]

                print(
                    f"  {prefix} nu{slot}: "
                    f"n={ratio.size}, "
                    f"median(calc_pt/pt)={median_or_nan(ratio):.6f}, "
                    f"median_px={median_or_nan(px[mask]):.4g}, "
                    f"median_py={median_or_nan(py[mask]):.4g}, "
                    f"median_pz={median_or_nan(pz[mask]):.4g}, "
                    f"median_pt={median_or_nan(pt[mask]):.4g}"
                )


if __name__ == "__main__":
    main()
