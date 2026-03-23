import uproot
import awkward as ak
import argparse


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("--infile",  required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--tree",    default="reco")
    ap.add_argument("--event",   type=int, required=True)
    args = ap.parse_args()

    with uproot.open(args.infile) as f:
        tree = f[args.tree]
        arrays = tree.arrays(library="ak")

    mask = arrays["eventNumber"] == args.event

    if ak.sum(mask) == 0:
        raise RuntimeError("Event not found")

    filtered = arrays[mask]

    # --- IMPORTANT PART ---
    branches = {name: filtered[name] for name in filtered.fields}

    with uproot.recreate(args.outfile) as fout:
        fout.mktree(args.tree, branches)
        fout[args.tree].extend(branches)

    print("Done")


if __name__ == "__main__":
    main()
