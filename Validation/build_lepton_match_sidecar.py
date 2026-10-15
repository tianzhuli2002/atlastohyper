#!/usr/bin/env python3
"""Copy only parent-W lepton-match scalars from raw ntuples into an identity sidecar."""

import argparse
import re
from pathlib import Path

import numpy as np
import uproot


MATCH_BRANCHES = tuple(
    [f"n_lep_matched_to_Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
     for side in ("t", "tbar") for daughter in (1, 2)]
    + [f"{flavour}_dRmatched_to_Ttbar_MC_Wdecay{daughter}_afterFSR_from_{side}"
       for side in ("t", "tbar") for daughter in (1, 2) for flavour in ("el", "mu")]
)
CAMPAIGN_BY_TAG = {"r13167": "mc20a", "r13144": "mc20d", "r13145": "mc20e"}


def build_sidecar(paths, output):
    chunks = {name: [] for name in ("mcChannelNumber", "eventNumber", "campaign", *MATCH_BRANCHES)}
    for path in paths:
        campaign = next((value for tag, value in CAMPAIGN_BY_TAG.items() if tag in str(path)), None)
        if campaign is None:
            match = re.search(r"mc20[ade]", str(path))
            campaign = match.group(0) if match else None
        if campaign is None:
            raise ValueError(f"{path}: cannot identify mc20 campaign from path")
        file = uproot.open(path)
        if "reco" not in file:
            raise KeyError(f"{path}: raw reco tree is absent")
        tree = file["reco"]
        required = ("mcChannelNumber", "eventNumber", *MATCH_BRANCHES)
        absent = sorted(set(required) - set(tree.keys()))
        if absent:
            raise KeyError(f"{path}: sidecar source branches absent: {absent}")
        arrays = tree.arrays(required, library="np")
        for name in required:
            chunks[name].append(np.asarray(arrays[name]))
        chunks["campaign"].append(np.full(tree.num_entries, campaign, dtype="U5"))
    values = {name: np.concatenate(parts) for name, parts in chunks.items()}
    keys = list(zip(values["mcChannelNumber"].astype(int),
                    values["eventNumber"].astype(int), values["campaign"]))
    if len(set(keys)) != len(keys):
        raise ValueError("Raw inputs contain duplicate (mcChannelNumber,eventNumber,campaign); sidecar join is unsafe")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **values)
    print(f"{output}: {len(keys)} unique raw event identities, {len(MATCH_BRANCHES)} match scalars")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_root", nargs="+")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    build_sidecar(args.raw_root, args.output)


if __name__ == "__main__":
    main()
