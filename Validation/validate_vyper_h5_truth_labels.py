#!/usr/bin/env python3
"""Check VyPER H5 neutrino labels against Run-2 parent-truth branches."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import uproot
import awkward as ak


def truth_neutrino(tree, entry: int, parent: str):
    """Return (PDG ID, px, py, pz) in GeV for the W neutrino from one parent."""
    for daughter in (1, 2):
        stem = f"Ttbar_MC_Wdecay{daughter}_afterFSR_from_{parent}"
        pdg = int(tree[f"{stem}_pdgId"].array(library="np", entry_start=entry, entry_stop=entry + 1)[0])
        if abs(pdg) not in (12, 14, 16):
            continue
        pt = float(tree[f"{stem}_pt"].array(library="np", entry_start=entry, entry_stop=entry + 1)[0])
        eta = float(tree[f"{stem}_eta"].array(library="np", entry_start=entry, entry_stop=entry + 1)[0])
        phi = float(tree[f"{stem}_phi"].array(library="np", entry_start=entry, entry_stop=entry + 1)[0])
        vector = np.asarray(
            [pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)], dtype=np.float64
        ) * 0.001
        return pdg, vector
    raise ValueError(f"entry {entry}: no W-decay neutrino found for parent {parent}")


def expected_slot_order(parent_neutrinos):
    """Match BuildNeutrinoInputVector: e, then mu, then remaining t/tbar parents."""
    ordered_parents = []
    for flavour in (12, 14):
        for parent in ("t", "tbar"):
            item = parent_neutrinos[parent]
            if abs(item[0]) == flavour:
                ordered_parents.append(parent)
    for parent in ("t", "tbar"):
        if parent not in ordered_parents:
            ordered_parents.append(parent)
    return [parent_neutrinos[parent] for parent in ordered_parents[:2]]


def load_h5(path: Path, expected_parity: int):
    with h5py.File(path, "r") as output:
        event_numbers = output["META/eventNumber"][:].astype(np.uint64, copy=False)
        mc_channels = output["META/mcChannelNumber"][:]
        labels = output["LABELS/NEUTRINO"][:]
        vectors = np.stack(
            [labels["px"], labels["py"], labels["pz"]], axis=-1
        ).astype(np.float64)
        electrons = output["INPUTS/ELECTRON"][:]
        muons = output["INPUTS/MUON"][:]
        electron_charge = electrons["charge"].astype(float)
        muon_charge = muons["charge"].astype(float)
        attrs = dict(output["META"].attrs)
    if np.any(event_numbers % 2 != expected_parity):
        raise AssertionError(f"{path}: event parity does not match split {expected_parity}")
    if len(np.unique(event_numbers)) != len(event_numbers):
        raise AssertionError(f"{path}: duplicate eventNumber values")
    return {
        int(event): {
            "mcChannelNumber": int(channel),
            "vectors": vector,
            "electron_charges": [float(value) for value in e_charge if np.isfinite(value)],
            "muon_charges": [float(value) for value in mu_charge if np.isfinite(value)],
        }
        for event, channel, vector, e_charge, mu_charge in zip(
            event_numbers, mc_channels, vectors, electron_charge, muon_charge
        )
    }, attrs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="VyPER FastFrames pilot ROOT")
    parser.add_argument("--raw", required=True, type=Path, help="Matching Run-2 raw input ROOT")
    parser.add_argument("--h5-even", required=True, type=Path)
    parser.add_argument("--h5-odd", required=True, type=Path)
    parser.add_argument("--report", type=Path, help="Optional JSON report destination")
    parser.add_argument("--tolerance-gev", type=float, default=0.003)
    args = parser.parse_args()

    even, attrs_even = load_h5(args.h5_even, 0)
    odd, attrs_odd = load_h5(args.h5_odd, 1)
    records = {**even, **odd}
    if len(records) != len(even) + len(odd):
        raise AssertionError("eventNumber appears in both parity split files")

    with uproot.open(args.root) as root_file, uproot.open(args.raw) as raw_file:
        output_tree = root_file["reco"]
        raw_tree = raw_file["reco"]
        root_entries = int(output_tree.num_entries)
        output_ids = output_tree["eventNumber"].array(library="np").astype(np.uint64)
        raw_ids = raw_tree["eventNumber"].array(library="np").astype(np.uint64)
        if len(np.unique(output_ids)) != len(output_ids):
            raise AssertionError("FastFrames pilot ROOT has duplicate eventNumber values")
        output_index = {int(event): index for index, event in enumerate(output_ids)}
        raw_index = {int(event): index for index, event in enumerate(raw_ids)}
        vectors = ak.to_numpy(output_tree["neutrino_vector_NOSYS"].array(library="ak"))
        vectors = np.asarray(vectors, dtype=np.float64).reshape(-1, 2, 3)

        order_matches = 0
        pair_matches = 0
        h5_root_matches = 0
        charge_checks = 0
        charge_matches = 0
        parent_flavours = Counter()
        slot_flavours = Counter()
        duplicates = 0
        missing_root = []
        missing_raw = []
        max_root_h5_difference = 0.0
        max_truth_difference = 0.0
        for event, record in records.items():
            if event not in output_index:
                missing_root.append(event)
                continue
            if event not in raw_index:
                missing_raw.append(event)
                continue
            root_row = output_index[event]
            raw_row = raw_index[event]
            root_vector = vectors[root_row]
            h5_vector = record["vectors"]
            difference = float(np.max(np.abs(root_vector - h5_vector)))
            max_root_h5_difference = max(max_root_h5_difference, difference)
            if difference <= 1.0e-6:
                h5_root_matches += 1

            parents = {
                parent: truth_neutrino(raw_tree, raw_row, parent)
                for parent in ("t", "tbar")
            }
            for pdg, _ in parents.values():
                parent_flavours[abs(pdg)] += 1
            expected = expected_slot_order(parents)
            expected_vectors = np.stack([item[1] for item in expected])
            truth_difference = float(np.max(np.abs(root_vector - expected_vectors)))
            max_truth_difference = max(max_truth_difference, truth_difference)
            if truth_difference <= args.tolerance_gev:
                order_matches += 1

            # Also count order-independent pair agreement to make an ordering
            # mismatch distinguishable from a wrong truth-neutrino source.
            direct = np.allclose(root_vector[0], expected_vectors[0], atol=args.tolerance_gev, rtol=0)
            reverse = np.allclose(root_vector[0], expected_vectors[1], atol=args.tolerance_gev, rtol=0)
            if (direct and np.allclose(root_vector[1], expected_vectors[1], atol=args.tolerance_gev, rtol=0)) or (
                reverse and np.allclose(root_vector[1], expected_vectors[0], atol=args.tolerance_gev, rtol=0)
            ):
                pair_matches += 1

            if np.allclose(root_vector[0], root_vector[1], atol=1.0e-6, rtol=0):
                duplicates += 1

            for slot, (pdg, _) in enumerate(expected):
                slot_flavours[abs(pdg)] += 1
                if abs(pdg) == 12:
                    charge_checks += 1
                    charge_matches += int(
                        any(np.sign(pdg) == np.sign(charge) for charge in record["electron_charges"])
                    )
                elif abs(pdg) == 14:
                    charge_checks += 1
                    charge_matches += int(
                        any(np.sign(pdg) == np.sign(charge) for charge in record["muon_charges"])
                    )

        if missing_root or missing_raw:
            raise AssertionError(
                f"Unmatched H5 events: missing from FastFrames ROOT={missing_root[:5]}, "
                f"missing from raw input={missing_raw[:5]}"
            )
        if h5_root_matches != len(records):
            raise AssertionError(
                f"H5 labels differ from the ROOT vector for {len(records) - h5_root_matches} events"
            )
        if order_matches != len(records) or pair_matches != len(records):
            raise AssertionError(
                f"ROOT neutrino vector disagrees with raw parent truth: ordered={order_matches}/"
                f"{len(records)}, unordered={pair_matches}/{len(records)}"
            )
        if duplicates:
            raise AssertionError(f"{duplicates} H5 events have duplicate neutrino target vectors")
        if charge_checks and charge_matches != charge_checks:
            raise AssertionError(
                f"Neutrino PDG sign disagrees with the same-flavour charged-lepton sign in "
                f"{charge_checks - charge_matches}/{charge_checks} checks"
            )

    result = {
        "root_entries": root_entries,
        "h5_events": len(records),
        "h5_even_events": len(even),
        "h5_odd_events": len(odd),
        "root_h5_event_identity_join": "eventNumber",
        "h5_neutrino_order_matches_raw_truth": order_matches,
        "h5_neutrino_pairs_match_raw_truth_unordered": pair_matches,
        "h5_targets_match_fastframes_vector": h5_root_matches,
        "max_h5_root_component_difference_GeV": max_root_h5_difference,
        "max_ordered_truth_component_difference_GeV": max_truth_difference,
        "neutrino_parent_flavour_counts_abs_pdg": dict(parent_flavours),
        "neutrino_slot_flavour_counts_abs_pdg": dict(slot_flavours),
        "direct_flavour_charge_checks": charge_checks,
        "direct_flavour_charge_matches": charge_matches,
        "duplicate_neutrino_target_events": duplicates,
        "meta_neutrino_slot_order": str(attrs_even.get("neutrino_slot_order", "")),
        "meta_neutrino_components": str(attrs_even.get("neutrino_components", "")),
        "note": "The fully-matched filter requires both truth b jets and exactly two opposite-charge selected e/mu leptons; same- and mixed-flavour events are retained.",
    }
    print(json.dumps(result, indent=2))
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
