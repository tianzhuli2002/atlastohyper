#!/usr/bin/env python3
"""ttbar dilepton ROOT -> H5: reconstructed inputs, partonic top targets,
regression truth extras, VyPER-style jet/lepton assignments and neutrino truth.

Run with --target-source partonic_after_fsr -fm -b -s -ds N.
-ds is the accepted event target PER PARITY OUTPUT when splitting is enabled.
Truth neutrinos are in e-flavour, mu-flavour order (not top/antitop order).
All stored momenta and masses are GeV; jet labels use the ORIGINAL jet indices.
Requires numpy, awkward, uproot and h5py on the machine reading the ROOT file.
"""
import argparse
import difflib
import json
from pathlib import Path
import awkward as ak
import h5py
import numpy as np
import uproot
TARGET_FIELDS = [
    "top_pt",
    "top_eta",
    "top_phi",
    "top_m",
    "antitop_pt",
    "antitop_eta",
    "antitop_phi",
    "antitop_m",
]
TARGET_BRANCHES_BY_SOURCE = {
    "reco_truth": {
        "top_pt": "top_reco_truth_pt_NOSYS",
        "top_eta": "top_reco_truth_eta_NOSYS",
        "top_phi": "top_reco_truth_phi_NOSYS",
        "top_m": "top_reco_truth_m_NOSYS",
        "antitop_pt": "antitop_reco_truth_pt_NOSYS",
        "antitop_eta": "antitop_reco_truth_eta_NOSYS",
        "antitop_phi": "antitop_reco_truth_phi_NOSYS",
        "antitop_m": "antitop_reco_truth_m_NOSYS",
    },
    "partonic_after_fsr": {
        "top_pt": "t_afterFSR_pt_NOSYS",
        "top_eta": "t_afterFSR_eta_NOSYS",
        "top_phi": "t_afterFSR_phi_NOSYS",
        "top_m": "t_afterFSR_m_NOSYS",
        "antitop_pt": "tbar_afterFSR_pt_NOSYS",
        "antitop_eta": "tbar_afterFSR_eta_NOSYS",
        "antitop_phi": "tbar_afterFSR_phi_NOSYS",
        "antitop_m": "tbar_afterFSR_m_NOSYS",
    },
}
OBJECT_DTYPE = np.dtype(
    [
        ("e", np.float32),
        ("eta", np.float32),
        ("phi", np.float32),
        ("cos_phi", np.float32),
        ("sin_phi", np.float32),
        ("pt", np.float32),
        ("btag", np.float32),
        ("charge", np.float32),
        ("id", np.float32),
    ]
)
GLOBAL_DTYPE = np.dtype([
    ("njet", np.float32),
    ("nelectron", np.float32),
    ("nmuon", np.float32),
    ("met_met", np.float32),
    ("met_phi", np.float32),
    ("met_cosphi", np.float32),
    ("met_sinphi", np.float32),
    ("nBJets_90", np.float32),
    ("nBJets_85", np.float32),
    ("nBJets_77", np.float32),
    ("nBJets_70", np.float32),
    ("nBJets_65", np.float32),
])
B_COUNT_BRANCHES = {
    "nBJets_90": "nBjets_GN2v01_90_NOSYS",
    "nBJets_85": "nBjets_GN2v01_85_NOSYS",
    "nBJets_77": "nBjets_GN2v01_77_NOSYS",
    "nBJets_70": "nBjets_GN2v01_70_NOSYS",
    "nBJets_65": "nBjets_GN2v01_65_NOSYS",
}
TARGET_DTYPE = np.dtype([(name, np.float32) for name in TARGET_FIELDS])
NEUTRINO_DTYPE = np.dtype([(name, np.float32) for name in ("px", "py", "pz")])
# VyPER label convention: 0 = background; 1 = top/b or positive lepton;
# 2 = antitop/bbar or negative lepton; NaN = padded object.

TRUTH_EXTRA_DTYPE = np.dtype(
    [
        ("truth_top_px", np.float32),
        ("truth_top_py", np.float32),
        ("truth_top_pz", np.float32),
        ("truth_top_E", np.float32),
        ("truth_top_m", np.float32),
        ("truth_antitop_px", np.float32),
        ("truth_antitop_py", np.float32),
        ("truth_antitop_pz", np.float32),
        ("truth_antitop_E", np.float32),
        ("truth_antitop_m", np.float32),
        ("truth_mttbar", np.float32),
        ("truth_ttbar_px", np.float32),
        ("truth_ttbar_py", np.float32),
        ("truth_ttbar_pz", np.float32),
        ("truth_ttbar_E", np.float32),
        ("truth_ttbar_pt", np.float32),
        ("truth_ttbar_y", np.float32),
        ("truth_ttbar_phi", np.float32),
        ("truth_log_mttbar", np.float32),
        ("truth_mttbar_bin", np.int16),
        ("truth_cos_theta_star", np.float32),
        ("truth_phi_star", np.float32),
    ]
)
TRUTH_MTTBAR_BIN_EDGES = np.asarray([250.0, 600.0, 1000.0, 1500.0, np.inf], dtype=np.float64)
def _tree_candidates(root_file):
    return [
        name
        for name, obj in root_file.items()
        if hasattr(obj, "num_entries") and hasattr(obj, "keys")
    ]
def _select_tree(root_file, tree_name):
    if tree_name:
        return tree_name, root_file[tree_name]
    candidates = _tree_candidates(root_file)
    if len(candidates) == 1:
        return candidates[0].split(";")[0], root_file[candidates[0]]
    clean = [name.split(";")[0] for name in candidates]
    if clean == ["reco"]:
        return "reco", root_file["reco"]
    if "reco" in clean:
        return "reco", root_file["reco"]
    raise RuntimeError(f"Ambiguous TTree selection. Candidates: {clean}")
def _closest(branches, patterns):
    matches = []
    for pattern in patterns:
        matches.extend([name for name in branches if pattern in name])
        matches.extend(difflib.get_close_matches(pattern, branches, n=8))
    return sorted(set(matches))
def _resolve_target_branches(tree, target_source):
    if target_source not in TARGET_BRANCHES_BY_SOURCE:
        raise ValueError(
f"Unsupported target_source={target_source!r}. "
f"Supported: {sorted(TARGET_BRANCHES_BY_SOURCE)}"
        )
    branches = set(tree.keys())
    target_branch_map = dict(TARGET_BRANCHES_BY_SOURCE[target_source])
    missing_logical = [
        field
        for field in TARGET_FIELDS
        if target_branch_map[field] not in branches
    ]
    missing_branches = [target_branch_map[field] for field in missing_logical]
    if missing_branches:
        print(f"Missing required target branches for target_source={target_source!r}:")
        for field, branch in zip(missing_logical, missing_branches):
            print(f"  {field}: {branch}")
        print(
            "Closest matches:",
            _closest(
list(branches),
                ["top", "antitop", "afterFSR", "after_FSR", "reco_truth", "tbar"],
            ),
        )
        raise SystemExit(2)
    return target_branch_map
def _target_branch_list(target_branch_map):
    return [target_branch_map[field] for field in TARGET_FIELDS]
def _first_existing(branches, names):
    for name in names:
        if name in branches:
            return name
    return None
def _str_to_bool(value):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected boolean value, got {value!r}")
def _pad_jagged(array, width, fill=np.nan):
    return ak.to_numpy(
        ak.fill_none(
            ak.pad_none(array, width, axis=1, clip=True),
            fill,
        )
    ).astype(np.float32)
def _energy_scale_from_branch(tree, branch_name, entry_stop):
    arr = tree[branch_name].array(entry_stop=entry_stop, library="ak")
    finite = ak.to_numpy(ak.flatten(arr, axis=None))
    finite = np.asarray(finite, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size and np.nanmedian(finite) > 10000.0:
        return 0.001, "MeV input magnitudes detected; writing GeV"
    return 1.0, "GeV input magnitudes detected; writing GeV"
def _make_object_array(e, eta, phi, pt, btag, charge, object_id, keep_mask):
    out = np.empty((int(keep_mask.sum()), e.shape[1]), dtype=OBJECT_DTYPE)
    for name in out.dtype.names:
        out[name] = np.nan
    kept_e = e[keep_mask]
    kept_eta = eta[keep_mask]
    kept_phi = phi[keep_mask]
    kept_pt = pt[keep_mask]
    kept_btag = btag[keep_mask]
    kept_charge = charge[keep_mask]
    valid = np.isfinite(kept_pt) & (kept_pt > 0.0)
    out["e"] = kept_e
    out["eta"] = kept_eta
    out["phi"] = kept_phi
    out["cos_phi"] = np.where(valid, np.cos(kept_phi), np.nan)
    out["sin_phi"] = np.where(valid, np.sin(kept_phi), np.nan)
    out["pt"] = kept_pt
    out["btag"] = kept_btag
    out["charge"] = kept_charge
    out["id"] = np.where(valid, float(object_id), np.nan)
    return out
def _stats(name, values):
    arr = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(arr)
    if not finite.any():
        return f"{name}: finite=0.000000 min=nan max=nan mean=nan std=nan"
    return (
f"{name}: finite={finite.mean():.6f} "
f"min={np.nanmin(arr[finite]):.6g} "
f"max={np.nanmax(arr[finite]):.6g} "
f"mean={np.nanmean(arr[finite]):.6g} "
f"std={np.nanstd(arr[finite]):.6g}"
    )
def _fourvectors_from_pt_eta_phi_m(pt, eta, phi, mass):
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    E = np.sqrt(np.maximum(px**2 + py**2 + pz**2 + np.maximum(mass, 0.0) ** 2, 0.0))
    return E, px, py, pz
def _mass(E, px, py, pz):
    return np.sqrt(np.maximum(E**2 - px**2 - py**2 - pz**2, 0.0))
def _rapidity(E, pz, eps=1.0e-6):
    return 0.5 * np.log(np.maximum(E + pz, eps) / np.maximum(E - pz, eps))
def _helicity_angles(top_E, top_px, top_py, top_pz, antitop_E, antitop_px, antitop_py, antitop_pz, eps=1.0e-9):
    """Top direction in the ttbar CM frame, using the ttbar lab flight direction as helicity axis."""
    top_p = np.stack([top_px, top_py, top_pz], axis=-1).astype(np.float64)
    tt_E = np.asarray(top_E + antitop_E, dtype=np.float64)
    tt_p = np.stack([top_px + antitop_px, top_py + antitop_py, top_pz + antitop_pz], axis=-1).astype(np.float64)
    beta = tt_p / np.maximum(tt_E[..., None], eps)
    beta2 = np.sum(beta * beta, axis=-1)
    gamma = 1.0 / np.sqrt(np.maximum(1.0 - beta2, eps))
    bp = np.sum(beta * top_p, axis=-1)
    factor = np.where(beta2 > eps, ((gamma - 1.0) * bp / np.maximum(beta2, eps) - gamma * top_E), 0.0)
    top_star = top_p + factor[..., None] * beta
    qhat = top_star / np.maximum(np.linalg.norm(top_star, axis=-1)[..., None], eps)
    axis_norm = np.linalg.norm(tt_p, axis=-1)
    axis = np.where(axis_norm[..., None] > eps, tt_p / np.maximum(axis_norm[..., None], eps), np.array([0.0, 0.0, 1.0]))
    cos_theta = np.clip(np.sum(qhat * axis, axis=-1), -1.0, 1.0)
    beam = np.broadcast_to(np.array([0.0, 0.0, 1.0]), axis.shape).copy()
    x_ref = beam - np.sum(beam * axis, axis=-1)[..., None] * axis
    x_norm = np.linalg.norm(x_ref, axis=-1)
    x_ref = np.where(x_norm[..., None] > eps, x_ref / np.maximum(x_norm[..., None], eps), np.array([1.0, 0.0, 0.0]))
    y_ref = np.cross(axis, x_ref)
    phi_star = np.arctan2(np.sum(qhat * y_ref, axis=-1), np.sum(qhat * x_ref, axis=-1))
    return cos_theta.astype(np.float32), phi_star.astype(np.float32)
def _build_truth_extras(target_matrix):
    top_E, top_px, top_py, top_pz = _fourvectors_from_pt_eta_phi_m(
        target_matrix[:, 0],
        target_matrix[:, 1],
        target_matrix[:, 2],
        target_matrix[:, 3],
    )
    at_E, at_px, at_py, at_pz = _fourvectors_from_pt_eta_phi_m(
        target_matrix[:, 4],
        target_matrix[:, 5],
        target_matrix[:, 6],
        target_matrix[:, 7],
    )
    tt_E = top_E + at_E
    tt_px = top_px + at_px
    tt_py = top_py + at_py
    tt_pz = top_pz + at_pz
    mtt = _mass(tt_E, tt_px, tt_py, tt_pz)
    tt_pt = np.sqrt(np.maximum(tt_px**2 + tt_py**2, 0.0))
    tt_phi = np.arctan2(tt_py, tt_px)
    cos_theta, phi_star = _helicity_angles(top_E, top_px, top_py, top_pz, at_E, at_px, at_py, at_pz)
    out = np.empty(target_matrix.shape[0], dtype=TRUTH_EXTRA_DTYPE)
    out["truth_top_px"] = top_px.astype(np.float32)
    out["truth_top_py"] = top_py.astype(np.float32)
    out["truth_top_pz"] = top_pz.astype(np.float32)
    out["truth_top_E"] = top_E.astype(np.float32)
    out["truth_top_m"] = target_matrix[:, 3].astype(np.float32)
    out["truth_antitop_px"] = at_px.astype(np.float32)
    out["truth_antitop_py"] = at_py.astype(np.float32)
    out["truth_antitop_pz"] = at_pz.astype(np.float32)
    out["truth_antitop_E"] = at_E.astype(np.float32)
    out["truth_antitop_m"] = target_matrix[:, 7].astype(np.float32)
    out["truth_mttbar"] = mtt.astype(np.float32)
    out["truth_ttbar_px"] = tt_px.astype(np.float32)
    out["truth_ttbar_py"] = tt_py.astype(np.float32)
    out["truth_ttbar_pz"] = tt_pz.astype(np.float32)
    out["truth_ttbar_E"] = tt_E.astype(np.float32)
    out["truth_ttbar_pt"] = tt_pt.astype(np.float32)
    out["truth_ttbar_y"] = _rapidity(tt_E, tt_pz).astype(np.float32)
    out["truth_ttbar_phi"] = tt_phi.astype(np.float32)
    out["truth_log_mttbar"] = np.log(np.maximum(mtt, 1.0e-6)).astype(np.float32)
    out["truth_mttbar_bin"] = np.digitize(mtt, TRUTH_MTTBAR_BIN_EDGES[1:-1], right=False).astype(np.int16)
    out["truth_cos_theta_star"] = cos_theta
    out["truth_phi_star"] = phi_star
    return out
def _append_dataset(ds, data):
    n_new = len(data)
    if n_new == 0:
        return
    start = ds.shape[0]
    stop = start + n_new
    ds.resize(stop, axis=0)
    ds[start:stop] = data
def _target_convention(target_source):
    if target_source == "reco_truth":
        return (
            "top and antitop from FastFrames top_reco_truth_* "
            "and antitop_reco_truth_* variables"
        )
    if target_source == "partonic_after_fsr":
        return (
            "top and antitop from partonic after-FSR branches "
            "t_afterFSR_* and tbar_afterFSR_*"
        )
    return f"unknown target source: {target_source}"
def _resolve_neutrino_branches(branches):
    """Resolve canonical parent-ordered neutrino truth: slot 0=top, slot 1=antitop."""
    vector = "truth_neutrino_top_antitop_vector_NOSYS"
    scalar = (
        "truth_nu1_px_NOSYS", "truth_nu1_py_NOSYS", "truth_nu1_pz_NOSYS",
        "truth_nu2_px_NOSYS", "truth_nu2_py_NOSYS", "truth_nu2_pz_NOSYS",
    )
    if vector in branches:
        return vector, [vector]
    missing = [name for name in scalar if name not in branches]
    if missing:
        raise RuntimeError(
            "Canonical parent-ordered neutrino truth is unavailable. Expected either "
            f"{vector!r} or all truth_nu1/2 px,py,pz branches. Missing: {missing}. "
            "This converter deliberately does not fall back to neutrino_vector_NOSYS, "
            "whose historical ordering is flavour-priority rather than guaranteed top/antitop."
        )
    return None, list(scalar)


def _scalar_neutrino_column(array, field):
    arr = ak.to_numpy(array)
    if arr.ndim == 2 and arr.shape[1] == 1:
        arr = arr[:, 0]
    if arr.ndim != 1:
        raise ValueError(f"Expected one value/event in {field}; got {arr.shape}")
    return np.asarray(arr, dtype=np.float32)


def _neutrino_matrix(arrays, vector_branch):
    """Return (N,2,3) in source units with exact order [top, antitop]."""
    if vector_branch:
        data = arrays[vector_branch]
        padded = ak.to_numpy(
            ak.fill_none(ak.pad_none(data, 6, axis=1, clip=True), np.nan)
        ).astype(np.float32)
        counts = np.asarray(ak.num(data, axis=1))
        if np.any(counts != 6):
            raise ValueError(
                f"{vector_branch}: expected six "
                "[top_px,top_py,top_pz,antitop_px,antitop_py,antitop_pz] "
                f"components/event; observed lengths {np.unique(counts).tolist()}"
            )
        return padded.reshape(-1, 2, 3)

    top = np.stack([
        _scalar_neutrino_column(arrays[f"truth_nu1_{c}_NOSYS"], f"truth_nu1_{c}_NOSYS")
        for c in ("px", "py", "pz")
    ], axis=-1)
    antitop = np.stack([
        _scalar_neutrino_column(arrays[f"truth_nu2_{c}_NOSYS"], f"truth_nu2_{c}_NOSYS")
        for c in ("px", "py", "pz")
    ], axis=-1)
    return np.stack([top, antitop], axis=1)


def _neutrino_scale(tree, branch_names, vector_branch, entry_stop, unit):
    if unit != "auto":
        return (0.001 if unit == "MeV" else 1.0), f"forced input unit: {unit}"
    probe = tree.arrays(branch_names, entry_stop=entry_stop, library="ak")
    values = _neutrino_matrix(probe, vector_branch)
    finite = np.abs(values[np.isfinite(values) & (values != 0.0)])
    if not finite.size:
        raise ValueError("Cannot infer neutrino units from empty/non-finite probe; "
                         "pass --neutrino-unit GeV or MeV explicitly")
    median = float(np.median(finite))
    if 500.0 <= median <= 2000.0:
        raise ValueError(f"Ambiguous neutrino component magnitude (median={median:g}); "
                         "set --neutrino-unit GeV or MeV explicitly")
    mev = median > 2000.0
    return (0.001 if mev else 1.0), f"auto: median |nonzero component|={median:.4g}, input {'MeV' if mev else 'GeV'}"

def _build_assignment_labels(arrays, jet_pt, el_pt, el_charge, mu_pt, mu_charge,
                             max_jets, max_electrons, max_muons):
    n = len(jet_pt)
    valid_jets = np.isfinite(jet_pt) & (jet_pt > 0)
    valid_el = np.isfinite(el_pt) & (el_pt > 0)
    valid_mu = np.isfinite(mu_pt) & (mu_pt > 0)
    jet_labels = np.full((n, max_jets), np.nan, dtype=np.float32)
    electron_labels = np.full((n, max_electrons), np.nan, dtype=np.float32)
    muon_labels = np.full((n, max_muons), np.nan, dtype=np.float32)
    jet_labels[valid_jets] = 0.0
    electron_labels[valid_el] = np.where(el_charge[valid_el] > 0, 1.0, 2.0)
    muon_labels[valid_mu] = np.where(mu_charge[valid_mu] > 0, 1.0, 2.0)
    def index_column(branch):
        raw = np.asarray(arrays[branch], dtype=np.float64)
        if raw.shape != (n,):
            raise ValueError(f"Expected one truth index per event in {branch}; got {raw.shape}")
        # Avoid casting NaN/sentinel values to arbitrary signed integers.
        return np.where(np.isfinite(raw), raw, -1).astype(np.int64)

    b_index = index_column("Truth_had_b_index")
    bbar_index = index_column("Truth_had_bbar_index")
    true_njets = np.asarray(ak.num(arrays["Jet_pt"], axis=1), dtype=np.int64)
    good_b = ((b_index >= 0) & (b_index < true_njets) & (b_index < max_jets))
    good_bbar = ((bbar_index >= 0) & (bbar_index < true_njets) & (bbar_index < max_jets))
    rows = np.arange(n)
    good_b &= valid_jets[rows, np.clip(b_index, 0, max_jets - 1)]
    good_bbar &= valid_jets[rows, np.clip(bbar_index, 0, max_jets - 1)]
    # Distinct b and bbar indices required; neither can overwrite the other.
    distinct = b_index != bbar_index
    matched_b = good_b & good_bbar & distinct
    # Preserve partial assignments when --fully-matched is not requested.
    assign_b = good_b & distinct
    assign_bbar = good_bbar & distinct
    jet_labels[rows[assign_b], b_index[assign_b]] = 1.0
    jet_labels[rows[assign_bbar], bbar_index[assign_bbar]] = 2.0
    all_charges = np.concatenate([
        np.where(valid_el, el_charge, 0.0),
        np.where(valid_mu, mu_charge, 0.0),
    ], axis=1)
    n_leptons = np.asarray(ak.num(arrays["Electron_pt"], axis=1)) + np.asarray(
        ak.num(arrays["Muon_pt"], axis=1))
    matched_leptons = ((n_leptons == 2)
                       & (np.sum(all_charges > 0, axis=1) == 1)
                       & (np.sum(all_charges < 0, axis=1) == 1)
                       & ((valid_el.sum(axis=1) + valid_mu.sum(axis=1)) == 2))
    return jet_labels, electron_labels, muon_labels, matched_b & matched_leptons


def _create_output_h5(
output,
args,
tree_name,
reco_scale,
truth_scale,
split_name,
target_per_split,
target_source,
target_branch_map,
):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    h5f = h5py.File(output, "w")
    inputs = h5f.create_group("INPUTS")
    targets = h5f.create_group("TARGETS")
    labels = h5f.create_group("LABELS")
    h5_chunk = int(args.h5_chunk_size)
    jet_chunk = (h5_chunk, args.max_jets)
    lep_chunk = (h5_chunk, args.max_electrons)
    mu_chunk = (h5_chunk, args.max_muons)
    inputs.create_dataset(
        "JET",
shape=(0, args.max_jets),
maxshape=(None, args.max_jets),
chunks=jet_chunk,
dtype=OBJECT_DTYPE,
compression=args.compression,
    )
    inputs.create_dataset(
        "ELECTRON",
shape=(0, args.max_electrons),
maxshape=(None, args.max_electrons),
chunks=lep_chunk,
dtype=OBJECT_DTYPE,
compression=args.compression,
    )
    inputs.create_dataset(
        "MUON",
shape=(0, args.max_muons),
maxshape=(None, args.max_muons),
chunks=mu_chunk,
dtype=OBJECT_DTYPE,
compression=args.compression,
    )
    inputs.create_dataset(
        "GLOBAL",
shape=(0,),
maxshape=(None,),
chunks=(h5_chunk,),
dtype=GLOBAL_DTYPE,
compression=args.compression,
    )
    top_ds = targets.create_dataset(
        "TOP",
shape=(0,),
maxshape=(None,),
chunks=(h5_chunk,),
dtype=TARGET_DTYPE,
compression=args.compression,
    )
    truth_ds = None
    if bool(args.write_top_regression_extras):
        truth_ds = targets.create_dataset(
            "TRUTH",
shape=(0,),
maxshape=(None,),
chunks=(h5_chunk,),
dtype=TRUTH_EXTRA_DTYPE,
compression=args.compression,
        )
    for name, width in (
        ("JET", args.max_jets),
        ("ELECTRON", args.max_electrons),
        ("MUON", args.max_muons),
    ):
        ds = labels.create_dataset(
            name, shape=(0, width), maxshape=(None, width),
            chunks=(h5_chunk, width), dtype=np.float32, compression=args.compression,
        )
        ds.attrs["label_convention"] = (
            "0=background; 1=top b or positive charged lepton; "
            "2=antitop bbar or negative charged lepton; NaN=padding"
        )

    nu_ds = targets.create_dataset(
        "NEUTRINO", shape=(0, 2), maxshape=(None, 2),
        chunks=(h5_chunk, 2), dtype=NEUTRINO_DTYPE, compression=args.compression,
    )
    nu_ds.attrs["units"] = "GeV"
    nu_ds.attrs["order"] = "top, antitop"
    nu_ds.attrs["slot_0"] = "neutrino from top (positive-lepton side)"
    nu_ds.attrs["slot_1"] = "antineutrino from antitop (negative-lepton side)"
    nu_ds.attrs["components"] = "px,py,pz"
    nu_ds.attrs["source_branches"] = json.dumps(args.neutrino_branch_names)
    if args.write_meta:
        meta = h5f.create_group("META")
        meta.create_dataset(
            "source_entry_index",
shape=(0,),
maxshape=(None,),
chunks=(h5_chunk,),
dtype=np.int64,
compression=args.compression,
        )
        meta.create_dataset(
            "eventNumber",
shape=(0,),
maxshape=(None,),
chunks=(h5_chunk,),
dtype=np.int64,
compression=args.compression,
        )
    h5f.attrs["compression"] = args.compression if args.compression else "none"
    h5f.attrs["input_feature_schema"] = "minimal_wp_v2"
    inputs["GLOBAL"].attrs["description"] = (
        "Original reconstructed multiplicities plus MET magnitude/direction; "
        "five GN2v01 WP b-jet counts; no engineered event-level invariant masses"
    )
    h5f.attrs["n_saved_events"] = 0
    h5f.attrs["tree_name"] = tree_name
    h5f.attrs["split_key"] = "eventNumber"
    h5f.attrs["split_name"] = split_name if split_name is not None else "none"
    h5f.attrs["target_per_split"] = -1 if target_per_split is None else int(target_per_split)
    h5f.attrs["target_source"] = target_source
    h5f.attrs["fully_matched_applied"] = bool(args.fully_matched)
    h5f.attrs["neutrino_energy_scale_applied"] = float(args.neutrino_scale)
    h5f.attrs["neutrino_source_branches"] = json.dumps(args.neutrino_branch_names)
    h5f.attrs["neutrino_order"] = "top,antitop"
    h5f.attrs["max_jets"] = int(args.max_jets)
    h5f.attrs["reco_energy_scale_applied"] = float(reco_scale)
    h5f.attrs["truth_energy_scale_applied"] = float(truth_scale)
    top_ds.attrs["target_name"] = "ttbar2l_top_antitop"
    top_ds.attrs["target_source"] = target_source
    top_ds.attrs["source_branch_map"] = json.dumps(target_branch_map, sort_keys=True)
    top_ds.attrs["source_branches"] = json.dumps(_target_branch_list(target_branch_map))
    top_ds.attrs["physical_parameterisation"] = "pt_eta_phi_m"
    top_ds.attrs["training_parameterisation"] = "logpt_eta_sincos_logm"
    top_ds.attrs["convention"] = _target_convention(target_source)
    top_ds.attrs["units"] = "GeV"
    top_ds.attrs["tree_name"] = tree_name
    top_ds.attrs["reco_energy_scale_applied"] = float(reco_scale)
    top_ds.attrs["truth_energy_scale_applied"] = float(truth_scale)
    if truth_ds is not None:
        truth_ds.attrs["target_name"] = "ttbar2l_top_regression_truth_extras"
        truth_ds.attrs["source_target_dataset"] = "/TARGETS/TOP"
        truth_ds.attrs["units"] = "GeV"
        truth_ds.attrs["mttbar_bin_edges"] = json.dumps(TRUTH_MTTBAR_BIN_EDGES.tolist())
        truth_ds.attrs["angular_convention"] = (
            "truth_cos_theta_star and truth_phi_star use the top direction after boosting "
            "to the ttbar centre-of-mass frame. The polar axis is the ttbar lab-frame "
            "flight direction; if degenerate, lab +z is used. phi_star is measured from "
            "lab +z projected perpendicular to that axis; if degenerate, lab +x is used."
        )
    return h5f
def _append_to_h5(
h5f,
jet_data,
el_data,
mu_data,
global_data,
target_data,
truth_extra_data=None,
source_entry_indices=None,
    event_numbers=None,
    jet_labels=None,
    electron_labels=None,
    muon_labels=None,
    neutrino_labels=None,
):
    _append_dataset(h5f["INPUTS/JET"], jet_data)
    _append_dataset(h5f["INPUTS/ELECTRON"], el_data)
    _append_dataset(h5f["INPUTS/MUON"], mu_data)
    _append_dataset(h5f["INPUTS/GLOBAL"], global_data)
    _append_dataset(h5f["TARGETS/TOP"], target_data)
    if truth_extra_data is not None and "TARGETS/TRUTH" in h5f:
        _append_dataset(h5f["TARGETS/TRUTH"], truth_extra_data)
    _append_dataset(h5f["LABELS/JET"], jet_labels)
    _append_dataset(h5f["LABELS/ELECTRON"], electron_labels)
    _append_dataset(h5f["LABELS/MUON"], muon_labels)
    _append_dataset(h5f["TARGETS/NEUTRINO"], neutrino_labels)
    if "META" in h5f:
        _append_dataset(h5f["META/source_entry_index"], source_entry_indices)
        _append_dataset(h5f["META/eventNumber"], event_numbers)
    h5f.attrs["n_saved_events"] = int(h5f["TARGETS/TOP"].shape[0])
def _build_chunk_data(
arrays,
args,
reco_scale,
truth_scale,
btag_branch,
chunk_entry_offset,
target_branch_map,
):
    first_target_branch = target_branch_map[TARGET_FIELDS[0]]
    n_input = len(arrays[first_target_branch])
    event_numbers = np.asarray(arrays["eventNumber"], dtype=np.int64)
    target_columns = [
        np.asarray(arrays[target_branch_map[field]], dtype=np.float32)
        for field in TARGET_FIELDS
    ]
    target_matrix_raw = np.stack(target_columns, axis=1).astype(np.float32)
    target_matrix = target_matrix_raw.copy()
    # Scale pt and mass fields to GeV if needed:
    # top_pt, top_m, antitop_pt, antitop_m
    target_matrix[:, [0, 3, 4, 7]] *= truth_scale
    finite = np.isfinite(target_matrix).all(axis=1)
    positive = (target_matrix[:, [0, 3, 4, 7]] > 0.0).all(axis=1)
    no_sentinel = ~np.isin(target_matrix, [-1.0, -99.0]).any(axis=1)
    keep = finite & positive & no_sentinel
    jet_e = _pad_jagged(arrays["Jet_e"] * reco_scale, args.max_jets)
    jet_eta = _pad_jagged(arrays["Jet_eta"], args.max_jets)
    jet_phi = _pad_jagged(arrays["Jet_phi"], args.max_jets)
    jet_pt = _pad_jagged(arrays["Jet_pt"] * reco_scale, args.max_jets)
    if btag_branch:
        jet_btag = _pad_jagged(arrays[btag_branch], args.max_jets)
    else:
        jet_btag = np.zeros_like(jet_pt, dtype=np.float32)
        jet_btag[~np.isfinite(jet_pt)] = np.nan
    jet_charge = np.zeros_like(jet_pt, dtype=np.float32)
    jet_charge[~np.isfinite(jet_pt)] = np.nan
    el_e = _pad_jagged(arrays["Electron_e"] * reco_scale, args.max_electrons)
    el_eta = _pad_jagged(arrays["Electron_eta"], args.max_electrons)
    el_phi = _pad_jagged(arrays["Electron_phi"], args.max_electrons)
    el_pt = _pad_jagged(arrays["Electron_pt"] * reco_scale, args.max_electrons)
    el_charge = _pad_jagged(arrays["Electron_q"], args.max_electrons)
    el_btag = np.zeros_like(el_pt, dtype=np.float32)
    el_btag[~np.isfinite(el_pt)] = np.nan
    mu_e = _pad_jagged(arrays["Muon_e"] * reco_scale, args.max_muons)
    mu_eta = _pad_jagged(arrays["Muon_eta"], args.max_muons)
    mu_phi = _pad_jagged(arrays["Muon_phi"], args.max_muons)
    mu_pt = _pad_jagged(arrays["Muon_pt"] * reco_scale, args.max_muons)
    mu_charge = _pad_jagged(arrays["Muon_q"], args.max_muons)
    mu_btag = np.zeros_like(mu_pt, dtype=np.float32)
    mu_btag[~np.isfinite(mu_pt)] = np.nan
    jet_labels, electron_labels, muon_labels, matched = _build_assignment_labels(
        arrays, jet_pt, el_pt, el_charge, mu_pt, mu_charge,
        args.max_jets, args.max_electrons, args.max_muons,
    )
    if args.fully_matched:
        keep &= matched
    nu_matrix = _neutrino_matrix(arrays, args.neutrino_branch) * args.neutrino_scale
    if args.require_neutrino_truth:
        keep &= np.isfinite(nu_matrix).all(axis=(1, 2))
    nu_data = np.empty((int(keep.sum()), 2), dtype=NEUTRINO_DTYPE)
    for j, component in enumerate(("px", "py", "pz")):
        nu_data[component] = nu_matrix[keep, :, j]
    jet_data = _make_object_array(
        jet_e,
        jet_eta,
        jet_phi,
        jet_pt,
        jet_btag,
        jet_charge,
        1,
        keep,
    )
    el_data = _make_object_array(
        el_e,
        el_eta,
        el_phi,
        el_pt,
        el_btag,
        el_charge,
        2,
        keep,
    )
    mu_data = _make_object_array(
        mu_e,
        mu_eta,
        mu_phi,
        mu_pt,
        mu_btag,
        mu_charge,
        3,
        keep,
    )
    met = np.asarray(arrays["Met_pt"], dtype=np.float32) * reco_scale
    met_phi = np.asarray(arrays["Met_phi"], dtype=np.float32)
    global_data = np.zeros(int(keep.sum()), dtype=GLOBAL_DTYPE)
    global_data["njet"] = ak.to_numpy(
        ak.num(arrays["Jet_pt"], axis=1)
    ).astype(np.float32)[keep]
    global_data["nelectron"] = ak.to_numpy(
        ak.num(arrays["Electron_pt"], axis=1)
    ).astype(np.float32)[keep]
    global_data["nmuon"] = ak.to_numpy(
        ak.num(arrays["Muon_pt"], axis=1)
    ).astype(np.float32)[keep]
    global_data["met_met"] = met[keep]
    global_data["met_phi"] = met_phi[keep]
    global_data["met_cosphi"] = np.cos(met_phi[keep])
    global_data["met_sinphi"] = np.sin(met_phi[keep])
    for field, branch in B_COUNT_BRANCHES.items():
        global_data[field] = np.asarray(arrays[branch], dtype=np.float32)[keep]
    target_data = np.empty(int(keep.sum()), dtype=TARGET_DTYPE)
    for i, field in enumerate(TARGET_FIELDS):
        target_data[field] = target_matrix[keep, i]
    truth_extra_data = _build_truth_extras(target_matrix[keep])
    local_entry_indices = np.arange(n_input, dtype=np.int64)
    source_entry_indices = chunk_entry_offset + local_entry_indices[keep]
    accepted_event_numbers = event_numbers[keep]
    return {
        "n_read": int(n_input),
        "n_accepted": int(keep.sum()),
        "n_rejected": int((~keep).sum()),
        "event_numbers": accepted_event_numbers,
        "source_entry_indices": source_entry_indices,
        "jet": jet_data,
        "electron": el_data,
        "muon": mu_data,
        "global": global_data,
        "target": target_data,
        "jet_labels": jet_labels[keep],
        "electron_labels": electron_labels[keep],
        "muon_labels": muon_labels[keep],
        "neutrino_labels": nu_data,
        "truth_extras": truth_extra_data,
        "target_matrix_all": target_matrix,
    }
def _resolve_outputs(args):
    output = Path(args.output) if args.output else Path("ttbar_dilepton_top_regression.h5")
    if args.parity_split:
        stem = output.with_suffix("")
        even_output = Path(args.output_even) if args.output_even else Path(f"{stem}_even.h5")
        odd_output = Path(args.output_odd) if args.output_odd else Path(f"{stem}_odd.h5")
        return output, [even_output, odd_output]
    return output, [output]
def _check_outputs(args):
    _, outputs_to_check = _resolve_outputs(args)
    for out in outputs_to_check:
        if out.exists() and not args.overwrite:
            raise FileExistsError(f"{out} exists. Use --overwrite to replace it.")
def _needed_branches(branches, btag_branch, target_branch_map, neutrino_branch_names):
    required = [
        "Jet_pt",
        "Jet_eta",
        "Jet_phi",
        "Jet_e",
        "Electron_pt",
        "Electron_eta",
        "Electron_phi",
        "Electron_e",
        "Electron_q",
        "Muon_pt",
        "Muon_eta",
        "Muon_phi",
        "Muon_e",
        "Muon_q",
        "Met_pt",
        "Met_phi",
        "eventNumber",
        "Truth_had_b_index", "Truth_had_bbar_index",
    ] + _target_branch_list(target_branch_map) + neutrino_branch_names
    required += list(B_COUNT_BRANCHES.values())
    optional = [btag_branch] if btag_branch else []
    missing = [name for name in required if name not in branches]
    if missing:
        raise RuntimeError(f"Missing required input branches: {missing}")
    # Use dict order to remove duplicates while preserving order.
    return list(dict.fromkeys(required + optional))
def convert(args):
    if args.dataset_size is not None and args.dataset_size <= 0:
        raise ValueError("--dataset-size must be a positive integer")
    if min(args.max_jets, args.max_electrons, args.max_muons, args.h5_chunk_size, args.chunk_size) < 1:
        raise ValueError("Object capacities and chunk sizes must be positive")
    if args.compression == "none":
        args.compression = None
    _check_outputs(args)
    root_file = uproot.open(args.input)
    tree_name, tree = _select_tree(root_file, args.tree_name)
    branches = set(tree.keys())
    target_branch_map = _resolve_target_branches(tree, args.target_source)
    print(f"Using tree: {tree_name}")
    print(f"Input entries: {tree.num_entries}")
    print(f"Target source: {args.target_source}")
    print("Resolved target branches:")
    for field in TARGET_FIELDS:
        print(f"  {field}: {target_branch_map[field]}")
    scale_probe_stop = min(int(args.scale_probe_events), tree.num_entries)
    reco_scale, reco_scale_note = _energy_scale_from_branch(tree, "Jet_pt", scale_probe_stop)
    truth_scale, truth_scale_note = _energy_scale_from_branch(
        tree,
        target_branch_map["top_m"],
        scale_probe_stop,
    )
    print(f"Reco scale: {reco_scale_note} (scale={reco_scale})")
    print(f"Truth scale: {truth_scale_note} (scale={truth_scale})")
    btag_branch = _first_existing(
        branches,
        [
            "Jet_GN2v01_PCBT_quantile_NOSYS",
            "Jet_GN2v01_score",
        ],
    )
    if btag_branch is None:
        print("WARNING: no Jet GN2 b-tag branch found; writing jet btag=0.")
    else:
        print(f"Using jet b-tag branch: {btag_branch}")
    args.neutrino_branch, neutrino_branch_names = _resolve_neutrino_branches(branches)
    args.neutrino_branch_names = neutrino_branch_names
    args.neutrino_scale, neutrino_note = _neutrino_scale(
        tree, neutrino_branch_names, args.neutrino_branch, scale_probe_stop,
        args.neutrino_unit,
    )
    print(f"Neutrino source: {args.neutrino_branch or neutrino_branch_names}")
    print(f"Neutrino scale: {neutrino_note} (scale={args.neutrino_scale})")
    needed_branches = _needed_branches(
        branches, btag_branch, target_branch_map, neutrino_branch_names,
    )
    output, outputs = _resolve_outputs(args)
    target_per_split = args.dataset_size
    if args.parity_split:
        even_output, odd_output = outputs
        print("Split key: eventNumber % 2")
        print(f"Even output: {even_output}")
        print(f"Odd output:  {odd_output}")
        if target_per_split is None:
            print("Target per split: all available accepted events")
        else:
            print(f"Target per split: {target_per_split}")
        even_h5 = _create_output_h5(
            even_output,
            args,
            tree_name,
            reco_scale,
            truth_scale,
split_name="even",
target_per_split=target_per_split,
target_source=args.target_source,
target_branch_map=target_branch_map,
        )
        odd_h5 = _create_output_h5(
            odd_output,
            args,
            tree_name,
            reco_scale,
            truth_scale,
split_name="odd",
target_per_split=target_per_split,
target_source=args.target_source,
target_branch_map=target_branch_map,
        )
        h5_files = [even_h5, odd_h5]
    else:
        print(f"Output: {output}")
        if target_per_split is None:
            print("Target events: all available accepted events")
        else:
            print(f"Target events: {target_per_split}")
        single_h5 = _create_output_h5(
            output,
            args,
            tree_name,
            reco_scale,
            truth_scale,
split_name=None,
target_per_split=target_per_split,
target_source=args.target_source,
target_branch_map=target_branch_map,
        )
        h5_files = [single_h5]
    total_read = 0
    total_accepted = 0
    total_rejected = 0
    chunk_index = 0
    target_stats_seen = False
    first_target_matrix_for_stats = None
    try:
        for arrays in tree.iterate(
            needed_branches,
library="ak",
step_size=int(args.chunk_size),
        ):
            if args.max_root_entries is not None and total_read >= args.max_root_entries:
                break
            if args.max_root_entries is not None:
                remaining_root = int(args.max_root_entries) - total_read
                if remaining_root <= 0:
                    break
                first_target_branch = target_branch_map[TARGET_FIELDS[0]]
                if len(arrays[first_target_branch]) > remaining_root:
                    arrays = arrays[:remaining_root]
            chunk = _build_chunk_data(
arrays=arrays,
args=args,
reco_scale=reco_scale,
truth_scale=truth_scale,
btag_branch=btag_branch,
chunk_entry_offset=total_read,
target_branch_map=target_branch_map,
            )
            chunk_index += 1
            total_read += chunk["n_read"]
            total_accepted += chunk["n_accepted"]
            total_rejected += chunk["n_rejected"]
            if not target_stats_seen and chunk["target_matrix_all"].size:
                first_target_matrix_for_stats = chunk["target_matrix_all"]
                target_stats_seen = True
            if chunk["n_accepted"] == 0:
                print(
f"Chunk {chunk_index}: read={chunk['n_read']} "
f"accepted=0 rejected={chunk['n_rejected']} "
f"saved_even={even_h5['TARGETS/TOP'].shape[0] if args.parity_split else '-'} "
f"saved_odd={odd_h5['TARGETS/TOP'].shape[0] if args.parity_split else '-'}"
                )
                continue
            if args.parity_split:
                event_numbers = chunk["event_numbers"]
                even_mask = (event_numbers % 2) == 0
                odd_mask = (event_numbers % 2) == 1
                even_saved = int(even_h5["TARGETS/TOP"].shape[0])
                odd_saved = int(odd_h5["TARGETS/TOP"].shape[0])
                if target_per_split is None:
                    even_take_idx = np.where(even_mask)[0]
                    odd_take_idx = np.where(odd_mask)[0]
                else:
                    even_remaining = max(0, int(target_per_split) - even_saved)
                    odd_remaining = max(0, int(target_per_split) - odd_saved)
                    even_take_idx = np.where(even_mask)[0][:even_remaining]
                    odd_take_idx = np.where(odd_mask)[0][:odd_remaining]
                if len(even_take_idx):
                    _append_to_h5(
                        even_h5,
                        chunk["jet"][even_take_idx],
                        chunk["electron"][even_take_idx],
                        chunk["muon"][even_take_idx],
                        chunk["global"][even_take_idx],
                        chunk["target"][even_take_idx],
                        chunk["truth_extras"][even_take_idx] if args.write_top_regression_extras else None,
                        chunk["source_entry_indices"][even_take_idx],
                        chunk["event_numbers"][even_take_idx],
                        jet_labels=chunk["jet_labels"][even_take_idx],
                        electron_labels=chunk["electron_labels"][even_take_idx],
                        muon_labels=chunk["muon_labels"][even_take_idx],
                        neutrino_labels=chunk["neutrino_labels"][even_take_idx],
                    )
                if len(odd_take_idx):
                    _append_to_h5(
                        odd_h5,
                        chunk["jet"][odd_take_idx],
                        chunk["electron"][odd_take_idx],
                        chunk["muon"][odd_take_idx],
                        chunk["global"][odd_take_idx],
                        chunk["target"][odd_take_idx],
                        chunk["truth_extras"][odd_take_idx] if args.write_top_regression_extras else None,
                        chunk["source_entry_indices"][odd_take_idx],
                        chunk["event_numbers"][odd_take_idx],
                        jet_labels=chunk["jet_labels"][odd_take_idx],
                        electron_labels=chunk["electron_labels"][odd_take_idx],
                        muon_labels=chunk["muon_labels"][odd_take_idx],
                        neutrino_labels=chunk["neutrino_labels"][odd_take_idx],
                    )
                even_saved = int(even_h5["TARGETS/TOP"].shape[0])
                odd_saved = int(odd_h5["TARGETS/TOP"].shape[0])
                print(
f"Chunk {chunk_index}: "
f"read_total={total_read} "
f"accepted_total={total_accepted} "
f"rejected_total={total_rejected} "
f"saved_even={even_saved} "
f"saved_odd={odd_saved}"
                )
                if (
                    target_per_split is not None
                    and even_saved >= int(target_per_split)
                    and odd_saved >= int(target_per_split)
                ):
                    print("Reached target per split.")
                    break
            else:
                saved = int(single_h5["TARGETS/TOP"].shape[0])
                if target_per_split is None:
                    take_idx = np.arange(chunk["n_accepted"])
                else:
                    remaining = max(0, int(target_per_split) - saved)
                    take_idx = np.arange(chunk["n_accepted"])[:remaining]
                if len(take_idx):
                    _append_to_h5(
                        single_h5,
                        chunk["jet"][take_idx],
                        chunk["electron"][take_idx],
                        chunk["muon"][take_idx],
                        chunk["global"][take_idx],
                        chunk["target"][take_idx],
                        chunk["truth_extras"][take_idx] if args.write_top_regression_extras else None,
                        chunk["source_entry_indices"][take_idx],
                        chunk["event_numbers"][take_idx],
                        jet_labels=chunk["jet_labels"][take_idx],
                        electron_labels=chunk["electron_labels"][take_idx],
                        muon_labels=chunk["muon_labels"][take_idx],
                        neutrino_labels=chunk["neutrino_labels"][take_idx],
                    )
                saved = int(single_h5["TARGETS/TOP"].shape[0])
                print(
f"Chunk {chunk_index}: "
f"read_total={total_read} "
f"accepted_total={total_accepted} "
f"rejected_total={total_rejected} "
f"saved={saved}"
                )
                if target_per_split is not None and saved >= int(target_per_split):
                    print("Reached target.")
                    break
        print(f"Read events: {total_read}")
        print(f"Accepted after target quality: {total_accepted}")
        print(f"Rejected target count: {total_rejected}")
        if first_target_matrix_for_stats is not None:
            print("Target stats from first processed chunk:")
            for field, values in zip(TARGET_FIELDS, first_target_matrix_for_stats.T):
                print(_stats(field, values))
        if args.parity_split:
            even_saved = int(even_h5["TARGETS/TOP"].shape[0])
            odd_saved = int(odd_h5["TARGETS/TOP"].shape[0])
            print(f"Wrote H5: {outputs[0]}")
            print(f"  saved events: {even_saved}")
            if args.write_meta and even_saved:
                evts = even_h5["META/eventNumber"][:]
                print(f"  eventNumber range: {int(np.min(evts))} -> {int(np.max(evts))}")
                print(f"  all even: {bool(np.all(evts % 2 == 0))}")
            print(f"Wrote H5: {outputs[1]}")
            print(f"  saved events: {odd_saved}")
            if args.write_meta and odd_saved:
                evts = odd_h5["META/eventNumber"][:]
                print(f"  eventNumber range: {int(np.min(evts))} -> {int(np.max(evts))}")
                print(f"  all odd: {bool(np.all(evts % 2 == 1))}")
        else:
            saved = int(single_h5["TARGETS/TOP"].shape[0])
            print(f"Wrote H5: {output}")
            print(f"  saved events: {saved}")
        print(
            "Output groups: "
            "/INPUTS/JET, /INPUTS/ELECTRON, /INPUTS/MUON, "
            "/INPUTS/GLOBAL, /TARGETS/TOP, /LABELS/JET, "
            "/LABELS/ELECTRON, /LABELS/MUON, /TARGETS/NEUTRINO"
        )
        if args.write_top_regression_extras:
            print("Truth extras dataset: /TARGETS/TRUTH")
        if args.write_meta:
            print("Metadata group: /META/source_entry_index, /META/eventNumber")
        print("Target dtype/order:", TARGET_DTYPE)
    finally:
        for h5f in h5_files:
            h5f.close()
        root_file.close()
def parse_args():
    parser = argparse.ArgumentParser(
description="Convert ttbar dilepton ntuple to joint regression + assignment + parent-ordered neutrino truth H5."
    )
    parser.add_argument(
        "-f",
        "--file",
        "--input",
dest="input",
required=True,
help="Input ROOT file.",
    )
    parser.add_argument(
        "-o",
        "--output",
default="ttbar_dilepton_combined_afterFSR.h5",
help=(
            "Output H5 file. With parity splitting, this is used as the basename. "
            "Default: ttbar_dilepton_combined_afterFSR.h5"
        ),
    )
    parser.add_argument(
        "--output-even",
default=None,
help="Explicit even-eventNumber output H5 path.",
    )
    parser.add_argument(
        "--output-odd",
default=None,
help="Explicit odd-eventNumber output H5 path.",
    )
    parser.add_argument(
        "--tree-name",
default=None,
help="TTree name. Default: auto-detect, preferring reco if present.",
    )
    parser.add_argument(
        "--target-source",
choices=sorted(TARGET_BRANCHES_BY_SOURCE),
default="partonic_after_fsr",
help=(
            "Target source for /TARGETS/TOP. "
            "reco_truth uses top_reco_truth_* and antitop_reco_truth_* branches. "
            "partonic_after_fsr uses t_afterFSR_* and tbar_afterFSR_* branches. "
            "Default: partonic_after_fsr."
        ),
    )
    parser.add_argument(
        "-ds",
        "--dataset-size",
dest="dataset_size",
type=int,
default=None,
help=(
            "Target number of accepted saved events. "
            "With parity splitting, this is the target per split, "
            "so -ds 2000000 writes up to 2M even and 2M odd."
        ),
    )
    parser.add_argument(
        "--max-root-entries",
type=int,
default=None,
help=(
            "Optional safety cap on total ROOT entries read. "
            "Usually leave unset when using -ds with parity splitting."
        ),
    )
    parser.add_argument(
        "--chunk-size",
type=int,
default=100000,
help="ROOT entries to process per chunk.",
    )
    parser.add_argument(
        "--compression", choices=("lzf", "gzip", "none"), default="lzf",
        help="HDF5 compression codec (default lzf for faster writing).",
    )
    parser.add_argument(
        "--h5-chunk-size",
type=int,
default=4096,
help="HDF5 chunk size along event axis.",
    )
    parser.add_argument(
        "--scale-probe-events",
type=int,
default=100000,
help="Number of initial entries to inspect for GeV/MeV scale detection.",
    )
    parser.add_argument(
        "--overwrite",
action="store_true",
help="Overwrite existing output H5 files.",
    )
    parser.add_argument("--max-jets", type=int, default=20)
    parser.add_argument("--max-electrons", type=int, default=2)
    parser.add_argument("--max-muons", type=int, default=2)
    parser.add_argument(
        "--neutrino-unit", choices=("auto", "GeV", "MeV"), default="auto",
        help="Units of neutrino input branches; output is GeV (default: auto).",
    )
    parser.add_argument(
        "--require-neutrino-truth", action="store_true",
        help="Reject events with non-finite neutrino truth components.",
    )
    parser.add_argument(
        "--parity-split",
action="store_true",
default=True,
help=(
            "Write separate even/odd H5 files using eventNumber parity. "
            "Default: enabled."
        ),
    )
    parser.add_argument(
        "--no-parity-split",
dest="parity_split",
action="store_false",
help="Disable even/odd splitting and write one H5 file.",
    )
    parser.add_argument(
        "--write-meta",
action="store_true",
default=True,
help="Write /META/eventNumber and /META/source_entry_index. Default: enabled.",
    )
    parser.add_argument(
        "--no-write-meta",
dest="write_meta",
action="store_false",
help="Do not write /META datasets.",
    )
    parser.add_argument(
        "--write-top-regression-extras",
nargs="?",
const=True,
default=True,
type=_str_to_bool,
help=(
            "Write /TARGETS/TRUTH with derived truth four-vector, ttbar-system, "
            "mttbar-bin, and helicity-angle quantities. Default: true."
        ),
    )
    # Compatibility flags matching your existing converter command style.
    parser.add_argument(
        "-s",
        "--save",
action="store_true",
help="Compatibility flag; output is always saved.",
    )
    parser.add_argument(
        "-b",
        "--btag",
action="store_true",
help=(
            "Compatibility flag. This converter always writes btag fields. "
            "If no btag branch exists, jet btag is filled with 0."
        ),
    )
    parser.add_argument(
        "-fm",
        "--fully-matched",
action="store_true",
default=True,
        help=(
            "Require two distinct in-range matched b/bbar jets and exactly "
            "two opposite-sign leptons (actual VyPER fully-matched selection)."
        ),
    )
    parser.add_argument("--no-fully-matched", dest="fully_matched", action="store_false")
    return parser.parse_args()
if __name__ == "__main__":
    convert(parse_args())
