"""Cheap gates for the truth-lepton sidecar and the identity join (numpy only).

Run with ``python3 -m pytest Validation/tests`` or ``python3 Validation/tests/test_truth_lepton_sidecar.py``.
"""

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import audit_truth_lepton_join as join_mod  # noqa: E402
import build_truth_lepton_sidecar as sidecar  # noqa: E402

E, MU, TAU, NUE, NUM, NUT = 11, 13, 15, 12, 14, 16
# rows: (t daughter1, t daughter2, tbar daughter1, tbar daughter2)
CASES = {
    "emu_ok":            (-E, NUE, MU, -NUM),     # e+ nu_e | mu- nubar_mu
    "swapped_order_ok":  (NUM, -MU, -NUE, E),     # nu first on both sides
    "tau_top":           (-TAU, NUT, E, -NUE),
    "hadronic_tbar":     (-E, NUE, 1, -2),
    "wrong_sign_top":    (E, NUE, MU, -NUM),      # e- from top
    "sibling_mismatch":  (-E, NUM, MU, -NUM),     # e+ with nu_mu
}


def _arrays(unit_scale=1000.0):
    names = list(CASES)
    n = len(names)
    arrays = {"eventNumber": np.arange(100, 100 + n), "mcChannelNumber": np.full(n, 410472)}
    rng = np.random.default_rng(1)
    for side_index, side in enumerate(("t", "tbar")):
        for d in (1, 2):
            pdg = np.array([CASES[k][2 * side_index + d - 1] for k in names])
            arrays[f"Ttbar_MC_Wdecay{d}_afterFSR_from_{side}_pdgId"] = pdg
            arrays[f"Ttbar_MC_Wdecay{d}_afterFSR_from_{side}_pt"] = rng.uniform(20, 120, n) * unit_scale
            arrays[f"Ttbar_MC_Wdecay{d}_afterFSR_from_{side}_eta"] = rng.uniform(-2, 2, n)
            arrays[f"Ttbar_MC_Wdecay{d}_afterFSR_from_{side}_phi"] = rng.uniform(-3, 3, n)
            mass = np.where(np.abs(pdg) == E, 0.000511, np.where(np.abs(pdg) == MU, 0.10566, 0.0))
            arrays[f"Ttbar_MC_Wdecay{d}_afterFSR_from_{side}_m"] = mass * unit_scale
    return names, arrays


def test_parent_ordering_signs_and_flags():
    names, arrays = _arrays()
    scale, _ = sidecar.detect_scale(arrays["Ttbar_MC_Wdecay1_afterFSR_from_t_pt"], "auto")
    assert scale == 1.0e-3
    out, counters = sidecar.extract(arrays, scale)
    row = {k: i for i, k in enumerate(names)}
    support = out["truth_spin_support"]
    assert support[row["emu_ok"]] and support[row["swapped_order_ok"]]
    for bad in ("tau_top", "hadronic_tbar", "wrong_sign_top", "sibling_mismatch"):
        assert not support[row[bad]], bad
    assert out["truth_lplus_pdgId"][row["emu_ok"]] == -E and out["truth_lminus_pdgId"][row["emu_ok"]] == MU
    assert out["truth_lplus_pdgId"][row["swapped_order_ok"]] == -MU and out["truth_lminus_pdgId"][row["swapped_order_ok"]] == E
    assert out["truth_tau_top"][row["tau_top"]] and out["truth_direct_emu_tbar"][row["tau_top"]]
    assert np.isnan(out["truth_lplus_E"][row["tau_top"]])
    assert np.isfinite(out["truth_lminus_E"][row["tau_top"]])
    assert np.isfinite(out["truth_lplus_E"][row["hadronic_tbar"]])
    assert np.isnan(out["truth_lminus_E"][row["hadronic_tbar"]])
    assert counters["n_wrong_sign"] == 1 and counters["n_sibling_neutrino_flavour_mismatch"] == 1
    assert counters["n_tau_top"] == 1 and counters["n_hadronic_tbar"] == 1
    assert counters["n_truth_spin_support"] == 2


def test_four_vector_units_and_choice_of_daughter():
    names, arrays = _arrays()
    out, _ = sidecar.extract(arrays, 1.0e-3)
    i = names.index("swapped_order_ok")  # the charged lepton is daughter 2 on the top side
    pt = arrays["Ttbar_MC_Wdecay2_afterFSR_from_t_pt"][i] * 1e-3
    eta = arrays["Ttbar_MC_Wdecay2_afterFSR_from_t_eta"][i]
    phi = arrays["Ttbar_MC_Wdecay2_afterFSR_from_t_phi"][i]
    expected = np.array([pt * np.cosh(eta), pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta)])
    got = np.array([out[f"truth_lplus_{c}"][i] for c in ("E", "px", "py", "pz")])
    assert np.allclose(got, expected, rtol=1e-5)
    assert np.isnan(out["truth_lplus_E"][names.index("tau_top")])


def test_campaign_from_run_number_and_path():
    assert list(sidecar.campaign_from("x", np.array([410000, 410000]))) == ["mc23a", "mc23a"]
    assert list(sidecar.campaign_from("/a/r13144/f.root", np.array([1]))) == ["mc20d"]


def _write_sidecar(path, eventnumbers, dsid, run, support):
    n = len(eventnumbers)
    arrays = {"mcChannelNumber": np.full(n, dsid), "runNumber": np.full(n, run), "eventNumber": np.asarray(eventnumbers),
              "campaign": np.full(n, "mc23a"), "source_file_index": np.zeros(n, np.int32), "source_entry": np.arange(n)}
    for name in join_mod.TRUTH_COLUMNS:
        arrays[name] = np.asarray(support) if name == "truth_spin_support" else np.ones(n)
    np.savez(path, **arrays)
    Path(path).with_suffix(".manifest.json").write_text(json.dumps({"output_sha256": "x"}))


def test_join_is_by_key_not_position_and_detects_ambiguity():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        keys = tmp / "keys.npz"
        np.savez(keys, event_id=np.array([7, 3, 5]), source_entry_index=np.array([0, 1, 2]),
                 mcChannelNumber=np.array([601230] * 3), runNumber=np.array([410000] * 3))
        side = tmp / "side.npz"
        _write_sidecar(side, [5, 9, 7, 3], 601230, 410000, [True, False, False, True])  # shuffled order
        report = join_mod.join(keys, [side], tmp / "aligned.npz", None)
        assert report["passed"] and report["matched"] == 3
        with np.load(tmp / "aligned.npz") as aligned:
            assert aligned["event_id"].tolist() == [7, 3, 5]
            assert aligned["truth_spin_support"].tolist() == [False, True, True]
        dup = tmp / "dup.npz"
        _write_sidecar(dup, [5, 5, 7, 3], 601230, 410000, [True, True, True, True])
        report = join_mod.join(keys, [dup], tmp / "aligned_dup.npz", None)
        assert not report["passed"] and report["ambiguous"] == 1
        missing = tmp / "missing.npz"
        _write_sidecar(missing, [5, 7], 601230, 410000, [True, True])
        report = join_mod.join(keys, [missing], tmp / "aligned_missing.npz", None)
        assert not report["passed"] and report["unmatched"] == 1


def test_validate_flags_wrong_sign_counter():
    names, arrays = _arrays()
    out, counters = sidecar.extract(arrays, 1.0e-3)
    n = len(names)
    ids = {"mcChannelNumber": np.full(n, 410472), "runNumber": np.full(n, 284500), "eventNumber": arrays["eventNumber"],
           "campaign": np.full(n, "mc20a"), "source_file_index": np.zeros(n, np.int32), "source_entry": np.arange(n)}
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "s.npz"
        np.savez(path, **ids, **out)
        path.with_suffix(".manifest.json").write_text(json.dumps({"rows": n, "counters": counters}))
        checks = sidecar.validate(path)
        assert checks["finite_on_support"] and checks["sibling_flavour_on_support"]
        assert checks["nan_off_support_lplus"] and checks["nan_off_support_lminus"]
        assert checks["max_abs_lepton_mass_gev"] < 0.2
        assert not checks["passed"]  # the synthetic sample deliberately contains one wrong-sign top


if __name__ == "__main__":
    for name, function in list(globals().items()):
        if name.startswith("test_"):
            function()
            print("PASS", name)
