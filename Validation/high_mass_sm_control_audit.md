# High-mass SM control audit

## MC20a nominal SM, DSID 410472

The exact local `r13167` directory contains 766 ROOT files. The raw `reco` tree scan found:

| Quantity | Count |
|---|---:|
| Raw tree rows | 7,246,983 |
| Finite positive after-FSR truth (m_{t\bar{t}}) | 7,246,982 |
| Direct after-FSR (W\to e/\mu) on both sides | 6,535,112 |
| Direct dilepton, (m_{t\bar{t}}>1.0\) TeV | 162,645 |
| Direct dilepton, (m_{t\bar{t}}>1.6\) TeV | 10,988 |
| Direct dilepton, (m_{t\bar{t}}>2.0\) TeV | 2,292 |
| Direct dilepton, (m_{t\bar{t}}>3.0\) TeV | 72 |
| Direct dilepton, (m_{t\bar{t}}>4.0\) TeV | 3 |
| Direct dilepton, (m_{t\bar{t}}>5.0\) TeV | 0 |

The largest direct-dilepton truth mass is 4.227 TeV. Direct-dilepton event counts in truth-mass bins [0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0] TeV are [3,592,632; 2,779,835; 146,047; 14,306; 2,220; 69; 3]. These are **raw truth counts only**. They are not selected, truth-reconstructible, or candidate-supported counts and must not be used as the final performance denominator.

This nominal MC20a sample can extend common response control beyond the current partial-slice limit near 1.6 TeV, likely into the 2–2.5 TeV range if selection and candidate support remain adequate. It does not provide useful nominal SM support through 5–6 TeV.

The high-mass cut and branch adapter were smoke-tested on the non-overlapping file `user.leet.52697341._000028.output.root`: 74 direct dilepton events above 1.6 TeV, 58 selected, 45 truth-reconstructible and candidate-supported. EM and NW pass the same W-mass and unit checks (median W masses 80.3999 and 80.4000 GeV), and assignment parent ordering validates. The event-level report is `Validation/high_mass_sm_control_pilot.json`. This is a one-file correctness check, not the full matched-control result.

Do not combine the local `r13144` and `r13145` DSID 410472 directories with `r13167` without an event-identity overlap audit. Exclude the already used SM presentation file from any high-mass add-on to avoid double counting.

## Catalogued alternate ttbar samples and productions

The local `metadata_kkgraviton/filelist.txt` classifies the method-bearing `r14861` DSID 601491 production as **MC20e fast simulation**, not MC23. The same DSID is a Powheg+Pythia8 ttbar pTHard variation family according to the [ATLAS generator classification](https://atlas-sw-doxygen.web.cern.ch/atlas-sw-doxygen/atlas_main--Doxygen/docs/html/d9/ddd/HadronOriginClassifier_8cxx_source.html). The `r14861` directory has 1,762 ROOT files and contains VyPER, EM, NW, and DA branches. Its raw support audit reports 24,283 direct-dilepton events above 1.6 TeV, 5,204 above 2 TeV, 149 above 3 TeV, 6 above 4 TeV, and none above 5 TeV; the maximum direct-dilepton truth mass is 4.534 TeV. This is an alternate fast-simulation pTHard sample, not the nominal full-simulation SM control. The one-file pilot has nontrivial `weight_mc_NOSYS` and `weight_pileup_NOSYS` distributions and does not establish a weighted high-mass response.

The same local catalogs list distinct MC23 full-simulation DSID 601491 TtNT productions: MC23a `r14622` (97 paths), MC23d `r15224` (4), and MC23e `r16083` (160) in the `metadata_ttbar2L_16_2` catalog. These are not the audited `r14861` sample. Their branch schemas, weights, and event-level high-mass support have not been validated. A separate `mc23ad_ttbar_dilepton_reco` directory inspected earlier lacks the required method branches, so it cannot stand in for those unverified TtNT productions.

The local method catalog also contains DSID 410470 MC20a full-simulation TtNT `r13167_p7266` (370 paths). It is an alternative DSID, but no high-mass support audit establishes that it is filtered or enriched in large truth mtt. The other catalogs list DSID 410472 `r13167_p6697` TtNT outputs with 164 or 168 paths; these are alternate productions of the same nominal DSID, not proven independent events. Do not combine these files with the `p7266` sample before schema and event-identity checks.

No repository configuration or local catalog entry identified a dedicated mtt-filtered or boosted SM ttbar sample. Metadata-catalogued presence alone does not prove branch compatibility, independence, or useful selected high-mass support.

## Interpretation and actual gap

There is no local validated **nominal** SM control at the same several-TeV support as the 3–6 TeV KK samples. The full nominal MC20a sample may support matched bins to roughly 2–2.5 TeV after selection. The audited MC20e fast-simulation pTHard variation has only 149 raw direct-dilepton events above 3 TeV and none above 5 TeV; it is not interchangeable with the nominal full-simulation denominator. Separate MC23 full-simulation DSID 601491 productions are present in the local metadata catalog but remain unvalidated and lack raw support counts here.

The minimum useful new control is nominal full-simulation SM ttbar processed with the same reconstruction branches, selection, and event weights, with at least 50 finite strict-population events per method and truth-mtt bin through as much of 5–6 TeV as feasible. Fifty is the current plotting support threshold, not a precision target; more events will be needed for stable uncertainty estimates and additional conditioning in truth top-pair pT. No such high-mass-filtered or enriched sample is validated locally.

The highest-information next action is high-mass-only method extraction from the remaining nominal DSID 410472 `r13167` files, using the same strict selection and assignment definitions and excluding the presentation file already used. That gives a same-campaign matched control into the supported roughly 2–2.5 TeV range. Then compare the separately labeled 601491 variation only after its campaign-specific event-weight treatment is resolved. Neither local SM sample provides statistically useful control through 5–6 TeV; that is the remaining genuine data gap. The complete catalog facts and validation states are in `Validation/high_mass_sm_catalogue_audit.json`.

## Filtered method extraction status (2026-10-08)

The filtered extraction path passed the one-file ROOT pilot, including the unchanged strict denominator, EM/NW units, W constraints, and parent ordering. The full command was launched with the 19-file presentation input list plus the exact nominal `r13167` directory and `m_{t\bar{t}}^{truth}>1.6` TeV direct-dilepton cut. Its initial `Path.resolve()` call blocked on a Ceph `lstat`; that unnecessary resolution was removed. A separate `*.root` directory-enumeration check still blocks in the Ceph metadata server, before any ROOT file is read. No full-control method metrics have been written. The pilot's 45 candidate-supported events are a correctness check only and must not be substituted for a high-stat matched result.

The smallest action that unblocks this control is restoring Ceph metadata/read access for the exact `r13167` directory, then rerunning:

```bash
python3 -m Plot.reco_plot.method_study \
  --summary-json plots/presentation/sm_bsm_bias/summary/presentation_metrics.json \
  --high-mass-sm-control /data/leet/ntuples/run2_ttbar_dilep/user.leet.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p7266.TtNT_2L_2Lsys_b71171e_output \
  --bootstrap-resamples 200 \
  --outdir plots/presentation/sm_bsm_bias
```

## Follow-up extraction status (2026-10-08)

The above raw support audit and one-file pilot remain valid. A full serial reconstruction pass was stopped after it had read 50/784 inputs in 13m39s and before it wrote any summary; that partial pass is not used. A bounded eight-worker rerun then blocked in `ceph_mdsc_wait_request` while enumerating the exact directory, before worker creation. Opening the already validated `_000028` pilot by its exact path also blocked in the same Ceph metadata wait. No additional high-mass SM method metrics or matched-response bins were produced. Existing result files and plots are unchanged.

`Plot/reco_plot/method_study.py` now accepts `--input-workers` and exact `--input-filelist` / `--high-mass-sm-control-filelist` options. The latter avoids repeated directory enumeration after an explicit file list has been captured. The serial run established that the method cut and adapter can read the first 50 files, but it was too slow to finish all 784; the parallel path has only passed unit and compile validation and still needs a real-input smoke run once Ceph responds.

The exact local metadata filelist identifies one MC20a `r13167` ROOT file for each of the 17 KK DSIDs. `Validation/kk_mass_grid_mc20a_representative_inputs.txt` combines those eleven previously unrepresented masses with the validated 19-file representative input list. This list has no duplicate paths and has one file per newly added DSID, but the eleven new files have not yet been opened for branch-schema or event-identity validation because Ceph reads are currently blocked.

The local metadata catalog contains all 766 exact paths for the nominal control. `Validation/high_mass_sm_control_inputs_r13167.txt` records the 765 paths remaining after excluding the SM presentation input `_000027`; the list was parsed by the method-study path-list reader and has no duplicate paths. This is a catalog-derived list, not proof that every file is currently readable. A fresh exact-path `stat` of the previously validated `_000028` ROOT file timed out after 4 seconds (exit 124), confirming the same Ceph metadata/read blocker. The nominal SM sample supports raw direct-dilepton truth through 4.227 TeV but has only 72 such events above 3 TeV and none above 5 TeV; even after recovery it cannot control KK response through 5–6 TeV.

After Ceph access returns, run the current representative-sample control extraction with the saved lists. Validate the eleven added mass-grid ROOT files before expanding the analysis to all 17 KK masses.

```bash
python3 -m Plot.reco_plot.method_study \
  --summary-json plots/presentation/sm_bsm_bias/summary/presentation_metrics.json \
  --high-mass-sm-control /data/leet/ntuples/run2_ttbar_dilep/user.leet.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p7266.TtNT_2L_2Lsys_b71171e_output \
  --high-mass-sm-control-filelist Validation/high_mass_sm_control_inputs_r13167.txt \
  --input-workers 8 \
  --bootstrap-resamples 200 \
  --outdir plots/presentation/sm_bsm_bias
```
