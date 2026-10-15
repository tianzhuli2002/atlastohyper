# VyPER Run-2 dilepton / KK-graviton study

Date: 2026-10-05

## Status

The Run-2 branch preflight, dedicated configs, 100-entry SM and KK pilot paths, H5 contract regression, native branch-copy check, validated `VYPER_TRUTHNU` and `VYPER_NUFLOWS` adapters, truth-assignment oracle, shared metrics, and canonical presentation plotting command are available. This is a validated pilot framework, not a full-statistics physics result. No full FastFrames production, external-checkpoint inference, or model training was run in this stage.

The final KK pilot uses DSID 574797 in mc20a, which the supplied DAOD list identifies as the 1 TeV member of `MGPy8EG_A14N23LO_RS_G_tt_2l_AutoWidth_c10`. The mapping for DSIDs 574792–574808 is recorded in `config/kkgraviton_samples.yml`. This mapping does not establish a fixed width, physical cross section, or filter efficiency. `M_G` labels the generated benchmark only; all reconstruction residuals use event-level truth `mtt`.

## Configs and input selection

- Production-oriented minimal H5 config: `FastFramesConfig/TtbarDileptonvyper_KKGraviton.yml`.
- Rich native reconstruction config: `FastFramesConfig/TtbarDileptonReco_VyPERStudy_KKGraviton.yml`.
- Both point to `/home/leet/FF/metadata_kkgraviton/filelist.txt` and `/home/leet/FF/metadata_kkgraviton/sum_of_weights.txt`.
- Both were parsed by `Validation/validate_fastframes_configs.py` under StatAnalysis 0.8.0; the report records one/ two regions and three campaign samples per config, and both metadata paths exist.
- These configs are shape-only (`normalise_to_luminosity: false`). Physics-yield normalization is not enabled.
- The supplied metadata list is mixed across 103 DSIDs. The exact 2L process/mass grid for DSIDs 574792–574808 is resolved by the supplied DAOD list; other DSIDs in the broader 574758–574808 range remain unclassified. AMI is still needed for width, physical cross section, filter efficiency, and complete normalization bookkeeping.
- The candidate filelist has 150 DSID/campaign combinations and only 64 sum-of-weights rows; 86 are missing. All 64 values are finite and positive, which does not establish bookkeeping correctness. AMI requests were blocked by expired/deactivated credentials or role errors. The complete audit is in `Validation/kkgraviton_metadata_audit.md`; the per-DSID candidate manifest is `config/kkgraviton_samples.yml`.

## Schema and native methods

`Validation/preflight_dilepton_inputs.py` checks the raw `reco` tree, required object and truth branches, DA/EM/NW/NuFlows/HyPER/VyPER families, spin branches, and branch types using Uproot rather than the generated MakeClass reader. The checked representatives were SM DSID 410472 mc20d (3,548 raw rows) and KK candidate DSID 574797 mc20a (34,010 raw rows). The report lists all 328 required branch checks and the 56 spin branches. Presence does not imply event-level validity; sentinels occur in native spin values.

Native method outputs are retained directly in the rich ROOT ntuple. The reproduced KK pilot joins by `eventNumber`: all 88 selected output rows match raw rows, and 172 method/spin branches (116 method values and 56 spin values) were copied with zero mismatches. This is a serialization and event-identity check, not an independent comparison of the algorithms. The rich frame now adds assignment/neutrino/reconstruction validity and failure codes for `VYPER_TRUTHNU` and `VYPER_NUFLOWS`; generic `has_native_*` adapters and corresponding flags for the other native methods are still absent.

The first constituent-unit audit found a consequential native EM/NW inconsistency in both representative raw samples. In the first 500 rows, EM and NW lepton constituent pT values are exactly `0.001` times the selected-lepton pT, while their b constituent pT values match selected jets at scale `1`. Their neutrino-pair transverse vector balance with MET is also substantially closer after multiplying the native neutrino four-vectors by 1000 (SM median residuals: EM 7.45 GeV and NW 8.62 GeV, versus 66.16 GeV unscaled; KK candidate: EM 26.99 GeV and NW 14.26 GeV, versus 113.90 GeV unscaled). The stored EM/NW top four-vectors close against the raw numerical constituent sums without conversion, yielding implausibly low median stored top masses (SM 8.10/10.57 GeV; KK candidate 6.06/13.46 GeV). Scaling lepton and neutrino constituents to MeV gives exploratory median masses of 171.00/204.12 GeV (SM) and 161.35/209.13 GeV (KK candidate), but this is not yet an accepted native-method adapter or a validated comparison result. Do not use the stored native EM/NW composite top branches for physics comparisons until the input-production convention is resolved and the explicit unit adapter is independently validated. The 500-row audit JSON files are under `Validation/pilots/` and are copied to `plots/diagnostics/` by the pilot script.

The current evidence does not identify the publication-facing expansion of `DA`; keep the stored label `DA`. EM and NW have native complete reconstruction branches, NuFlows supplies neutrino outputs, and HyPER has native assignments and stored HyPER+EM/NW/NuFlows combinations. Per analysis-owner confirmation, `VYPER_NATIVE` is the existing SM-trained VyPER prediction and is the zero-shot SM→KK result for this initial study. It is not a newly inferred checkpoint identity; future alternative checkpoints must retain separate labels and provenance.

The standalone NuFlows ordering audit establishes that `nuflows_nu_out_NOSYS` is ordered by the selected electron collection in input order followed by the selected muon collection in input order. Mapping positive lepton charge to the top neutrino and negative charge to the antitop antineutrino agrees with the native constituent branches in 88/88 rich-pilot events, including ee/em/mm. The C++ adapter uses the same selected object masks and ordering. Do not substitute pT-sorted charge arrays for this native slot order.

The separate native HyPER+NuFlows top branch still needs an independently validated unit adapter before use as a baseline; its current constituent audit found mixed-scale components.

## Truth, units, matching, and H5 convention

- The pilot truth mass uses partonic `Ttbar_MC_t_afterFSR` and `Ttbar_MC_tbar_afterFSR`, combined as `ttbar_afterFSR`. ROOT four-vector and mass branches are in MeV; plotting converts mass values to GeV. Legacy H5 object four-vectors remain in GeV.
- Legacy H5 jet targets continue to use `Truth_had_b_index` and `Truth_had_bbar_index`, derived from the existing after-FSR parton-match indices and converted to the existing pT-sorted jet order. No new truth-matching definition was substituted.
- The H5 neutrino target remains the existing six Cartesian values in `neutrino_vector_NOSYS`, split into two `(px, py, pz)` entries. Components are in GeV. The validated order is the existing electron-flavour candidate first, then muon-flavour candidate, with remaining slots filled from the parent top then antitop. No H5 target reordering was performed.
- `-fm` preserves the existing filter: two truth b labels plus exactly two selected opposite-charge electrons/muons. It retains ee, mumu, and e-mu channels. The FastFrames H5 region also requires the existing dilepton selection, at least two jets, at least one b-tagged jet, and valid nonzero neutrino inputs.
- The stricter W-daughter lepton match is an evaluation denominator only. It does not change the established H5 fully-matched filtering, labels, ordering, or split behavior.
- Legacy split mode remains `eventNumber % 2`: even IDs go to `_even.h5`, odd IDs to `_odd.h5`. Optional `META` stores event identity and convention attributes without changing `INPUTS` or `LABELS`.
- H5 object four-vectors are stored in GeV; eta, phi, charge, b-tag category, object labels, field order, fixed widths, NaN padding, and float dtypes match the existing converter reference.

## Pilot execution and validation

Reproduce the small run with:

```bash
bash Validation/run_kk_candidate_pilots.sh
```

The script first parses the two full KK configs without running events, then runs the one-file pilots.

Before any pilot event processing, the script checks one direct SM and KK input with `preflight_dilepton_inputs.py` and audits the native constituent units over the first 500 entries. The repeated preflight uses `--shape-only` and skips opening every path in the 49,868-row production filelist; it still parses the full filelist/SOW tables and opens the two selected representative ROOT files. Its report is `Validation/pilots/kk_candidate_preflight_reproduced.json`.

The script uses one representative input file per config and `--max_event 100`. Event-key inspection confirms the FastFrames ROOT outputs contain only selected events from the first 100 input rows. It does not process the remaining entries from those files. Required external known-good H5 reference files are under `/data/leet/FFTutorial/output_h5_*_pilot/`.

The presentation study requires the raw selected-lepton-to-W-daughter match branches, which were not copied by the original 100-event rich pilot. A separate 100-input-event config adds only those match branches and writes `Validation/pilots/artifacts/root/kk574797_reco_leptonmatch_pilot.root`. Reproduce that small branch-complete input (with the repository's configured StatAnalysis/build environment) using:

```bash
python3 /home/leet/FF/fastframes/python/FastFrames.py \
  -c Validation/pilots/TtbarDileptonReco_VyPERStudy_KKGraviton_LeptonMatch_Pilot.yml \
  --step n --max_event 100 --split_n_jobs 1 --job_index 0
cp /data/leet/FFTutorial/output_ntuples_kk574797_leptonmatch_pilot_reco/kk574797_mc20a_pilot_reco_574797_mc20a_fullsim_Njobs_1_jobIndex_0.root \
  Validation/pilots/artifacts/root/kk574797_reco_leptonmatch_pilot.root
```

The resulting 88 selected rows have the 12 required scalar lepton-match index/count branches. The original `kk574797_reco_pilot.root` remains available for legacy smoke tests, but cannot establish the strict selected-W-lepton reconstructibility denominator.

Observed output counts:

| Pilot | Input rows available | First rows processed | Output ROOT rows | H5 rows |
|---|---:|---:|---:|---:|
| SM VyPER, DSID 410472 mc20d | 3,548 | 100 | 72 | 72 (29 even, 43 odd) |
| KK candidate VyPER, DSID 574797 mc20a | 34,010 | 100 | 74 | 74 (36 even, 38 odd) |
| KK candidate rich comparison | 34,010 | 100 | 88 | n/a |

Both SM and KK H5 split files match known-good upstream-converter `INPUTS`/`LABELS` schemas, field names/order, event counts, labels, and NaN masks. Across the 8 datasets in each split, the maximum numerical difference was `5.96e-8` (trigonometric float rounding); labels and neutrino targets match exactly. Categorical jet/electron/muon labels are `float32`. The H5 regression unit test guards these dtype and field contracts.

The independent neutrino check joins by event identity. It confirms H5 labels equal `neutrino_vector_NOSYS` for all 72 SM and 74 KK events. They also match the raw after-FSR W neutrinos in the established order: maximum component differences are `3.24e-5 GeV` (SM) and `8.09e-5 GeV` (KK), with 139/139 and 137/137 lepton-flavour/charge checks passing, respectively. No duplicate target vectors were found.

The artifact ROOT/H5 files and JSON validation reports are under `Validation/pilots/artifacts/` and `Validation/pilots/`. The 1 TeV pilot is the only pilot used for the current presentation suite; no 500 GeV or 5 TeV pilot was run because the existing pilot was sufficient to validate the framework and the 500k productions are still in progress.

## Initial KK pilot observations

An earlier exploratory assignment report used the legacy 74-event b-pair-valid denominator. The present study adds parent-W lepton matching and uses the stricter 63-event truth-reconstructible denominator; the old custom-method numbers are not carried forward as directly comparable results.

For native VyPER, 88/88 selected rows have finite reconstructed and partonic-after-FSR truth `mtt`. Reconstructed median `mtt` is 922.38 GeV versus truth median 979.44 GeV. The median residual is -45.28 GeV; the 16th/84th percentiles are -238.22/+115.91 GeV; the central-68% width is 354.12 GeV. The relative residual median is -4.51%, with 28.4% of these events outside `|relative residual| < 0.2`. DSID 574797 is verified as the 1 TeV `AutoWidth_c10` benchmark, and the native output is the existing SM-trained VyPER zero-shot prediction. These are descriptive pilot observations only; the event-level truth `mtt`, rather than the benchmark mass, defines residuals.

The rich frame now writes `VYPER_TRUTHNU`: it keeps the native VyPER `b`, `bbar`, positive-lepton, and negative-lepton constituents and replaces only the two neutrinos with the existing parent-ordered truth pair (`nu_from_t`, `nu_from_tbar`). On this pilot all 88/88 assignments and truth-neutrino pairs pass the finite/nonzero checks; both b constituent vectors match distinct selected jets, both lepton vectors match the expected charge slots, and the native top/antitop constituent sums close within 0.09 MeV. An independent four-vector sum reproduces the stored oracle `mtt` within 0.031 MeV.

The rich frame also writes `VYPER_NUFLOWS`: it holds the native VyPER b/lepton assignment fixed and substitutes charge-mapped NuFlows neutrinos. All 88/88 pilot rows pass assignment, neutrino, and reconstruction validity with failure code zero. The maximum charge-mapped neutrino kinematic difference from an independent calculation is 0.52 MeV, and independent `mtt` differs by at most 0.18 MeV. The scalar pT, pz, and mass outputs are in MeV; the NuFlows source vector is in GeV.

On all 88 selected rows, native VyPER and VyPER+NuFlows are the deployable/all-selected comparison. Oracle metrics are restricted to the 63 strict truth-reconstructible rows; older all-row `VYPER_TRUTHNU` summaries under `plots/mttbar/` are superseded for the current denominator. The paired component comparison is reported below and is a low-statistics pilot check, not a model-domain conclusion.

The canonical pilot plots under `plots/presentation/pilot_1TeV/` were generated from the branch-complete 88-row ROOT file and visually inspected. Truth-binned metrics retain empty/low-denominator bins; differential display points with fewer than five pairs and spectrum-ratio bins with fewer than ten truth events are masked. The mass response remains statistically coarse. Native spin plots are reconstructed distributions only; this pilot has no truth-spin reference for response or preservation claims. External prediction joins were separately exercised with synthetic CSVs whose rows were reversed; matching uses `(mcChannelNumber, eventNumber)`, and duplicate/missing keys are reported without using row order. There are no real checkpoint CSVs to join.

## Metrics-first support and reconstruction decomposition

The current plotting/metrics implementation lives in `Plot/reco_plot/metrics.py`, `observables.py`, `style.py`, and `study.py`. It uses a strict event-paired contract: reco/truth arrays, support masks, and paired-bootstrap arrays with unequal lengths raise an error instead of being silently truncated.

The selected population is all rows in the `ttll_analysis` tree, whose FastFrames selection is two-lepton region, trigger pass, at least two baseline-JVT jets, and at least one b-tag. No truth-matched filter is applied to this denominator. The stages are:

1. **Truth-reconstructible:** valid distinct unique matched reco jets for the after-FSR top-b and antitop-bbar; exactly two tight selected leptons with one positive and one negative charge, each uniquely matched to a W daughter from its corresponding parent top with the expected charge; finite parent-ordered truth neutrinos and positive finite partonic after-FSR `mtt`.
2. **Candidate-supported:** both H5-aligned, pT-ranked truth jet indices are among the actual baseline-JVT VyPER input candidates after the 20-jet converter cap. The candidate count is measured before truncation. The input is the pT-ordered baseline-JVT jet collection, without a b-tag cut; b-tag is a feature. The converter clips this collection to 20 entries. This is the converter input contract; it does not claim a distinct smaller internal model top-k.
3. **Reconstruction-valid:** `VYPER_NATIVE_reco_valid_NOSYS == 1` and finite positive native `mtt`.
4. **Assignment-correct:** both first native VyPER top/antitop assignment indices match the truth b/bbar ranks, evaluated on the candidate-supported denominator.

Reported fractions are `N_truth_reconstructible/N_selected`, `N_supported/N_truth_reconstructible`, `N_valid/N_selected` (and also valid/supported), and `N_correct/N_supported`. In the metrics table, deployable native outputs are reported on all selected events; oracle methods `VYPER_TRUTHNU` and `TRUTHASSIGN_TRUTHNU` use only the strict truth-reconstructible population. A second `decomposition/` suite masks every method to the same truth-reconstructible events, so component differences are paired on identical events. Other methods' raw DA/EM/NW/HyPER assignment comparisons use the strict truth-reconstructible denominator and raw matched-jet ordering; native EM/NW composite four-vectors are excluded because their constituent scales are unresolved.

`TRUTHASSIGN_TRUTHNU` is the visible-object/reconstructed-object oracle: distinct correctly matched reconstructed b and bbar jets, selected positive and negative leptons, and the parent-top/parent-antitop truth neutrinos. It is separate from `VYPER_TRUTHNU`, which retains VyPER's assignment. The validator independently sums raw constituents and compares with stored `mttbar_TRUTH_NOSYS`; 74 rows pass the b-oracle branch/closure checks (mean difference -0.0032 MeV; maximum absolute difference 0.0619 MeV). Of these, 63 also pass the strict selected-W-lepton pairing definition used for reconstruction response. The residual truth reference remains the partonic after-FSR top pair, so the visible oracle and partonic truth answer different questions.

The mechanism-variable registry also carries event-level truth `pT(t)`/`pT(tbar)` after FSR in GeV, reconstructed physical `N_jets` before model truncation, and the baseline-JVT candidate count before the converter's 20-entry cap. The local training-input converter confirms the actual jet tensor has 20 pT-ordered baseline-JVT candidates, no b-tag-based candidate removal, and truncates only beyond entry 20. No separate smaller top-k is applied in the available converter/model-input path; the native inference checkpoint/config itself is external to this checkout. `min DeltaR(b,l)` is defined as the smaller of the reconstructed matched top-b to selected positive-lepton and matched antitop-bbar to selected negative-lepton distances, using raw matched jet indices and tight-selected object four-vectors. It is defined only for truth-reconstructible events; its efficiency denominator is restricted accordingly. The b-jet/lepton pairing is therefore explicit and does not substitute a nearest-jet choice.

The mass residual is `(mtt_reco - mtt_truth)/mtt_truth`. `sigma68` always means the full central width `q84 - q16`. Relative tail fractions use strict thresholds `< -0.2`, `> +0.2`, and `abs(residual) > 0.2`. Angular `DeltaPhi(t,tbar)` is absolute wrapped separation in `[0, pi]`; periodic angular residuals use wrapping to `[-pi, pi)`. The signed production-angle proxy is `cos(theta*) = tanh((y_t-y_tbar)/2)`, preserving the top direction relative to the positive beam axis; it is not folded to absolute value. Other listed core observables use central registry definitions in `observables.py`.

Algorithm-performance metrics are unweighted. Spectrum plots explicitly choose `--weighting unweighted` or `--weighting physics`; either is normalized over all finite events/weights, including values outside the displayed histogram range. Under/overflow counts and the visible integral are recorded, so a cropped axis cannot silently renormalize away tails. This command does not produce absolute yield normalization. Cross section/filter-efficiency normalization remains unavailable and is never inferred. Metrics and JSON record selection, denominator, truth definition, unit, and weighting. Efficiencies use 68% Wilson intervals. Migration and continuous-response maps divide by all valid reco/truth pairs in each displayed truth bin; reco/response underflow and overflow stay in the denominator and are reported. Residual summaries include deterministic bootstrap intervals, and method differences use event-paired bootstrap with the same event resampling indices. The default is 200 replicas for quick pilot validation; increase this for final full-stat results. Bootstrap intervals from the 88-event pilot are software checks, not useful final uncertainties.

## Canonical plotting command and outputs

From the repository root:

```bash
python3 Plot/compare_reco.py study \
  Validation/pilots/artifacts/root/kk574797_reco_leptonmatch_pilot.root \
  --sample-manifest config/kkgraviton_samples.yml \
  --outdir plots/presentation/pilot_1TeV \
  --weighting unweighted \
  --bootstrap-resamples 200
```

The command produces per-sample outputs under `plots/presentation/pilot_1TeV/kk_574797_mg1000/` and the aggregate summary under `plots/presentation/pilot_1TeV/summary/`. Plot families include assignment and candidate support vs physical pre-truncation `N_jets`; reconstructibility/support/validity/assignment rates vs truth `mtt`, top and antitop `pT`, min `DeltaR(b,l)`, physical `N_jets`, and candidate-jet count; `mtt` bias and resolution vs topology; score-vs-correctness diagnostics; spectra with ratio panels; residuals; differential median and central-width response with bootstrap intervals; per-truth-bin migration matrices; continuous `mtt_reco/mtt_truth` maps; resonance-mass scan/linearity scaffolding; and native spin distributions. The additional `decomposition/` directory repeats the central `mtt` spectra, residual, differential response, migration, and continuous response on the same strict truth-reconstructible event population for every method. The 1 TeV mass scan contains exactly one point; missing mass points are not fabricated. The mechanism and topology figures reserve a dedicated legend area, and the differential plots consume the shared differential metrics rather than reimplementing residual definitions.

The strict pilot ladder is `88 selected -> 63 truth-reconstructible -> 63 candidate-supported -> 88 finite-valid native reconstructions`; both b jets are in candidate support for all 63 strict events, and 57/63 supported native assignments are correct. The component checks give 74/88 unique b-pair matches and 75/88 selected W-lepton pairs before requiring both components together. No event reaches the 20-jet cap (maximum candidate count: 7), so this pilot does not validate high-multiplicity truncation behavior. On the same 63-event population, relative `mtt` median / `sigma68` are -3.74% / 36.48% for native VyPER, -3.12% / 32.81% for VyPER+NuFlows, -1.95% / 7.07% for VyPER+truth neutrinos, and -1.81% / 6.66% for truth assignment+truth neutrinos. The last comparison exposes a small residual assignment/oracle difference; these 63 events are only a software/pilot check, not a method ranking. Inclusive all-selected native VyPER has 88 pairs, -4.51% median relative residual, 36.54% `sigma68`, and low/high/absolute tail fractions 19.3%/9.1%/28.4%. The score plot combines the minimum of the two scalar top/antitop score branches for a diagnostic only. It is not a calibrated probability; no candidate-score margin is reported because the candidate score ordering semantics are not established in the current output contract.

The regenerated JSON contains all metric rows and bootstrap intervals. CSV outputs are `summary/metrics.csv`, `summary/differential_metrics.csv`, `summary/assignment_metrics.csv`, `summary/assignment_differential_metrics.csv`, `summary/topology_response_metrics.csv`, and `summary/mass_scan_metrics.csv`; `summary/study_metrics.json` contains provenance, definitions, and plot paths. `summary/truthassign_truthnu_validation.json` is produced by the separate closure command below rather than by the plotting command.

The native spin figure covers all eight copied native spin quantities and the available DA/EM/NW/HyPER/VyPER branches, but the pilot has no truth-spin branches. It is therefore a reconstructed-distribution diagnostic only; no angular/spin preservation claim is possible yet. Native EM/NW and HyPER+NuFlows constituent-unit problems also prevent full-reconstruction comparisons with those branches. The validated central decomposition presently comprises `VYPER_NATIVE`, `VYPER_NUFLOWS`, `VYPER_TRUTHNU`, `TRUTHASSIGN_TRUTHNU`, and partonic truth.

## Current limitations and next full-stat commands

- No full FastFrames production, external checkpoint inference, or BSM retraining was started in this stage. The 500k files are still being produced.
- No conclusion about SM→BSM degradation, mass dependence, matched-SM response, or generator dependence is supported by this low-stat pilot. The pilot is the 1 TeV KK sample only.
- `VYPER_EM`/`TRUTHASSIGN_EM` and assignment-conditioned NW are not implemented. They are not required to run the present native VyPER decomposition.
- EM/NW primitive and composite units remain unresolved for physics reconstruction. Their assignment-only results are available with explicitly stated denominators.
- Native spin distributions are copied and plotted, but truth spin values are absent, so truth/reco response and angular bias are blocked. A validated common hybrid-spin calculator is not present.
- `DA` expansion remains unknown; keep the native branch label `DA`.
- The alternative-generator SM, high-mass SM support control, and SM/KK overlap/reweighting test are not run. The plotting metrics can consume the same observables from verified inputs once provided. Extremely high KK masses with negligible SM overlap must be reported as extrapolation, not reweighting evidence.
- AMI-level physical cross sections/filter efficiencies and complete campaign sums of weights remain unresolved. Shape-normalized comparisons are available; physics yields are not.

When the 500k outputs finish, first analyze the representative grid points (400/500, 1000, 2000, 3000, and 5000 GeV) by passing their rich ROOT files together:

```bash
python3 Plot/compare_reco.py study \
  <rich-ROOT-400-or-500-GeV.root> \
  <rich-ROOT-1000-GeV.root> \
  <rich-ROOT-2000-GeV.root> \
  <rich-ROOT-3000-GeV.root> \
  <rich-ROOT-5000-GeV.root> \
  --sample-manifest config/kkgraviton_samples.yml \
  --outdir plots/presentation/kk_priority_mass_scan \
  --weighting unweighted \
  --bootstrap-resamples 1000
```

For a DSID represented by multiple campaign files, pass all campaign outputs for that DSID to the same command; identity checks report duplicate `(mcChannelNumber,eventNumber,campaign)` keys and reject duplicated inputs. To add SM or generator-control ROOT files, first add verified sample entries to a manifest (do not assign process/domain labels from filenames alone). The command then groups by DSID and records input campaigns. After inspecting the priority grid, rerun over all 17 listed KK masses; the mass-scan plots summarize the trend without producing 17 copies of each presentation spectrum.

Run the independent oracle closure on any new rich ntuple before analysis:

```bash
python3 Validation/validate_truthassign_truthnu.py \
  <rich-ROOT.root> \
  --json <output-dir>/truthassign_truthnu_validation.json
```

No production-scale or retraining action is currently required to validate the plotting framework. Before 500k analysis, retain the existing required FastFrames pilot gate, inspect each new ROOT schema and event identity, and use the independent oracle closure and generated plot suite as above.

## SM-vs-KK presentation suite

The focused presentation command is separate from the broad diagnostic tree and writes PNG, PDF, JSON and CSV outputs under `plots/presentation/sm_bsm_bias/`. It uses the strict direct after-FSR W-to-e/mu truth-reconstructible population with VyPER candidate support; the event-level after-FSR partonic top-pair mass is the residual reference. The reco-object oracle remains a separate simulation-only reference built from truth-matched reconstructed b jets and selected leptons with parent truth neutrinos. A sample summary can be rendered again without reopening ROOT files.

The checked-in pilot fixtures still do not contain a fair SM/KK pair: `kk574797_reco_leptonmatch_pilot.root` has the complete rich schema (88 selected rows, 63 supported truth-reconstructible rows), while `sm410472_vyper_pilot.root` is reduced. The separate 500k rich-output directory remains absent. However, a narrow raw-ntuple adapter now provides the event-level after-FSR partonic `mtt`, selected-object and truth-match information, native SM-trained VyPER top vectors, and the inputs needed to build the reco-object oracle. It replays the exact `ttll_analysis` selection and is closed against the first 100 KK 1 TeV raw event identities: 88 selected, 63/63 strict events supported, native `mtt` agrees at floating-point precision, and the oracle closes within 0.062 MeV. The raw presentation adapter therefore supports a real, though incomplete, SM/KK comparison without another FastFrames production.

The current representative input list is `Validation/representative_raw_topology_filelist.txt`: one SM DSID 410472 file and three event-key-disjoint MC20a derivation-tag files for each of KK 400, 500, 1000, 2000, 3000 and 5000 GeV. These are partial raw file slices, not complete samples or campaigns. Within each DSID, the plotting adapter rejects duplicate `(mcChannelNumber,eventNumber)` identities before combining files. The raw results are unweighted event-level reconstruction metrics; the displayed spectra are unit-area shapes. They should be read as representative slice comparisons, not final full-statistics performance.

Reproduce the current visual/functional pilot validation:

```bash
python3 Plot/compare_reco.py presentation \
  Validation/pilots/artifacts/root/kk574797_reco_leptonmatch_pilot.root \
  --sample-manifest config/kkgraviton_samples.yml \
  --outdir plots/presentation/sm_bsm_bias \
  --bootstrap-resamples 200
```

Run the representative raw-slice comparison, including all files in the checked-in list:

```bash
mapfile -t inputs < <(rg '^[^#[:space:]]' Validation/representative_raw_topology_filelist.txt)
python3 Plot/compare_reco.py presentation "${inputs[@]}" \
  --raw-input \
  --sample-manifest config/kkgraviton_samples.yml \
  --outdir plots/presentation/sm_bsm_bias \
  --bootstrap-resamples 200
```

Re-render the figures from their saved plot-ready summary without event-level processing:

```bash
python3 Plot/compare_reco.py presentation \
  --summary-json plots/presentation/sm_bsm_bias/summary/presentation_metrics.json \
  --outdir plots/presentation/sm_bsm_bias
```

For a full-campaign update, use audited rich outputs after they appear, schema-audit their truth/native/oracle branches, and use the validated identity-keyed lepton-match sidecar only if those match scalars are absent. The command requires identical truth-mtt bin edges and masks domain-difference points unless both samples have at least 50 finite response pairs in that bin. It records exact SM/KK counts and the support range for each plotted summary point. At high truth mass, bins without populated SM support remain masked.

Current representative raw-slice outputs are under `plots/presentation/sm_bsm_bias/`. Figures 01–03 compare SM and KK response only in common truth-`mtt` bins; figures 04–05 summarize inclusive sample-slice response and its separate reco-object-oracle floor; figures 06–08 show unit-area `mtt` spectra, row-normalized response matrices, and residual shapes for selected representative samples. `figure_manifest.md` and `summary/` record population, truth/reference definitions, binning, support counts, histogram counts, migration probabilities, and uncertainties. The raw scope and per-sample file/count information are preserved in `summary/presentation_metrics.json`.

### First SM-vs-KK response result (2026-10-06)

The matched-domain presentation uses the 19-file representative raw list (506,859 raw entries total), grouped into one SM slice and KK 400, 500, 1000, 2000, 3000 and 5000 GeV slices. These are partial MC20a inputs, not the full campaign or the absent 500k rich production. The performance population remains direct after-FSR W-to-e/mu, selected, truth-reconstructible and VyPER-supported. Selected / strict-supported counts are: SM 41,120 / 28,608; KK400 75,156 / 49,461; KK500 94,762 / 62,088; KK1000 116,365 / 75,147; KK2000 52,327 / 25,896; KK3000 21,430 / 8,833; KK5000 11,036 / 4,854. In these inputs every strict reconstructible event is within candidate support.

The inclusive SM-trained VyPER median bias / sigma68 / absolute `|Delta m/m| > 0.2` tail are respectively: SM `-3.82% / 31.05% / 23.58%`; KK400 `-0.66% / 30.10% / 19.38%`; KK500 `-2.72% / 30.89% / 21.92%`; KK1000 `-5.79% / 34.14% / 28.65%`; KK2000 `-10.69% / 40.17% / 39.78%`; KK3000 `-11.26% / 43.03% / 43.45%`; KK5000 `-9.67% / 40.85% / 39.35%`. The reco-object oracle sigma68 stays between 9.35% and 12.55% over the same samples, while the VyPER-minus-oracle width gap grows from about 18 pp for SM/KK400–500 to about 31 pp for KK2000–5000. Thus the inclusive high-mass loss is substantially above the reconstructed-object floor in these slices.

At common event-level truth-mtt, using bins with at least 50 finite pairs in both domains, the Nmin-weighted KK-minus-SM median-bias differences are `+2.47`, `+1.28`, `+3.85`, `+4.11`, `+4.95`, and `+2.83` percentage points for KK400, 500, 1000, 2000, 3000, and 5000, respectively. Their conservative 68% envelopes are above zero; matched sigma68 differences range from `-0.15` to `+1.42` points and all envelopes include zero. Common support is only 350–1600 GeV truth-mtt for KK400–2000 and 500–1600 GeV for KK3000/5000; the several-TeV resonance regions therefore have no SM comparison here. The low-mass control shows a small, positive median-response difference but no resolved conditional broadening. Since the comparison conditions on truth-mtt only and uses partial file slices, it is evidence of a response difference at matched truth-mtt, not by itself proof of a BSM-specific learned bias.

Plots 01–08 were regenerated from the saved event-level summary and visually checked. The log-scale spectrum ratio exposes the strong low-mass migration without flattening ratios near unity; differential domain points below the explicit 50-pair common-support cut are masked. Machine-readable values and support are in `plots/presentation/sm_bsm_bias/summary/`; the concise figure definitions are in `plots/presentation/sm_bsm_bias/figure_manifest.md`.

## Stage validation snapshot (2026-10-05)

| Area | Status | Evidence / scope |
|---|---|---|
| 2L DSID→mass map | PASS | 574792–574808 mapped; 574797 is 1000 GeV; width and XS fields remain unresolved. |
| Existing Run-2 FastFrames/H5 contract | PASS | Prior first-100 pilots and H5 regression remain intact; this stage did not alter that interface. |
| Candidate support contract | PASS | Baseline-JVT pT-ordered H5 input, no b-tag cut, converter cap 20 established from converter/frame code; pilot candidate count inspected. |
| Reconstructibility/support/valid/correct ladder | PASS | Strict 88→63→63→88 populations; 57/63 VyPER assignments correct. Component counts and denominators are exported to JSON/CSV. |
| `TRUTHASSIGN_TRUTHNU` | PASS | Independent raw-branch closure on 74 b-oracle-valid rows; strict response population is 63 after selected W-lepton matching; max closure difference <0.062 MeV. |
| Shared metrics and uncertainty utilities | PASS | The repository unit/regression suite passes; it covers paired residuals, angular wrap, relative masks, sigma68, signed tails, Wilson intervals, support denominators, deterministic/paired bootstrap, population masking, and shape mismatch failures. |
| Canonical pilot suite | PASS | `study` CLI regenerated all-selected and paired truth-reconstructible plots/JSON/CSV for the actual 1 TeV rich ROOT pilot; spectra, differential response, migration, support and topology figures inspected. |
| Candidate-score calibration | NOT CLAIMED | Raw score diagnostic only; scalar score not calibrated and candidate score margin semantics are unknown. |
| Truth/reco spin response | BLOCKED | No truth-spin branches in the available pilot. Native spin spectra only. |
| Full-stat SM/KK conclusions | NOT RUN | 500k outputs remain in progress; pilot statistics are not used for conclusions. |
| Full FastFrames or BSM training | NOT RUN | Explicitly out of scope for this stage. |

## Recovery and next physics diagnostics (2026-10-05)

After the interrupted session, the branch-complete pilot and source changes were intact. Python compile checks passed, the existing 29 Python tests passed before this extension, and the canonical pilot regenerated 88 selected, 74 unique b-pair matches, 75 W-lepton pair matches, 63 strict truth-reconstructible, 63 supported, and 57 correct native assignments. Both the all-selected and paired truth-reconstructible plot suites remained present. No partial source edit needed repair. The existing working tree also contains older C++/FastFrames edits and generated build files; they were preserved.

The new mutually exclusive assignment categories are evaluated on supported strict events: `CORRECT` has the ordered truth b/bbar pair, `SWAPPED` has the same two truth jets with reversed parents, `ONE_WRONG` has exactly one correct parent-side jet, `BOTH_WRONG` has neither correct parent-side jet after removing swaps, and `INVALID` has no valid distinct candidate pair. `UNSUPPORTED` is counted separately against the truth-reconstructible denominator; selected events outside that population are not mislabeled as unsupported. The 1 TeV pilot has 57 correct, 1 swapped, 5 one-wrong, 0 both-wrong, 0 invalid among 63 supported, and 0 unsupported among 63 strict reconstructible events.

Individual top and antitop `pT`, rapidity, mass, and parent-paired `m(b,l+)`/`m(bbar,l-)` now use the same strict event set across native VyPER, VyPER+NuFlows, VyPER+truth neutrinos, and the visible-object truth-assignment oracle. The raw visible oracle closes against stored `mttbar_TRUTH_NOSYS` to better than 1 MeV in the strict population. Unlike total `mtt`, these observables retain sensitivity to parent assignment. On 63 pilot events, native top-`pT` relative sigma68 is 64.6%, falling to 8.9% with truth neutrinos; native top-mass sigma68 is 130.2%, falling to 11.6%. Top-side `m_bl` sigma68 stays 32.7% across the three VyPER neutrino routes because their visible constituents are identical, then becomes 31.5% for the truth-assignment oracle. These numbers are pilot diagnostics only.

Parent-ordered neutrino diagnostics now record `Delta px`, `Delta py`, `Delta pz`, `Delta pT`, `Delta eta`, wrapped `Delta phi`, `Delta E`, `Delta |pz|`, and wrong-`pz`-sign fractions. For each neutrino route and parent, `Delta pz` is binned against truth `pz`, event-level truth `mtt`, and truth top `pT`. Pilot native `Delta pz` sigma68 is 254/157 GeV on top/antitop sides versus 144/127 GeV for NuFlows; wrong-sign fractions are 30.2%/30.2% native and 19.0%/19.0% NuFlows. At fixed VyPER visible assignment, the pilot correlation of the neutrino-induced `mtt` shift with summed `Delta E_nu` is 0.924 native and 0.936 NuFlows; with summed `Delta pz` it is only 0.218/0.339. This suggests that neutrino-energy error is a more direct total-mass mechanism than signed summed `pz` alone. The correlation is descriptive and does not resolve nonlinear dependence or establish a full-stat ranking.

The production rich-output directory configured in `FastFramesConfig/TtbarDileptonReco_VyPERStudy_KKGraviton.yml` (`/data/leet/FFTutorial/output_ntuples_ttbarDilepton_reco_VyPERStudy_KKGraviton`) does not exist as of this audit. No 500k rich output can yet be assigned a branch-completeness verdict, and no full-stat SM/KK reconstruction or matched-SM domain comparison was run. The current production YAML includes the 12 strict lepton-match branch patterns. The older 88-row rich pilot lacks them; the current branch-complete pilot has them. `Validation/audit_rich_output_schema.py` confirms that both pilots have the core identity/truth/object/native/mechanism groups, and that only the older pilot needs a sidecar. The sample manifest now includes an explicit SM DSID 410472 baseline label, without assigning it a resonance mass or inferring unaudited generator details.

The raw input has all 12 required lepton-match scalars and `(mcChannelNumber,eventNumber)`. `Validation/build_lepton_match_sidecar.py` writes only these scalars plus identity/campaign to a compressed NPZ. The study joins by exact identity, rejects duplicate or missing keys, and requires campaign disambiguation when the same channel/event pair occurs across campaigns. A sidecar built from the 34,010-row raw 1 TeV file reproduces the branch-complete pilot's 88 identities and all strict reconstructibility/support/category masks from the older pilot: 88/63/63/57. This establishes a safe lightweight route if completed 500k rich outputs were made before the lepton-match pattern was added. Do not rerun full FastFrames solely for these scalar branches.

The raw pre-selection ladder uses direct after-FSR W daughters of flavours e or mu on both parents as its truth denominator; tau decays are outside that denominator. It replays the frame's exactly two tight `pT>25 GeV` opposite-charge leptons, at least two baseline-JVT jets, at least one baseline 77%-tagged jet with `pT>30 GeV`, and the mc20 single-lepton trigger OR. Its selected event keys match the branch-complete 100-input-event pilot exactly (88/88, no missing or extra keys). It then counts unique raw matched b jets and selected W-lepton matches. Raw combined matching lacks the rich oracle/neutrino-validity checks, so it is a proxy before full rich processing; on the 100-row pilot it coincides with the strict 63 events.

One existing mc20a raw file per DSID gives this **diagnostic slice**, not a 500k or campaign result:

| Sample | Raw rows | Direct truth 2L | Analysis selected / direct truth 2L | Unique reco b pair / all selected | Combined raw match / all selected |
|---|---:|---:|---:|---:|---:|
| SM 410472 | 49,386 | 44,621 | 83.3% | 77.1% | 69.6% |
| KK 400 | 23,420 | 20,914 | 82.1% | 74.2% | 66.5% |
| KK 500 | 28,737 | 25,330 | 83.0% | 74.9% | 66.2% |
| KK 1000 | 34,010 | 28,814 | 85.9% | 76.6% | 64.7% |
| KK 2000 | 16,784 | 14,038 | 77.7% | 59.3% | 49.2% |
| KK 3000 | 7,965 | 6,586 | 67.2% | 49.8% | 40.4% |
| KK 5000 | 4,108 | 3,440 | 67.7% | 53.7% | 44.5% |

The truth-parton minimum `Delta R` between each after-FSR b and its charged W daughter supplies a topology diagnostic before reconstruction. In these slices, the share of direct-truth-2L events with min `Delta R<0.4` rises from 1.1% at 500 GeV to 34.4% at 2 TeV and 50.3% at 3 TeV (46.1% at 5 TeV). Among selected events in that close-pair region, unique reconstructed b-pair matching is 32.1%, 10.7%, 9.0%, and 10.7% at 500, 2000, 3000, and 5000 GeV. Truth-`mtt`, top-`pT`, and min-`Delta R` differential acceptance and b-pair counts are saved with explicit bin denominators. This is evidence for an early resolved-object/reconstructibility limitation in the available slices; it does not yet measure VyPER assignment at those masses.

Reproduce the current diagnostics from the repository root:

```bash
python3 -m unittest discover -s Tests -p 'test_*.py' -q
python3 Plot/compare_reco.py study Validation/pilots/artifacts/root/kk574797_reco_leptonmatch_pilot.root --sample-manifest config/kkgraviton_samples.yml --outdir plots/presentation/pilot_1TeV --weighting unweighted --bootstrap-resamples 200
python3 Validation/audit_raw_acceptance.py --filelist Validation/representative_raw_mc20a_filelist.txt --output Validation/raw_acceptance_representative.json
python3 Validation/summarize_raw_acceptance.py Validation/raw_acceptance_representative.json --outdir plots/presentation/raw_acceptance
python3 Validation/audit_rich_output_schema.py Validation/pilots/artifacts/root/kk574797_reco_pilot.root Validation/pilots/artifacts/root/kk574797_reco_leptonmatch_pilot.root --output Validation/pilots/kk574797_rich_schema_audit.json
python3 Validation/build_lepton_match_sidecar.py /data/leet/ntuples/run2_ttbar_dilep/user.leet.574797.MGPy8EG.DAOD_PHYS.e8588_s3797_r13167_p7266.TtNT_2L_2Lsys_b71171e_output/user.leet.52697131._000002.output.root --output Validation/pilots/kk574797_lepton_match_sidecar.npz
```

When the configured rich-output directory appears, audit the exact output paths before the study:

```bash
python3 Validation/audit_rich_output_schema.py /data/leet/FFTutorial/output_ntuples_ttbarDilepton_reco_VyPERStudy_KKGraviton/*.root --output Validation/full_rich_schema_audit.json
```

For the representative full-stat study, pass the exact audited SM 410472 and KK 574794/574797/574801/574803/574807 ROOT paths as positional arguments to the canonical `study` command above, change `--outdir` to `plots/presentation/representative_full`, and use `--bootstrap-resamples 1000`. Supply `--lepton-match-sidecar <matching-raw-sidecar.npz>` only if the schema audit says the strict lepton-match group is absent. Inspect representative results and input identity before supplying all 17 KK DSIDs. The exact full-stat filenames cannot be written yet because the configured directory is absent; guessing production paths would risk analyzing wrong or incomplete files. A matched-SM control likewise requires those rich SM and low-mass KK outputs. No model-training intervention is justified by the pilot or raw slices alone.

New machine-readable outputs are under `plots/presentation/pilot_1TeV/summary/`: `assignment_categories.csv`, `assignment_category_differential.csv`, `topside_metrics.csv`, `neutrino_metrics.csv`, `neutrino_differential.csv`, and `neutrino_mass_relation.csv`, alongside the prior JSON/CSV. New compact figures are in the sample's `physics/` directory and the aggregate `mass_scan/` directory. Raw stage reports and sidecar evidence are in `Validation/`; their aggregate CSV and figures are in `plots/presentation/raw_acceptance/`.

## Causal neutrino diagnostics and multi-tag raw topology study (2026-10-05)

The canonical pilot suite now also evaluates the two neutrino routes on both the strict `truth_reconstructible` population and the paired `truth_reconstructible_and_assignment_CORRECT` population. It adds two simulation-only component oracles with the native visible assignment: `VYPER_TRUTHNU_PT` substitutes both truth-neutrino transverse components `(px,py)` while keeping VyPER `pz`; `VYPER_TRUTHNU_PZ` substitutes truth `pz` while keeping VyPER `(px,py)`. Neutrino energies are recomputed as massless four-vectors. `VYPER_METPROJECTED` is a diagnostic-only equal-share projection of the two predicted transverse vectors onto reconstructed MET; it leaves both `pz` values unchanged and also recomputes energy. These are causal probes, not deployable reconstruction routes.

On the 1 TeV pilot, the strict population remains 63/88 selected events and contains all 63 supported events; 57 are correctly assigned. On those 57 correct-assignment events, relative `mtt` `sigma68` is 40.6% for native VyPER and 33.9% for NuFlows. Replacing truth `(px,py)` gives 11.3%, replacing truth `pz` gives 33.0%, and replacing both gives 7.0%; the MET-projected diagnostic gives 38.5%. The paired width recovered by full truth-neutrino replacement is 33.6 percentage points for native and 27.0 points for NuFlows. The broad native and NuFlows neutrino errors therefore persist after removing all six assignment failures. Within this low-stat pilot, transverse-neutrino prediction/partition has the largest recoverable effect on total `mtt`; the pz-only oracle helps less, and imposing pair-sum MET closure alone changes the width only slightly. Do not promote these pilot values to full-stat rankings.

The direct parent-neutrino residuals support that diagnosis. On correct-assignment events, native `Delta pz` `sigma68` is 240 GeV for the top neutrino and 165 GeV for the antitop neutrino; wrong-sign fractions are 31.6% (18/57) and 26.3% (15/57). NuFlows gives 130 and 131 GeV, with 19.3% (11/57) and 14.0% (8/57) wrong signs. Native median pair-to-MET mismatch is 40.0 GeV (relative mismatch 0.44), compared with 29.0 GeV (0.31) for NuFlows. The equal-share projection closes the pair sum to floating-point precision, but does not recover the individual-neutrino partition and yields only a small `mtt` width change. The truth-neutrino pair itself does not close exactly to reconstructed MET because detector-level MET is not the truth-neutrino pair.

The correct-assignment neutrino-energy response ratios are broad: native median `E_reco/E_truth` is 1.00 (top) and 0.94 (antitop), with central-68% widths 2.27 and 1.66. In truth-energy bins, native median response falls below one above 100 GeV (about 0.88 in 100–200 GeV, with 17 top-side and 21 antitop-side entries); the 200–400 GeV medians are about 0.34 and 0.46 (14 and 17 entries), and higher bins have too few events to support a scale trend. This is a pilot hint of high-scale compression, not evidence for an SM-specific bias.

Local parent-mass consistency is also recorded on the same 57 events. For native VyPER, `m(l nu)` residual `sigma68` is 140 GeV on the top side and 233 GeV on the antitop side; `m(b l nu)` widths are 217 and 252 GeV. NuFlows reduces these to about 11 GeV for `m(l nu)` and 35–39 GeV for `m(b l nu)`. Full truth-neutrino replacement gives 1.4–1.5 GeV and 18–19 GeV, respectively. Individual-top observables remain assignment-sensitive: the native top and antitop `pT` relative widths are about 63% and 62%, falling to 8.8% and 8.5% with truth neutrinos. The fixed-visible-assignment neutrino oracles therefore identify large conditional neutrino error even when the two b jets are paired correctly.

The raw cutflow was extended from one-file slices to 19 explicitly listed MC20a raw inputs: one SM input and three non-overlapping input files for each of KK 400, 500, 1000, 2000, 3000, and 5000 GeV. The inputs contain 506,859 rows. The three KK derivation tags are `r13144`, `r13145`, and `r13167`; they are input derivation tags, not separate MC campaigns. The `(mcChannelNumber,eventNumber)` audit found no within-file duplicate keys and no cross-file overlap within any DSID (`Validation/representative_raw_topology_identity_audit.json`). These are not the missing rich VyPER outputs; the presentation response metrics are produced separately by the raw-branch adapter documented above.

The truth denominator is explicitly direct/prompt after-FSR `W -> e/mu` on both top-parent sides; it excludes `W -> tau` events. The selected analysis can contain events with a W tau daughter, counted separately as a proxy category because the inspected W-daughter branches do not identify the tau decay's final e/mu descendant. This category is 9.6% of all selected SM events, rising to 15.2–16.8% at 1–5 TeV, so it is reported separately and is not folded into direct-e/mu efficiency. The raw replay uses exactly two tight selected opposite-charge leptons with `pT > 25 GeV`; the configured single-lepton trigger OR; `n_jets_NOSYS >= 2` baseline-JVT jets (the frame count applies no extra jet-pT threshold); and at least one baseline-JVT, 77%-tagged jet with `pT > 30 GeV`. Both the direct-truth-2L selected count and all selected events are saved. Lepton-match and b-pair efficiencies use the direct-truth-2L selected denominator. Counts and nominal `weight_mc_NOSYS` efficiencies are both saved; for these samples, weighted and unweighted selected/direct-truth efficiencies differ by at most about 0.12 percentage points.

The updated cutflow was smoke-tested on the first 100 raw 1 TeV input rows against the branch-complete rich pilot: 85 are direct truth 2L, 75 pass selection within that denominator, and 88 pass the complete selection across truth categories. All 88 selected event keys match the rich pilot exactly; 13 selected events have a W tau daughter, and 63 pass the raw unique-b-pair plus selected-lepton-match stages. This validates that the direct-e/mu denominator and all-selected tau category remain separate.

| Sample | Direct truth 2L | All selected | Tight OS pair / truth 2L | Trigger / pair | >=2 JVT jets / trigger | >=1 b-tag / prior stage | Direct-2L selected / truth 2L | Unique b pair / direct-2L selected | Strict object match / direct-2L selected |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SM | 44,621 | 41,120 | 99.6% | 99.94% | 100% | 83.7% | 83.3% | 77.2% | 77.0% |
| KK 400 | 82,330 | 75,156 | 99.6% | 99.99% | 100% | 82.1% | 81.7% | 73.7% | 73.5% |
| KK 500 | 100,902 | 94,762 | 99.5% | 99.99% | 100% | 83.3% | 82.8% | 74.4% | 74.3% |
| KK 1000 | 114,696 | 116,365 | 99.1% | 100% | 100% | 86.8% | 86.1% | 76.2% | 76.1% |
| KK 2000 | 56,595 | 52,327 | 98.3% | 100% | 100% | 78.7% | 77.4% | 59.2% | 59.1% |
| KK 3000 | 26,449 | 21,430 | 97.0% | 100% | 100% | 69.5% | 67.4% | 49.7% | 49.5% |
| KK 5000 | 13,675 | 11,036 | 96.3% | 99.97% | 100% | 69.9% | 67.3% | 53.0% | 52.7% |

The explicit event selection and trigger do not cause the largest high-mass loss in these slices: trigger and the `>=2` baseline-JVT multiplicity stage are nearly fully efficient once the previous stage passes. Tight-lepton-pair efficiency falls from 99.5–99.6% at low mass to 96.3% at 5 TeV, while the conditional b-tag stage falls to about 69–70% at 3–5 TeV and unique truth b-pair matching among selected events falls from 74–77% at 400–1000 GeV to about 50–53% at 3–5 TeV. The b-tag and unique-b-pair stages therefore carry most of the measured loss; the raw cutflow alone cannot identify whether merging, overlap removal, or another jet response effect causes each lost match.

Truth `min DeltaR(b,l)` makes the resolved-topology effect visible. The fraction of direct truth-2L events with `min DeltaR < 0.4` is 1.1% at 400–500 GeV, 6.3% at 1 TeV, 34.1% at 2 TeV, 50.3% at 3 TeV, and 46.5% at 5 TeV. Within this close-pair region, the direct-2L selected fraction is 61% at 2 TeV, 51% at 3 TeV, and 48% at 5 TeV; unique b-pair matching among those selected events is only 10.4%, 9.1%, and 10.1%. The curves do not collapse at fixed separation: in the `0.7 < min DeltaR < 1.0` bin, selected unique-b-pair fractions are 81.3% for the SM slice, 75.0% at 1 TeV, 64.7% at 2 TeV, 62.1% at 3 TeV, and 66.9% at 5 TeV. DeltaR strongly tracks the loss but does not explain all of the mass dependence. Among selected events with identifiable matched b/lepton counterparts, reconstructed-minus-truth minimum DeltaR has median shifts within about 0.002 across the KK samples and `sigma68` about 0.039–0.055; this resolution statement applies only to matched survivors.

The raw-topology and pilot-neutrino conclusions above were reached before the representative response plots. The subsequent partial-slice presentation comparison now supplies provisional conditional VyPER results as described above, but it still does not constitute a full-campaign result or justify a training intervention.

The configured rich production directory was checked again on 2026-10-06 and remains absent: `/data/leet/FFTutorial/output_ntuples_ttbarDilepton_reco_VyPERStudy_KKGraviton` has zero ROOT files. The check is recorded in `Validation/full_rich_output_status.json`; therefore no full-rich schema audit, representative full-stat reconstruction, or matched-SM control could be run. Once those files appear, audit their exact paths with `Validation/audit_rich_output_schema.py`; attach the validated lepton-match sidecar only if the audit reports those scalar branches missing.

Reproduce the current outputs from the repository root:

```bash
python3 -m unittest discover -s Tests -p 'test_*.py' -q
python3 Plot/compare_reco.py study Validation/pilots/artifacts/root/kk574797_reco_leptonmatch_pilot.root --sample-manifest config/kkgraviton_samples.yml --outdir plots/presentation/pilot_1TeV --weighting unweighted --bootstrap-resamples 200
python3 Validation/audit_raw_acceptance.py --filelist Validation/representative_raw_topology_filelist.txt --output Validation/raw_acceptance_representative_topology.json
python3 Validation/summarize_raw_acceptance.py Validation/raw_acceptance_representative_topology.json --outdir plots/presentation/raw_acceptance_multitag
```

The principal new pilot figure is `plots/presentation/pilot_1TeV/kk_574797_mg1000/physics/neutrino_causal_diagnostics.png`; the parent-side response figure is `.../physics/topside_response.png`. Machine-readable neutrino, oracle, W/top, and differential response rows are in `plots/presentation/pilot_1TeV/summary/neutrino_metrics.csv` and `neutrino_differential.csv`, with the full study record in `summary/study_metrics.json`. The assertion-backed pilot check is `Validation/pilots/kk574797_causal_validation.json`. Raw cutflow tables and figures are in `plots/presentation/raw_acceptance_multitag/`; the full input report, identity audit, and current rich-output status are `Validation/raw_acceptance_representative_topology.json`, `Validation/representative_raw_topology_identity_audit.json`, and `Validation/full_rich_output_status.json`.

## Branch-method comparison on representative MC20a slices (2026-10-08)

This section supersedes earlier statements that representative full-method and SM-control response comparisons were unavailable. It uses the exact 19 non-overlapping raw ROOT inputs in `plots/presentation/sm_bsm_bias/summary/presentation_metrics.json` (506,859 raw rows), with the established strict population unchanged: direct after-FSR W-to-e/mu on both parent sides, analysis selected, truth reconstructible, and VyPER candidate-supported. The performance metrics are unweighted. The seven inclusive samples remain SM plus KK 400, 500, 1000, 2000, 3000, and 5000 GeV; no filtered control sample is included in those inclusive points.

The audited route definitions and limitations are in `Validation/reconstruction_method_inventory.md`. Native EM and NW were rebuilt from their stored constituents using the validated mixed-unit conversion (b p4 in MeV; lepton and neutrino p4 in GeV multiplied by 1000). W-mass, MET-balance, parent ordering, and assignment-index checks pass. DA indices are usable for assignment-only metrics, but DA complete top p4 is excluded because its W-mass constraint fails. No local EM or NW solver accepts caller-supplied b indices, so `VYPER_EM` and `VYPER_NW` could not be implemented honestly. The fixed-assignment neutrino comparison therefore contains the validated stored VyPER, NuFlows, and truth-neutrino routes; it does not substitute native EM/NW neutrinos for a controlled hybrid.

### Assignment efficiency

For each method, the efficiency is `CORRECT / method-valid` among truth-reconstructible events; method validity/support is shown separately. The strict VyPER-supported population is retained as the event input. Parent swaps remain a separate category because they can leave total mtt nearly unchanged.

| Sample | VyPER | HyPER | Ellipse assignment | NW assignment | DA assignment |
|---|---:|---:|---:|---:|---:|
| SM | 84.0% (100% valid) | 83.9% (100%) | 77.1% (72.4%) | 47.9% (95.2%) | 83.8% (100%) |
| KK 400 | 77.0% (100%) | 76.8% (100%) | 70.7% (72.5%) | 51.5% (95.5%) | 76.5% (100%) |
| KK 500 | 89.2% (100%) | 89.2% (100%) | 84.9% (70.1%) | 48.9% (97.5%) | 88.6% (100%) |
| KK 1000 | 96.2% (100%) | 96.1% (100%) | 96.6% (65.6%) | 59.1% (97.8%) | 95.4% (100%) |
| KK 2000 | 96.0% (100%) | 95.9% (100%) | 97.7% (53.7%) | 62.8% (96.3%) | 93.4% (100%) |
| KK 3000 | 94.8% (100%) | 95.2% (100%) | 97.8% (46.2%) | 61.0% (94.8%) | 89.9% (100%) |
| KK 5000 | 95.0% (100%) | 94.9% (100%) | 97.4% (50.7%) | 61.0% (95.4%) | 89.9% (100%) |

VyPER's conditional assignment accuracy is high in the supported high-mass population, while the pre-reconstruction candidate/reconstructibility stages lose events as the topology becomes more boosted. This distinction prevents the assignment percentage from being misread as an inclusive selection efficiency. At 3 TeV the VyPER truth-reconstructible category fractions are 94.8% correct, 0.1% swapped, 5.0% one-wrong, 0.2% both-wrong, and 0.0% invalid; unsupported is reported separately in the machine-readable rows. Full categories and Wilson intervals are in `summary/method_study_metrics.json` and `summary/assignment_metrics.csv`.

### Complete methods and controlled neutrino routes

The complete-method figure compares methods on the same valid events within each sample. At 5 TeV that shared subset has 2,440 events (50.3% of the strict candidate-supported denominator):

| Route | Median relative bias | Relative sigma68 | `P(|r|>0.2)` |
|---|---:|---:|---:|
| SM-trained VyPER | -10.8% | 38.8% | 38.0% |
| Native Ellipse Method | +2.7% | 41.9% | 33.0% |
| Native Neutrino Weighting | +13.9% | 79.3% | 60.2% |
| VyPER + NuFlows | -13.9% | 37.9% | 41.9% |
| Reco-object oracle | -1.8% | 8.2% | 1.6% |

The matched-valid population is reported per sample in the JSON/CSV; this table is not a comparison of unequal method-valid subsets. The native Ellipse Method has a smaller high-mass bias than VyPER but a similar sigma68. Thus these branches do not show that VyPER uniquely broadens relative to the conventional ellipse route. The native NW route performs substantially worse in this slice and has different assignment semantics/support. It is not evidence about `VYPER_NW`.

The fixed-VyPER-assignment oracle ladder uses the same visible assignment and strict event set. At KK 5000 GeV, native VyPER has bias -9.67% and sigma68 40.85%; replacing only the neutrinos with parent truth neutrinos gives -1.59% and 10.19%; additionally supplying the truth-matched reco b assignment gives -1.57% and 9.91%. The 30.7 percentage-point native-to-truth-neutrino sigma68 reduction is much larger than the 0.3-point truth-neutrino-to-reco-object-oracle gap. The same ordering is present at SM and the other representative KK masses. This is the high-stat causal result: the dominant conditional loss is in neutrino reconstruction, not the visible assignment among supported events. It does not yet distinguish the learned VyPER neutrino head from a physics-constrained ellipse solution with the assignment fixed; that requires the missing callable ellipse solver.

### Matched domains and additional observables (pre-grid representative-only results; superseded below)

Using only the representative SM slice, matched mtt response requires at least 50 finite events in each domain and currently ends at truth mtt 1.6 TeV. Across the common bins, VyPER KK-minus-SM median-bias differences are about +1 to +7 percentage points; conditional sigma68 differences are compatible with zero for VyPER, while native EM/NW have method-dependent shifts. The fixed-assignment truth-neutrino and reco-object-oracle response differences remain around 1–3 percentage points in median bias and near zero in sigma68. These are conditional responses inside common support, not evidence for extrapolated SM behavior above 1.6 TeV.

The covariate comparison and minimum added `(mtt, truth top pT)` response are saved separately. Inside common mtt support, the largest residual phase-space difference is truth ttbar pT: the standardized median shift reaches 4.4 SM sigma68 units in a 3 TeV KK truth-mtt bin. Truth top pT reaches 1.2 units, while min-DeltaR(b,l) reaches about 0.5. The first two-dimensional response control uses `(mtt, truth top pT)`; it has 22 common cells for KK 1 TeV, 5 for 3 TeV, and 7 for 5 TeV under the 50-event-per-domain requirement. Median-bias offsets remain in some cells, but the high-mass cells are sparse and do not establish a coherent conditional broadening. This control is not a reweighting. Truth ttbar pT is the strongest candidate for the next conditioning variable once the raw arrays can be reread. The first reconstructed observable promoted beyond mtt is top pT, alongside antitop pT and charge-paired m(b,l). All use shared method populations and the same method colors as the headline plots.

The matched assignment plot does show a residual KK-versus-SM difference inside some truth-mtt bins. In 1100–1300 GeV, VyPER is correct on 98.8% of SM events (260 valid) versus 92.1% at KK 3 TeV (529 valid), a -6.8 percentage-point difference with the plotted independent Wilson 68% envelope at about [-8.6,-4.8] points. The corresponding differences are -5.0 points at KK 2 TeV and -4.0 points at KK 5 TeV. These are assignment efficiencies conditional on the unchanged strict supported population and matched only in mtt. Because truth ttbar pT and other topology distributions still differ within mtt bins, the result is not evidence by itself for a process-specific assignment failure; assignment comparison after a second kinematic condition needs the event-level arrays and stronger SM support.

The nominal MC20a DSID 410472 audit counts 10,988 raw direct-dilepton events above 1.6 TeV, 2,292 above 2 TeV, 72 above 3 TeV, and none above 5 TeV. A filtered method extraction from the non-overlapping high-mass part of that sample was attempted. The initial `Path.resolve()` metadata lookup was removed, but ROOT-file enumeration and an exact-file `stat` both time out in the mounted Ceph metadata server; no high-mass-control method metrics were produced. A 765-path control list has now been derived from the local catalog after removing the presentation input overlap (`Validation/high_mass_sm_control_inputs_r13167.txt`), so recovery no longer requires directory enumeration; the list does not validate current file readability. Until Ceph read access recovers and the extraction is rerun, the measured matched-response endpoint remains 1.6 TeV. This sample is the right next same-campaign control into approximately 2–2.5 TeV, but it cannot control the 3–6 TeV KK region. The method-bearing DSID 601491 `r14861` production is catalogued as MC20e fast simulation; its pTHard variation has nontrivial weights and is not substituted for the nominal denominator. Separate MC23 full-simulation DSID 601491 TtNT paths exist in the local catalog but lack branch and support validation. The full 17-mass KK sample expansion also depends on ROOT files under the same mounted path; only the six representative KK masses are measured in this report.

The initial plots are under `plots/presentation/sm_bsm_bias/{headline,assignment,neutrino_methods,matched_response,observables,spectra,matrices,diagnostics}/`. Single-sample event uncertainties use 200 bootstrap resamples; the final matched-domain comparisons use 1,000. Regenerate figures without reopening ROOT files with:

```bash
python3 -m Plot.reco_plot.method_study \
  --render-summary-json plots/presentation/sm_bsm_bias/summary/method_study_metrics.json \
  --outdir plots/presentation/sm_bsm_bias
```

## Follow-up execution status (2026-10-08, superseded by the 2026-10-09 execution below)

The representative results above remain the current measured method comparison. The nominal SM raw support audit establishes that the available MC20a sample can only be a partial matched control: raw direct-dilepton support falls from 2,292 events above 2 TeV to 72 above 3 TeV and zero above 5 TeV. The strict selected/reconstructible/candidate-supported counts are lower. The full filtered method extraction did not complete: a serial pass reached 50 of 784 inputs in 13m39s and was stopped before summary writing; the bounded-worker retry blocked in Ceph metadata before opening any inputs. The validated one-file SM pilot remains a branch/unit/ordering check and is not promoted to a high-stat result. The existing matched-response endpoint remains 1.6 TeV, and no BSM-specific high-energy response conclusion is supported.

The reconstruction-method inventory shows native VyPER, complete native Ellipse and Neutrino Weighting branches, NuFlows neutrinos, HyPER and DA assignments, truth assignment, parent truth neutrinos, partonic truth, and the reconstructed-object oracle. No available callable ellipse or NW solver accepts externally supplied VyPER jet indices. `VYPER_EM` and `VYPER_NW` remain unavailable without new solver implementation; no unvalidated solver or assignment fallback was introduced. The fixed-assignment truth-neutrino comparison remains the causal result: at 5 TeV, native VyPER has relative `mtt` bias/sigma68 of about -9.7%/40.9%, VyPER plus truth neutrinos -1.6%/10.2%, and the reco-object oracle -1.6%/9.9%. This supports neutrino reconstruction as the dominant high-mass contribution on supported events, while the assignment-only oracle gap is small. It does not distinguish a VyPER-specific neutrino extrapolation from generic high-energy extrapolation.

An explicit 30-file list, `Validation/kk_mass_grid_mc20a_representative_inputs.txt`, adds one exact `r13167` file for each of the eleven KK masses absent from the current six-mass presentation study. The local metadata filelist confirms the DSID-to-mass inputs; the new ROOT files still need schema, event-count, and identity validation when Ceph read access returns. The study reader now supports exact path-list files and bounded input workers, preserving the existing strict denominator and deterministic input order. Compile checks and all 60 tests pass; the parallel reader has not had a successful real-input run.

Current interpretation: **observation** — representative samples show VyPER broadening at high KK mass and the same-assignment truth-neutrino substitution removes most of the 5 TeV width; **interpretation** — neutrino reconstruction dominates over visible b assignment on the supported population; **hypothesis** — the effect may be generic high-energy phase-space extrapolation; **conclusion not established** — nominal SM control lacks supported common bins above 1.6 TeV, and no nominal SM sample has useful support through 5–6 TeV. The highest-information next action is restoring Ceph metadata/read access, validating the eleven new mass-grid ROOT files, then extracting the remaining nominal SM `r13167` control and extending the matched comparison into its actual selected support.

## Controlled study execution (2026-10-09)

This section supersedes the earlier access, coverage, matched-domain, and mass-grid status. On the healthy node, one established ROOT input opened normally and its `reco` tree and required branches read with ordinary latency. The full analysis then read all 30 paths in `Validation/kk_mass_grid_mc20a_representative_inputs.txt` and all 765 paths in `Validation/high_mass_sm_control_inputs_r13167.txt` (795 files total). Existing representative files were reread because the new cache needed event identity, parent neutrino four-vectors, weights, method-validity flags, assignments, and reconstructed four-vectors for the controlled comparisons. The extraction completed and wrote the summary and caches.

Reproduce the extraction from the exact path lists with:

```bash
python3 -u -m Plot.reco_plot.method_study \
  --input-filelist Validation/kk_mass_grid_mc20a_representative_inputs.txt \
  --high-mass-sm-control /data/leet/ntuples/run2_ttbar_dilep/user.leet.410472.PhPy8EG.DAOD_PHYS.e6348_s3681_r13167_p7266.TtNT_2L_2Lsys_b71171e_output \
  --high-mass-sm-control-filelist Validation/high_mass_sm_control_inputs_r13167.txt \
  --input-workers 6 --bootstrap-resamples 200 \
  --event-cache-dir plots/presentation/sm_bsm_bias/cache \
  --outdir plots/presentation/sm_bsm_bias
```

The strict population remains selected `ttll_analysis`, direct after-FSR W→e/μ on both sides, truth reconstructible, and candidate supported. Reported event-level performance is unweighted; available nominal weight branches are preserved in the cache but are not used for these shape/performance comparisons. The 17 KK inputs are representative file slices, not complete campaigns. Each added mass input passed the reader's required-branch, identity, validity, and assignment checks; the first added 750 GeV file was separately inspected as the pilot. The cache contains 19 strict-population NPZs (SM representative, high-mass SM control, and 17 KK points), with JSON sidecars mapping event identities to source files and documenting the assignment-category encoding. Cache creation code is in `Plot/reco_plot/method_study.py` and `Plot/reco_plot/bias_diagnostics.py`; current repository base is `5be71643f0293ee6aabd6088f113d52694661abd`, with the study changes in the working tree. Later matching and plotting use the saved event arrays and do not reopen ROOT files.

Repeat only the matched-domain analysis and figures from these caches, with no ROOT read:

```bash
python3 -m Validation.recompute_cached_domain_matches \
  --cache-dir plots/presentation/sm_bsm_bias/cache \
  --summary-json plots/presentation/sm_bsm_bias/summary/method_study_metrics.json \
  --outdir plots/presentation/sm_bsm_bias \
  --bootstrap-resamples 1000
```

### Coverage and effective performance

Coverage is (C=N(\mathrm{valid}\land\mathrm{strict})/N(\mathrm{strict})). Conditional assignment accuracy is (A=N(\mathrm{valid}\land\mathrm{correct})/N(\mathrm{valid}\land\mathrm{strict})). Effective correct assignment is (E_{\rm corr}=N(\mathrm{valid}\land\mathrm{correct})/N(\mathrm{strict})). Complete-method (E_{20}) is the strict-population fraction with a valid result and |Δmtt/mtt|<0.2. These denominators are stored separately in `summary/coverage_effective_metrics.csv`.

At 5 TeV (strict (N=4{,}854)):

| Route | Coverage | Conditional assignment accuracy | Effective correct assignment | Effective (E_{20}) |
|---|---:|---:|---:|---:|
| VyPER assignment / native | 100.0% | 95.0% | 95.0% | 60.7% |
| Ellipse assignment / complete | 50.7% | 97.4% | 49.4% | 34.0% |
| NW assignment / complete | 95.4% | 61.0% | 58.2% | 37.3% |
| VyPER + NuFlows | 100.0% | 95.0% | 95.0% | 56.4% |
| Reco-object oracle | 100.0% | — | — | 95.0% |

The Ellipse conditional assignment accuracy is high, but only about half the strict events have a valid Ellipse solution. This validity loss is a material end-to-end cost. Across the full grid, Ellipse coverage decreases from about 70% at 350–500 GeV to 46% at 3–3.5 TeV, then rises to 54% at 6 TeV. VyPER validity remains 100%; NW validity is about 95–99%. Assignment-only HyPER and DA remain assignment routes and are not assigned complete-method (E_{20}) values.

### Ellipse-validity fairness and fixed-assignment diagnostic

The Ellipse-validity comparison uses identical event subsets for both methods. At 5 TeV, VyPER mtt bias/sigma68 is −9.67%/40.85% on all 4,854 strict VyPER-valid events, −10.81%/38.73% on the 2,462 Ellipse-valid events (50.7% of strict), and −10.80%/38.85% on the 2,440 common-valid complete-method events. Ellipse on those latter subsets is +2.75%/41.68% and +2.75%/41.95%, respectively. At 3 TeV, restricting to Ellipse-valid events shifts VyPER from −11.26% on all strict events to −13.89% on the 4,084-event subset, while Ellipse is +1.06%. Ellipse selection therefore does not explain its centring relative to VyPER; on these same events VyPER is slightly more negative.

The same-assignment diagnostic requires exact equality of the ordered top/antitop jet indices. At 5 TeV, 2,418/4,854 strict events (49.8%) share the assignment; 2,388/4,854 (49.2%) have both methods assignment-correct. On the same-assignment subset, native VyPER gives mtt bias/sigma68 −10.73%/38.53%, while native Ellipse gives +2.81%/41.33%. Replacing VyPER neutrinos with parent truth neutrinos gives −1.77%/8.28%; the reco-object oracle gives −1.78%/8.20%. The both-correct subset gives essentially the same comparison. This is a diagnostic subset, not a full-population `VYPER_EM` result.

The same exact subset also shows top-pT bias/sigma68 of −21.24%/62.27% for VyPER, +2.78%/43.92% for Ellipse, −1.30%/9.65% for VyPER with truth neutrinos, and −1.29%/9.62% for the oracle. For reconstructed top mass, the corresponding values are +2.54%/224.1%, −0.47%/5.02%, −1.62%/14.83%, and −1.62%/14.69%. The top-mass width is particularly broad for native VyPER. The sample-by-sample bias, width, tail fraction, count, and strict-population fraction for mtt, top pT, and top mass are stored in `summary/same_assignment_diagnostic.csv`; the new top-pT/top-mass comparison is in `neutrino_methods/same_assignment_toppt_topmass.png`.

The charge-paired visible-side control remains narrower on the common-valid complete-method population: at 5 TeV the top-side m(b,ℓ) sigma68 is 32.5% for VyPER, 32.7% for Ellipse, and 32.1% for the reco-object oracle. This contrasts with the large same-assignment top-pT and mtt response differences and is consistent with a neutrino-sensitive failure rather than a comparable visible-side resolution loss.

No callable Ellipse solver accepting externally supplied VyPER assignments was recovered in the checked repository, adjacent analysis exports, or installed solver dependencies. Stored Ellipse branches are native-assignment outputs and cannot reproduce `VYPER_EM`. The production solver's units, ordering, ranking, and no-solution behavior therefore could not be reproduction-tested; `VYPER_EM` remains unavailable. The same-assignment result independently shows a neutrino/global-solution difference, but does not replace that hybrid.

### Neutrino response and transverse decomposition

On the 5 TeV same-assignment subset, VyPER median reconstructed/truth responses are 0.687 (top-neutrino energy), 0.663 (top-neutrino pT), 0.758 (antitop-neutrino energy), 0.748 (antitop-neutrino pT), and 0.736 for the summed pair energy. Ellipse's corresponding energy/pT responses are about 1.15/1.14 on both sides, with pair-energy response 1.140. The high-scale VyPER deficit is consistent with neutrino-kinematic compression, but this is an observed method-response pattern; without `VYPER_EM` and the exact model/training provenance it does not isolate a neural-head causal mechanism.

The (S/D) decomposition uses the 4,611 strict events with a correct VyPER assignment at 5 TeV. For mtt, native gives −9.55%/40.31% bias/sigma68; truth (S) with native (D) gives −8.45%/39.14%; native (S) with truth (D) gives −2.19%/21.86%; truth transverse pair with native pz gives −2.25%/20.86%; full truth neutrinos give −1.59%/9.72%. Top-pT bias/sigma68 changes from −21.04%/64.10% (native), to −20.90%/63.60% (truth (S) only), to −2.04%/18.56% (truth (D) only), then −1.16%/10.80% with the truth transverse pair. The transverse partition (D) is the dominant missing component in this oracle decomposition; truth pz further improves mtt width.

### Assignment independence and full KK grid

VyPER and HyPER have the same ordered assignment in 4,764/4,854 strict 5 TeV events (98.15%). Of the 90 disagreements, VyPER is correct while HyPER is wrong in 42 events, HyPER is correct while VyPER is wrong in 39, and both are wrong differently in 9. A further 195 events are both wrong with the same ordered assignment; they are included in the agreement count. This high agreement means HyPER is not an independent headline baseline in this slice.

The inclusive strict-population VyPER mtt trend and end-to-end coverage summary are:

| (M_G) [GeV] | Strict (N) | VyPER bias [%] | σ₆₈ [%] | |Δm/m|>20% [%] | Ellipse coverage [%] | VyPER (E_{20}) [%] | Oracle (E_{20}) [%] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 350 | 14,047 | −2.23 | 32.29 | 23.56 | 70.35 | 76.44 | 96.55 |
| 400 | 49,461 | −0.66 | 30.10 | 19.38 | 72.52 | 80.62 | 96.66 |
| 500 | 62,088 | −2.72 | 30.89 | 21.92 | 70.11 | 78.08 | 96.41 |
| 750 | 19,538 | −4.48 | 32.58 | 25.83 | 67.93 | 74.17 | 96.79 |
| 850 | 19,576 | −4.96 | 33.13 | 26.99 | 67.49 | 73.01 | 96.70 |
| 1,000 | 75,147 | −5.79 | 34.14 | 28.65 | 65.58 | 71.35 | 97.01 |
| 1,250 | 15,942 | −6.83 | 35.26 | 31.39 | 62.75 | 68.61 | 97.10 |
| 1,500 | 12,162 | −8.12 | 37.16 | 34.25 | 60.55 | 65.75 | 97.09 |
| 1,750 | 8,853 | −9.13 | 39.10 | 37.19 | 57.31 | 62.81 | 96.93 |
| 2,000 | 25,896 | −10.69 | 40.17 | 39.78 | 53.74 | 60.22 | 96.80 |
| 2,500 | 3,673 | −11.05 | 42.36 | 41.44 | 50.29 | 58.56 | 96.95 |
| 3,000 | 8,833 | −11.26 | 43.03 | 43.45 | 46.24 | 56.55 | 95.71 |
| 3,500 | 1,436 | −11.50 | 42.89 | 43.38 | 46.38 | 56.62 | 94.29 |
| 4,000 | 1,165 | −10.51 | 43.69 | 41.46 | 47.90 | 58.54 | 95.45 |
| 4,500 | 1,111 | −9.66 | 41.85 | 38.34 | 48.78 | 61.66 | 94.42 |
| 5,000 | 4,854 | −9.67 | 40.85 | 39.35 | 50.72 | 60.65 | 94.99 |
| 6,000 | 2,007 | −9.38 | 39.85 | 37.42 | 53.71 | 62.58 | 96.01 |

The inclusive mass trend is an **inclusive high-energy degradation**: bias and width worsen through roughly 2.5–3.5 TeV generated benchmarks and remain broad at 4–6 TeV. These are representative slices and do not carry physical-yield normalization. The all-mass values, method-specific denominators, and Wilson intervals are in the summary CSV/JSON.

### Nominal SM support and matched-domain controls

The 765-file nominal MC20a DSID 410472 high-mass control contributes 6,221 strict events above truth mtt 1.6 TeV, 1,239 above 2 TeV, 202 above 2.5 TeV, and 39 above 3 TeV; there are none above 4 TeV, and the strict maximum is 3.64 TeV. The representative SM and added control identities are disjoint. This extends the same-campaign reference into the 2–2.5 TeV range but does not provide a robust nominal SM denominator in the 3–6 TeV KK resonance region.

The mtt-only assignment comparison in truth mtt 1.1–1.3 TeV gives SM 98.85% (257/260), KK 2 TeV 93.93% (1,796/1,912), KK 3 TeV 92.06% (487/529), and KK 5 TeV 94.82% (531/560). Under one-dimensional truth-mtt matching, the effective correct-assignment differences at 2/3/5/6 TeV are −1.50/−3.21/−2.60/−2.44 percentage points. After common ((mtt,p_T^{tt})) cell matching they shrink to −0.27 [−0.51,−0.05], −1.28 [−1.79,−0.80], −0.40 [−0.86,+0.06], and −0.89 [−1.63,−0.20] pp (68% bootstrap intervals). Thus the larger mtt-only assignment gap is mostly explained by differing truth-pTtt composition; residual differences are at most about 1.3 pp in this supported comparison. At 6 TeV the effective KK support is about 861 events.

All final matched-domain intervals below use 1,000 independent weighted event-bootstrap resamples. Weights equalize each common bin/cell to the smaller domain count, have no clipping, and preserve the strict population. The one-dimensional match has 9 common mtt bins at 2/3/5 TeV and 8 at 6 TeV. Its KK−SM mtt bias differences are +0.90 [+0.65,+1.25], +3.46 [+2.91,+3.86], +2.59 [+2.00,+3.04], and +1.61 [+0.93,+2.19] pp. After adding truth pTtt matching, the mtt bias differences are −0.69 [−1.06,−0.36], −0.25 [−0.98,+0.48], +0.42 [−0.19,+0.97], and +0.05 [−0.96,+1.05] pp. The corresponding 2D sigma68 shifts are +0.82 [+0.06,+1.35], +2.59 [+1.50,+3.33], +0.47 [−0.54,+1.52], and −1.63 [−2.83,+0.33] pp. This indicates little residual mtt centring difference after 2D matching, with a small width increase resolved at 3 TeV.

The absolute 2D matched response is negative in both domains: at the 5 TeV benchmark, matched mtt median bias is −11.0% for SM and −10.6% for KK, while top-pT bias is −18.7% and −24.0%; at 3 TeV the corresponding values are −12.1%/−12.4% in mtt and −19.4%/−27.2% in top pT. Thus the high-mtt negative response is already present in nominal SM, while the additional KK top-pT shift is larger. The 2D top-pT bias differences are −5.81 [−6.45,−5.30], −7.78 [−8.90,−6.65], −5.30 [−6.29,−3.96], and −2.42 [−4.10,−0.78] pp at 2/3/5/6 TeV. Adding coarse truth neutrino-pair-energy cells [0,500,1000,2000,5000,10000] GeV reduces these to −2.57 [−3.24,−1.93], −3.70 [−4.97,−2.32], and −1.70 [−3.11,−0.92] pp at 2/3/5 TeV. At the same time the 3D mtt bias differences are +1.96 [+1.59,+2.30], +2.64 [+1.98,+3.34], and +2.20 [+1.31,+2.89] pp. Neutrino-pair-energy composition therefore accounts for much of the 2D top-pT gap, although a smaller response difference remains in the 3D cells. The 3D control retains 24/17/17 common cells at 2/3/5 TeV, with effective support SM/KK of about 6,546/10,347, 2,326/1,506, and 2,669/1,305. At 6 TeV only four 3D cells remain (effective support 541/245), so the apparent width differences are unstable and not interpreted. The 1D, 2D, and 3D cells all lie in common truth-mtt support below 3.5 TeV; this does not test the 5–6 TeV resonance tail. Raw support, effective counts, cell boundaries, masked regions, weights, and bootstrap intervals are in the matched CSVs and JSON.

The optional MC23 full-simulation pilot was attempted on one catalogued DSID 601491 mc23a file:
`/data/leet/ttbar2L/user.leet.601491.PhPy8EG.DAOD_PHYS.e8549_s4162_r14622_p6697.TtNT_2L_2Lsys_5b9ac5d_output/user.leet.48537881._000001.output.root`.
That exact metadata path returned `FileNotFoundError`; no schema, weight, method, or high-mtt claims are made for the candidate, and no wider catalog search was run. The cached MC20e fast-simulation pTHard variation remains a secondary control and is not substituted for the nominal SM denominator.

### Training provenance and interpretation

The checked repository/configuration and adjacent analysis exports do not identify the exact SM-trained VyPER checkpoint, config, training split, or training ntuple/H5 source. Therefore no numeric claim about what fraction of KK events lies outside the actual training support is made. A high-energy neutrino-response compression hypothesis is supported by the observed ratios and the (S/D) oracle decomposition, but the data do not isolate the learned neutrino head as the sole causal mechanism. The same-assignment comparison shows a sizeable neutrino/global-solution difference independent of visible assignment; Ellipse validity selection does not explain away that difference; and the (D) partition is the leading transverse component in the oracle intervention.

The strongest supported conclusion is: **VyPER's negative high-mass response is dominated by neutrino-side reconstruction, especially the transverse partition between the two neutrinos; the same negative mtt response is also present in matched nominal SM events, supporting a generic high-energy reconstruction component within common support.** The extra KK top-pT shift shrinks after controlling neutrino-pair energy but remains a few percentage points, so a residual domain effect is plausible below 3.5 TeV. A learned-neutrino-head causal claim is not established without `VYPER_EM` and training provenance. A genuine SM-to-BSM difference in the 3–6 TeV KK resonance region remains inconclusive because nominal-like SM support ends at 3.64 TeV. The single highest-information next input is a nominal-like SM full-simulation sample with the same method branches and useful strict support through truth mtt of at least 5 TeV; it would test whether the KK-tail response is generic or process-dependent.

### Final validation and rerendering

The cache-based domain comparison can be repeated without reading ROOT files:

```bash
python3 -m Validation.recompute_cached_domain_matches \
  --cache-dir plots/presentation/sm_bsm_bias/cache \
  --summary-json plots/presentation/sm_bsm_bias/summary/method_study_metrics.json \
  --outdir plots/presentation/sm_bsm_bias \
  --bootstrap-resamples 1000
```

The presentation figures can be rerendered from the saved summary with:

```bash
python3 -m Plot.reco_plot.method_study \
  --render-summary-json plots/presentation/sm_bsm_bias/summary/method_study_metrics.json \
  --outdir plots/presentation/sm_bsm_bias
```

Final validation passed: `python3 -m py_compile Plot/reco_plot/raw_presentation.py Plot/reco_plot/method_study.py Plot/reco_plot/bias_diagnostics.py Validation/recompute_cached_domain_matches.py Tests/test_reco_plot_bias_diagnostics.py` and `python3 -m unittest discover -s Tests -p 'test_*.py' -q` (72 tests). The cache-matching CLI completed a one-resample integration smoke run over all 19 caches. A cache audit confirmed unique event identities within every cache, no overlap between the SM representative and high-mass control, strict masks/provenance alignment, finite truth and valid-method response arrays, and the expected 17 KK masses. Final matched-domain tables use 1,000 weighted event-bootstrap resamples.
