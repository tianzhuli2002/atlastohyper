#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

export ATLAS_LOCAL_ROOT_BASE="${ATLAS_LOCAL_ROOT_BASE:-/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase}"
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh"
asetup StatAnalysis,0.8.0
source build_ff08/setup.sh
source build_custom/setup.sh

cmake --build build_ff08 -j2
cmake --install build_ff08 --prefix "$repo_root/install"

python3 Validation/validate_fastframes_configs.py \
  FastFramesConfig/TtbarDileptonvyper_KKGraviton.yml \
  FastFramesConfig/TtbarDileptonReco_VyPERStudy_KKGraviton.yml \
  --report Validation/pilots/kkgraviton_fastframes_config_parse.json

converter_python="${VYPER_CONVERTER_PYTHON:-/home/leet/miniforge3/bin/python3}"
validation_python="${VYPER_VALIDATION_PYTHON:-$converter_python}"
sm_raw="$(awk 'NR == 1 {print $4}' Validation/pilots/sm410472_mc20d_filelist.txt)"
kk_raw="$(awk 'NR == 1 {print $4}' Validation/pilots/kk574797_mc20a_filelist.txt)"
"$validation_python" Validation/preflight_dilepton_inputs.py \
  --metadata-dir /home/leet/FF/metadata_kkgraviton \
  --sm-dsid 410472 --sm-campaign mc20d --sm-file "$sm_raw" \
  --kk-dsid 574797 --kk-campaign mc20a --kk-file "$kk_raw" \
  --sample-events 500 --shape-only --skip-path-check \
  --json-report Validation/pilots/kk_candidate_preflight_reproduced.json
"$validation_python" Validation/audit_native_reco_units.py "$sm_raw" \
  --max-events 500 --report Validation/pilots/sm410472_native_unit_audit.json
"$validation_python" Validation/audit_native_reco_units.py "$kk_raw" \
  --max-events 500 --report Validation/pilots/kk574797_native_unit_audit.json

python3 /home/leet/FF/fastframes/python/FastFrames.py \
  -c Validation/pilots/TtbarDileptonvyper_SM_Pilot.yml \
  --step n --max_event 100 \
  2>&1 | tee Validation/pilots/sm410472_fastframes_reproduced.log

python3 /home/leet/FF/fastframes/python/FastFrames.py \
  -c Validation/pilots/TtbarDileptonvyper_KKGraviton_Pilot.yml \
  --step n --max_event 100 --split_n_jobs 1 --job_index 0 \
  2>&1 | tee Validation/pilots/kk574797_fastframes_vyper_reproduced.log

python3 /home/leet/FF/fastframes/python/FastFrames.py \
  -c Validation/pilots/TtbarDileptonReco_VyPERStudy_KKGraviton_Pilot.yml \
  --step n --max_event 100 --split_n_jobs 1 --job_index 0 \
  2>&1 | tee Validation/pilots/kk574797_fastframes_reco_reproduced.log

mkdir -p Validation/pilots/artifacts/root Validation/pilots/artifacts/h5

sm_root=/data/leet/FFTutorial/output_ntuples_sm_vyper_pilot/sm_pilot_410472_mc20d_410472_mc20d_fullsim.root
kk_vyper_root=/data/leet/FFTutorial/output_ntuples_kk574797_pilot_reproduced_vyper/kk574797_mc20a_pilot_vyper_574797_mc20a_fullsim_Njobs_1_jobIndex_0.root
kk_reco_root=/data/leet/FFTutorial/output_ntuples_kk574797_pilot_reproduced_reco/kk574797_mc20a_pilot_reco_574797_mc20a_fullsim_Njobs_1_jobIndex_0.root
cp "$sm_root" Validation/pilots/artifacts/root/sm410472_vyper_pilot.root
cp "$kk_vyper_root" Validation/pilots/artifacts/root/kk574797_vyper_pilot.root
cp "$kk_reco_root" Validation/pilots/artifacts/root/kk574797_reco_pilot.root

"$converter_python" H5Converters/TtbarDiLeptonVyPER.py \
  -f "$sm_root" -o Validation/pilots/artifacts/h5/sm410472_vyper_pilot.h5 -b -fm -s \
  2>&1 | tee Validation/pilots/sm410472_h5_convert_reproduced.log
"$converter_python" H5Converters/TtbarDiLeptonVyPER.py \
  -f "$kk_vyper_root" -o Validation/pilots/artifacts/h5/kk574797_vyper_pilot.h5 -b -fm -s \
  2>&1 | tee Validation/pilots/kk574797_h5_convert_reproduced.log

for split in even odd; do
  "$validation_python" Validation/compare_vyper_h5_contract.py \
    --current "Validation/pilots/artifacts/h5/sm410472_vyper_pilot_${split}.h5" \
    --reference "/data/leet/FFTutorial/output_h5_sm_vyper_pilot/sm410472_reference_${split}.h5" \
    --report "Validation/pilots/sm410472_artifact_h5_contract_${split}.json"
  "$validation_python" Validation/compare_vyper_h5_contract.py \
    --current "Validation/pilots/artifacts/h5/kk574797_vyper_pilot_${split}.h5" \
    --reference "/data/leet/FFTutorial/output_h5_kk_graviton_pilot/kkgrav_574797_reference_${split}.h5" \
    --report "Validation/pilots/kk574797_reproduced_h5_contract_${split}.json"
done

"$validation_python" Validation/validate_vyper_h5_truth_labels.py \
  --root "$sm_root" --raw "$sm_raw" \
  --h5-even Validation/pilots/artifacts/h5/sm410472_vyper_pilot_even.h5 \
  --h5-odd Validation/pilots/artifacts/h5/sm410472_vyper_pilot_odd.h5 \
  --report Validation/pilots/sm410472_artifact_neutrino_truth_validation.json
"$validation_python" Validation/validate_vyper_h5_truth_labels.py \
  --root "$kk_vyper_root" --raw "$kk_raw" \
  --h5-even Validation/pilots/artifacts/h5/kk574797_vyper_pilot_even.h5 \
  --h5-odd Validation/pilots/artifacts/h5/kk574797_vyper_pilot_odd.h5 \
  --report Validation/pilots/kk574797_reproduced_neutrino_truth_validation.json
"$validation_python" Validation/validate_native_branch_copy.py \
  --raw "$kk_raw" --fastframes "$kk_reco_root" \
  --report Validation/pilots/kk574797_reproduced_native_branch_copy_validation.json
"$validation_python" Validation/validate_vyper_truthnu.py "$kk_reco_root" \
  --report Validation/pilots/kk574797_vyper_truthnu_validation.json
"$validation_python" Validation/validate_nuflows_ordering.py "$kk_reco_root" \
  --max-events 100 --report Validation/pilots/kk574797_nuflows_ordering_validation.json
"$validation_python" Validation/validate_vyper_nuflows.py "$kk_reco_root" \
  --report Validation/pilots/kk574797_vyper_nuflows_validation.json
"$validation_python" Validation/summarize_reco_pilot.py "$kk_reco_root" \
  --json Validation/pilots/kk574797_reproduced_pilot_metrics.json \
  --csv Validation/pilots/kk574797_reproduced_pilot_metrics.csv
"$validation_python" Validation/inspect_vyper_h5.py \
  --report Validation/pilots/kk574797_reproduced_h5_inspection.json \
  Validation/pilots/artifacts/h5/kk574797_vyper_pilot_even.h5 \
  Validation/pilots/artifacts/h5/kk574797_vyper_pilot_odd.h5

python3 Plot/compare_reco.py vyper-ablation "$kk_reco_root" \
  --outdir plots --region ttll_analysis --truth-source partonic_after_fsr

cp Validation/pilots/kk574797_vyper_truthnu_validation.json \
  plots/diagnostics/kk574797_vyper_truthnu_validation.json
cp Validation/pilots/kk574797_nuflows_ordering_validation.json \
  plots/diagnostics/kk574797_nuflows_ordering_validation.json
cp Validation/pilots/kk574797_vyper_nuflows_validation.json \
  plots/diagnostics/kk574797_vyper_nuflows_validation.json
cp Validation/pilots/sm410472_native_unit_audit.json \
  plots/diagnostics/sm410472_native_unit_audit.json
cp Validation/pilots/kk574797_native_unit_audit.json \
  plots/diagnostics/kk574797_native_unit_audit.json
cp Validation/pilots/kk_candidate_preflight_reproduced.json \
  plots/diagnostics/kk_candidate_preflight_reproduced.json
