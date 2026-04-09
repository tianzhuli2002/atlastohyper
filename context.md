**Purpose**: Comprehensive run & implementation notes for `TtbarDileptonFrameReco` (DiLepton HyPER integration).

- **Location**: /data/leet/FFTutorial/atlastohyper/Root/TtbarDileptonFrameReco.cc
- **Config**: /data/leet/FFTutorial/atlastohyper/FastFramesConfig/TtbarDileptonReco.yml

**Summary of recent edits**:
- Added HyPER performance exposure (indices, IDs, convergence, per-method mttbar outputs and correctness flags) inside `TtbarDileptonFrameReco` only — MainFrame was NOT modified.
- Implemented local, input-only fallbacks for optional HyPER / nuflows input branches so Snapshot/ntuple selection does not abort when inputs are missing.
- Removed creation of placeholder output-named columns in fallback stage to prevent duplicate `RDataFrame::Define` errors.
- Added presence-check guards before calling `systematicDefine` for `mttbar_*` and correctness flags to avoid re-defining already-existing columns.
- Exposed non-systematic ntuple aliases (e.g., `mttbar_HYPER` from `mttbar_HYPER_NOSYS`) using `Define` on NOSYS columns.

**Canonical input branch names used** (match `check/601230.h`):
- `hyper_TtbarDiLepton_Top1_Indices_NOSYS`
- `hyper_TtbarDiLepton_Top1_IDs_NOSYS`
- `hyper_TtbarDiLepton_Top2_Indices_NOSYS`
- `hyper_TtbarDiLepton_Top2_IDs_NOSYS`
- `nuflows_nu_out_NOSYS`
- `nuflows_loglik_NOSYS`

If your input files use different names, either rename or update the frame to match.

**What the code now does (behavioral overview)**:
- On entry to `defineVariables(...)` the frame inspects `mainNode.GetColumnNames()`.
- If HyPER/nuflows input branches are missing, the code will NOT run HyPER-dependent `systematicDefine` blocks — it will skip them and continue safely.
- The frame only defines fallbacks for OPTIONAL INPUT branches (vector-int/float placeholders). It does NOT create placeholders for output-named branches (e.g., `mttbar_HYPER_NOSYS`) during fallback.
- Before each `systematicDefine` call for `mttbar_*` and correctness flags, the code re-checks `GetColumnNames()` and only defines when the column is not already present.

**Key functions & helpers referenced**:
- `m_partonMatcher->systematicDefine(node, name, functor, {dependencies})` — systematics-aware wrapper used throughout.
- `TtbarDileptonVars::buildNuFromNuFlows(...)` — constructs neutrino TLVs from `nuflows_nu_out_NOSYS`.
- `TtbarDileptonVars::getAssignment(Method, ...)` — returns assignment + indices for a given method (HYPER/CHI2/PT/DELTAR).
- `TtbarDileptonVars::compute_mttbar(...)` — computes m(ttbar) from TLVs and neutrinos.
- `HyPERPerformanceDiLepton` / `HyPERPerformanceSingleLepton` helpers — truth-match indices, convergence, correctness.

**Build & install (recommended)**:
1. From repository root (or where your `build_custom` was configured):

```bash
setupATLAS 
asetup StatAnalysis,0.7.2
source build/setup.sh 
cmake --build build_custom -j4 --target install
source build_custom/setup.sh
```

If setupATLAS doesn't work, you can try the more manual approach:
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
If this works please add it to the documentation and remove the note about setupATLAS not working.

To simply recompile the custom frame 
```bash
cmake --build build_custom -j4 --target install
```

if already built/compiled then you can simply run 
```bash
source build_custom/setup.sh
```
1. Confirm install artifacts exist:

- `/data/leet/FFTutorial/install/lib/libFastFrames.so`
- `/data/leet/FFTutorial/install/bin/fast-frames.exe`

**Environment for running FastFrames (important)**:
- Ensure the installed libs are on `LD_LIBRARY_PATH` for Python wrapper or when running compiled binary:

```bash
export LD_LIBRARY_PATH=/data/leet/FFTutorial/install/lib:${LD_LIBRARY_PATH}
```

- Python wrapper may need same `LD_LIBRARY_PATH` exported in the shell used to run `python3`.
- If you see GLIBCXX symbol mismatch errors when importing compiled extensions in Python, that indicates a libstdc++ mismatch between the build toolchain and runtime interpreter. Build and run using the same toolchain/environment.
- If ONNX inference is required, ensure `libonnxruntime.so` (matching built-on runtime version) is available on `LD_LIBRARY_PATH`.

**Run commands**:


- Full-split 
```bash
export LD_LIBRARY_PATH=/data/leet/FFTutorial/install/lib:$LD_LIBRARY_PATH
python3 python/FastFrames.py -c /data/leet/FFTutorial/atlastohyper/FastFramesConfig/TtbarDileptonReco.yml --step n --split_n_jobs 500 --job_index 0
```

**Typical failure modes & diagnostics**:
- Duplicate Define error (RDataFrame::Define: cannot define column "X")
  - Cause: column named `X` was already created earlier (e.g., placeholder fallback or other frame).
  - How to debug: inspect logs for "Defining fallback <name>" messages; search repo for where `Define("X"` or `systematicDefine(... "X"` occurs.
  - Fix: prefer Redefine when you must replace an existing column, or ensure only one place defines that output. The code in `TtbarDileptonFrameReco` now guards defines with existence checks.

- Missing HyPER/nuflows input columns lead to early Snapshot abort — resolved by local fallback & presence-guards. Confirm by checking logs for "Skipping HyPER-dependent defines" message.

- GLIBCXX symbol mismatch when importing extension into Python
  - Symptom: Python import failure with undefined symbol from libstdc++.
  - Fix: run Python with the same GCC/libstdc++ used to build the C++ modules; or build in the environment matching the cluster's Python.

- Missing `libonnxruntime` when running compiled binary
  - Symptom: runtime loader error complaining `libonnxruntime.so.X` not found.
  - Fix: install or provide `libonnxruntime` on `LD_LIBRARY_PATH` or disable ONNX inference codepaths.

**How to verify produced ntuple contents**:
- The ntuple YAML is at `/data/leet/FFTutorial/atlastohyper/FastFramesConfig/TtbarDileptonReco.yml` — inspect the `Variables`/`Snapshot` branch list for requested outputs such as `mttbar_HYPER_NOSYS`.
- After a successful run, open the output ROOT file and list tree branches:

```bash
root -l output_file.root
root[0] T->Print();
```

- Or use a small python snippet with uproot/awkward to list branches.

**Where to look in logs for HyPER activity**:
- Look for these log lines (messages present in the frame):
  - "Defining fallback <branch> (empty) in reco frame"
  - "Skipping HyPER-dependent defines: HyPER/nuflows branches not present in input"
  - "Adding variable: <name>_NOSYS" for each `systematicDefine`
  - "Adding ntuple branch: mttbar_HYPER" for the non-systematic alias defines

**Notes on design decisions taken**:
- All changes were confined to `TtbarDileptonFrameReco.cc` per user request — no edits to `MainFrame.cc`.
- Fallbacks are intentionally conservative: only provide placeholders for optional INPUT branches so Snapshot's column existence checks succeed without pre-creating output-named columns that will be defined properly later.
- For safety, the frame does a runtime check of `mainNode.GetColumnNames()` before defining each mttbar or correctness column — this prevents duplicate Define when another module already provided the column.

**Next steps & suggestions**:
- Confirm single-split run finishes without RDataFrame duplicate-Define error and inspect output ROOT for `mttbar_*` branches.
- If duplicate-Define still occurs, search the repo for any other `Define("mttbar_..."` or `systematicDefine(..., "mttbar_..."` lines and convert the earlier one to `Redefine` or add a guard.
- If you plan to run on many jobs (100 splits), ensure `LD_LIBRARY_PATH` and ONNX runtime availability on the compute nodes.

**Quick reference — important file locations**:
- Frame: /data/leet/FFTutorial/atlastohyper/Root/TtbarDileptonFrameReco.cc
- YAML: /data/leet/FFTutorial/atlastohyper/FastFramesConfig/TtbarDileptonReco.yml
- Installed libs: /data/leet/FFTutorial/install/lib/
- Binaries: /data/leet/FFTutorial/install/bin/fast-frames.exe

**Contact / testing notes**:
- If you share a failing log, include the first ~200 lines around the failure and the `GetColumnNames()` listing (the frame prints fallback defines). That will make diagnosing duplicate-Define vs. missing input trivial.

---
End of context.md
