# KK-graviton sample metadata audit

Date: 2026-10-05

## Scope and evidence

Audited `/home/leet/FF/metadata_kkgraviton/filelist.txt` and `sum_of_weights.txt`, local repository references, the supplied task/collection label, and the official ATLAS pyAMI client. The file list has four whitespace-separated fields: DSID, campaign, simulation type, and a local ROOT-file path. The sum-of-weights file has five fields: DSID, campaign, simulation type, systematic label, and sum of weights.

The supplied file list is mixed: it contains 49,868 ROOT paths across 103 DSIDs, including Powheg/Herwig/Pythia, Sherpa, MGPy8EG, and diboson-tagged entries. It is not a KK-only list. The 51 contiguous DSIDs 574758–574808 form an MGPy8EG-tagged, full-simulation group and are treated only as candidate KK-graviton signal points because the user-provided task and collection identify the samples as KK-graviton. This classification has not been independently confirmed from AMI. All 51 IDs occur in the file list; their 162 candidate ROOT paths exist on disk, and no duplicate path was found among those candidates. No mass is inferred from DSID order.

## Campaign and weight coverage

The candidate group has 150 DSID/campaign combinations in the list (48 DSIDs have mc20a/d/e; DSIDs 574759, 574761, and 574781 each have only two of the three campaigns). There are 64 matching sum-of-weights rows and 86 missing rows for campaign combinations present in the list. Therefore the provided weight file is insufficient for normalization of the candidate group as listed. Across the entire mixed list, 304 DSID/campaign/simulation combinations are represented but only 140 have sum-of-weights rows, leaving 164 missing. There are no duplicate sum-of-weights keys, malformed rows, non-finite values, or non-positive values. Positive finite values alone do not establish correctness or completeness against production bookkeeping.

`filelist.txt` has 49,868 parseable rows, no malformed rows, and no duplicate ROOT paths. Candidate file paths use local `/data/leet/ntuples/run2_ttbar_dilep/...` locations; all 162 candidate paths currently exist. This is a path-existence check only, not a ROOT-open or schema check.

## Authoritative metadata blocker

Local AMI access could not retrieve dataset/production metadata. After `setupATLAS` and `lsetup pyami`, the default client hit an expired TLS client certificate. The no-proxy query on `atlas-v2` and `atlas-replica-v2` returned `wrong role for command SearchQuery, expected role(s) [AMI_USER], found role(s) []`; `atlas-replica` reported the `atlas` account deactivated. No local authoritative record for these DSIDs was found in the repository metadata. The supplied files also do not provide cross sections or filter efficiencies.

Consequently, the manifest records the process/domain only as user-context labels and leaves `mass_GeV`, `width`, and `coupling_or_benchmark` null. The generator is recorded as the `MGPy8EG` path tag, while tune/settings remain unverified. Exact decay configuration, production benchmark, cross section, and filter efficiency remain unresolved. Do not make yield-normalised signal claims from this metadata alone.

## Reproduction and required follow-up

The local count checks can be rerun with the Python standard library by parsing the two paths above. With a valid ATLAS AMI client role/certificate, query every candidate DSID and each available campaign using pyAMI, e.g.:

```bash
source "$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh"
lsetup pyami
ami list datasets 'mc16_13TeV.574758.%' --fields physics_short,generator_name,physics_comment,cross_section,generator_filter_efficiency,events
ami show dataset info <exact-mc16_13TeV-LDN>
```

Repeat for DSIDs 574758–574808 and validate the associated EVNT/production request metadata before filling any mass, width, coupling, or benchmark values. Reconcile the 86 missing sum-of-weights rows and three campaign omissions before processing the full candidate list. A dedicated signal input list should be selected explicitly; the supplied mixed `filelist.txt` must not be treated as KK-only.
