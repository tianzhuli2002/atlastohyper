import h5py
import numpy as np
from tqdm import tqdm
import uproot
import awkward as ak

from Common import argparser, main
"""
python TtbarDilepton.py -f /eos/user/l/leet/Rebuild/output_ntuple_new/ttll_601230_mc23d_fullsim_Njobs_4_jobIndex_0.root -o tll_601230_mc23d_fullsim.h5 -ds 2000000 -b -s 
python TtbarDilepton.py -f /data/leet/FFTutorial/output_ntuples_run_3_comparison/ttll_601230_mc23d_fullsim_Njobs_100_jobIndex_10.root -o llbb_compare.h5 -b -fm -ic

"""

def MakeDatasetTtbarDilepton(args, loop_storage, output: str, split: bool = False, saveEven: bool = True):

    if loop_storage is not None:
        raise ValueError('LoopStorage is not used in this converter.')

    includeBTag  = args.btag
    fullyMatch   = args.fully_matched
    dataset_size = args.dataset_size
    includerandom = args.include_random
    includeClassification = args.include_classification

    pad_to_jet = 20
    pad_to_lepton = 3

    input = args.input
    root_file = uproot.open(input)
    tree = root_file['reco']
    num_entries = tree.num_entries

    # -------------------------
    # Load data from ROOT
    # -------------------------

    jet_bTag  = tree["Jet_GN2v01_PCBT_quantile_NOSYS"].array(library='np')
    jet_eta   = tree["Jet_eta"].array(library='np')
    jet_phi   = tree["Jet_phi"].array(library='np')
    jet_e     = tree["Jet_e"].array(library='np')
    jet_pt    = tree["Jet_pt"].array(library='np')

    el_eta = tree["Electron_eta"].array(library='np')
    el_phi = tree["Electron_phi"].array(library='np')
    el_e   = tree["Electron_e"].array(library='np')
    el_pt  = tree["Electron_pt"].array(library='np')
    el_q   = tree["Electron_q"].array(library='np')

    mu_eta = tree["Muon_eta"].array(library='np')
    mu_phi = tree["Muon_phi"].array(library='np')
    mu_e   = tree["Muon_e"].array(library='np')
    mu_pt  = tree["Muon_pt"].array(library='np')
    mu_q   = tree["Muon_q"].array(library='np')

    # ---- MET (only for GLOBAL) ----
    met_phi = tree["Met_phi"].array(library='np')
    met_pt  = tree["Met_pt"].array(library='np')

    lep_eta = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_eta, mu_eta)]
    lep_phi = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_phi, mu_phi)]
    lep_e   = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_e, mu_e)]
    lep_pt  = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_pt, mu_pt)]
    lep_q   = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_q, mu_q)]

    isTtbarDilepton = tree["isTtbarDilepton_NOSYS"].array(library='np').astype(np.int32)
#THIS IS HARD CODING IT TO BE ONE PLEASE REMOVE LATER
    if args.force_ttbar_label:
        isTtbarDilepton = np.ones(num_entries, dtype=np.int32)
    else:
        isTtbarDilepton = tree["isTtbarDilepton_NOSYS"].array(library='np').astype(np.int32)
    
    event_number = tree["eventNumber"].array(library='np')
    njet = (ak.count(tree["Jet_pt"].array(), axis=1).to_numpy()).reshape(-1,1)
    nlep = (
        ak.count(tree["Electron_pt"].array(),axis=1).to_numpy().reshape(-1,1) +
        ak.count(tree["Muon_pt"].array(),axis=1).to_numpy().reshape(-1,1)
    )

    nbTagged_90 = tree["nBjets_GN2v01_90_NOSYS"].array(library='np')
    nbTagged_85 = tree["nBjets_GN2v01_85_NOSYS"].array(library='np')
    nbTagged_77 = tree["nBjets_GN2v01_77_NOSYS"].array(library='np')
    nbTagged_70 = tree["nBjets_GN2v01_70_NOSYS"].array(library='np')
    nbTagged_65 = tree["nBjets_GN2v01_65_NOSYS"].array(library='np')

    bhad_index    = tree["Truth_had_b_index"].array(library='np')
    bbarhad_index = tree["Truth_had_bbar_index"].array(library='np')

    # -------------------------
    # Count events
    # -------------------------

    num_events_to_save = num_entries if not split else sum(
        (event_number[i] % 2 == 0) == saveEven for i in range(num_entries)
    )

    if dataset_size is not None:
        num_events_to_save = min(num_events_to_save, dataset_size)

    # -------------------------
    # Define datatypes
    # -------------------------

    if includeBTag:
        node_dt = np.dtype([
            ('e', np.float32), ('eta', np.float32),
            ('phi', np.float32), ('pt', np.float32),
            ('btag', np.float32), ('charge', np.float32),
            ('id', np.float32)
        ])
    else:
        node_dt = np.dtype([
            ('e', np.float32), ('eta', np.float32),
            ('phi', np.float32), ('pt', np.float32),
            ('charge', np.float32), ('id', np.float32)
        ])

    jet_data = np.full((num_events_to_save, pad_to_jet), np.nan, dtype=node_dt)
    lep_data = np.full((num_events_to_save, pad_to_lepton), np.nan, dtype=node_dt)

    if includeBTag:
        global_dt = np.dtype([
            ('njet', np.float32),
            ('met_pt', np.float32),
            ('met_phi', np.float32),
            ('nBJets_90', np.float32),
            ('nBJets_85', np.float32),
            ('nBJets_77', np.float32),
            ('nBJets_70', np.float32),
            ('nBJets_65', np.float32)
        ])
    else:
        global_dt = np.dtype([
            ('njet', np.float32),
            ('met_pt', np.float32),
            ('met_phi', np.float32)
        ])

    global_data = np.zeros((num_events_to_save, 1), dtype=global_dt)

    jet_labels    = np.full((num_events_to_save, pad_to_jet), np.nan)
    lepton_labels = np.full((num_events_to_save, pad_to_lepton), np.nan)
    global_labels = np.full((num_events_to_save, 1), np.nan)

    # -------------------------
    # HDF5
    # -------------------------

    h5_file = h5py.File(output, 'w')
    inputs_group = h5_file.create_group('INPUTS')
    labels_group = h5_file.create_group('LABELS')

    save_index = 0

    for i in tqdm(range(num_entries), desc='Creating HDF5 dataset', unit='events'):

        if split and ((event_number[i] % 2 == 0) != saveEven):
            continue

        if save_index >= num_events_to_save:
            break

        num_jets = njet[i][0]
        for j in range(num_jets):
            if includeBTag:
                jet_data[save_index][j] = (
                    jet_e[i][j], jet_eta[i][j], jet_phi[i][j],
                    jet_pt[i][j], jet_bTag[i][j], 0, 1
                )
            else:
                jet_data[save_index][j] = (
                    jet_e[i][j], jet_eta[i][j], jet_phi[i][j],
                    jet_pt[i][j], 0, 1
                )

        num_leps = nlep[i][0]
        for j in range(num_leps):
            if includeBTag:
                lep_data[save_index][j] = (
                    lep_e[i][j], lep_eta[i][j], lep_phi[i][j],
                    lep_pt[i][j], 0, lep_q[i][j], 2
                )
            else:
                lep_data[save_index][j] = (
                    lep_e[i][j], lep_eta[i][j], lep_phi[i][j],
                    lep_pt[i][j], lep_q[i][j], 2
                )

        # ---- GLOBAL (MET kept here only) ----
        if includeBTag:
            global_data[save_index] = (
                njet[i], met_pt[i], met_phi[i],
                nbTagged_90[i], nbTagged_85[i],
                nbTagged_77[i], nbTagged_70[i], nbTagged_65[i]
            )
        else:
            global_data[save_index] = (
                njet[i], met_pt[i], met_phi[i]
            )

        existent_jet_labels = np.array(
            [0]*num_jets + [np.nan]*(pad_to_jet-num_jets)
        )
        existent_lep_labels = np.array(
            [0]*num_leps + [np.nan]*(pad_to_lepton-num_leps)
        )

        if 0 <= bhad_index[i] < num_jets:
            existent_jet_labels[bhad_index[i]] = 1
        if 0 <= bbarhad_index[i] < num_jets:
            existent_jet_labels[bbarhad_index[i]] = 2
                #     existent_jet_labels[bbarhad_index[i]] = 2
        #existent_lep_labels[i] = 1
        for j in range(num_leps):
            charge = lep_q[i][j]       # Usually +1 or -1
        # Label the leptons
            if charge > 0:
                    # Positively charged lepton (ℓ⁺)
                existent_lep_labels[j] = 1
            else:
                    # Negatively charged lepton (ℓ⁻)
                existent_lep_labels[j] = 2
        
        
        jet_labels[save_index] = existent_jet_labels
        lepton_labels[save_index] = existent_lep_labels

        if includeClassification:
            global_labels[save_index][0] = isTtbarDilepton[i]

        save_index += 1

    if fullyMatch:
        required_jet_ids = {1,2}
        keep_indices = [
            i for i in range(len(jet_labels))
            if required_jet_ids.issubset(
                set(np.nan_to_num(jet_labels[i], nan=-1))
            )
        ]

        jet_data = jet_data[keep_indices]
        lep_data = lep_data[keep_indices]
        global_data = global_data[keep_indices]
        jet_labels = jet_labels[keep_indices]
        lepton_labels = lepton_labels[keep_indices]
        global_labels = global_labels[keep_indices]

    inputs_group.create_dataset("JET", data=jet_data)
    inputs_group.create_dataset("LEPTON", data=lep_data)
    inputs_group.create_dataset("GLOBAL", data=global_data)

    labels_group.create_dataset("JET", data=jet_labels)
    labels_group.create_dataset("LEPTON", data=lepton_labels)
    labels_group.create_dataset("GLOBAL", data=global_labels)

    h5_file.close()
    root_file.close()

if __name__ == "__main__":

    parser = argparser()

    parser.add_argument('-ds', '--dataset-size', type=int, help='Number of events to save', required=False)
    parser.add_argument('-fm', '--fully-matched', action='store_true', help='Keep only fully matched events', required=False)
    parser.add_argument('-b', '--btag', action='store_true', help='Have btag information in the h5 file', required=False)
    parser.add_argument('-ir', '--include-random', action='store_true', help='Include random node in the h5 file', required=False)
    parser.add_argument('-ic', '--include-classification', action='store_true', help='Include classification labels in the h5 file', required=False)
    parser.add_argument(
    '--force-ttbar-label',
    action='store_true',
    help='Hardcode isTtbarDilepton label to 1 for all events',
    required=False
)
    args = parser.parse_args()    
   
    main(MakeDatasetTtbarDilepton, args)