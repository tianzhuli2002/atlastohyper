import h5py
import numpy as np
from tqdm import tqdm
import uproot
import awkward as ak

from Common import argparser, main

def MakeDatasetTtbarSingleLeptonVyper(args, loop_storage, output: str, split: bool = False, saveEven: bool = True):
    r"""Function to construct HyPER dataset for ttbar single-lepton.
    Electrons and Muons are loaded and saved as separate datasets.

    Args:
        input (str): input ROOT file.
        output (str): output HyPER dataset in .h5 format.
        split (bool): save only events with odd/even event number.
        saveEven (bool): if True, will save events with even event number
    """

    if loop_storage != None:
        raise ValueError('LoopStorage is not used in this converter.')

    includeBTag  = args.btag
    fullyMatch   = args.fully_matched
    dataset_size = args.dataset_size

    # Max number of objects (padding)
    pad_to_jet = 20
    pad_to_lepton = 3
    pad_to_nu = 2  # Fixed: was missing in original

    # Load the ROOT file
    input = args.input
    root_file = uproot.open(input)
    tree = root_file['reco']
    num_entries = tree.num_entries

    # ---------------------------------------------------
    #                 Load data from ROOT
    # ---------------------------------------------------
    # Jets features
    jet_bTag  = tree["Jet_GN2v01_PCBT_quantile_NOSYS"].array(library='np')
    jet_eta   = tree["Jet_eta"].array(library='np')
    jet_phi   = tree["Jet_phi"].array(library='np')
    jet_e     = tree["Jet_e"].array(library='np')
    jet_pt    = tree["Jet_pt"].array(library='np')

    # Leptons features (KEPT SEPARATE)
    el_eta   = tree["Electron_eta"].array(library='np')
    el_phi   = tree["Electron_phi"].array(library='np')
    el_e     = tree["Electron_e"].array(library='np')
    el_pt    = tree["Electron_pt"].array(library='np')
    el_q     = tree["Electron_q"].array(library='np')

    mu_eta   = tree["Muon_eta"].array(library='np')
    mu_phi   = tree["Muon_phi"].array(library='np')
    mu_e     = tree["Muon_e"].array(library='np')
    mu_pt    = tree["Muon_pt"].array(library='np')
    mu_q     = tree["Muon_q"].array(library='np')

    # Neutrino features
    nu_px  = tree["nu_px"].array(library='np')
    nu_py  = tree["nu_py"].array(library='np')
    nu_pz  = tree["nu_pz"].array(library='np")

    # Met features
    met_phi   = tree["Met_phi"].array(library='np')
    met_pt    = tree["Met_pt"].array(library='np')

    # Global features
    event_number = tree["eventNumber"].array(library='np')
    njet         = (ak.count(tree["Jet_pt"].array(), axis=1).to_numpy()).reshape(-1, 1)
    nbTagged_90  = tree["nBjets_GN2v01_90_NOSYS"].array(library='np')
    nbTagged_85  = tree["nBjets_GN2v01_85_NOSYS"].array(library='np')
    nbTagged_77  = tree["nBjets_GN2v01_77_NOSYS"].array(library='np')
    nbTagged_70  = tree["nBjets_GN2v01_70_NOSYS"].array(library='np')
    nbTagged_65  = tree["nBjets_GN2v01_65_NOSYS"].array(library='np')

    # Labels (from parton truth)
    up_index     = tree["Truth_up_index"].array(library='np')
    down_index   = tree["Truth_down_index"].array(library='np')
    bhad_index   = tree["Truth_had_b_index"].array(library='np')
    blep_index   = tree["Truth_lep_b_index"].array(library='np')

    # ---------------------------------------------------
    #               Count number of events to save
    # ---------------------------------------------------
    num_events_to_save = num_entries if not split else sum(
        (event_number[i] % 2 == 0) == saveEven for i in range(num_entries)
    )

    if dataset_size is not None:
        if num_events_to_save > dataset_size:
            num_events_to_save = dataset_size

    # ---------------------------------------------------
    #               Define Input Datatype
    # ---------------------------------------------------
    # Removed 'id' field to avoid np.nan -> int32 type errors.
    # Dataset names (ELECTRON/MUON) now distinguish the lepton type.
    if includeBTag:
        node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), 
                            ('pt', np.float32), ('btag', np.float32), ('charge', np.float32)])
    else:
        node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), 
                            ('pt', np.float32), ('charge', np.float32)])

    jet_data = np.full((num_events_to_save, pad_to_jet), np.nan, dtype=node_dt)
    el_data  = np.full((num_events_to_save, pad_to_lepton), np.nan, dtype=node_dt)
    mu_data  = np.full((num_events_to_save, pad_to_lepton), np.nan, dtype=node_dt)
    met_data = np.full((num_events_to_save, 1), np.nan, dtype=node_dt)

    nu_dt = np.dtype([('px', np.float32), ('py', np.float32), ('pz', np.float32)])
    nu_data  = np.full((num_events_to_save, pad_to_nu), np.nan, dtype=nu_dt)
    nu_labels = np.full((num_events_to_save, pad_to_nu), np.nan, dtype=nu_dt)

    # Global dtype
    if includeBTag:
        global_dt = np.dtype([('njet', np.float32), ('nBJets_90', np.float32), ('nBJets_85', np.float32), 
                              ('nBJets_77', np.float32), ('nBJets_70', np.float32), ('nBJets_65', np.float32)])
    else:
        global_dt = np.dtype([('njet', np.float32)])
    
    global_data = np.zeros((num_events_to_save, 1), dtype=global_dt)

    # Labels
    jet_labels      = np.full((num_events_to_save, pad_to_jet), np.nan)
    el_labels       = np.full((num_events_to_save, pad_to_lepton), np.nan)
    mu_labels       = np.full((num_events_to_save, pad_to_lepton), np.nan)
    met_labels      = np.full((num_events_to_save, 1), np.nan)
    global_labels   = np.zeros((num_events_to_save, 1), dtype=np.float32)

    # ---------------------------------------------------
    #               Define HDF5 Structure
    # ---------------------------------------------------
    h5_file = h5py.File(output, 'w')
    inputs_group = h5_file.create_group('INPUTS')
    labels_group = h5_file.create_group('LABELS')

    print(f'Saving {num_events_to_save} events')
    save_index = 0
    for i in tqdm(range(num_entries), desc='Creating HDF5 dataset', unit='events'):
        # Split logic
        if split and ((event_number[i] % 2 == 0) != saveEven):
            continue
        if save_index >= num_events_to_save:
            break

        # --- Fill Jets ---
        num_jets = int(njet[i][0])
        for j in range(min(num_jets, pad_to_jet)):
            if includeBTag:
                jet_data[save_index][j] = (jet_e[i][j], jet_eta[i][j], jet_phi[i][j], jet_pt[i][j], jet_bTag[i][j], 0)
            else:
                jet_data[save_index][j] = (jet_e[i][j], jet_eta[i][j], jet_phi[i][j], jet_pt[i][j], 0)

        # --- Fill Electrons ---
        num_el = el_pt[i].size
        for j in range(min(num_el, pad_to_lepton)):
            if includeBTag:
                el_data[save_index][j] = (el_e[i][j], el_eta[i][j], el_phi[i][j], el_pt[i][j], 0, el_q[i][j])
            else:
                el_data[save_index][j] = (el_e[i][j], el_eta[i][j], el_phi[i][j], el_pt[i][j], el_q[i][j])
            el_labels[save_index][j] = 1 if el_q[i][j] > 0 else 2

        # --- Fill Muons ---
        num_mu = mu_pt[i].size
        for j in range(min(num_mu, pad_to_lepton)):
            if includeBTag:
                mu_data[save_index][j] = (mu_e[i][j], mu_eta[i][j], mu_phi[i][j], mu_pt[i][j], 0, mu_q[i][j])
            else:
                mu_data[save_index][j] = (mu_e[i][j], mu_eta[i][j], mu_phi[i][j], mu_pt[i][j], mu_q[i][j])
            mu_labels[save_index][j] = 1 if mu_q[i][j] > 0 else 2

        # --- Fill Neutrino ---
        if nu_px[i].size > 0:
            nu_labels[save_index][0] = (nu_px[i][0], nu_py[i][0], nu_pz[i][0])

        # --- Fill MET ---
        if includeBTag:
            met_data[save_index][0] = (met_pt[i], 0, met_phi[i], met_pt[i], 0, 0)
        else:
            met_data[save_index][0] = (met_pt[i], 0, met_phi[i], met_pt[i], 0)
        met_labels[save_index][0] = 1

        # --- Fill Global ---
        if includeBTag:
            global_data[save_index] = (njet[i], nbTagged_90[i], nbTagged_85[i], nbTagged_77[i], nbTagged_70[i], nbTagged_65[i])
        else:
            global_data[save_index] = (njet[i],)

        # --- Truth Labels ---
        existent_jet_labels = np.array([0] * num_jets + [np.nan] * (pad_to_jet - num_jets))
        if 0 <= up_index[i] < num_jets:     existent_jet_labels[up_index[i]] = 3
        if 0 <= down_index[i] < num_jets:   existent_jet_labels[down_index[i]] = 4
        if 0 <= bhad_index[i] < num_jets:   existent_jet_labels[bhad_index[i]] = 2
        if 0 <= blep_index[i] < num_jets:   existent_jet_labels[blep_index[i]] = 1

        jet_labels[save_index]      = existent_jet_labels
        met_labels[save_index][0]   = 1
        # Lepton labels are already filled above during data insertion

        save_index += 1

    # ---------------------------------------------------
    #               Filter Fully Matched Events
    # ---------------------------------------------------
    if fullyMatch:
        required_jet_ids = {1, 2, 3, 4}
        required_met_ids = {1}

        keep_indices = []
        for i in range(save_index):
            jet_ids = set(np.nan_to_num(jet_labels[i], nan=-1))
            met_ids = set(np.nan_to_num(met_labels[i], nan=-1))
            # Check if either electron or muon is present with label 1
            el_has_lep = 1 in set(np.nan_to_num(el_labels[i], nan=-1))
            mu_has_lep = 1 in set(np.nan_to_num(mu_labels[i], nan=-1))
            
            if required_jet_ids.issubset(jet_ids) & required_met_ids.issubset(met_ids) & (el_has_lep or mu_has_lep):
                keep_indices.append(i)

        keep_indices = np.array(keep_indices)

        jet_labels    = jet_labels[keep_indices]
        el_labels     = el_labels[keep_indices]
        mu_labels     = mu_labels[keep_indices]
        met_labels    = met_labels[keep_indices]
        jet_data      = jet_data[keep_indices]
        el_data       = el_data[keep_indices]
        mu_data       = mu_data[keep_indices]
        met_data      = met_data[keep_indices]
        global_data   = global_data[keep_indices]
        global_labels = global_labels[keep_indices]
        nu_labels     = nu_labels[keep_indices]
        nu_data       = nu_data[keep_indices]

    # ---------------------------------------------------
    #               Save to HDF5
    # ---------------------------------------------------
    inputs_group.create_dataset("JET",      data=jet_data)
    inputs_group.create_dataset("ELECTRON", data=el_data)
    inputs_group.create_dataset("MUON",     data=mu_data)
    inputs_group.create_dataset("MET",      data=met_data)
    inputs_group.create_dataset("GLOBAL",   data=global_data)

    labels_group.create_dataset("JET",      data=jet_labels)
    labels_group.create_dataset("ELECTRON", data=el_labels)
    labels_group.create_dataset("MUON",     data=mu_labels)
    labels_group.create_dataset("MET",      data=met_labels)
    labels_group.create_dataset("NEUTRINO", data=nu_labels)

    h5_file.close()
    root_file.close()
    print(f"Successfully saved {save_index if not fullyMatch else len(keep_indices)} events to {output}")


if __name__ == "__main__":
    parser = argparser()
    parser.add_argument('-ds', '--dataset-size', type=int, help='Number of events to save', required=False)
    parser.add_argument('-fm', '--fully-matched', action='store_true', help='Keep only fully matched events', required=False)
    parser.add_argument('-b', '--btag', action='store_true', help='Have btag information in the h5 file', required=False)
    args = parser.parse_args()    
   
    main(MakeDatasetTtbarSingleLeptonVyper, args)