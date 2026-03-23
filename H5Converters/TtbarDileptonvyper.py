import h5py
import numpy as np
from tqdm import tqdm
import uproot
import awkward as ak

from Common import argparser, main
"""
python TtbarDileptonvyper.py -f /eos/user/l/leet/Rebuild/output_ntuple_new/ttll_601230_mc23d_fullsim_Njobs_4_jobIndex_0.root -o tll_601230_mc23d_fullsim.h5 -ds 2000000 -b -s 
"""


def MakeDatasetTtbarDileptonvyper(args, loop_storage, output: str, split: bool = False, saveEven: bool = True):
    r"""Function to construct HyPER dataset.

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
    #includerandom = args.include_random
    #includeClassification = args.include_classification

    # Max number of objects (padding)
    pad_to_jet = 20
    pad_to_lepton = 2

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

    # Leptons features
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
    

    

    nu_e_px  = tree["nu_e_px"].array(library='np')
    nu_e_py = tree["nu_e_py"].array(library='np')
    nu_e_pz = tree["nu_e_pz"].array(library='np')
    
    nu_mu_px  = tree["nu_mu_px"].array(library='np')
    nu_mu_py = tree["nu_mu_py"].array(library='np')
    nu_mu_pz = tree["nu_mu_pz"].array(library='np')
    
    nu_px = [np.concatenate([np.atleast_1d(nu_e), np.atleast_1d(nu_mu)], axis=0)
          for nu_e, nu_mu in zip(nu_e_px, nu_mu_px)]
    nu_py = [np.concatenate([np.atleast_1d(nu_e), np.atleast_1d(nu_mu)], axis=0)
          for nu_e, nu_mu in zip(nu_e_py, nu_mu_py)]
    nu_pz = [np.concatenate([np.atleast_1d(nu_e), np.atleast_1d(nu_mu)], axis=0)
          for nu_e, nu_mu in zip(nu_e_pz, nu_mu_pz)]
    
    
    # nu_from_t_eta  = tree["nu_from_t_eta"].array(library='np')
    # nu_from_t_phi  = tree["nu_from_t_phi"].array(library='np')
    # nu_from_t_pt   = tree["nu_from_t_pt"].array(library='np')
    # nu_from_t_e   = tree["nu_from_t_e"].array(library='np')
    
    # nu_from_tbar_eta  = tree["nu_from_tbar_eta"].array(library='np')
    # nu_from_tbar_phi  = tree["nu_from_tbar_phi"].array(library='np')
    # nu_from_tbar_pt   = tree["nu_from_tbar_pt"].array(library='np')
    # nu_from_tbar_e   = tree["nu_from_tbar_e"].array(library='np')

    # nu_eta = [np.concatenate([np.atleast_1d(nu_t), np.atleast_1d(nu_tbar)], axis=0)
    #       for nu_t, nu_tbar in zip(nu_from_t_eta, nu_from_tbar_eta)]
    # nu_phi = [np.concatenate([np.atleast_1d(nu_t), np.atleast_1d(nu_tbar)], axis=0)
    #         for nu_t, nu_tbar in zip(nu_from_t_phi, nu_from_tbar_phi)]
    # nu_pt  = [np.concatenate([np.atleast_1d(nu_t), np.atleast_1d(nu_tbar)], axis=0)
    #           for nu_t, nu_tbar in zip(nu_from_t_pt, nu_from_tbar_pt)]
    # nu_e   = [np.concatenate([np.atleast_1d(nu_t), np.atleast_1d(nu_tbar)], axis=0)
    #          for nu_t, nu_tbar in zip(nu_from_t_e, nu_from_tbar_e)]
    # lep_phi = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_phi, mu_phi)]
    # lep_e   = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_e, mu_e)]
    # lep_pt  = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_pt, mu_pt)]
    # lep_q   = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_q, mu_q)]
    
    #isTtbarDilepton = tree["is_fully_matched_NOSYS"].array(library='np').astype(np.int32)

    # Build per-event lep_id array: electrons first, then muons (since we do concat [el, mu])
    # el_id = np.full((num_entries, pad_to_lepton), 0, dtype=np.int32)  # 0 = empty
    # mu_id = np.full((num_entries, pad_to_lepton), 0, dtype=np.int32) 
    # get per-event counts robustly
    def count_nonzero_length(arrays):
        counts = []
        for arr in arrays:
            counts.append(0 if arr.size == 0 else arr.size)
        return np.array(counts, dtype=int)

    num_electrons_per_event = count_nonzero_length(el_pt)
    num_muons_per_event     = count_nonzero_length(mu_pt)

    #fill lep_id for each event: electrons -> id 2, muons -> id 3
    # for i_e in range(num_entries):
    #     ne = num_electrons_per_event[i_e]
    #     nm = num_muons_per_event[i_e]
    #     # electrons occupy indices [0:ne]
    #     if ne > 0:
    #         lep_id[i_e, :min(ne, pad_to_lepton)] = 2
    #     # muons occupy next slots
    #     start = ne
    #     end = ne + nm
    #     if start < pad_to_lepton and nm > 0:
    #         lep_id[i_e, start:min(end, pad_to_lepton)] = 3

    # Met features
    met_phi   = tree["Met_phi"].array(library='np')
    met_pt    = tree["Met_pt"].array(library='np')

    # Global features (eventNumber only used for the splitting)
    N = len(tree["Jet_pt"].array(library='np'))
    event_number = np.arange(N)
    #event_number = tree["eventNumber"].array(library='np')
    njet         = (ak.count(tree["Jet_pt"].array(),axis=1).to_numpy()).reshape(-1,1)
    nlep         = (ak.count(tree["Electron_pt"].array(),axis=1).to_numpy()).reshape(-1,1) + (ak.count(tree["Muon_pt"].array(),axis=1).to_numpy()).reshape(-1,1) #np.count_nonzero(~np.isnan(lep_pt), axis=1)
    #nbTagged_100 = tree["nBjets_GN2v01_100_NOSYS"].array(library='np')
    nbTagged_90  = tree["nBjets_GN2v01_90_NOSYS"].array(library='np')
    nbTagged_85  = tree["nBjets_GN2v01_85_NOSYS"].array(library='np')
    nbTagged_77  = tree["nBjets_GN2v01_77_NOSYS"].array(library='np')
    nbTagged_70  = tree["nBjets_GN2v01_70_NOSYS"].array(library='np')
    nbTagged_65  = tree["nBjets_GN2v01_65_NOSYS"].array(library='np')

    # Labels (from parton truth)
    #up_index = tree["Truth_up_index"].array(library='np')
    #down_index = tree["Truth_down_index"].array(library='np')
    bhad_index = tree["Truth_had_b_index"].array(library='np')
    bbarhad_index = tree["Truth_had_bbar_index"].array(library='np')

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

    if includeBTag:
        node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), ('pt', np.float32), ('btag', np.float32), ('charge', np.float32)])
    else:
        node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), ('pt', np.float32), ('charge', np.float32)])

    jet_data = np.full((num_events_to_save, pad_to_jet), np.nan, dtype=node_dt)
    el_data = np.full((num_events_to_save, pad_to_lepton), np.nan, dtype=node_dt)
    mu_data = np.full((num_events_to_save, pad_to_lepton), np.nan, dtype=node_dt)
    met_data = np.full((num_events_to_save, 1), np.nan, dtype=node_dt)
    
    pad_to_nu = 2

    # nu_dt = np.dtype([
    #     ('eta', np.float32),
    #     ('phi', np.float32),
    #     ('pt', np.float32),
    # ])
    nu_dt = np.dtype([('px',np.float32), ('py', np.float32), ('pz', np.float32)])

    nu_data = np.full((num_events_to_save, pad_to_nu), np.nan, dtype=node_dt)
    nu_labels = np.full((num_events_to_save, pad_to_nu), np.nan, dtype = nu_dt)


    #Define global var
    if includeBTag:
        global_dt = np.dtype([('njet', np.float32), ('met_pt', np.float32), ('met_phi', np.float32),('nBJets_90', np.float32), ('nBJets_85', np.float32), ('nBJets_77', np.float32), ('nBJets_70', np.float32), ('nBJets_65', np.float32)])
    else:
        global_dt = np.dtype([('njet', np.float32), ('met_pt', np.float32), ('met_phi', np.float32)])
    
    global_data = np.zeros((num_events_to_save, 1), dtype=global_dt)


    # We need `VertexID` to label truth edges and hyperedges
    jet_labels      = np.full((num_events_to_save, pad_to_jet), np.nan)  # use -9 for none filled values
    el_labels   = np.full((num_events_to_save, pad_to_lepton), np.nan)
    mu_labels   = np.full((num_events_to_save, pad_to_lepton), np.nan)
    met_labels      = np.full((num_events_to_save, 1), np.nan)
    global_labels   = np.full((num_events_to_save, 1), np.nan)

    # ---------------------------------------------------
    #               Define HDF5 Structure
    # ---------------------------------------------------
    h5_file = h5py.File(output, 'w')
    inputs_group = h5_file.create_group('INPUTS')
    labels_group = h5_file.create_group('LABELS')

    print(f'Saving {num_events_to_save} events')
    save_index = 0
    for i in tqdm(range(num_entries), desc='Creating HDF5 dataset', unit='events'):

        # Skip if splitting and event doesn't match
        if split and ((event_number[i] % 2 == 0) != saveEven):
            continue
        if save_index >= num_events_to_save:
            break

        # Fill jets
        num_jets = njet[i][0]
        for j in range(num_jets):
            if includeBTag:
                jet_data[save_index][j] = (jet_e[i][j], jet_eta[i][j], jet_phi[i][j], jet_pt[i][j], jet_bTag[i][j], 0)
            else:
                jet_data[save_index][j] = (jet_e[i][j], jet_eta[i][j], jet_phi[i][j], jet_pt[i][j], 0)

        # Fill leptons
        num_el = el_pt[i].size
        num_mu = mu_pt[i].size
        #num_nu = nu_pt[i].size

        existent_el_labels = np.full(pad_to_lepton, np.nan)
        existent_mu_labels = np.full(pad_to_lepton, np.nan)
        existent_nu_labels = np.full(pad_to_nu, np.nan)

        for j in range(min(num_el, pad_to_lepton)):
            if includeBTag:
                el_data[save_index][j] = (el_e[i][j], el_eta[i][j], el_phi[i][j], el_pt[i][j], 0, el_q[i][j])
            else:
                el_data[save_index][j] = (el_e[i][j], el_eta[i][j], el_phi[i][j], el_pt[i][j], el_q[i][j])
                
            existent_el_labels[j] = 1 if el_q[i][j] > 0 else 2

        for j in range(min(num_mu, pad_to_lepton)):
            if includeBTag:
                mu_data[save_index][j] = (mu_e[i][j], mu_eta[i][j], mu_phi[i][j], mu_pt[i][j], 0, mu_q[i][j])
            else:
                mu_data[save_index][j] = (mu_e[i][j], mu_eta[i][j], mu_phi[i][j], mu_pt[i][j], mu_q[i][j])
            existent_mu_labels[j] = 1 if mu_q[i][j] > 0 else 2


        # for j in range(min(num_nu, pad_to_nu)):
        #     if includeBTag:
        #         nu_data[save_index][j] = (nu_e[i][j], nu_eta[i][j], nu_phi[i][j], nu_pt[i][j], 0, 0)
        #     else:
        #         nu_data[save_index][j] = (nu_e[i][j], nu_eta[i][j], nu_phi[i][j], nu_pt[i][j], 0)
                
        # if num_nu >= 1:
        #     nu_labels[save_index][0] = 1
        # if num_nu >= 2:
        #     nu_labels[save_index][1] = 2

        nu_labels[save_index][0] = (nu_px[i][0], nu_py[i][0], nu_pz[i][0])
        nu_labels[save_index][1] = (nu_px[i][1], nu_py[i][1], nu_pz[i][1])
        # Fill MET
        
        if includeBTag:
            met_data[save_index][0] = (met_pt[i], 0, met_phi[i], met_pt[i], 0, 0)
        else:
            met_data[save_index][0] = (met_pt[i], 0, met_phi[i], met_pt[i], 0)
        met_labels[save_index][0] = 1

        # Global data
        if includeBTag:
            global_data[save_index] = (njet[i], met_pt[i], met_phi[i], nbTagged_90[i], nbTagged_85[i], nbTagged_77[i], nbTagged_70[i], nbTagged_65[i])
        else:
            global_data[save_index] = (njet[i], met_pt[i], met_phi[i])

        # Jet labels
        existent_jet_labels = np.array([0] * num_jets + [np.nan] * (pad_to_jet - num_jets))
        if 0 <= bhad_index[i] < num_jets:
            existent_jet_labels[bhad_index[i]] = 1
        if 0 <= bbarhad_index[i] < num_jets:
            existent_jet_labels[bbarhad_index[i]] = 2

        jet_labels[save_index] = existent_jet_labels
        el_labels[save_index] = existent_el_labels
        mu_labels[save_index] = existent_mu_labels
        #nu_labels[save_index] = existent_nu_labels

        save_index += 1

    # Save datasets

    #Eventually filter the dataset to keep only the fully matched events (preferably done with FastFrames)
    if fullyMatch:
        required_jet_ids = {1, 2} #DILEPTON NTUPLES ARE ALREADY FULLY MATCHED IN FASTFRAMES CODE
        #required_lep_ids = {1, 2}
        #required_met_ids = {1}

        keep_indices = [
            i for i in range(len(jet_labels)) 
            if required_jet_ids.issubset( #Check indices in union of the sets
                set(np.nan_to_num(jet_labels[i], nan=-1))         
            ) and ((el_labels[i][0] == 1 and mu_labels[i][0] == 2) or (el_labels[i][0] == 2 and mu_labels[i][0] == 1))
            # ) & required_met_ids.issubset(
            #     set(np.nan_to_num(met_labels[i], nan=-1))
            # )
        ]

        jet_labels    = jet_labels[keep_indices]
        #lepton_labels = lepton_labels[keep_indices]
        el_labels     = el_labels[keep_indices]
        mu_labels     = mu_labels[keep_indices]
        met_labels    = met_labels[keep_indices]
        jet_data      = jet_data[keep_indices]
        #lep_data      = lep_data[keep_indices]
        el_data       = el_data[keep_indices]
        mu_data       = mu_data[keep_indices]
        met_data      = met_data[keep_indices]
        global_data   = global_data[keep_indices]
        global_labels = global_labels[keep_indices]
        nu_labels     = nu_labels[keep_indices]
        nu_data       = nu_data[keep_indices]
    # random node (placeholder)

    #Inputs h5 datasets
    inputs_group.create_dataset("JET", data=jet_data)
    inputs_group.create_dataset("ELECTRON", data=el_data)
    inputs_group.create_dataset("MUON", data=mu_data)
    inputs_group.create_dataset("MET", data=met_data)
    inputs_group.create_dataset("GLOBAL", data=global_data)
    #inputs_group.create_dataset("NEUTRINO", data=nu_data)

    labels_group.create_dataset("JET", data=jet_labels)
    labels_group.create_dataset("ELECTRON", data=el_labels)
    labels_group.create_dataset("MUON", data=mu_labels)
    labels_group.create_dataset("MET", data=met_labels)
    labels_group.create_dataset("NEUTRINO", data=nu_labels)
    #labels_group.create_dataset("GLOBAL", data=global_labels)
    # Close the HDF5 file
    h5_file.close()
    root_file.close()


if __name__ == "__main__":

    parser = argparser()

    parser.add_argument('-ds', '--dataset-size', type=int, help='Number of events to save', required=False)
    parser.add_argument('-fm', '--fully-matched', action='store_true', help='Keep only fully matched events', required=False)
    parser.add_argument('-b', '--btag', action='store_true', help='Have btag information in the h5 file', required=False)
    #parser.add_argument('-ir', '--include-random', action='store_true', help='Include random node in the h5 file', required=False)
    #parser.add_argument('-ic', '--include-classification', action='store_true', help='Include classification labels in the h5 file', required=False)
    args = parser.parse_args()    
   
    main(MakeDatasetTtbarDileptonvyper, args)