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
    includerandom = args.include_random
    includeClassification = args.include_classification

    # Max number of objects (padding)
    pad_to_jet = 20
    pad_to_lepton = 3

    # Load the ROOT file
    input = args.input
    root_file = uproot.open(input)
    tree = root_file['reco']
    num_entries = tree.num_entries

    # ---------------------------------------------------
    #                 Load data from ROOT
    # ---------------------------------------------------
    # jets features

    #jet_bTag  = tree["jet_GN2v01_PCBT_quantile_NOSYS"].array(library='np')
    jet_bTag = tree["jet_GN2v01_Continuous_quantile"].array(library='np')
    jet_eta   = tree["jet_eta"].array(library='np')
    jet_phi   = tree["jet_phi"].array(library='np')
    jet_e     = tree["jet_e_NOSYS"].array(library='np')/1000.0
    jet_pt    = tree["jet_pt_NOSYS"].array(library='np')/1000.0  #convert to GeV

    # Leptons features
    el_eta   = tree["el_eta"].array(library='np')
    el_phi   = tree["el_phi"].array(library='np')
    el_e     = tree["el_e_NOSYS"].array(library='np')/1000.0
    el_pt    = tree["el_pt_NOSYS"].array(library='np')/1000.0
    el_charge     = tree["el_charge"].array(library='np')

    mu_eta   = tree["mu_eta"].array(library='np')
    mu_phi   = tree["mu_phi"].array(library='np')
    mu_e     = tree["mu_e_NOSYS"].array(library='np')/1000.0
    mu_pt    = tree["mu_pt_NOSYS"].array(library='np')/1000.0
    mu_charge     = tree["mu_charge"].array(library='np')

    lep_eta = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_eta, mu_eta)]
    lep_phi = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_phi, mu_phi)]
    lep_e   = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_e, mu_e)]#/1000.0
    lep_pt  = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_pt, mu_pt)]#/1000.0
    lep_q   = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_charge, mu_charge)]
    
    #isTtbarDilepton = tree["is_fully_matched_NOSYS"].array(library='np').astype(np.int32)
    #isTtbarDilepton = tree["signal_NOSYS"].array(library='np').astype(np.int32)
    #isTtbarDilepton = tree["isTtbarDilepton_NOSYS"].array(library='np').astype(np.int32)
    #n_events = tree.num_entries
    isTtbarDilepton = np.ones(num_entries, dtype=np.int32)
    # Build per-event lep_id array: els first, then mus (since we do concat [el, mu])
    lep_id = np.full((num_entries, pad_to_lepton), 0, dtype=np.int32)  # 0 = empty
    # get per-event counts robustly
    def count_nonzero_length(arrays):
        counts = []
        for arr in arrays:
            counts.append(0 if arr.size == 0 else arr.size)
        return np.array(counts, dtype=int)

    num_els_per_event = count_nonzero_length(el_pt)
    num_mus_per_event     = count_nonzero_length(mu_pt)

    #fill lep_id for each event: els -> id 2, mus -> id 3
    for i_e in range(num_entries):
        ne = num_els_per_event[i_e]
        nm = num_mus_per_event[i_e]
        # els occupy indices [0:ne]
        if ne > 0:
            lep_id[i_e, :min(ne, pad_to_lepton)] = 2
        # mus occupy next slots
        start = ne
        end = ne + nm
        if start < pad_to_lepton and nm > 0:
            lep_id[i_e, start:min(end, pad_to_lepton)] = 3

    # Met features
    met_phi   = tree["met_phi_NOSYS"].array(library='np')
    met_pt    = tree["met_met_NOSYS"].array(library='np')/1000.0

    # Global features (eventNumber only used for the splitting)
    N = len(tree["jet_pt_NOSYS"].array(library='np'))
    #event_number = np.arange(N)
    event_number = tree["eventNumber"].array(library='np')
    njet         = (ak.count(tree["jet_pt_NOSYS"].array(),axis=1).to_numpy()).reshape(-1,1)
    nlep         = (ak.count(tree["el_pt_NOSYS"].array(),axis=1).to_numpy()).reshape(-1,1) + (ak.count(tree["mu_pt_NOSYS"].array(),axis=1).to_numpy()).reshape(-1,1) #np.count_nonzero(~np.isnan(lep_pt), axis=1)
    #nbTagged_100 = tree["nBjets_GN2v01_100_NOSYS"].array(library='np')
   

    # Load fixed-cut flags (jagged arrays: events × jets)
    FCBEff_90 = tree["jet_GN2v01_FixedCutBEff_90_select"].array()
    FCBEff_85 = tree["jet_GN2v01_FixedCutBEff_85_select"].array()
    FCBEff_77 = tree["jet_GN2v01_FixedCutBEff_77_select"].array()
    FCBEff_70 = tree["jet_GN2v01_FixedCutBEff_70_select"].array()
    FCBEff_65 = tree["jet_GN2v01_FixedCutBEff_65_select"].array()

    nbTagged_90 = ak.sum(FCBEff_90, axis=1).to_numpy()
    nbTagged_85 = ak.sum(FCBEff_85, axis=1).to_numpy()
    nbTagged_77 = ak.sum(FCBEff_77, axis=1).to_numpy()
    nbTagged_70 = ak.sum(FCBEff_70, axis=1).to_numpy()
    nbTagged_65 = ak.sum(FCBEff_65, axis=1).to_numpy()


    # nbTagged_90  = tree["nBjets_GN2v01_90_NOSYS"].array(library='np')
    # nbTagged_85  = tree["nBjets_GN2v01_85_NOSYS"].array(library='np')
    # nbTagged_77  = tree["nBjets_GN2v01_77_NOSYS"].array(library='np')
    # nbTagged_70  = tree["nBjets_GN2v01_70_NOSYS"].array(library='np')
    # nbTagged_65  = tree["nBjets_GN2v01_65_NOSYS"].array(library='np')

    # Labels (from parton truth)
    #up_index = tree["Truth_up_index"].array(library='np')
    #down_index = tree["Truth_down_index"].array(library='np')
    bhad_index = np.ones(num_entries, dtype=np.int32)
#tree["Truth_had_b_index"].array(library='np')
    bbarhad_index = np.ones(num_entries, dtype=np.int32) 
    #tree["Truth_had_bbar_index"].array(library='np')

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
        node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), ('pt', np.float32), ('btag', np.float32), ('charge', np.float32), ('id', np.float32)])
    else:
        node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), ('pt', np.float32), ('charge', np.float32), ('id', np.float32)])

    jet_data = np.full((num_events_to_save, pad_to_jet), np.nan, dtype=node_dt)
    lep_data = np.full((num_events_to_save, pad_to_lepton), np.nan, dtype=node_dt)
    met_data = np.full((num_events_to_save, 1), np.nan, dtype=node_dt)

    #Define global var
    if includeBTag:
        global_dt = np.dtype([('njet', np.float32), ('met_pt', np.float32), ('met_phi', np.float32),('nBjets_90', np.float32), ('nBjets_85', np.float32), ('nBjets_77', np.float32), ('nBjets_70', np.float32), ('nBjets_65', np.float32)])
    else:
        global_dt = np.dtype([('njet', np.float32), ('met_pt', np.float32), ('met_phi', np.float32)])
    
    global_data = np.zeros((num_events_to_save, 1), dtype=global_dt)


    # We need `VertexID` to label truth edges and hyperedges
    jet_labels      = np.full((num_events_to_save, pad_to_jet), np.nan)  # use -9 for none filled values
    lepton_labels   = np.full((num_events_to_save, pad_to_lepton), np.nan)
    met_labels      = np.full((num_events_to_save, 1), np.nan)
    global_labels   = np.full((num_events_to_save, 1), np.nan)

    # ---------------------------------------------------
    #               Define HDF5 Structure
    # ---------------------------------------------------
    h5_file = h5py.File(output, 'w')
    inputs_group = h5_file.create_group('INPUTS')
    labels_group = h5_file.create_group('LABELS')

    print('Saving %d events' % num_events_to_save)
    save_index = 0
    for i in tqdm(range(num_entries), desc='Creating HDF5 dataset', unit='events'):
        
        #Check if event should be saved
        if split and ((event_number[i] % 2 == 0) != saveEven):
            continue

        #If there is a limit on the number of events to save, check it
        if save_index >= num_events_to_save:
            break
                
        num_jets = njet[i][0]
        for j in range(num_jets):
            if includeBTag:
                jet_data[save_index][j] = (
                    jet_e[i][j], jet_eta[i][j], jet_phi[i][j], jet_pt[i][j],
                    jet_bTag[i][j], 0, 1
                )
            else:
                jet_data[save_index][j] = (
                    jet_e[i][j], jet_eta[i][j], jet_phi[i][j], jet_pt[i][j], 0, 1
                )

        #Assuming two leptons per event
        num_leps = nlep[i][0]
        for j in range(num_leps):
            if includeBTag:
                lep_data[save_index][j] = (
                    lep_e[i][j], lep_eta[i][j], lep_phi[i][j], lep_pt[i][j],
                    0, lep_q[i][j], lep_id[i][j]  # Filling in 0 for the btag , lep_id[i][j] 
                )
            else:
                lep_data[save_index][j] = (
                    lep_e[i][j], lep_eta[i][j], lep_phi[i][j], lep_pt[i][j],
                    lep_q[i][j], lep_id[i][j]  #, lep_id[i][j] 
                )

        if includeBTag:
            met_data[save_index][0] = (
                met_pt[i], 0, met_phi[i], met_pt[i], 0, 0, 4
            )
        else:
            met_data[save_index][0] = (
                met_pt[i], 0, met_phi[i], met_pt[i], 0, 4
                )
            
        if includeBTag:
            global_data[save_index] = (njet[i], met_pt[i], met_phi[i], nbTagged_90[i], nbTagged_85[i], nbTagged_77[i], nbTagged_70[i], nbTagged_65[i])
        else:
            global_data[save_index] = (njet[i], met_pt[i], met_phi[i])

        existent_jet_labels = np.array([0] * num_jets + [np.nan] * (pad_to_jet - num_jets))
        existent_lep_labels = np.array([0] * num_leps  + [np.nan] * (pad_to_lepton - num_leps))
        existent_met_labels = np.array([0] * 1)
        existent_global_labels = np.array([0] * 1)
        # if up_index[i] != -1:
        #     existent_jet_labels[up_index[i]] = 3
        # if down_index[i] != -1:
        #     existent_jet_labels[down_index[i]] = 4
        #existent_jet_labels = np.array([0] * num_jets + [np.nan] * (pad_to_jet - num_jets))

        if 0 <= bhad_index[i] < num_jets:
            existent_jet_labels[bhad_index[i]] = 1
        if 0 <= bbarhad_index[i] < num_jets:
            existent_jet_labels[bbarhad_index[i]] = 2

        # if bhad_index[i] != -1:
        #     existent_jet_labels[bhad_index[i]] = 1
        # if bbarhad_index[i] != -1:
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
        existent_met_labels[0] = 1

        jet_labels[save_index]      = existent_jet_labels
        lepton_labels[save_index]   = existent_lep_labels
        met_labels[save_index]      = existent_met_labels

        if includeClassification:
            existent_global_labels[0] = isTtbarDilepton[i]
            global_labels[save_index] = existent_global_labels
            
        save_index += 1

    #Eventually filter the dataset to keep only the fully matched events (preferably done with FastFrames)
    if fullyMatch:
        required_jet_ids = {1, 2} #DILEPTON NTUPLES ARE ALREADY FULLY MATCHED IN FASTFRAMES CODE
        required_lep_ids = {1, 2}
        #required_met_ids = {1}

        keep_indices = [
            i for i in range(len(jet_labels)) 
            if required_jet_ids.issubset( #Check indices in union of the sets
                set(np.nan_to_num(jet_labels[i], nan=-1))         
            ) & required_lep_ids.issubset(
                set(np.nan_to_num(lepton_labels[i], nan=-1)))
            # ) & required_met_ids.issubset(
            #     set(np.nan_to_num(met_labels[i], nan=-1))
            # )
        ]

        jet_labels    = jet_labels[keep_indices]
        lepton_labels = lepton_labels[keep_indices]
        met_labels    = met_labels[keep_indices]
        jet_data      = jet_data[keep_indices]
        lep_data      = lep_data[keep_indices]
        met_data      = met_data[keep_indices]
        global_data   = global_data[keep_indices]
        global_labels = global_labels[keep_indices]
    
    # random node (placeholder)

    #Inputs h5 datasets
    inputs_group.create_dataset("JET",     data=jet_data)
    inputs_group.create_dataset("LEPTON",  data=lep_data)
    inputs_group.create_dataset("MET",  data=met_data)
    inputs_group.create_dataset("GLOBAL", data=global_data)

    #Labels h5 datasets
    labels_group.create_dataset("JET", data=jet_labels) 
    labels_group.create_dataset("LEPTON", data=lepton_labels)
    labels_group.create_dataset("MET", data=met_labels)
    labels_group.create_dataset("GLOBAL", data=global_labels)
    
    if includerandom:
    # Create structured array of zeros first
        random_data = np.zeros((num_events_to_save, 1), dtype=node_dt)

        # Fill fields with random values in realistic ranges
        random_data['e'] = np.random.uniform(0, 500, size=(num_events_to_save, 1))
        random_data['pt'] = np.random.uniform(0, 300, size=(num_events_to_save, 1))
        random_data['eta'] = np.random.uniform(-2.5, 2.5, size=(num_events_to_save, 1))
        random_data['phi'] = np.random.uniform(-np.pi, np.pi, size=(num_events_to_save, 1))
        random_data['charge'] = np.random.choice([-1, 0, 1], size=(num_events_to_save, 1))
        
        if includeBTag:
            random_data['btag'] = np.random.uniform(0, 1, size=(num_events_to_save, 1))

        # Labels (placeholder)
        random_labels = np.ones((num_events_to_save, 1), dtype=np.float32)

        # Save to HDF5
        inputs_group.create_dataset("RANDOM", data=random_data)
        labels_group.create_dataset("RANDOM", data=random_labels)
    # Close the HDF5 file
    h5_file.close()
    root_file.close()


if __name__ == "__main__":

    parser = argparser()

    parser.add_argument('-ds', '--dataset-size', type=int, help='Number of events to save', required=False)
    parser.add_argument('-fm', '--fully-matched', action='store_true', help='Keep only fully matched events', required=False)
    parser.add_argument('-b', '--btag', action='store_true', help='Have btag information in the h5 file', required=False)
    parser.add_argument('-ir', '--include-random', action='store_true', help='Include random node in the h5 file', required=False)
    parser.add_argument('-ic', '--include-classification', action='store_true', help='Include classification labels in the h5 file', required=False)
    args = parser.parse_args()    
   
    main(MakeDatasetTtbarDilepton, args)