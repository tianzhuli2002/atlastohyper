import h5py
import numpy as np
from tqdm import tqdm
import uproot
import awkward as ak

from Common import argparser, main


def MakeDatasetTtbarSingleLepton(args, loop_storage, output: str, split: bool = False, saveEven: bool = True):
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

    lep_eta = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_eta, mu_eta)]
    lep_phi = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_phi, mu_phi)]
    lep_e   = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_e, mu_e)]
    lep_pt  = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_pt, mu_pt)]
    lep_q   = [np.concatenate([el, mu], axis=0) for el, mu in zip(el_q, mu_q)]

    lep_id = np.full((num_entries, pad_to_lepton), 3) #id is 3 for muons, 2 for electrons
    el_pt_forCount = [] #used only to count the number of electrons per event
    for arr in el_pt: #fill empty with Nan to use isnan
        if arr.size == 0:
            el_pt_forCount.append([np.nan]) 
        else:
            el_pt_forCount.append(arr)
    el_pt_forCount = np.array(el_pt_forCount, dtype=np.float32)
    num_electrons_per_event = np.count_nonzero(~np.isnan(el_pt_forCount), axis=1)
    for i, num_electrons in enumerate(num_electrons_per_event): #replace 3 with 2 in the presence of electrons
        lep_id[i, :num_electrons] = 2

    lep_id[lep_e == 0] = 0 #empty events are assigned an id of 0

    # Met features
    met_phi   = tree["Met_phi"].array(library='np')
    met_pt    = tree["Met_pt"].array(library='np')

    # Global features (eventNumber only used for the splitting)
    event_number = tree["eventNumber"].array(library='np')
    njet         = (ak.count(tree["Jet_pt"].array(),axis=1).to_numpy()).reshape(-1,1)
    nlep         = np.count_nonzero(~np.isnan(lep_pt), axis=1)
    nbTagged_100 = tree["nBjets_GN2v01_100_NOSYS"].array(library='np')
    nbTagged_90  = tree["nBjets_GN2v01_90_NOSYS"].array(library='np')
    nbTagged_85  = tree["nBjets_GN2v01_85_NOSYS"].array(library='np')
    nbTagged_77  = tree["nBjets_GN2v01_77_NOSYS"].array(library='np')
    nbTagged_70  = tree["nBjets_GN2v01_70_NOSYS"].array(library='np')
    nbTagged_65  = tree["nBjets_GN2v01_65_NOSYS"].array(library='np')

    # Labels (from parton truth)
    up_index = tree["Truth_up_index"].array(library='np')
    down_index = tree["Truth_down_index"].array(library='np')
    bhad_index = tree["Truth_had_b_index"].array(library='np')
    blep_index = tree["Truth_lep_b_index"].array(library='np')

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
        node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), ('pt', np.float32), ('charge', np.float32), ('id', np.int32)])

    jet_data = np.full((num_events_to_save, pad_to_jet), np.nan, dtype=node_dt)
    lep_data = np.full((num_events_to_save, pad_to_lepton), np.nan, dtype=node_dt)
    met_data = np.full((num_events_to_save, 1), np.nan, dtype=node_dt)

    #Define global var
    if includeBTag:
        global_dt = np.dtype([('njet', np.float32), ('nBJets_90', np.float32), ('nBJets_85', np.float32), ('nBJets_77', np.float32), ('nBJets_70', np.float32), ('nBJets_65', np.float32)])
    else:
        global_dt = np.dtype([('njet', np.float32)])
    
    global_data = np.zeros((num_events_to_save, 1), dtype=global_dt)


    # We need `VertexID` to label truth edges and hyperedges
    jet_labels      = np.full((num_events_to_save, pad_to_jet), np.nan)  # use -9 for none filled values
    lepton_labels   = np.full((num_events_to_save, pad_to_lepton), np.nan)
    met_labels      = np.full((num_events_to_save, 1), np.nan)

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

        #Assuming only one lepton
        num_lep = 1
        if includeBTag:
            lep_data[save_index][0] = (
                lep_e[i][0], lep_eta[i][0], lep_phi[i][0], lep_pt[i][0],
                0, lep_q[i][0], lep_id[i][0]   # Filling in 0 for the btag
            )
        else:
            lep_data[save_index][0] = (
                lep_e[i][0], lep_eta[i][0], lep_phi[i][0], lep_pt[i][0],
                lep_q[i][0], lep_id[i][0] 
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
            global_data[save_index] = (njet[i], nbTagged_90[i], nbTagged_85[i], nbTagged_77[i], nbTagged_70[i], nbTagged_65[i])
        else:
            global_data[save_index] = (njet[i])

        existent_jet_labels = np.array([0] * num_jets + [np.nan] * (pad_to_jet - num_jets))
        existent_lep_labels = np.array([0] * num_lep  + [np.nan] * (pad_to_lepton - num_lep))
        existent_met_labels = np.array([0] * 1)
        
        if up_index[i] != -1:
            existent_jet_labels[up_index[i]] = 3
        if down_index[i] != -1:
            existent_jet_labels[down_index[i]] = 4
        if bhad_index[i] != -1:
            existent_jet_labels[bhad_index[i]] = 2
        if blep_index[i] != -1:
            existent_jet_labels[blep_index[i]] = 1
        existent_lep_labels[0] = 1
        existent_met_labels[0] = 1

        jet_labels[save_index]      = existent_jet_labels
        lepton_labels[save_index]   = existent_lep_labels
        met_labels[save_index]      = existent_met_labels

        save_index += 1

    #Eventually filter the dataset to keep only the fully matched events (preferably done with FastFrames)
    if fullyMatch:
        required_jet_ids = {1, 2, 3, 4}
        required_lep_ids = {1}
        required_met_ids = {1}

        keep_indices = [
            i for i in range(len(jet_labels)) 
            if required_jet_ids.issubset( #Check indices in union of the sets
                set(np.nan_to_num(jet_labels[i], nan=-1))         
            ) & required_lep_ids.issubset(
                set(np.nan_to_num(lepton_labels[i], nan=-1))
            ) & required_met_ids.issubset(
                set(np.nan_to_num(met_labels[i], nan=-1))
            )
        ]

        jet_labels    = jet_labels[keep_indices]
        lepton_labels = lepton_labels[keep_indices]
        met_labels    = met_labels[keep_indices]
        jet_data      = jet_data[keep_indices]
        lep_data      = lep_data[keep_indices]
        met_data      = met_data[keep_indices]
        global_data   = global_data[keep_indices]

    #Inputs h5 datasets
    inputs_group.create_dataset("JET",     data=jet_data)
    inputs_group.create_dataset("LEPTON",  data=lep_data)
    inputs_group.create_dataset("MET",  data=met_data)
    inputs_group.create_dataset("GLOBAL", data=global_data)

    #Labels h5 datasets
    labels_group.create_dataset("JET", data=jet_labels) 
    labels_group.create_dataset("LEPTON", data=lepton_labels)
    labels_group.create_dataset("MET", data=met_labels)

    # Close the HDF5 file
    h5_file.close()
    root_file.close()


if __name__ == "__main__":

    parser = argparser()

    parser.add_argument('-ds', '--dataset-size', type=int, help='Number of events to save', required=False)
    parser.add_argument('-fm', '--fully-matched', action='store_true', help='Keep only fully matched events', required=False)
    parser.add_argument('-b', '--btag', action='store_true', help='Have btag information in the h5 file', required=False)

    args = parser.parse_args()    
   
    main(MakeDatasetTtbarSingleLepton, args)