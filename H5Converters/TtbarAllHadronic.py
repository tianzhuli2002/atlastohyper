import h5py
import numpy as np
from tqdm import tqdm
import uproot

from Common import argparser, main

def MakeDatasetTtbarAllHad(input: str, output: str, split: bool = False, saveEven = True):
    r"""Function to construct HyPER dataset.

    Args:
        input (str): input ROOT file.
        output (str): output HyPER dataset in .h5 format.
        split (bool): split dataset in two according to (eventNumber mod 2).
    """
    # Max number of objects (padding)
    pad_to_jet = 20

    # Load the ROOT file
    root_file = uproot.open(input)
    tree = root_file['reco']
    num_entries = tree.num_entries
    #num_entries = 10000

    # ---------------------------------------------------
    #                 Load data from ROOT
    # ---------------------------------------------------
    jet_bTag  = tree["GN2v01_PCBT_quantile_NOSYS"].array(library='np')
    jet_eta   = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fEta"].array(library='np')
    jet_phi   = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPhi"].array(library='np')
    jet_e     = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fE"].array(library='np')
    jet_pt    = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPt"].array(library='np')
    event_number = tree["eventNumber"].array(library='np')
    # ......
    njet      = tree["nJets_good_NOSYS"].array(library='np')
    nbTagged  = tree["nBjets_GN2v01_77_NOSYS"].array(library='np')
    # ...... labels
    q1_index = tree["q1_index_NOSYS"].array(library='np')
    q2_index = tree["q2_index_NOSYS"].array(library='np')
    q1bar_index = tree["q1bar_index_NOSYS"].array(library='np')
    q2bar_index = tree["q2bar_index_NOSYS"].array(library='np')
    b_index = tree["b_index_NOSYS"].array(library='np')
    bbar_index = tree["bbar_index_NOSYS"].array(library='np')

    # ---------------------------------------------------
    #               Count number of events to save
    # ---------------------------------------------------
    num_events_to_save = num_entries
    if split:
        num_events_to_save = 0
        for i in range(num_entries):
            if saveEven:
                if event_number[i] % 2 == 0:
                    num_events_to_save += 1
            else:
                if event_number[i] % 2 != 0:
                    num_events_to_save += 1


    # ---------------------------------------------------
    #               Define Input Datatype
    # ---------------------------------------------------
    node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), ('pt', np.float32), ('btag', np.float32), ('id', np.float32)])
    jet_data = np.zeros((num_events_to_save, pad_to_jet), dtype=node_dt)
    # ......
    global_dt = np.dtype([('njet', np.float32), ('nbTagged', np.float32)])
    global_data = np.zeros((num_events_to_save, 1), dtype=global_dt)

    # We need `VertexID` to label truth egdes and hyperedges
    VertexID_data = np.full((num_events_to_save, pad_to_jet), -9)  # use -9 for none filled values
    IndexSelect_data = np.zeros(num_events_to_save)

    # ---------------------------------------------------
    #               Define HDF5 Structure
    # ---------------------------------------------------
    h5_file = h5py.File(output, 'w')
    inputs_group = h5_file.create_group('INPUTS')
    labels_group = h5_file.create_group('LABELS')

    print('Saving %d events' % num_events_to_save)
    save_index = 0
    for i in tqdm(range(num_entries), desc='Creating HDF5 dataset', unit='events'):
        # ---------------------------------------------------
        #               Decide if event will be saved
        # ---------------------------------------------------
        if split:
            if saveEven:
                if event_number[i] % 2 != 0:
                    continue
            else:
                if event_number[i] % 2 == 0:
                    continue
        # ---------------------------------------------------

        num_jets = njet[i]
        for j in range(num_jets):
            jet_data[save_index][j] = (
                jet_e[i][j], jet_eta[i][j], jet_phi[i][j], jet_pt[i][j],
                (jet_bTag[i][j] >= 4 ), 1
            )

        global_data[save_index] = (njet[i], nbTagged[i])

        existent_jet_labels = np.array([0]*num_jets)
        existent_jet_labels[b_index[i]] = 1
        existent_jet_labels[q1_index[i]] = 2
        existent_jet_labels[q2_index[i]] = 3
        existent_jet_labels[bbar_index[i]] = 4
        existent_jet_labels[q1bar_index[i]] = 5
        existent_jet_labels[q2bar_index[i]] = 6

        VertexID_data[save_index, :num_jets] = existent_jet_labels
        IndexSelect_data[save_index] = True

        save_index += 1

    inputs_group.create_dataset("jet", data=jet_data)
    inputs_group.create_dataset("global", data=global_data)
    labels_group.create_dataset("VertexID", data=VertexID_data)
    labels_group.create_dataset("IndexSelect", data=IndexSelect_data)


    # Close the HDF5 file
    h5_file.close()
    root_file.close()


if __name__ == "__main__":
    main(MakeDatasetTtbarAllHad)