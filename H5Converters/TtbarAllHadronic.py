import h5py
import numpy as np
from tqdm import tqdm
import uproot

from Common import argparser, main

def MakeDatasetTtbarAllHad(args, input_file: str, output: str, split: bool = False, saveEven = True):
    r"""Function to construct HyPER dataset.

    Args:
        input (str): input ROOT file.
        output (str): output HyPER dataset in .h5 format.
        split (bool): split dataset in two according to (eventNumber mod 2).
    """
    # Ask user if they want to include the KLF info.
    include_klf = False
    answer_klf = input("Do you want to include the KLF info? (y/n): ")
    while answer_klf not in ['y', 'n']:
        answer_klf = input("Invalid input. Please enter 'y' or 'n': ")
    include_klf = answer_klf == 'y'

    # Max number of objects (padding)
    pad_to_jet = 20

    # Load the ROOT file
    branches = ["GN2v01_PCBT_quantile_NOSYS",
                "jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fEta",
                "jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPhi",
                "jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fE",
                "jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPt",
                "eventNumber",
                "nJets_good_NOSYS",
                "q1_index_NOSYS",
                "q2_index_NOSYS",
                "q1bar_index_NOSYS",
                "q2bar_index_NOSYS",
                "b_index_NOSYS",
                "bbar_index_NOSYS"]  # Choose relevant branches
    if include_klf:
        branches += ["KLF_b_from_top1_jetIndex_NOSYS",
                     "KLF_b_from_top2_jetIndex_NOSYS",
                     "KLF_lj1_from_top1_jetIndex_NOSYS",
                     "KLF_lj1_from_top2_jetIndex_NOSYS",
                     "KLF_lj2_from_top1_jetIndex_NOSYS",
                     "KLF_lj2_from_top2_jetIndex_NOSYS"]


    root_file = uproot.open(input_file)
    num_entries = root_file["reco"].num_entries
    tree = root_file["reco"].arrays(branches, library="np")

    # If you want to manually limit the number of events to save.
    limit_entries = True
    total_num_entries = num_entries
    if limit_entries:
        max_num_entries = 250000
        if max_num_entries > num_entries:
            print('You are requesting more events than available in the file. Exiting...')
            print('Number of events in the file: %d' % num_entries)
            exit(1)
        num_entries = max_num_entries

    # ---------------------------------------------------
    #                 Load data from ROOT
    # ---------------------------------------------------
    jet_bTag  = tree["GN2v01_PCBT_quantile_NOSYS"]
    jet_eta   = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fEta"]
    jet_phi   = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPhi"]
    jet_e     = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fE"]
    jet_pt    = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPt"]
    event_number = tree["eventNumber"]
    # ......
    njet      = tree["nJets_good_NOSYS"]
    # ...... labels
    q1_index = tree["q1_index_NOSYS"]
    q2_index = tree["q2_index_NOSYS"]
    q1bar_index = tree["q1bar_index_NOSYS"]
    q2bar_index = tree["q2bar_index_NOSYS"]
    b_index = tree["b_index_NOSYS"]
    bbar_index = tree["bbar_index_NOSYS"]
    if include_klf:
        klf_b_from_top1 = tree["KLF_b_from_top1_jetIndex_NOSYS"]
        klf_b_from_top2 = tree["KLF_b_from_top2_jetIndex_NOSYS"]
        klf_lj1_from_top1 = tree["KLF_lj1_from_top1_jetIndex_NOSYS"]
        klf_lj1_from_top2 = tree["KLF_lj1_from_top2_jetIndex_NOSYS"]
        klf_lj2_from_top1 = tree["KLF_lj2_from_top1_jetIndex_NOSYS"]
        klf_lj2_from_top2 = tree["KLF_lj2_from_top2_jetIndex_NOSYS"]

    # ---------------------------------------------------
    #               Count number of events to save
    # ---------------------------------------------------
    num_events_to_save = num_entries
    num_events_to_loop = num_entries
    if split:
        num_events_to_loop = 0
        num_events_to_save = 0
        for i in range(total_num_entries):
            if saveEven:
                if event_number[i] % 2 == 0:
                    num_events_to_save += 1
                num_events_to_loop += 1
            else:
                if event_number[i] % 2 != 0:
                    num_events_to_save += 1
                num_events_to_loop += 1
            if num_events_to_save >= num_entries:
                break

    # ---------------------------------------------------
    #               Define Input Datatype
    # ---------------------------------------------------
    node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), ('pt', np.float32), ('btag', np.float32)])
    jet_data = np.full((num_events_to_save, pad_to_jet), np.nan, dtype=node_dt)
    # ......
    if include_klf:
        global_dt = np.dtype([('njet', np.float32), ('nbTagged90', np.float32), ('nbTagged85', np.float32), ('nbTagged77', np.float32), ('nbTagged70', np.float32), ('nbTagged65', np.float32),
                              ('klf_b_from_top1', np.float32), ('klf_b_from_top2', np.float32),
                              ('klf_lj1_from_top1', np.float32), ('klf_lj1_from_top2', np.float32),
                              ('klf_lj2_from_top1', np.float32), ('klf_lj2_from_top2', np.float32)])
    else:
        global_dt = np.dtype([('njet', np.float32), ('nbTagged90', np.float32), ('nbTagged85', np.float32), ('nbTagged77', np.float32), ('nbTagged70', np.float32), ('nbTagged65', np.float32)])
    global_data = np.full((num_events_to_save, 1), np.nan, dtype=global_dt)

    # We need `VertexID` to label truth egdes and hyperedges
    jet_indices = np.full((num_events_to_save, pad_to_jet), np.nan)  # use -9 for none filled values
    fully_matched = np.full(num_events_to_save, np.nan)

    # ---------------------------------------------------
    #               Define HDF5 Structure
    # ---------------------------------------------------
    h5_file = h5py.File(output, 'w')
    inputs_group = h5_file.create_group('INPUTS')
    labels_group = h5_file.create_group('LABELS')
    metadata_group = h5_file.create_group("METADATA")

    print('Saving %d events' % num_events_to_save)
    save_index = 0
    for i in tqdm(range(num_events_to_loop), desc='Creating HDF5 dataset', unit='events'):
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
        num_btags_90 = 0
        num_btags_85 = 0
        num_btags_77 = 0
        num_btags_70 = 0
        num_btags_65 = 0
        for j in range(num_jets):
            jet_data[save_index][j] = (
                jet_e[i][j]/1000, # MeV to GeV
                jet_eta[i][j],
                jet_phi[i][j],
                jet_pt[i][j]/1000, # MeV to GeV
                (jet_bTag[i][j])
            )
            if jet_bTag[i][j] > 1:
                num_btags_90 += 1
                if jet_bTag[i][j] > 2:
                    num_btags_85 += 1
                    if jet_bTag[i][j] > 3:
                        num_btags_77 += 1
                        if jet_bTag[i][j] > 4:
                            num_btags_70 += 1
                            if jet_bTag[i][j] > 5:
                                num_btags_65 += 1

        if include_klf:
            global_data[save_index] = (njet[i], num_btags_90, num_btags_85, num_btags_77, num_btags_70, num_btags_65,
                                       klf_b_from_top1[i], klf_b_from_top2[i],
                                       klf_lj1_from_top1[i], klf_lj1_from_top2[i],
                                       klf_lj2_from_top1[i], klf_lj2_from_top2[i])
        else:
            global_data[save_index] = (njet[i], num_btags_90, num_btags_85, num_btags_77, num_btags_70, num_btags_65)

        existent_jet_labels = np.array([0]*num_jets)
        existent_jet_labels[b_index[i]] = 1
        existent_jet_labels[q1_index[i]] = 2
        existent_jet_labels[q2_index[i]] = 3
        existent_jet_labels[bbar_index[i]] = 4
        existent_jet_labels[q1bar_index[i]] = 5
        existent_jet_labels[q2bar_index[i]] = 6

        jet_indices[save_index, :num_jets] = existent_jet_labels
        fully_matched[save_index] = True
        save_index += 1
        

        if save_index >= num_events_to_save:
            print('Saved %d events' % save_index)
            break

    inputs_group.create_dataset("JET", data=jet_data)
    inputs_group.create_dataset("GLOBAL", data=global_data)
    labels_group.create_dataset("JET", data=jet_indices)
    metadata_group.create_dataset("FullyMatched", data=fully_matched)


    # Close the HDF5 file
    h5_file.close()
    root_file.close()


if __name__ == "__main__":

    parser = argparser()

    args = parser.parse_args()

    main(MakeDatasetTtbarAllHad, args)