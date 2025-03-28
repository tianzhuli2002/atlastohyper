import h5py
import numpy as np
from tqdm import tqdm
import uproot

from Common import argparser, main, LoopStorage

def load_tree(input_file, loop_storage):
    # Ask user if they want to include the KLF info.
    include_klf = False
    answer_klf = input("Do you want to include the KLF info? (y/n): ")
    while answer_klf not in ['y', 'n']:
        answer_klf = input("Invalid input. Please enter 'y' or 'n': ")
    include_klf = answer_klf == 'y'


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
    loop_storage.n_events_in_tree = root_file["reco"].num_entries
    tree = root_file["reco"].arrays(branches, library="np")
    loop_storage.include_klf = include_klf
    loop_storage.tree = tree

    root_file.close()

    return loop_storage

def fill_storage(loop_storage, args, num_events_to_save = -1, range_to_save = [0,0]):
    # ---------------------------------------------------
    #               Count number of events to save
    # ---------------------------------------------------
    total_num_entries = loop_storage.n_events_in_tree
    event_number = loop_storage.tree["eventNumber"]
    n_even_events = 0
    n_odd_events = 0
    even_array = []
    odd_array = []
    whole_array = []
    num_events_to_save = total_num_entries if num_events_to_save == -1 else num_events_to_save
    
    for i in tqdm(range(total_num_entries), desc='Counting events in HDF5 dataset', unit='events'):
        if event_number[i] % 2 == 0:
            n_even_events += 1
            even_array.append(i)
        if event_number[i] % 2 != 0:
            n_odd_events += 1
            odd_array.append(i)

        whole_array.append(i)
        
    # Check that there are enough events to save.
    if args.split:
        if (num_events_to_save > n_even_events) or (num_events_to_save > n_odd_events):
            print('Error:: You are requesting more events than available in the file.')
            print('Info:: N even events = %d' % n_even_events)
            print('Info:: N odd events = %d' % n_odd_events)
            exit(1)
    else:
        if num_events_to_save > total_num_entries:
            print('Error:: You are requesting more events than available in the file.')
            print('Info:: N events = %d' % total_num_entries)
            exit(1)
    
    # Check that the range is valid if provided.
    if range_to_save != [0,0]:
        if range_to_save[1] - range_to_save[0] != num_events_to_save:
            print('Error:: The range provided does not match the number of events to save.')
            print('Info:: Range provided = [%d, %d]' % (range_to_save[0], range_to_save[1]))
            exit(1)
        # Check that upper limit is not greater than the total number of events.
        if args.split:
            if range_to_save[1] > n_even_events:
                print('Error:: The upper limit of the range is greater than the number of even events.')
                print('Info:: Upper limit = %d' % range_to_save[1])
                print('Info:: N even events = %d' % n_even_events)
                exit(1)
            if range_to_save[1] > n_odd_events:
                print('Error:: The upper limit of the range is greater than the number of odd events.')
                print('Info:: Upper limit = %d' % range_to_save[1])
                print('Info:: N odd events = %d' % n_odd_events)
                exit(1)
        else: 
            if range_to_save[1] > total_num_entries:
                print('Error:: The upper limit of the range is greater than the total number of events.')
                print('Info:: Upper limit = %d' % range_to_save[1])
                print('Info:: N events = %d' % total_num_entries)
                exit(1)

        

    else:
        range_to_save = [0, num_events_to_save]

    loop_storage.n_events_to_save = num_events_to_save
    loop_storage.whole_array = whole_array[range_to_save[0]:range_to_save[1]]
    if args.split:
        loop_storage.even_array = even_array[range_to_save[0]:range_to_save[1]]
        loop_storage.odd_array = odd_array[range_to_save[0]:range_to_save[1]]

    return loop_storage

def MakeDatasetTtbarAllHad(args, loop_storage, output: str, split: bool = False, saveEven = True):
    r"""Function to construct HyPER dataset.

    Args:
        input (str): input ROOT file.
        output (str): output HyPER dataset in .h5 format.
        split (bool): split dataset in two according to (eventNumber mod 2).
    """
    # Max number of objects (padding)
    pad_to_jet = 20

    if loop_storage == None:
        print('Error:: Loop storage not provided.')
        exit(1)

    include_klf = loop_storage.include_klf
    tree = loop_storage.tree
    num_events_to_save = loop_storage.n_events_to_save
    # Set the appropriate vector to loop over.
    if split:
        if saveEven:
            loop_array = loop_storage.even_array
        else:
            loop_array = loop_storage.odd_array
    else:
        loop_array = loop_storage.whole_array


    # ---------------------------------------------------
    #                 Load data from ROOT
    # ---------------------------------------------------
    jet_bTag  = tree["GN2v01_PCBT_quantile_NOSYS"]
    jet_eta   = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fEta"]
    jet_phi   = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPhi"]
    jet_e     = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fE"]
    jet_pt    = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPt"]
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

    eventNumber = tree["eventNumber"]

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

    eventNumber_dt = np.dtype(np.int32)
    eventNumber_data = np.full((num_events_to_save, 1), np.nan, dtype=eventNumber_dt)

    # ---------------------------------------------------
    #               Define HDF5 Structure
    # ---------------------------------------------------
    h5_file = h5py.File(output, 'w')
    inputs_group = h5_file.create_group('INPUTS')
    labels_group = h5_file.create_group('LABELS')
    metadata_group = h5_file.create_group("METADATA")

    print('Saving %d events' % num_events_to_save)
    save_index = 0
    for i in tqdm(loop_array, desc='Creating HDF5 dataset', unit='events'):
        
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
        eventNumber_data[save_index] = eventNumber[i]
        save_index += 1
        

        if save_index >= num_events_to_save:
            print('Saved %d events' % save_index)
            break

    inputs_group.create_dataset("JET", data=jet_data)
    inputs_group.create_dataset("GLOBAL", data=global_data)
    labels_group.create_dataset("JET", data=jet_indices)
    metadata_group.create_dataset("FullyMatched", data=fully_matched)
    metadata_group.create_dataset("eventNumber", data=eventNumber_data)


    # Close the HDF5 file
    h5_file.close()

# Inherits from LoopStorage
class LoopStorageTtbarAllHad(LoopStorage):
    def __init__(self):
        super().__init__()
        self.whole_array = []
        self.even_array = []
        self.odd_array = []
        self.include_klf = False
        self.include_hyper = False


if __name__ == "__main__":

    parser = argparser()

    args = parser.parse_args()

    loop_storage = LoopStorageTtbarAllHad()

    loop_storage = load_tree(args.input, loop_storage)

    # INPUT!!!!! Define the number of events you want to save and the range.
    num_events_to_save = -1
    range_to_save = [0,0]
    loop_storage = fill_storage(loop_storage, args, num_events_to_save, range_to_save)

    main(MakeDatasetTtbarAllHad, args, loop_storage)