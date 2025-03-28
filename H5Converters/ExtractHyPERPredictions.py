import h5py
import numpy as np
from tqdm import tqdm
import uproot

from Common import argparser, main, LoopStorage

def load_tree(input_file, loop_storage):

    # Load the ROOT file
    branches = ["eventNumber",
                "HyPER_b_from_top1_jetIndex_NOSYS",
                "HyPER_b_from_top2_jetIndex_NOSYS",
                "HyPER_lj1_from_top1_jetIndex_NOSYS",
                "HyPER_lj1_from_top2_jetIndex_NOSYS",
                "HyPER_lj2_from_top1_jetIndex_NOSYS",
                "HyPER_lj2_from_top2_jetIndex_NOSYS"]  # Choose relevant branches

    root_file = uproot.open(input_file)
    loop_storage.n_events_in_tree = root_file["reco"].num_entries
    tree = root_file["reco"].arrays(branches, library="np")
    loop_storage.tree = tree

    root_file.close()

    return loop_storage

def fill_storage(loop_storage, args, num_events_to_save = -1, range_to_save = [0,0]):
    # ---------------------------------------------------
    #               Count number of events to save
    # ---------------------------------------------------
    total_num_entries = loop_storage.n_events_in_tree
    whole_array = []
    num_events_to_save = total_num_entries if num_events_to_save == -1 else num_events_to_save
    
    for i in tqdm(range(total_num_entries), desc='Counting events in HDF5 dataset', unit='events'):
        whole_array.append(i)
        
    # Check that there are enough events to save.
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
        if range_to_save[1] > total_num_entries:
            print('Error:: The upper limit of the range is greater than the total number of events.')
            print('Info:: Upper limit = %d' % range_to_save[1])
            print('Info:: N events = %d' % total_num_entries)
            exit(1)

    range_to_save = [0, num_events_to_save]

    loop_storage.n_events_to_save = num_events_to_save
    loop_storage.whole_array = whole_array[range_to_save[0]:range_to_save[1]]

    return loop_storage

def MakeDatasetExtractHyPERPredictions(args, loop_storage, output: str, split: bool = False, saveEven = True):
    r"""Function to construct HyPER dataset.

    Args:
        input (str): input ROOT file.
        output (str): output HyPER dataset in .h5 format.
        split (bool): split dataset in two according to (eventNumber mod 2).
    """
    if loop_storage == None:
        print('Error:: Loop storage not provided.')
        exit(1)

    tree = loop_storage.tree
    num_events_to_save = loop_storage.n_events_to_save
    # Set the appropriate vector to loop over.
    loop_array = loop_storage.whole_array


    # ---------------------------------------------------
    #                 Load data from ROOT
    # ---------------------------------------------------
    eventNumber = tree["eventNumber"]
    hyper_b_from_top1 = tree["HyPER_b_from_top1_jetIndex_NOSYS"]
    hyper_b_from_top2 = tree["HyPER_b_from_top2_jetIndex_NOSYS"]
    hyper_lj1_from_top1 = tree["HyPER_lj1_from_top1_jetIndex_NOSYS"]
    hyper_lj1_from_top2 = tree["HyPER_lj1_from_top2_jetIndex_NOSYS"]
    hyper_lj2_from_top1 = tree["HyPER_lj2_from_top1_jetIndex_NOSYS"]
    hyper_lj2_from_top2 = tree["HyPER_lj2_from_top2_jetIndex_NOSYS"]

    # ---------------------------------------------------
    #               Define Input Datatype
    # ---------------------------------------------------

    top1_dt = np.dtype([('p1', np.float32), ('p2', np.float32), ('p3', np.float32)])
    top2_dt = np.dtype([('p1', np.float32), ('p2', np.float32), ('p3', np.float32)])
    w1_dt = np.dtype([('p1', np.float32), ('p2', np.float32)])
    w2_dt = np.dtype([('p1', np.float32), ('p2', np.float32)])
    eventNumber_dt = np.dtype(np.int32)

    top1_data = np.full((num_events_to_save, 1), np.nan, dtype=top1_dt)
    top2_data = np.full((num_events_to_save, 1), np.nan, dtype=top2_dt)
    w1_data = np.full((num_events_to_save, 1), np.nan, dtype=w1_dt)
    w2_data = np.full((num_events_to_save, 1), np.nan, dtype=w2_dt)
    eventNumber_data = np.full((num_events_to_save, 1), np.nan, dtype=eventNumber_dt)

    # ---------------------------------------------------
    #               Define HDF5 Structure
    # ---------------------------------------------------
    h5_file = h5py.File(output, 'w')
    tct_group = h5_file.create_group('TCT_OUTPUT')

    print('Saving %d events' % num_events_to_save)
    save_index = 0
    for i in tqdm(loop_array, desc='Creating HDF5 dataset', unit='events'):
        top1_data[save_index] = (hyper_b_from_top1[i], hyper_lj1_from_top1[i], hyper_lj2_from_top1[i])
        top2_data[save_index] = (hyper_b_from_top2[i], hyper_lj1_from_top2[i], hyper_lj2_from_top2[i])
        w1_data[save_index] = (hyper_lj1_from_top1[i], hyper_lj2_from_top1[i])
        w2_data[save_index] = (hyper_lj1_from_top2[i], hyper_lj2_from_top2[i])
        eventNumber_data[save_index] = int(eventNumber[i])
        save_index += 1

    tct_group.create_dataset("HyPER_best_top1", data=top1_data)
    tct_group.create_dataset("HyPER_best_top2", data=top2_data)
    tct_group.create_dataset("HyPER_best_w1", data=w1_data)
    tct_group.create_dataset("HyPER_best_w2", data=w2_data)
    tct_group.create_dataset("eventNumber", data=eventNumber_data)

    # Close the HDF5 file
    h5_file.close()

if __name__ == "__main__":

    parser = argparser()

    args = parser.parse_args()

    if args.split:
        print('Info:: Splitting dataset not supported for this conversion.')
        exit(1)

    loop_storage = LoopStorage()

    loop_storage = load_tree(args.input, loop_storage)

    # INPUT!!!!! Define the number of events you want to save and the range.
    num_events_to_save = -1
    range_to_save = [0,0]
    loop_storage = fill_storage(loop_storage, args, num_events_to_save, range_to_save)

    main(MakeDatasetExtractHyPERPredictions, args, loop_storage)