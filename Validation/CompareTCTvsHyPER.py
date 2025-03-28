import h5py
import numpy as np
from tqdm import tqdm

def compare_tct_hyper():
    # Open the two HDF5 files
    tct_file = h5py.File('../H5Converters/tct_definitive_test.h5', 'r')
    hyper_file = h5py.File('../H5Converters/predicted_definitive_test.h5', 'r')

    # Get the data from the HDF5 files
    event_number = tct_file['TCT_OUTPUT']['eventNumber'][:]
    tct_top1 = tct_file['TCT_OUTPUT']['HyPER_best_top1']
    
    hyper_top1 = hyper_file['OUTPUT']['HyPER_best_top1']

    different_events = []
    # Compare the two datasets
    for i in tqdm(range(len(event_number)), desc='Comparing datasets', unit='events'):
        if i < 10:
            print(i, event_number[i])
        tct = np.array([ int(j) for j in tct_top1[i][0]])
        hyper = np.array([ int(j) for j in hyper_top1[i][0]])
        tct.sort()
        hyper.sort()
        if not np.array_equal(tct, hyper):
            different_events.append(event_number[i])
            if len(different_events) < 10:
                print(tct, hyper, event_number[i])
    
    print('Number of events: %d' % len(event_number))
    print('Number of different events: %d' % len(different_events))

    print(different_events[0])

    return different_events


#This is skeleton code to print HyPER inference data using pickle files.
def print_pickle_prediction():
    import pandas as pd
    import numpy as np

    predict_path = 'predicted_mini_test.pkl'
    predict = pd.read_pickle(predict_path)

    columns_to_display = ['HyPER_HE_RAW', 'HyPER_GE_RAW','HyPER_best_top1','HyPER_best_top2'] 
    print(predict[columns_to_display].head(20))
    a = predict['HyPER_HE_RAW'][0]
    b = predict['HyPER_GE_RAW'][0]
    a.sort()
    b.sort()
    print("HyPER HE scores:")
    print(a)
    print("HyPER E scores:")
    print(b)

if __name__ == "__main__":
    diff_list = compare_tct_hyper()
    #print_pickle_prediction()
