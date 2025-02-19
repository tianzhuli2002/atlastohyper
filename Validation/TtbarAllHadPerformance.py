import h5py
import numpy as np
from tqdm import tqdm
import argparse
import matplotlib.pyplot as plt

def argparser():
    parser = argparse.ArgumentParser(description='Calculated the performance of HyPER on ttbar all-hadronic events given two files: one containing the labels and the other the predictions.')
    parser.add_argument('-f', '--inputs',  type=str, required=True, help='H5 input file containing the labels.')
    parser.add_argument('-o', '--predictions', type=str, required=True, help='H5 input file containing the predictions.')
    parser.add_argument('-k', '--klf-input', type=str, help='Input file with KLF info to measure KLF Performance.',default='')
    parser.add_argument('-p', '--plot', action='store_true', help='Plot the results.')
    args = parser.parse_args()

    return args

# This function takes a top represented by indices. It compares with a list of available tops (represented by indices) and returns True if the top was found in the list.
# It also returns the list of available tops without the top that was found.
def check_if_top_was_correct(top, available_tops):
    if len(available_tops) == 0:
        print('You tried to check if a top was correct, but there are no available tops.')
        exit(1)

    top_correctly_predicted = False
    top_found = []
    for available_top in available_tops:
        if np.array_equal(top, available_top):
            top_correctly_predicted = True
            top_found = available_top
            break
    
    if top_correctly_predicted:
        new_available_tops = []
        for available_top in available_tops:
            if not np.array_equal(available_top, top_found):
                new_available_tops.append(available_top)
        return True, new_available_tops
    
    return False, available_tops

# This function stores the results of the evaluation in dictionaries.
# It does the splitting of the results by the number of jets in the event.
def store_results(total_events_dict, two_correct_dict, one_correct_dict, n_jets_per_event, top1_correct, top2_correct):
    if n_jets_per_event > 10:
        n_jets_per_event_str = '10+'
    else:
        n_jets_per_event_str = str(int(n_jets_per_event))

    # Record the total number of events
    total_events_dict[n_jets_per_event_str] += 1

    # Record the number of events where both top 1 and top 2 were correct
    if top1_correct and top2_correct:
        two_correct_dict[n_jets_per_event_str] += 1
    
    # Record the number of events with only one correct top
    if (top1_correct and not top2_correct) or (not top1_correct and top2_correct):
        one_correct_dict[n_jets_per_event_str] += 1

def measure_hyper_performance(args):
    # Open the files.
    labels = h5py.File(args.inputs, 'r')
    predictions = h5py.File(args.predictions, 'r')

    # Create the dictionaries that will store the results.
    top_results_dict = {}
    top_results_dict['total_events'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    top_results_dict['two_correct'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    top_results_dict['one_correct'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}

    w_results_dict = {}
    w_results_dict['total_events'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    w_results_dict['two_correct'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    w_results_dict['one_correct'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}

    # Extract the data that we need.
    global_features = labels['INPUTS']['GLOBAL']
    n_input_events = len(labels['LABELS']['JET'])
    n_label_events = len(predictions['OUTPUT']['HyPER_best_top1'])

    # Sanity check. Number of events in the labels and predictions files should match.
    if n_input_events != n_label_events:
        print('The number of events in the labels and predictions files do not match.')
        exit(1)

    print('Evaluating the perfomance of %d events.' % n_input_events)
    for i_event in tqdm(range(n_input_events), desc='Evaluating performance...', unit='events'):
        # Top 1 and top 2 predictions are easy to get!
        top1_pred = list(predictions['OUTPUT']['HyPER_best_top1'][i_event][0])
        top2_pred = list(predictions['OUTPUT']['HyPER_best_top2'][i_event][0])
        w1_pred = list(predictions['OUTPUT']['HyPER_best_w1'][i_event][0])
        w2_pred = list(predictions['OUTPUT']['HyPER_best_w2'][i_event][0])

        # Sanity check. The number of jets in the event should be 6 or more.
        n_jets_per_event = global_features[i_event][0][0]
        if n_jets_per_event < 6:
            print('The number of jets in the event is not 6. What...')
            exit(1)

        # Build top 1 and top2 from input labels
        event_labels = labels['LABELS']['JET'][i_event]
        top1 = np.array([
            event_labels[top1_pred[0]],
            event_labels[top1_pred[1]],
            event_labels[top1_pred[2]],
            ])
        top1.sort()
        w1 = np.array([
            event_labels[w1_pred[0]],
            event_labels[w1_pred[1]],
            ])
        w1.sort()

        top2 = np.array([
            event_labels[top2_pred[0]],
            event_labels[top2_pred[1]],
            event_labels[top2_pred[2]],
            ])
        top2.sort()
        w2 = np.array([
            event_labels[w2_pred[0]],
            event_labels[w2_pred[1]],
            ])
        w2.sort()

        # Sanity check. Top 1 and top 2 should be different.
        if np.array_equal(top1, top2):
            print('Top 1 and top 2 are the same!')
            print('What?')
            exit(1)

        # Check if the top 1 and top 2 predictions are correct
        available_tops = [np.array([1,2,3]), np.array([4,5,6])]
        top1_correct, available_tops = check_if_top_was_correct(top1, available_tops)
        top2_correct, available_tops = check_if_top_was_correct(top2, available_tops)

        # Check if the W1 and W2 predictions are correct
        available_ws = [np.array([2,3]), np.array([5,6])]
        w1_correct, available_ws = check_if_top_was_correct(w1, available_ws)
        w2_correct, available_ws = check_if_top_was_correct(w2, available_ws)

        # Record results
        store_results(top_results_dict['total_events'], top_results_dict['two_correct'], top_results_dict['one_correct'], n_jets_per_event, top1_correct, top2_correct)
        store_results(w_results_dict['total_events'], w_results_dict['two_correct'], w_results_dict['one_correct'], n_jets_per_event, w1_correct, w2_correct)


    # Process the results.
    # For tops
    full_reco_eff = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    at_least_one_eff = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    # Calculate the efficiencies.
    for category in full_reco_eff:
        full_reco_eff[category] = 100*top_results_dict['two_correct'][category] / top_results_dict['total_events'][category]
        at_least_one_eff[category] = 100*(top_results_dict['two_correct'][category] + top_results_dict['one_correct'][category]) / top_results_dict['total_events'][category]
    # For Ws
    w_full_reco_eff = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    w_at_least_one_eff = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    # Calculate the efficiencies.
    for category in w_full_reco_eff:
        w_full_reco_eff[category] = 100*w_results_dict['two_correct'][category] / w_results_dict['total_events'][category]
        w_at_least_one_eff[category] = 100*(w_results_dict['two_correct'][category] + w_results_dict['one_correct'][category]) / w_results_dict['total_events'][category]


    # Print the results.
    print('Results:')
    print('Full reconstruction efficiency:')
    for category in full_reco_eff:
        print('For events with %s jets: %.3f' % (category, full_reco_eff[category]))

    print('At least one correct top efficiency:')
    for category in at_least_one_eff:
        print('For events with %s jets: %.3f' % (category, at_least_one_eff[category]))

    print('W Full reconstruction efficiency:')
    for category in w_full_reco_eff:
        print('For events with %s jets: %.3f' % (category, w_full_reco_eff[category]))
    
    print('W At least one correct top efficiency:')
    for category in w_at_least_one_eff:
        print('For events with %s jets: %.3f' % (category, w_at_least_one_eff[category]))

    # Plot the results in a histogram if requested
    if args.plot:
        x = np.arange(len(full_reco_eff))
        width = 0.35

        fig, ax = plt.subplots()
        rects1 = ax.bar(x , list(full_reco_eff.values()), width/4, label='Two tops reconstructed.')
        rects2 = ax.bar(x - width/2, list(at_least_one_eff.values()), width/4, label='At least one top reconstructed.')
        rects3 = ax.bar(x - width/4, list(w_full_reco_eff.values()), width/4, label='Two Ws reconstructed.')
        rects4 = ax.bar(x - 3*width/4 , list(w_at_least_one_eff.values()), width/4, label='At least one W reconstructed.')

        ax.set_ylabel('Efficiency')
        ax.set_xlabel('N jets')
        ax.set_title('HyPER reconstruction')
        ax.set_xticks(x)
        ax.set_xticklabels(full_reco_eff.keys())
        ax.legend()

        post_fix = 'even' if 'even' in args.inputs else'odd'
        plt.savefig('PerformanceTtbarAllHad_%s.png' % post_fix)

# The KLF predictions and labels are stored in the same file. The predictions are stored in the globals array.
def measure_klf_performance(args):
    # Open the files.
    h5_file = h5py.File(args.klf_input, 'r')

    # Create the dictionaries that will store the results.
    top_results_dict = {}
    top_results_dict['total_events'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    top_results_dict['two_correct'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    top_results_dict['one_correct'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}

    w_results_dict = {}
    w_results_dict['total_events'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    w_results_dict['two_correct'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    w_results_dict['one_correct'] = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}

    # Extract the data that we need.
    global_features = h5_file['INPUTS']['GLOBAL']
    n_input_events = len(h5_file['LABELS']['JET'])

    print('Evaluating the perfomance of %d events.' % n_input_events)
    for i_event in tqdm(range(n_input_events), desc='Evaluating KLF performance...', unit='events'):
        # Top 1 and top 2 predictions are easy to get!
        top1_pred = [int(global_features[i_event][0][6]), int(global_features[i_event][0][8]), int(global_features[i_event][0][10])]
        top2_pred = [int(global_features[i_event][0][7]), int(global_features[i_event][0][9]), int(global_features[i_event][0][11])]

        bad_prediction = False
        for i in range(3):
            if top1_pred[i] >= 20 or top1_pred[i] < 0:
                print('Invalid top 1 prediction!')
                print(top1_pred)
                bad_prediction = True
                break
            if top2_pred[i] >= 20 or top2_pred[i] < 0:
                print('Invalid top 2 prediction!')
                print(top2_pred)
                bad_prediction = True
                break
        if bad_prediction:
            continue

        # Sanity check. The number of jets in the event should be 6 or more.
        n_jets_per_event = global_features[i_event][0][0]
        if n_jets_per_event < 6:
            print('The number of jets in the event is not 6. What...')
            exit(1)

        # Build top 1 and top2 from input labels
        event_labels = h5_file['LABELS']['JET'][i_event]
        top1 = np.array([
            event_labels[top1_pred[0]],
            event_labels[top1_pred[1]],
            event_labels[top1_pred[2]],
            ])
        top1.sort()

        top2 = np.array([
            event_labels[top2_pred[0]],
            event_labels[top2_pred[1]],
            event_labels[top2_pred[2]],
            ])
        top2.sort()

        w1 = top1[1:]
        w2 = top2[1:]

        # Sanity check. Top 1 and top 2 should be different.
        if np.array_equal(top1, top2):
            print('Top 1 and top 2 are the same!')
            continue

        # Check if the top 1 and top 2 predictions are correct
        available_tops = [np.array([1,2,3]), np.array([4,5,6])]
        top1_correct, available_tops = check_if_top_was_correct(top1, available_tops)
        top2_correct, available_tops = check_if_top_was_correct(top2, available_tops)

        # Check if the W1 and W2 predictions are correct
        available_ws = [np.array([2,3]), np.array([5,6])]
        w1_correct, available_ws = check_if_top_was_correct(w1, available_ws)
        w2_correct, available_ws = check_if_top_was_correct(w2, available_ws)
        
        # Record results
        store_results(top_results_dict['total_events'], top_results_dict['two_correct'], top_results_dict['one_correct'], n_jets_per_event, top1_correct, top2_correct)
        store_results(w_results_dict['total_events'], w_results_dict['two_correct'], w_results_dict['one_correct'], n_jets_per_event, w1_correct, w2_correct)

    # Process the results.
    # For tops
    full_reco_eff = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    at_least_one_eff = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    # Calculate the efficiencies.
    for category in full_reco_eff:
        full_reco_eff[category] = 100*top_results_dict['two_correct'][category] / top_results_dict['total_events'][category]
        at_least_one_eff[category] = 100*(top_results_dict['two_correct'][category] + top_results_dict['one_correct'][category]) / top_results_dict['total_events'][category]
    # For Ws
    w_full_reco_eff = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    w_at_least_one_eff = {'6': 0, '7': 0, '8': 0, '9': 0, '10': 0, '10+': 0}
    # Calculate the efficiencies.
    for category in w_full_reco_eff:
        w_full_reco_eff[category] = 100*w_results_dict['two_correct'][category] / w_results_dict['total_events'][category]
        w_at_least_one_eff[category] = 100*(w_results_dict['two_correct'][category] + w_results_dict['one_correct'][category]) / w_results_dict['total_events'][category]


    # Print the results.
    print('Results:')
    print('Full reconstruction efficiency:')
    for category in full_reco_eff:
        print('For events with %s jets: %.3f' % (category, full_reco_eff[category]))

    print('At least one correct top efficiency:')
    for category in at_least_one_eff:
        print('For events with %s jets: %.3f' % (category, at_least_one_eff[category]))

    print('W Full reconstruction efficiency:')
    for category in w_full_reco_eff:
        print('For events with %s jets: %.3f' % (category, w_full_reco_eff[category]))
    
    print('W At least one correct top efficiency:')
    for category in w_at_least_one_eff:
        print('For events with %s jets: %.3f' % (category, w_at_least_one_eff[category]))

    # Plot the results in a histogram if requested
    if args.plot:
        x = np.arange(len(full_reco_eff))
        width = 0.35

        fig, ax = plt.subplots()
        rects1 = ax.bar(x , list(full_reco_eff.values()), width/4, label='Two tops reconstructed.')
        rects2 = ax.bar(x - width/2, list(at_least_one_eff.values()), width/4, label='At least one top reconstructed.')
        rects3 = ax.bar(x - width/4, list(w_full_reco_eff.values()), width/4, label='Two Ws reconstructed.')
        rects4 = ax.bar(x - 3*width/4 , list(w_at_least_one_eff.values()), width/4, label='At least one W reconstructed.')

        ax.set_ylabel('Efficiency')
        ax.set_xlabel('N jets')
        ax.set_title('KLFitter reconstruction')
        ax.set_xticks(x)
        ax.set_xticklabels(full_reco_eff.keys())
        ax.legend()

        post_fix = 'even' if 'even' in args.inputs else'odd'
        plt.savefig('KLF_PerformanceTtbarAllHad_%s.png' % post_fix)



def main():
    # Parse arguments.
    args = argparser()

    if args.klf_input != '':
        # Measure the performance of KLF.
        measure_klf_performance(args)
    else:
        # Measure the performance of HyPER.
        measure_hyper_performance(args)

if __name__ == '__main__':
    main()