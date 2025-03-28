import argparse
import uproot

class LoopStorage:
    def __init__(self):
        self.tree = None
        self.n_events_in_tree = 0
        self.n_events_to_save = 0
        self.whole_array = []
        self.even_array = []
        self.odd_array = []

def argparser():

    parser = argparse.ArgumentParser(description='Convert ROOT dataset to HyPER format')
    parser.add_argument('-f', '--input',  type=str, required=True, help='ROOT input file.')
    parser.add_argument('-o', '--output', type=str, help='Output .h5 file', required=False)
    parser.add_argument('-p', '--print-branches', action='store_true', help='Print branches in the input ROOT file.', required=False)
    parser.add_argument('-s', '--split', action='store_true', help='Split dataset in two according to (eventNumber mod 2)', required=False)

    
    return parser

def print_branches(root_file_name):
    tree_name = 'reco'
    root_file = uproot.open(root_file_name)
    tree = root_file[tree_name]
    branches = tree.keys()
    for branch in branches:
        print(branch)

def MakeDataset(make_dataset_function, args, loop_storage):
    output = args.output
    if not args.split:
        make_dataset_function(args, loop_storage, output)
    else:
        make_dataset_function(args, loop_storage, output.replace('.h5', '_even.h5'), split=True, saveEven=True)
        make_dataset_function(args, loop_storage, output.replace('.h5', '_odd.h5'), split=True, saveEven=False)

def main(make_dataset_function, args, loop_storage = None):

    # Check that when using -p only the input file is provided
    if args.print_branches and (args.output or args.split):
        raise ValueError('Print branches option is not compatible with output file or split option.')

    if args.print_branches:
        print_branches(args.input)
    else:
        # Validate output name
        if args.output is None:
            raise ValueError('Output file is required.')
            exit(1)
        if args.output.split('.')[-1] != 'h5':
            raise ValueError('Output file must be in .h5 format.')
            exit(1)

        MakeDataset(make_dataset_function, args, loop_storage)