import argparse
import uproot

def argparser():
    parser = argparse.ArgumentParser(description='Convert Subgrpah dataset to SPANet format')
    parser.add_argument('-f', '--input',  type=str, required=True, help='ROOT input file.')
    parser.add_argument('-o', '--output', type=str, help='Output .h5 file', required=False)
    parser.add_argument('-p', '--print-branches', action='store_true', help='Print branches in the input ROOT file.', required=False)
    parser.add_argument('-s', '--split', action='store_true', help='Split dataset in two according to (eventNumber mod 2)', required=False)

    args = parser.parse_args()

    # Check that when using -p only the input file is provided
    if args.print_branches and (args.output or args.split):
        raise ValueError('Print branches option is not compatible with output file or split option.')

    return args

def print_branches(root_file_name):
    tree_name = 'reco'
    root_file = uproot.open(root_file_name)
    tree = root_file[tree_name]
    branches = tree.keys()
    for branch in branches:
        print(branch)

def MakeDataset(make_dataset_function , input: str, output: str, split: bool = False):
    if not split:
        make_dataset_function(input, output)
    else:
        make_dataset_function(input, output.replace('.h5', '_even.h5'), split=True, saveEven=True)
        make_dataset_function(input, output.replace('.h5', '_odd.h5'), split=True, saveEven=False)

def main(make_dataset_function):
    args = argparser()

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
        
        if args.split:
            MakeDataset(make_dataset_function, args.input, args.output, True)
        else:
            MakeDataset(make_dataset_function, args.input, args.output)