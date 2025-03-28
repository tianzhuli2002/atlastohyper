import yaml
import os
import argparse
import ROOT
from Logger import INFO, WARNING, ERROR, DEBUG, Logger

from HistogramHelpers import stackPlot, templatesDict

def get_template(dict_name : str):
    if dict_name not in templatesDict:
        ERROR.log('Dictionary %s not found... do you need to define it in HistogramHelpers.py ? ' % dict_name)
        exit(1)
    return templatesDict[dict_name]

def check_file_and_move_dir(file : str):
    if not os.path.isfile(file):
        print('File %s not found' % file)
        exit(1)
    os.chdir(os.path.dirname(file) if os.path.dirname(file) != '' else '.')

def load_config(config_file : str):
    check_file_and_move_dir(config_file) # Check the file exists and move to the directory containing it.
    config_file_name = os.path.basename(config_file)
    with open(config_file_name) as file:
        INFO.log('Loading configuration file ----> %s' % config_file)
        config = yaml.load(file, Loader=yaml.FullLoader)
    return config

def print_config(config, indent = ''):
    if type(config) == dict:
        for key, value in config.items():
            if type(value) == list:
                INFO.log(key+' ----> ')
                for i in value:
                    print_config(i,indent+'  ')
                    INFO.log('-------------------------------')
            else:
                INFO.log(indent+str(key) + '  ----> ' + str(value))
    else:
        INFO.log(indent+str(config))

def verify():
    print('')
    INFO.log('Run this configuration? [y/n]')
    answer = input().lower()

    while answer != 'y' and answer != 'n':
        WARNING.log('Please enter y or n')
        answer = input().lower()

    return answer == 'y'

def calculate_color_operations(color_operations : str):
    operators = []
    numbers = []
    for i in color_operations:
        if i in ['+','-']:
            operators.append(i)
        elif i.isdigit():
            numbers.append(int(i))
        else:
            ERROR.log('Invalid string %s in color.' % i)
            exit(1)

    # Validate before doing calculation.
    if len(operators) != len(numbers):
        ERROR.log('Different number of operators and numbers. Make sure %s is a valid operation.' % color_operations)
        exit(1)

    # Calculate the number to add to the color.
    number_to_add = 0
    for i in range(len(operators)):
        if operators[i] == '+':
            number_to_add += numbers[i]
        elif operators[i] == '-':
            number_to_add -= numbers[i]

    return number_to_add


def color_to_ROOT(color_input : str):
    color_dict = { 'kAzure' : ROOT.kAzure,
                  'kGreen' : ROOT.kGreen,
                    'kRed' : ROOT.kRed,
                    'kYellow' : ROOT.kYellow,
                    'kCyan' : ROOT.kCyan,
                    'kViolet' : ROOT.kViolet,
                    'kTeal' : ROOT.kTeal,}
    color_found = False
    for predef_color in color_dict:
        if predef_color in color_input:
            root_color = color_dict[predef_color]
            color_found = True
            break
    if not color_found:
        ERROR.log('Color %s not found in predefined dict! Add it!' % color_input)
        exit(1)
    
    color_operations = color_input.replace(predef_color,'')

    # Just a color inputed, return it.
    if color_operations == '':
        return root_color
    
    # Calculate the color.
    return root_color + calculate_color_operations(color_operations)

def build_sample_pack(samples : list):    
    sample_pack = {}
    for sample_dict in samples:
        try:
            sample_name = sample_dict['name']
            sample_file = sample_dict['file']
            sample_color = color_to_ROOT(sample_dict['color'])
        except KeyError as e:
            if e == 'name':
                ERROR.log('Name not given to one of the background samples. Check!')
            else:
                ERROR.log('Key % s not found in background sample: %s' % (e, sample_name))
            exit(1)

        sample_pack[sample_name] = [sample_file, sample_color, 0]

    return sample_pack

def build_zprime_pack(zprime_pack : dict):
    if zprime_pack is None:
        return None
    
    zp_pack = build_background(zprime_pack, 'dummy')
    zp_pack.pop('QCDjj') # Remove the dummy QCDjj

    return zp_pack

def generate_watermark(after_fit : bool):
    if after_fit:
        return '_AfterFit'
    return ''

def get_args():
    parser = argparse.ArgumentParser(description='Plot the histograms.')
    parser.add_argument('config_file', type=str, help='The configuration file to use.')
    parser.add_argument('--debug', action='store_true', help='Set the log level to debug.')
    return parser.parse_args()

def Plot():
    # Get arguments and set the log level.
    args = get_args()
    Logger.LOGLEVEL = 3
    if args.debug:
        Logger.LOGLEVEL = 4
    # Get the configuration file.
    config = load_config(args.config_file)
    print_config(config)
    run = verify()
    if not run:
        INFO.log('Exiting... without running!')
        exit(0)
    INFO.log('Running...')
    # Set the configuration.
    # First the always needed parameters.
    try:
        SIGNAL = build_sample_pack(config['signal'])
        BACKGROUND = build_sample_pack(config['background'])
        HISTOGRAMS = get_template(config['histograms'])
        APPLY_SF = True if bool(config['is_post_fit']) else False
        FINAL_STATE = config['final_state']
        REGION = config['region']
        WATERMARK = generate_watermark(config['is_post_fit'])
        PRINT_VERSION = config['print_version']
    except KeyError as e:
        ERROR.log('Key %s not found in configuration file.' % e)
        exit(1)
    # Now the optional parameters.
    optional_defaults = {'data':'NO_DATA',
                         'blind':True,
                         'unblind_purity_limit':0.0,
                         'print_overflows':False,
                         'purity_multiplier':1.0}
    for key in optional_defaults:
        if key not in config:
            config[key] = optional_defaults[key]
    
    DATA = {"Data":[config['data'],ROOT.kBlack,0]} if config['data'] != 'NO_DATA' else 'NO_DATA'
    BLIND = config['blind']
    UNBLIND_PURITY_LIMIT = config['unblind_purity_limit']
    PRINT_OVERFLOWS = config['print_overflows']
    PURITY_MULTIPLIER = config['purity_multiplier']

    # Print the configuration for debugging.
    DEBUG.log('DATA ----> %s' % DATA)
    DEBUG.log('SIGNAL ----> %s' % SIGNAL)
    DEBUG.log('BACKGROUND ----> %s' % BACKGROUND)
    DEBUG.log('HISTOGRAMS ----> %s' % HISTOGRAMS)
    DEBUG.log('FINAL_STATE ----> %s' % FINAL_STATE)
    DEBUG.log('REGION ----> %s' % REGION)
    DEBUG.log('WATERMARK ----> %s' % WATERMARK)
    DEBUG.log('PRINT_VERSION ----> %s' % PRINT_VERSION)
    DEBUG.log('BLIND ----> %s' % BLIND)
    DEBUG.log('UNBLIND_PURITY_LIMIT ----> %s' % UNBLIND_PURITY_LIMIT)
    DEBUG.log('PRINT_OVERFLOWS ----> %s' % PRINT_OVERFLOWS)
    DEBUG.log('PURITY_MULTIPLIER ----> %s' % PURITY_MULTIPLIER)

    # Run the plotting.
    
    stackPlot(data = DATA,
              signal = SIGNAL,
              background = BACKGROUND,
              histograms = HISTOGRAMS,
              watermark = WATERMARK,
              scales_dict = None,
              final_state = FINAL_STATE,
              regionLabel = REGION,
              blind = BLIND,
              unblindPurityLimit = UNBLIND_PURITY_LIMIT,
              printVersion = PRINT_VERSION,
              printOverflows = PRINT_OVERFLOWS,
              purityMultiplier = PURITY_MULTIPLIER
              )
    
    
if __name__ == '__main__':
    Plot()