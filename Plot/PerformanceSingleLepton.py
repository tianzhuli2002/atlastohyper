import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import argparse
import copy
import ROOT as r

def load_results(root_file , categories : list[str], histogram_names : list[str], algorith_name: str = 'HyPER'):

    # First, create the dictionary templates
    base_dictionary = {'Bad': 0.0, 'Good': 0.0, 'BadError': 0.0, 'GoodError': 0.0}
    results_dictionary_one_category = {'TopLep': copy.deepcopy(base_dictionary),
                                       'TopHad' : copy.deepcopy(base_dictionary),
                                       'WHad': copy.deepcopy(base_dictionary),
                                       'TwoTop': copy.deepcopy(base_dictionary)}
    results_dictionary = {category: copy.deepcopy(results_dictionary_one_category) for category in categories}

    # Adaptor dictionaries from plotting script to ROOT data.
    categories_adaptor_dictionary = {'4 jets': '4Jets', '5 jets': '5Jets', '6 jets': '6Jets', '7 jets': '7Jets', '8 jets': '8Jets', '9+ jets': '9pJets'}
    histograms_adaptor_dictionary = {'TopLep_Correct': 'TopLep', 'TopHad_Correct': 'TopHad', 'WHad_Correct': 'WHad', 'TwoTop_Correct': 'TwoTop'}

    # Load the results from the ROOT fil
    
    file = r.TFile(root_file)
    
    def get_histogram_name(histo_category, category, algorithm_name):
        root_file_dir = 'NOSYS/%s/' % (categories_adaptor_dictionary[category])
        return '%s%s_%s_%s' % (root_file_dir,histo_category, algorithm_name, categories_adaptor_dictionary[category])
    
    for cat in categories:
        for histo in histogram_names:
            # Get the histogram
            histo_name = get_histogram_name(histo, cat, algorith_name)
            histo_obj = file.Get(histo_name)
            
            # Check if the histogram exists
            if not histo_obj:
                print(f"Histogram {histo_name} not found in {root_file}")
                continue
            
            # Get the values and errors from the histogram
            good_bin_value = histo_obj.GetBinContent(2)
            bad_bin_value = histo_obj.GetBinContent(1)
            good_bin_error = histo_obj.GetBinError(2)
            bad_bin_error = histo_obj.GetBinError(1)
            
            

            # Store the values in the dictionary
            histo_key = histograms_adaptor_dictionary[histo]

            
            results_dictionary[cat][histo_key]['Good'] = good_bin_value
            results_dictionary[cat][histo_key]['Bad'] = bad_bin_value
            results_dictionary[cat][histo_key]['GoodError'] = good_bin_error
            results_dictionary[cat][histo_key]['BadError'] = bad_bin_error

    file.Close()

    return results_dictionary

def calculate_efficiencies(results_dictionary : dict[str, dict[str, dict[str,float]]]):
    eff_dict = {}
    for n_jets, reco_cat in results_dictionary.items():
        eff_dict[n_jets] = {i: {} for i in reco_cat.keys()}
        for reco_type, values in reco_cat.items():
            good = values['Good']
            bad = values['Bad']
            # good_error = values['GoodError']
            # bad_error = values['BadError']

            # Calculate the overall efficiency
            overall_efficiency = good / (good + bad)
            overall_efficiency_error = 100 * np.sqrt(overall_efficiency*(1-overall_efficiency)/(good+bad))
            overall_efficiency = 100 * overall_efficiency

            # Store the results in a dictionary
            eff_dict[n_jets][reco_type] = {
                'Efficiency': overall_efficiency,
                'EfficiencyError': overall_efficiency_error
            }
            
    return eff_dict

def plot_eff(efficiencies_dict_hyper, efficiencies_dict_klf, run_nb, btagged, categories, key_to_plot = 'TopLep'):
    if btagged:
        tag = "btag"
    else:
        tag = "nobtag"
    n_categories = len(categories)

    HyPER_overall_eff = [round(efficiencies_dict_hyper[n_jets][key_to_plot]['Efficiency'],2) for n_jets in categories]
    KLFitter_overall_eff = [round(efficiencies_dict_klf[n_jets][key_to_plot]['Efficiency'],2) for n_jets in categories]
    HyPER_overall_eff_err = [round(efficiencies_dict_hyper[n_jets][key_to_plot]['EfficiencyError'],2) for n_jets in categories]
    KLFitter_overall_eff_err = [round(efficiencies_dict_klf[n_jets][key_to_plot]['EfficiencyError'],2) for n_jets in categories]

    myAlpha = 0.3 #transparency of the error bars

    xerrors = [0.25] * n_categories

    x = np.arange(n_categories) + 0.5

    fig, ax = plt.subplots(figsize=(9, 6))

    for i in range(len(x)):
        # Shade the area between the x-error limits
        ax.fill_betweenx([HyPER_overall_eff[i] - HyPER_overall_eff_err[i], HyPER_overall_eff[i] + HyPER_overall_eff_err[i]], x[i] - 0.25 - xerrors[i], x[i] - 0.25 + xerrors[i], color='blue', alpha=myAlpha)
        ax.fill_betweenx([KLFitter_overall_eff[i] - KLFitter_overall_eff_err[i], KLFitter_overall_eff[i] + KLFitter_overall_eff_err[i]], x[i] + 0.25 - xerrors[i], x[i] + 0.25 + xerrors[i], color='red', alpha=myAlpha)
        #pass

    HyPER_plot    = ax.errorbar(x - 0.25, HyPER_overall_eff, xerr=xerrors, fmt=' ', color='blue', ecolor='blue', elinewidth=2, capsize=0, label='HyPER')
    KLFitter_plot = ax.errorbar(x + 0.25, KLFitter_overall_eff, xerr=xerrors, fmt=' ', color='red', ecolor='red', elinewidth=2, capsize=0, label='KLFitter')

    ax.set_ylabel('Reconstruction efficiency [%]', fontsize=14)
    plt.yticks(fontsize=12)
    ax.set_ylim(0, 130)

    ax.set_xlim(0., n_categories)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=14)

    if btagged:
        ax.legend(loc='lower left', bbox_to_anchor=(0., 0.8), fontsize=14, frameon=True, edgecolor='black',)
    else:
        ax.legend(loc='lower left', fontsize=14, frameon=True, edgecolor='black',)

    # Add value above each point
    for i, v in enumerate(HyPER_overall_eff):
        ax.text(x[i] - 0.25, v + 2, f"{v:.2f}", ha='center', color='black', fontsize=12)
    for i, v in enumerate(KLFitter_overall_eff):
        ax.text(x[i] + 0.25, v + 2, f"{v:.2f}", ha='center', color='black', fontsize=12)

    ax.grid(True, which='both', linestyle='--', linewidth=0.6, axis='y')  # Dashed lines for y-axis
    for pos in np.arange(0., n_categories, 1.0):
        ax.axvline(pos, color='gray', alpha=0.7, linestyle='-', linewidth=0.5)

    # Place error bar and dataset information within the plot, above the x-axis
    ax.text(3., 120., 'Training dataset:\nEvaluated dataset:', 
            ha='center', va='top', fontsize=14, fontweight='bold')
    ax.text(5., 120, 'All fully matched events\nAll fully matched events', 
            ha='center', va='top', fontsize=14)
    ax.text(4., 105, 'error bars multiplied by 10', ha='center', fontsize=12, fontstyle='italic')

    ax.set_title(f'Run {run_nb}, {key_to_plot}, ttbar SL', fontsize=16, fontweight='bold')

    # Save the plot as an image file
    plt.tight_layout()
    plt.savefig('Run%d_%s.png'% (run_nb, key_to_plot) )
    plt.close()

#Plot HyPER and KLFitter overall efficiencies as a function of number of jets
def main(args):

    run_nb  = args.run
    btagged = args.wasbtagged
    eff_to_plot = args.eff_to_plot
    
    categories = ['4 jets', '5 jets', '6 jets', '7 jets', '8 jets', '9+ jets']
    
    # Load the results...
    base_path = '/home/dbaronmo/Work/resultsTopRecoSingleLepton/'
    root_file = 'ttbar_nonallhad_run%d.root' % (run_nb)
    histograms = ['TopLep_Correct', 'TopHad_Correct','WHad_Correct','TwoTop_Correct']

    # HyPER results
    hyper_results = load_results(base_path + root_file, categories, histograms)
    hyper_efficiencies = calculate_efficiencies(hyper_results)
    # KLFitter results
    klf_results = load_results(base_path + root_file, categories, histograms, 'KLF')
    klf_efficiencies = calculate_efficiencies(klf_results)

    plot_eff(hyper_efficiencies, klf_efficiencies, run_nb, btagged, categories, eff_to_plot)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot HyPER and KLFitter comparison')

    parser.add_argument('run', type=int, help='Run number (2 or 3)')
    parser.add_argument('wasbtagged', type=lambda x: x.lower() == 'true', help='Did HyPER have access to btagging info during the training')
    parser.add_argument('eff_to_plot', type=str, help='Efficiency to plot (TopLep, TopHad, WHad)')

    args = parser.parse_args()

    main(args)