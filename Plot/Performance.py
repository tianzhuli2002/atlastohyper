import matplotlib.pyplot as plt
import numpy as np
import ROOT as r

# Results storage
N_Correct_Tops = {"KLF" : {}, "HyPER": {}}
N_Correct_Ws = {"KLF" : {}, "HyPER": {}}
categories = ['6Jets', '7Jets', '8Jets', '9Jets', '10pJets']
categories = [i+'_Even' for i in categories]

# External parameters
file_name = "../../resultsTopRecoAllHadronic/ttbar_allhad_run3.root"
correct_tops_histo_name = "NOSYS/N_corect_Top"
correct_ws_histo_name = "NOSYS/N_corect_W"

def fill_results(file_name, histogram_to_scan, results_dict, algorithm_name):
    file = r.TFile(file_name)
    for category in categories:
        hist = file.Get("%s_%s_%s" % (histogram_to_scan, algorithm_name, category))
        results_dict[algorithm_name][category] = {-1 : round(hist.GetBinContent(2),0), 0 : round(hist.GetBinContent(3),0), 1 : round(hist.GetBinContent(4),0), 2 : round(hist.GetBinContent(5),0)}
    file.Close()

    return results_dict

N_Correct_Tops = fill_results(file_name, correct_tops_histo_name, N_Correct_Tops, "KLF")
N_Correct_Tops = fill_results(file_name, correct_tops_histo_name, N_Correct_Tops, "HyPER")
N_Correct_Ws = fill_results(file_name, correct_ws_histo_name, N_Correct_Ws, "KLF")
N_Correct_Ws = fill_results(file_name, correct_ws_histo_name, N_Correct_Ws, "HyPER")

def sum_indices(dictionary, indices):
    s = 0
    for index in indices:
        s += dictionary[index]
    return s

# Create datasets you want to plot.
Two_Tops_Reconstructed = {"KLF" : [N_Correct_Tops["KLF"][category][2] for category in categories],
                          "HyPER" : [N_Correct_Tops["HyPER"][category][2] for category in categories]}
Reco_Performance = {"KLF" : [100*(N_Correct_Tops["KLF"][category][2]/(sum_indices(N_Correct_Tops["KLF"][category],[-1,0,1,2]))) for category in categories],
                    "HyPER" : [100*(N_Correct_Tops["HyPER"][category][2]/(sum_indices(N_Correct_Tops["HyPER"][category],[-1,0,1,2]))) for category in categories]}
Two_Ws_Reconstructed = {"KLF" : [N_Correct_Ws["KLF"][category][2] for category in categories],
                          "HyPER" : [N_Correct_Ws["HyPER"][category][2] for category in categories]}

print(Reco_Performance)

def make_plot(data_dict, plot_name, plot_file_name, y_label):
    species = ("6 jets", "7 jets", "8 jets","9 jets","10+ jets")
    categories = data_dict

    x = np.arange(len(species))  # the label locations
    width = 0.25  # the width of the bars
    multiplier = 0

    fig, ax = plt.subplots(layout='constrained')

    for attribute, measurement in categories.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute)
        ax.bar_label(rects, padding=3)
        multiplier += 1

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(y_label)
    ax.set_title(plot_name)
    ax.set_xticks(x + width, species)
    ax.legend(loc='upper right', ncols=3)
    #ax.set_ylim(0, 250)

    plt.savefig("%s.png" % plot_file_name)


make_plot(Two_Tops_Reconstructed, "Two tops correctly reconstructed" ,"Two_Tops_Reconstructed_Run3_Even", "N Events")
make_plot(Two_Ws_Reconstructed, "Two Ws correctly reconstructed" ,"Two_Ws_Reconstructed_Run3_Even", "N Events")
make_plot(Reco_Performance, "Full reconstruction performance" ,"Reco_Performance_Run3_Even", "Percentage of events")
