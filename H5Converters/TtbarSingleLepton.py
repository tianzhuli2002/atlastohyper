import h5py
import numpy as np
from tqdm import tqdm
import uproot

from Common import argparser, main


def MakeDatasetTtbarSingleLepton(input: str, output: str, split: bool = False, split_value: int = 2):
    r"""Function to construct HyPER dataset.

    Args:
        input (str): input ROOT file.
        output (str): output HyPER dataset in .h5 format.
    """
    # Max number of objects (padding)
    pad_to_jet = 20
    # pad_to_el  = 2
    # pad_to_mu  = 2

    # Load the ROOT file
    root_file = uproot.open(input)
    tree = root_file['reco']
    num_entries = tree.num_entries

    # ---------------------------------------------------
    #                 Load data from ROOT
    # ---------------------------------------------------
    jet_bTag  = tree["GN2v01_PCBT_quantile_NOSYS"].array(library='np')
    jet_eta   = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fEta"].array(library='np')
    jet_phi   = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPhi"].array(library='np')
    jet_e     = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fE"].array(library='np')
    jet_pt    = tree["jet_Ordered_TLV_NOSYS/jet_Ordered_TLV_NOSYS.fCoordinates.fPt"].array(library='np')
    event_number = tree["eventNumber"].array(library='np')
    #jet_match = tree["jet_truthmatch"].array(library='np')
    # ......
    njet      = tree["nJets_good_NOSYS"].array(library='np')
    nbTagged  = tree["nBjets_GN2v01_77_NOSYS"].array(library='np')
    # ...... labels
    q1_index = tree["q1_index_NOSYS"].array(library='np')
    q2_index = tree["q2_index_NOSYS"].array(library='np')
    q1bar_index = tree["q1bar_index_NOSYS"].array(library='np')
    q2bar_index = tree["q2bar_index_NOSYS"].array(library='np')
    b_index = tree["b_index_NOSYS"].array(library='np')
    bbar_index = tree["bbar_index_NOSYS"].array(library='np')

    # ---------------------------------------------------

    num_entries = 200
    # ---------------------------------------------------
    #               Define Input Datatype
    # ---------------------------------------------------
    node_dt = np.dtype([('e', np.float32), ('eta', np.float32), ('phi', np.float32), ('pt', np.float32), ('btag', np.float32), ('id', np.float32)])
    jet_data = np.zeros((num_entries, pad_to_jet), dtype=node_dt)
    # el_data = np.zeros((num_entries, pad_to_el), dtype=node_dt)
    # mu_data = np.zeros((num_entries, pad_to_mu), dtype=node_dt)
    # ......
    global_dt = np.dtype([('njet', np.float32), ('nbTagged', np.float32)])
    global_data = np.zeros((num_entries, 1), dtype=global_dt)

    # We need `VertexID` to label truth egdes and hyperedges
    VertexID_data = np.full((num_entries, pad_to_jet), -9)  # use -9 for none filled values
    # VertexID_data = np.zeros((num_entries, pad_to_jet+pad_to_el+pad_to_mu))
    IndexSelect_data = np.zeros(num_entries)

    # ---------------------------------------------------
    #               Define HDF5 Structure
    # ---------------------------------------------------
    h5_file = h5py.File(output, 'w')
    inputs_group = h5_file.create_group('INPUTS')
    labels_group = h5_file.create_group('LABELS')

    for i in tqdm(range(num_entries), desc='Creating HDF5 dataset', unit='events'):
        num_jets = njet[i]
        for j in range(num_jets):
            jet_data[i][j] = (
                jet_e[i][j], jet_eta[i][j], jet_phi[i][j], jet_pt[i][j],
                (jet_bTag[i][j] >= 4 ), 1
            )
        # num_el = nel[i]
        # for j in range(num_el):
        #     el_data[i][j] = (
        #         el_e[i][j], el_eta[i][j], el_phi[i][j], el_pt[i][j],
        #         0, el_charge[i][j], -1*el_charge[i][j]   # Filling in 0 for the btag
        #     )
        # num_mu = nmu[i]
        # for j in range(num_mu):
        #     mu_data[i][j] = (
        #         mu_e[i][j], mu_eta[i][j], mu_phi[i][j], mu_pt[i][j],
        #         0, mu_charge[i][j], -2*mu_charge[i][j]   # Filling in 0 for the btag
        #     )

        # for j in range(num_nu):
        #     nu_data[i][j] = (
        #         nu_e[i][j], nu_eta[i][j], nu_phi[i][j], nu_pt[i][j]#, nu_pid[i][j]
        #     )
        global_data[i] = (njet[i], nbTagged[i])

        existent_jet_labels = np.array([0]*num_jets)
        existent_jet_labels[b_index[i]] = 1
        existent_jet_labels[q1_index[i]] = 2
        existent_jet_labels[q2_index[i]] = 3
        existent_jet_labels[bbar_index[i]] = 4
        existent_jet_labels[q1bar_index[i]] = 5
        existent_jet_labels[q2bar_index[i]] = 6
        VertexID_data[i, :num_jets] = existent_jet_labels
        # VertexID_data[i, :num_jets+num_el+num_mu] = np.concatenate((jet_match[i],el_match[i],mu_match[i]))
        IndexSelect_data[i] = True


    inputs_group.create_dataset("jet", data=jet_data)
    # inputs_group.create_dataset("electron", data=el_data)
    # inputs_group.create_dataset("muon", data=mu_data)
    inputs_group.create_dataset("global", data=global_data)
    labels_group.create_dataset("VertexID", data=VertexID_data)
    labels_group.create_dataset("IndexSelect", data=IndexSelect_data)


    # Close the HDF5 file
    h5_file.close()
    root_file.close()


if __name__ == "__main__":
    main(MakeDatasetTtbarSingleLepton)