#include "PartonMatcher/TtbarSingleLeptonFrameVyper.h"
#include "FastFrames/UniqueSampleID.h"
#include "PartonMatcher/Variables.h"

ROOT::RDF::RNode TtbarSingleLeptonFrameVyper::defineVariables(ROOT::RDF::RNode mainNode,
                                              const std::shared_ptr<Sample>& /*sample*/,
                                              const UniqueSampleID& /*id*/) {

  // You can also use the UniqueSampleID object to apply a custom define
  // based on the sample:
  //   id.dsid() returns sample DSID
  //   id.campaign() returns sample campaign
  //   id.simulation() return simulation flavour
  // You can use it in your functions to apply only per sample define

  // Add all the variables needed for the analysisx
  LOG(INFO) << "You are running custom class: TtbarSingleLeptonFrameVyper" << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Char to bool ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: globalTriggerMatch_bool_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "globalTriggerMatch_bool_NOSYS",
                                         Variables::charToBool,
                                         {"globalTriggerMatch_NOSYS"});  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Object multiplicity ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: nElectrons_tight_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nElectrons_tight_NOSYS",
                                         TtbarSingleLeptonVars::nGoodLeptons,
                                         {"el_pt_NOSYS", "el_select_tight_NOSYS", "el_select_tight_NOSYS"});

  LOG(INFO) << "Adding variable: nMuons_tight_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nMuons_tight_NOSYS",
                                         TtbarSingleLeptonVars::nGoodLeptons,
                                         {"mu_pt_NOSYS", "mu_select_tight_NOSYS", "mu_select_tight_NOSYS"});

  LOG(INFO) << "Adding variable: nJets_good_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nJets_good_NOSYS",
                                         TtbarSingleLeptonVars::nGoodJets,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS", "jet_select_baselineJvt_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Pt indices ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  LOG(INFO) << "Adding variable: jets_ptOrderIndex_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "jets_ptOrderIndex_NOSYS",
                                         Variables::indexPtOrderedGoodOneSelec,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS"});  
  
  LOG(INFO) << "Adding variable: el_ptOrderIndex_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "el_ptOrderIndex_NOSYS",
                                         Variables::indexPtOrderedGoodOneSelec,
                                         {"el_pt_NOSYS", "el_select_tight_NOSYS"});

  LOG(INFO) << "Adding variable: mu_ptOrderIndex_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "mu_ptOrderIndex_NOSYS",
                                         Variables::indexPtOrderedGoodOneSelec,
                                         {"mu_pt_NOSYS", "mu_select_tight_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Jet KLFitter indices ////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // LOG(INFO) << "Adding variable: bHad_indexKLF_NOSYS\n";
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "bHad_indexKLF_NOSYS",
  //                                        TtbarSingleLeptonVars::fromKLFToNtupleIndex,
  //                                        {"klfitter_bhad_jetIndex_NOSYS"});
    
  // LOG(INFO) << "Adding variable: bLep_indexKLF_NOSYS\n";
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "bLep_indexKLF_NOSYS",
  //                                        TtbarSingleLeptonVars::fromKLFToNtupleIndex,
  //                                        {"klfitter_blep_jetIndex_NOSYS"});

  // LOG(INFO) << "Adding variable: q1_indexKLF_NOSYS\n";
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "q1_indexKLF_NOSYS",
  //                                        TtbarSingleLeptonVars::fromKLFToNtupleIndex,
  //                                        {"klfitter_lq1_jetIndex_NOSYS"});

  // LOG(INFO) << "Adding variable: q2_indexKLF_NOSYS\n";
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "q2_indexKLF_NOSYS",
  //                                        TtbarSingleLeptonVars::fromKLFToNtupleIndex,
  //                                        {"klfitter_lq2_jetIndex_NOSYS"});                                                                                              
        
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Define the truth indices ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: parton_truth_down_index_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "parton_truth_down_index_NOSYS",
                                         TtbarSingleLeptonVars::getPartonTruthDown,
                                         {"event_jet_truth_idx"});
  
  LOG(INFO) << "Adding variable: parton_truth_up_index_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "parton_truth_up_index_NOSYS",
                                         TtbarSingleLeptonVars::getPartonTruthUp,
                                         {"event_jet_truth_idx"});

  LOG(INFO) << "Adding variable: parton_truth_had_b_index_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "parton_truth_had_b_index_NOSYS",
                                         TtbarSingleLeptonVars::getPartonTruthBHad,
                                         {"event_jet_truth_idx"});
  
  LOG(INFO) << "Adding variable: parton_truth_lep_b_index_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "parton_truth_lep_b_index_NOSYS",
                                         TtbarSingleLeptonVars::getPartonTruthBLep,
                                         {"event_jet_truth_idx"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Define the electron selection cuts ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Since datasets were produced in different runs and 2023 data has not yet good trigger SFs, etc we define this variable for this datasets.
  // Lambda funtion to return just true
  auto colNames = mainNode.GetColumnNames();
  if (std::find(colNames.begin(), colNames.end(), "globalTriggerMatch_NOSYS") == colNames.end()){
    auto returnTrue = [](const char& dummy){(void)(dummy); return true;};
    LOG(INFO) << "Adding variable: globalTriggerMatch_bool_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "globalTriggerMatch_bool_NOSYS",
                                         returnTrue,
                                         {"pass_ejets_NOSYS"});

    auto dummySF = [](const char& dummy){(void)(dummy); return 1.0f;};                             
    LOG(INFO) << "Adding variable: globalTriggerEffSF_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "globalTriggerEffSF_NOSYS",
                                         dummySF,
                                         {"pass_ejets_NOSYS"}); 
  }

  LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "Electron_pass_triggers_NOSYS",
                                         TtbarSingleLeptonVars::ElectronSelection::triggers,
                                         {"trigPassed_HLT_e60_lhmedium_L1EM22VHI", "trigPassed_HLT_e140_lhloose_L1EM22VHI","trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI",
                                         "trigPassed_HLT_e26_lhtight_nod0_ivarloose","trigPassed_HLT_e140_lhloose_nod0","trigPassed_HLT_e60_lhmedium_nod0"});

  LOG(INFO) << "Adding variable: Electron_pass_basic_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "Electron_pass_basic_NOSYS",
                                         TtbarSingleLeptonVars::ElectronSelection::basic,
                                         {"pass_ejets_NOSYS", "nElectrons_tight_NOSYS","nJets_good_NOSYS","Electron_pass_triggers_NOSYS","globalTriggerMatch_bool_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Define the muon selection cuts ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "Muon_pass_triggers_NOSYS",
                                         TtbarSingleLeptonVars::MuonSelection::triggers,
                                         {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH", "trigPassed_HLT_mu50_L1MU14FCH",
                                         "trigPassed_HLT_mu50","trigPassed_HLT_mu26_ivarmedium"});

  LOG(INFO) << "Adding variable: Muon_pass_basic_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "Muon_pass_basic_NOSYS",
                                         TtbarSingleLeptonVars::MuonSelection::basic,
                                         {"pass_mujets_NOSYS", "nMuons_tight_NOSYS","nJets_good_NOSYS","Muon_pass_triggers_NOSYS","globalTriggerMatch_bool_NOSYS"});

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Define the jets selection cuts ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: Jets_fullyMatched_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "Jets_fullyMatched_NOSYS",
                                         TtbarSingleLeptonVars::JetsSelection::fullyMatched,
                                         {"parton_truth_down_index_NOSYS", "parton_truth_up_index_NOSYS",
                                         "parton_truth_had_b_index_NOSYS","parton_truth_lep_b_index_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////// btag //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  LOG(INFO) << "Adding variable: jet_GN2v01_Continuous_quantile_NOSYS\n";   
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "jet_GN2v01_Continuous_quantile_NOSYS",
                                         Variables::getPCBT,
                                         {"jet_GN2v01_FixedCutBEff_65_select", "jet_GN2v01_FixedCutBEff_70_select", "jet_GN2v01_FixedCutBEff_77_select", "jet_GN2v01_FixedCutBEff_85_select", "jet_GN2v01_FixedCutBEff_90_select"});

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Truth Neutrinos ////////////////////////////////////////////////////////  

LOG(INFO) << "Adding variable: TruthNeutrino_NOSYS" << std::endl;

mainNode = m_partonMatcher->systematicDefine(
    mainNode,
    "TruthNeutrino_NOSYS",
    TtbarSingleLeptonVyper::ExtractTruthNeutrino,
    {
        // from top
        "Ttbar_MC_Wdecay1_afterFSR_from_t_pt",
        "Ttbar_MC_Wdecay1_afterFSR_from_t_eta",
        "Ttbar_MC_Wdecay1_afterFSR_from_t_phi",
        "Ttbar_MC_Wdecay1_afterFSR_from_t_m",
        "Ttbar_MC_Wdecay1_afterFSR_from_t_pdgId",

        "Ttbar_MC_Wdecay2_afterFSR_from_t_pt",
        "Ttbar_MC_Wdecay2_afterFSR_from_t_eta",
        "Ttbar_MC_Wdecay2_afterFSR_from_t_phi",
        "Ttbar_MC_Wdecay2_afterFSR_from_t_m",
        "Ttbar_MC_Wdecay2_afterFSR_from_t_pdgId",

        // from tbar
        "Ttbar_MC_Wdecay1_afterFSR_from_tbar_pt",
        "Ttbar_MC_Wdecay1_afterFSR_from_tbar_eta",
        "Ttbar_MC_Wdecay1_afterFSR_from_tbar_phi",
        "Ttbar_MC_Wdecay1_afterFSR_from_tbar_m",
        "Ttbar_MC_Wdecay1_afterFSR_from_tbar_pdgId",

        "Ttbar_MC_Wdecay2_afterFSR_from_tbar_pt",
        "Ttbar_MC_Wdecay2_afterFSR_from_tbar_eta",
        "Ttbar_MC_Wdecay2_afterFSR_from_tbar_phi",
        "Ttbar_MC_Wdecay2_afterFSR_from_tbar_m",
        "Ttbar_MC_Wdecay2_afterFSR_from_tbar_pdgId"
    }
);

using TtbarSingleLeptonVyper::TruthNuInfo;

mainNode = m_partonMatcher->systematicDefine(
    mainNode,
    "nu_pt_NOSYS",
    [](const TruthNuInfo& nu){ return static_cast<float>(nu.nu.Pt()); },
    {"TruthNeutrino_NOSYS"}
);

mainNode = m_partonMatcher->systematicDefine(
    mainNode,
    "nu_eta_NOSYS",
    [](const TruthNuInfo& nu){ return static_cast<float>(nu.nu.Eta()); },
    {"TruthNeutrino_NOSYS"}
);

mainNode = m_partonMatcher->systematicDefine(
    mainNode,
    "nu_phi_NOSYS",
    [](const TruthNuInfo& nu){ return static_cast<float>(nu.nu.Phi()); },
    {"TruthNeutrino_NOSYS"}
);

mainNode = m_partonMatcher->systematicDefine(
    mainNode,
    "nu_m_NOSYS",
    [](const TruthNuInfo& nu){ return static_cast<float>(nu.nu.M()); },
    {"TruthNeutrino_NOSYS"}
);

mainNode = m_partonMatcher->systematicDefine(
    mainNode,
    "nu_e_NOSYS",
    [](const TruthNuInfo& nu){ return static_cast<float>(nu.nu.E()); },
    {"TruthNeutrino_NOSYS"}
);

mainNode = m_partonMatcher->systematicDefine(
    mainNode,
    "nu_px_NOSYS",
    [](const TruthNuInfo& nu){ return static_cast<float>(nu.nu.Px()); },
    {"TruthNeutrino_NOSYS"}
);

mainNode = m_partonMatcher->systematicDefine(
    mainNode,
    "nu_py_NOSYS",
    [](const TruthNuInfo& nu){ return static_cast<float>(nu.nu.Py()); },
    {"TruthNeutrino_NOSYS"}
);
mainNode= m_partonMatcher->systematicDefine(
    mainNode,
    "nu_pz_NOSYS",
    [](const TruthNuInfo& nu){ return static_cast<float>(nu.nu.Pz()); },
    {"TruthNeutrino_NOSYS"}
);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////// HyPER results ////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  /*LOG(INFO) << "Adding variable: HyPER_TopHad_Indices_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_TopHad_Indices_NOSYS",
                              Variables::getSelecHyPERIndices,
                              { "hyper_TtbarLJets_TopHad_Indices_NOSYS",
                                "hyper_TtbarLJets_TopHad_IDs_NOSYS",
                                "jets_ptOrderIndex_NOSYS",
                                "el_ptOrderIndex_NOSYS",
                                "mu_ptOrderIndex_NOSYS",
                                "mu_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: HyPER_TopLep_Indices_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_TopLep_Indices_NOSYS",
                              Variables::getSelecHyPERIndices,
                              { "hyper_TtbarLJets_TopLep_Indices_NOSYS",
                                "hyper_TtbarLJets_TopLep_IDs_NOSYS",
                                "jets_ptOrderIndex_NOSYS",
                                "el_ptOrderIndex_NOSYS",
                                "mu_ptOrderIndex_NOSYS",
                                "mu_ptOrderIndex_NOSYS"}); //muon as placeholder for taus

  LOG(INFO) << "Adding variable: HyPER_WHad_Indices_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_WHad_Indices_NOSYS",
                              Variables::getSelecHyPERIndices,
                              { "hyper_TtbarLJets_WHad_Indices_NOSYS",
                                "hyper_TtbarLJets_TopHad_IDs_NOSYS",
                                "jets_ptOrderIndex_NOSYS",
                                "el_ptOrderIndex_NOSYS",
                                "mu_ptOrderIndex_NOSYS",
                                "mu_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: HyPER_WLep_Indices_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_WLep_Indices_NOSYS",
                              Variables::getSelecHyPERIndices,
                              { "hyper_TtbarLJets_WLep_Indices_NOSYS",
                                "hyper_TtbarLJets_TopLep_IDs_NOSYS",
                                "jets_ptOrderIndex_NOSYS",
                                "el_ptOrderIndex_NOSYS",
                                "mu_ptOrderIndex_NOSYS",
                                "mu_ptOrderIndex_NOSYS"});*/

  return mainNode;
}

ROOT::RDF::RNode TtbarSingleLeptonFrameVyper::defineVariablesNtuple(ROOT::RDF::RNode mainNode,
                                                    const std::shared_ptr<Sample>& sample,
                                                    const UniqueSampleID& id) {
  mainNode = this->defineVariables(mainNode, sample, id);

  // Define the collections to store.
  // Electrons
  LOG(INFO) << "Adding variable: Electron_pt\n";
  mainNode = mainNode.Define("Electron_pt",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"el_pt_NOSYS","el_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Electron_eta\n";
  mainNode = mainNode.Define("Electron_eta",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"el_eta","el_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Electron_phi\n";
  mainNode = mainNode.Define("Electron_phi",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"el_phi","el_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Electron_e\n";
  mainNode = mainNode.Define("Electron_e",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"el_e_NOSYS","el_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Electron_q\n";
  mainNode = mainNode.Define("Electron_q",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"el_charge","el_ptOrderIndex_NOSYS"});
  
  // Muons
  LOG(INFO) << "Adding variable: Muon_pt\n";
  mainNode = mainNode.Define("Muon_pt",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"mu_pt_NOSYS","mu_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Muon_eta\n";
  mainNode = mainNode.Define("Muon_eta",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"mu_eta","mu_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Muon_phi\n";
  mainNode = mainNode.Define("Muon_phi",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"mu_phi","mu_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Muon_e\n";
  mainNode = mainNode.Define("Muon_e",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"mu_e_NOSYS","mu_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Muon_q\n";
  mainNode = mainNode.Define("Muon_q",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"mu_charge","mu_ptOrderIndex_NOSYS"});

  // Jets
  LOG(INFO) << "Adding variable: Jet_pt\n";
  mainNode = mainNode.Define("Jet_pt",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"jet_pt_NOSYS","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Jet_eta\n";
  mainNode = mainNode.Define("Jet_eta",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_eta","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Jet_phi\n";
  mainNode = mainNode.Define("Jet_phi",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_phi","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Jet_e\n";
  mainNode = mainNode.Define("Jet_e",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"jet_e_NOSYS","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Jet_SVmass\n";
  LOG(INFO) << "Jet_SVmass ntuple branch disabled for this config" << std::endl;

  //btag  
  LOG(INFO) << "Adding variable: Jet_GN2v01_score\n";
  mainNode = mainNode.Define("Jet_GN2v01_score",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_GN2v01","jets_ptOrderIndex_NOSYS"});

  LOG(INFO) << "Adding variable: Jet_GN2v01_PCBT_quantile_NOSYS\n";   
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "Jet_GN2v01_PCBT_quantile_NOSYS",
                                         Variables::pTSortedContainer<int>,
                                         {"jet_GN2v01_Continuous_quantile_NOSYS", "jets_ptOrderIndex_NOSYS"});

  LOG(INFO) << "Adding variable: nBjets_GN2v01_100_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nBjets_GN2v01_100_NOSYS",
                                         Variables::nBJets(100),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile_NOSYS"});

  LOG(INFO) << "Adding variable: nBjets_GN2v01_90_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nBjets_GN2v01_90_NOSYS",
                                         Variables::nBJets(90),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile_NOSYS"});

  LOG(INFO) << "Adding variable: nBjets_GN2v01_85_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nBjets_GN2v01_85_NOSYS",
                                         Variables::nBJets(85),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile_NOSYS"});

  LOG(INFO) << "Adding variable: nBjets_GN2v01_77_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nBjets_GN2v01_77_NOSYS",
                                         Variables::nBJets(77),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile_NOSYS"});

  LOG(INFO) << "Adding variable: nBjets_GN2v01_70_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nBjets_GN2v01_70_NOSYS",
                                         Variables::nBJets(70),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile_NOSYS"});
  
  LOG(INFO) << "Adding variable: nBjets_GN2v01_65_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nBjets_GN2v01_65_NOSYS",
                                         Variables::nBJets(65),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile_NOSYS"});
  

  // MET
  LOG(INFO) << "Adding variable: Met_pt\n";
  mainNode = mainNode.Define("Met_pt",
                              [](float met_pt) { return met_pt/1000.f; },
                              {"met_met_NOSYS"});
  LOG(INFO) << "Adding variable: Met_phi\n";
  mainNode = mainNode.Define("Met_phi",
                              [](float met_phi) { return met_phi/1.f; },
                              {"met_phi_NOSYS"});

  // Parton labels
  LOG(INFO) << "Adding variable: Truth_down_index\n";
  mainNode = mainNode.Define("Truth_down_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_down_index_NOSYS","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Truth_up_index\n";
  mainNode = mainNode.Define("Truth_up_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_up_index_NOSYS","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Truth_had_b_index\n";
  mainNode = mainNode.Define("Truth_had_b_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_had_b_index_NOSYS","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Truth_lep_b_index\n";
  mainNode = mainNode.Define("Truth_lep_b_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_lep_b_index_NOSYS","jets_ptOrderIndex_NOSYS"});

  // // KLFitter results
  // LOG(INFO) << "Adding variable: KLF_lq1_index\n";
  // mainNode = mainNode.Define("KLF_lq1_index",
  //                             Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
  //                             {"q1_indexKLF_NOSYS","jets_ptOrderIndex_NOSYS"});                      
  // LOG(INFO) << "Adding variable: KLF_lq2_index\n";
  // mainNode = mainNode.Define("KLF_lq2_index",
  //                             Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
  //                             {"q2_indexKLF_NOSYS","jets_ptOrderIndex_NOSYS"}); 
  // LOG(INFO) << "Adding variable: KLF_had_b_index\n";
  // mainNode = mainNode.Define("KLF_had_b_index",
  //                             Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
  //                             {"bHad_indexKLF_NOSYS","jets_ptOrderIndex_NOSYS"}); 
  // LOG(INFO) << "Adding variable: KLF_lep_b_index\n";
  // mainNode = mainNode.Define("KLF_lep_b_index",
  //                             Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
  //                             {"bLep_indexKLF_NOSYS","jets_ptOrderIndex_NOSYS"});


  LOG(INFO) << "Adding variable: nu_px/\n";
  mainNode = mainNode.Define("nu_px",
                              [](float nu_px) { return nu_px/1000.f; },
                              {"nu_px_NOSYS"});
                              
  LOG(INFO) << "Adding variable: nu_py/\n";
  mainNode = mainNode.Define("nu_py",
                              [](float nu_py) { return nu_py/1000.f; },
                              {"nu_py_NOSYS"});

  LOG(INFO) << "Adding variable: nu_pz/\n";
  mainNode = mainNode.Define("nu_pz",
                              [](float nu_pz) { return nu_pz/1000.f; },
                              {"nu_pz_NOSYS"});
    
  // HyPER results
  /*LOG(INFO) << "Adding variable: HyPER_TopHad_IDs\n";
  mainNode = mainNode.Define("HyPER_TopHad_IDs",
                              [](const std::vector<long>& topHadIDs) { return topHadIDs; },
                              {"hyper_TtbarLJets_TopHad_IDs_NOSYS"});
  LOG(INFO) << "Adding variable: HyPER_TopLep_IDs\n";
  mainNode = mainNode.Define("HyPER_TopLep_IDs",
                              [](const std::vector<long>& topLepIDs) { return topLepIDs; },
                              {"hyper_TtbarLJets_TopLep_IDs_NOSYS"});

  LOG(INFO) << "Adding variable: HyPER_TopHad_Indices\n";
  mainNode = mainNode.Define("HyPER_TopHad_Indices",
                              [](const std::vector<long>& topHadIndices) { return topHadIndices; },
                              {"HyPER_TopHad_Indices_NOSYS"});
  LOG(INFO) << "Adding variable: HyPER_TopLep_Indices\n";
  mainNode = mainNode.Define("HyPER_TopLep_Indices",
                              [](const std::vector<long>& topLepIndices) { return topLepIndices; },
                              {"HyPER_TopLep_Indices_NOSYS"});
  LOG(INFO) << "Adding variable: HyPER_WHad_Indices\n";
  mainNode = mainNode.Define("HyPER_WHad_Indices",
                              [](const std::vector<long>& WHadIndices) { return WHadIndices; },
                              {"HyPER_WHad_Indices_NOSYS"});
  LOG(INFO) << "Adding variable: HyPER_WLep_Indices\n";
  mainNode = mainNode.Define("HyPER_WLep_Indices",
                              [](const std::vector<long>& WLepIndices) { return WLepIndices; },
                              {"HyPER_WLep_Indices_NOSYS"});*/

  return mainNode;
}

ROOT::RDF::RNode TtbarSingleLeptonFrameVyper::defineVariablesTruth(ROOT::RDF::RNode node,
                                                   const std::string& treeName,
                                                   const std::shared_ptr<Sample>& /*sample*/,
                                                   const UniqueSampleID& /*sampleID*/) {

  if (treeName != "parton") return node; // do not do anything for other trees ("truth" = parton truth)

  return node;
}

ROOT::RDF::RNode TtbarSingleLeptonFrameVyper::defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                           const std::string& /*treeName*/,
                                                           const std::shared_ptr<Sample>& /*sample*/,
                                                           const UniqueSampleID& /*sampleID*/) {
  return node;
}