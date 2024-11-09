#include "PartonMatcher/PartonMatcher.h"
#include "FastFrames/UniqueSampleID.h"
#include "PartonMatcher/Variables.h"

ROOT::RDF::RNode TtbarSingleLeptonFrame::defineVariables(ROOT::RDF::RNode mainNode,
                                              const std::shared_ptr<Sample>& /*sample*/,
                                              const UniqueSampleID& /*id*/) {

  // You can also use the UniqueSampleID object to apply a custom define
  // based on the sample:
  //   id.dsid() returns sample DSID
  //   id.campaign() returns sample campaign
  //   id.simulation() return simulation flavour
  // You can use it in your functions to apply only per sample define

  // Add all the variables needed for the analysisx
  LOG(INFO) << "You are running custom class: TtbarSingleLeptonFrame" << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Object multiplicity ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::string el_OR_flag = "el_select_passesOR_NOSYS"; // "el_select_or_NOSYS"
  std::string mu_OR_flag = "mu_select_passesOR_NOSYS"; // "mu_select_or_NOSYS"
  std::string jet_OR_flag = "jet_select_passesOR_NOSYS"; // "jet_select_or_NOSYS"

  LOG(INFO) << "Adding variable: nElectrons_tight_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nElectrons_tight_NOSYS",
                                         TtbarSingleLeptonVars::nGoodLeptons,
                                         {"el_pt_NOSYS", "el_select_tight_NOSYS",el_OR_flag});

  LOG(INFO) << "Adding variable: nMuons_tight_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nMuons_tight_NOSYS",
                                         TtbarSingleLeptonVars::nGoodLeptons,
                                         {"mu_pt_NOSYS", "mu_select_tight_NOSYS",mu_OR_flag});

  LOG(INFO) << "Adding variable: nJets_good_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nJets_good_NOSYS",
                                         TtbarSingleLeptonVars::nGoodJets,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS",jet_OR_flag});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Pt indices ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  LOG(INFO) << "Adding variable: jets_ptOrderIndex_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "jets_ptOrderIndex_NOSYS",
                                         TtbarSingleLeptonVars::indexPtOrderedGoodTwoSelec,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS",jet_OR_flag});  
  
  LOG(INFO) << "Adding variable: el_ptOrderIndex_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "el_ptOrderIndex_NOSYS",
                                         TtbarSingleLeptonVars::indexPtOrderedGoodTwoSelec,
                                         {"el_pt_NOSYS", "el_select_tight_NOSYS",el_OR_flag});

  LOG(INFO) << "Adding variable: mu_ptOrderIndex_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "mu_ptOrderIndex_NOSYS",
                                         TtbarSingleLeptonVars::indexPtOrderedGoodTwoSelec,
                                         {"mu_pt_NOSYS", "mu_select_tight_NOSYS",mu_OR_flag});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Jet KLFitter indices ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: topHad_indexKLF_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "topHad_indexKLF_NOSYS",
                                         TtbarSingleLeptonVars::fromKLFToNtupleIndex,
                                         {"b_from_tophad_jetIndex_NOSYS"});
    
  LOG(INFO) << "Adding variable: topLep_indexKLF_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "topLep_indexKLF_NOSYS",
                                         TtbarSingleLeptonVars::fromKLFToNtupleIndex,
                                         {"b_from_toplep_jetIndex_NOSYS"});

  LOG(INFO) << "Adding variable: W1_indexKLF_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "W1_indexKLF_NOSYS",
                                         TtbarSingleLeptonVars::fromKLFToNtupleIndex,
                                         {"jet1_from_W_jetIndex_NOSYS"});

  LOG(INFO) << "Adding variable: W2_indexKLF_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "W2_indexKLF_NOSYS",
                                         TtbarSingleLeptonVars::fromKLFToNtupleIndex,
                                         {"jet2_from_W_jetIndex_NOSYS"});                                                                                              
        
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Define the electron selection cuts ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Since datasets were produced in different runs and 2023 data has not yet good trigger SFs, etc we define this variable for this datasets.
  // Lambda funtion to return just true
  auto colNames = mainNode.GetColumnNames();
  if (std::find(colNames.begin(), colNames.end(), "globalTriggerMatch_NOSYS") == colNames.end()){
    auto returnTrue = [](const char& dummy){(void)(dummy); return true;};
    LOG(INFO) << "Adding variable: globalTriggerMatch_NOSYS\n";
    mainNode = MainFrame::systematicDefine(mainNode,
                                         "globalTriggerMatch_NOSYS",
                                         returnTrue,
                                         {"pass_ejets_NOSYS"});

    auto dummySF = [](const char& dummy){(void)(dummy); return 1.0f;};                             
    LOG(INFO) << "Adding variable: globalTriggerEffSF_NOSYS\n";
    mainNode = MainFrame::systematicDefine(mainNode,
                                         "globalTriggerEffSF_NOSYS",
                                         dummySF,
                                         {"pass_ejets_NOSYS"}); 
  }

  LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "Electron_pass_triggers_NOSYS",
                                         TtbarSingleLeptonVars::ElectronSelection::triggers,
                                         {"trigPassed_HLT_e60_lhmedium_L1EM22VHI", "trigPassed_HLT_e140_lhloose_L1EM22VHI","trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI",
                                         "trigPassed_HLT_e26_lhtight_nod0_ivarloose","trigPassed_HLT_e140_lhloose_nod0","trigPassed_HLT_e60_lhmedium_nod0"});

  LOG(INFO) << "Adding variable: Electron_pass_basic_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "Electron_pass_basic_NOSYS",
                                         TtbarSingleLeptonVars::ElectronSelection::basic,
                                         {"pass_ejets_NOSYS", "nElectrons_tight_NOSYS","nJets_good_NOSYS","Electron_pass_triggers_NOSYS","globalTriggerMatch_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Define the muon selection cuts ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "Muon_pass_triggers_NOSYS",
                                         TtbarSingleLeptonVars::MuonSelection::triggers,
                                         {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH", "trigPassed_HLT_mu50_L1MU14FCH",
                                         "trigPassed_HLT_mu50","trigPassed_HLT_mu26_ivarmedium"});

  LOG(INFO) << "Adding variable: Muon_pass_basic_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "Muon_pass_basic_NOSYS",
                                         TtbarSingleLeptonVars::MuonSelection::basic,
                                         {"pass_mujets_NOSYS", "nMuons_tight_NOSYS","nJets_good_NOSYS","Muon_pass_triggers_NOSYS","globalTriggerMatch_NOSYS"});

  return mainNode;
}

ROOT::RDF::RNode TtbarSingleLeptonFrame::defineVariablesNtuple(ROOT::RDF::RNode mainNode,
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
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
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
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
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
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_e_NOSYS","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Jet_GN2v01_score\n";
  mainNode = mainNode.Define("Jet_GN2v01_score",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_GN2v01","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Jet_SVmass\n";
  mainNode = mainNode.Define("Jet_SVmass",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_SV1_masssvx","jets_ptOrderIndex_NOSYS"});

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
                              TtbarSingleLeptonVars::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_down_index","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Truth_up_index\n";
  mainNode = mainNode.Define("Truth_up_index",
                              TtbarSingleLeptonVars::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_up_index","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Truth_had_b_index\n";
  mainNode = mainNode.Define("Truth_had_b_index",
                              TtbarSingleLeptonVars::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_had_b_index","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Truth_lep_b_index\n";
  mainNode = mainNode.Define("Truth_lep_b_index",
                              TtbarSingleLeptonVars::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_lep_b_index","jets_ptOrderIndex_NOSYS"});
                          

  return mainNode;
}

ROOT::RDF::RNode TtbarSingleLeptonFrame::defineVariablesTruth(ROOT::RDF::RNode node,
                                                   const std::string& treeName,
                                                   const std::shared_ptr<Sample>& /*sample*/,
                                                   const UniqueSampleID& /*sampleID*/) {

  if (treeName != "parton") return node; // do not do anything for other trees ("truth" = parton truth)

  return node;
}

ROOT::RDF::RNode TtbarSingleLeptonFrame::defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                           const std::string& /*treeName*/,
                                                           const std::shared_ptr<Sample>& /*sample*/,
                                                           const UniqueSampleID& /*sampleID*/) {
  return node;
}