#include "PartonMatcher/ZllCalibrationFrame.h"
#include "FastFrames/UniqueSampleID.h"
#include "PartonMatcher/Variables.h"
#include "FastFrames/DefineHelpers.h"
#include <algorithm>
#include <vector>

ROOT::RDF::RNode ZllCalibrationFrame::defineVariables(ROOT::RDF::RNode mainNode,
                                              const std::shared_ptr<Sample>& sample,
                                              const UniqueSampleID& id) {

  // Add all the variables needed for the analysisx
  LOG(INFO) << "You are running custom class: ZllCalibrationFrame" << std::endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Char to bool ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // LOG(INFO) << "Adding variable: globalTriggerMatch_bool_NOSYS\n";
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "globalTriggerMatch_bool_NOSYS",
  //                                        Variables::charToBool,
  //                                        {"globalTriggerMatch_NOSYS"});  
                   
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Object multiplicity ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto numberOfJets30 = [](const ROOT::VecOps::RVec<float>& ptV,
    const ROOT::VecOps::RVec<char>& selection) {
      return DefineHelpers::numberOfObjects(ptV, 0, selection);
  };
            
auto numberOfBJets30 = []( const ROOT::VecOps::RVec<TLV>& tlv){
      std::size_t nBJets = 0;
      for (const auto& jet : tlv) {
          if (jet.Pt() > 30000) nBJets++;
      }
      return nBJets;
  };
auto numberOfLeptons25 = [](const ROOT::VecOps::RVec<float>& ptV,
    const ROOT::VecOps::RVec<char>& selection) {
      return DefineHelpers::numberOfObjects(ptV, 25000, selection);
  };


  // Muons

LOG(INFO) << "Adding variable: n_muons_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "n_muons_NOSYS",
                                         numberOfLeptons25,
                                         {"mu_pt_NOSYS", "mu_select_tight_NOSYS"});

  // Electrons
  LOG(INFO) << "Adding variable: n_electrons_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "n_electrons_NOSYS",
                                         numberOfLeptons25,
                                         {"el_pt_NOSYS", "el_select_tight_NOSYS"});
  LOG(INFO) << "Adding variable: n_jets_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "n_jets_NOSYS",
                                         numberOfJets30,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS"});  

  LOG(INFO) << "Adding variable: nElectrons_tight_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nElectrons_tight_NOSYS",
                                         numberOfLeptons25,
                                         {"el_pt_NOSYS", "el_select_tight_NOSYS"});

  LOG(INFO) << "Adding variable: nMuons_tight_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nMuons_tight_NOSYS",
                                         numberOfLeptons25,
                                         {"mu_pt_NOSYS", "mu_select_tight_NOSYS"});

  LOG(INFO) << "Adding variable: nJets_good_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "nJets_good_NOSYS",
                                         numberOfJets30,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS"});

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
  ////////////////////////////////////////////////////////////// pT sorted objects /////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Specify the type for the pT sorting functions
  using sorted_particle_1sel = ROOT::VecOps::RVec<TLV>(*)(const ROOT::VecOps::RVec<TLV>&,const ROOT::VecOps::RVec<char>&);
  // Jets
  LOG(INFO) << "Adding variable: sorted_jet_TLV_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "sorted_jet_TLV_NOSYS",
                                          static_cast<sorted_particle_1sel>(DefineHelpers::sortedPassedVector),
                                          {"jet_TLV_NOSYS", "jet_select_baselineJvt_NOSYS"});

  // Muons
  LOG(INFO) << "Adding variable: sorted_mu_TLV_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "sorted_mu_TLV_NOSYS",
                                          static_cast<sorted_particle_1sel>(DefineHelpers::sortedPassedVector),
                                          {"mu_TLV_NOSYS", "mu_select_tight_NOSYS"});

  // Electrons
  LOG(INFO) << "Adding variable: sorted_el_TLV_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "sorted_el_TLV_NOSYS",
                                          static_cast<sorted_particle_1sel>(DefineHelpers::sortedPassedVector),
                                          {"el_TLV_NOSYS", "el_select_tight_NOSYS"});

  // b-tagged jets.
  // using sorted_particle_2sel = ROOT::VecOps::RVec<TLV>(*)(const ROOT::VecOps::RVec<TLV>&,const ROOT::VecOps::RVec<char>&, const ROOT::VecOps::RVec<char>&);
  // LOG(INFO) << "Adding variable: sorted_bjet_TLV_NOSYS" << std::endl;
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                         "sorted_bjet_TLV_NOSYS",
  //                                         static_cast<sorted_particle_2sel>(DefineHelpers::sortedPassedVector),
  //                                         {"jet_TLV_NOSYS", "jet_select_baselineJvt_NOSYS","jet_GN2v01_FixedCutBEff_77_select"});

  // Lepton charges
  LOG(INFO) << "Adding variable: sorted_mu_charge_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "sorted_mu_charge_NOSYS",
                                         TtbarDileptonVars::sortedPassedChargeVector25,
                                         {"mu_charge", "mu_pt_NOSYS", "mu_select_tight_NOSYS"});

  LOG(INFO) << "Adding variable: sorted_el_charge_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "sorted_el_charge_NOSYS",
                                         TtbarDileptonVars::sortedPassedChargeVector25,
                                         {"el_charge", "el_pt_NOSYS", "el_select_tight_NOSYS"});

  LOG(INFO) << "Adding variable: region_name_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "region_name_NOSYS",
                                          ZllCalibrationVars::regionName,
                                          {"n_electrons_NOSYS", "n_muons_NOSYS", "sorted_el_TLV_NOSYS", "sorted_mu_TLV_NOSYS"});
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Analysis specific variables ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // LOG(INFO) << "Adding variable: n_bjets_NOSYS" << std::endl;
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                         "n_bjets_NOSYS",
  //                                         numberOfBJets30,
  //                                         {"sorted_bjet_TLV_NOSYS"});
  //LOG(INFO) << "Adding variable: PT diff"                                        

  // LOG(INFO) << "Adding variable: parton_truth_had_b_index_NOSYS" << std::endl;
  // mainNode = m_partonMatcher->systematicDefine(mainNode, "parton_truth_had_b_index_NOSYS",
  //                                         TtbarDileptonVars::getTruthIndexBHad,
  //                                         {"event_jet_truth_idx",
  //                                          "event_jet_truth_candidates"
  //                                          });
  // LOG(INFO) << "Adding variable: parton_truth_had_bbar_index_NOSYS" << std::endl;
  // mainNode = m_partonMatcher->systematicDefine(mainNode, "parton_truth_had_bbar_index_NOSYS",
  //                                         TtbarDileptonVars::getTruthIndexBbarHad,
  //                                         {"event_jet_truth_idx",
  //                                          "event_jet_truth_candidates"
  //                                          });                                         
  // // Fully matched mlb pairing masses added

 
  
  LOG(INFO) << "Adding variable: jet_GN2v01_pb_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "jet_GN2v01_pb_NOSYS",
                                         TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                                         {"jet_GN2v01_pb", "jets_ptOrderIndex_NOSYS"});

  LOG(INFO) << "Adding variable: jet_GN2v01_pc_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "jet_GN2v01_pc_NOSYS",
                                         TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                                         {"jet_GN2v01_pc", "jets_ptOrderIndex_NOSYS"});

  LOG(INFO) << "Adding variable: jet_GN2v01_pu_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "jet_GN2v01_pu_NOSYS",
                                         TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                                         {"jet_GN2v01_pu", "jets_ptOrderIndex_NOSYS"});
                                       
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Jet KLFitter indices ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  ////////////////////////////////////////////////////////////// Define the single lepton trigger ////////////////////////////////////////////////////
  // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //THESE ARE DIFFERENT DEPENDING ON WHICH SAMPLE
//   std::string campaign = id.campaign();
//   //std::cout << "Campaign is: " << campaign << std::endl;                                    
//   if (campaign == "mc23d") {
//     LOG(INFO) << "Sample is mc23d, using mc23d trigger definitions" << std::endl;
//     LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
//     mainNode = m_partonMatcher->systematicDefine(mainNode,
//                                           "Electron_pass_triggers_NOSYS",
//                                           TtbarSingleLeptonVars::ElectronSelection::triggers,
//                                           {"trigPassed_HLT_e60_lhmedium_L1EM22VHI", "trigPassed_HLT_e140_lhloose_L1EM22VHI","trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI",
//                                           "trigPassed_HLT_e26_lhtight_nod0_ivarloose","trigPassed_HLT_e140_lhloose_nod0","trigPassed_HLT_e60_lhmedium_nod0"});


//     LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
//     mainNode = m_partonMatcher->systematicDefine(mainNode,
//                                           "Muon_pass_triggers_NOSYS",
//                                           TtbarSingleLeptonVars::MuonSelection::triggers,
//                                           {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH", "trigPassed_HLT_mu50_L1MU14FCH",
//                                           "trigPassed_HLT_mu50","trigPassed_HLT_mu26_ivarmedium"});
//   } else if (campaign == "mc23e") {
//     LOG(INFO) << "Sample is mc23e, using mc23e trigger definitions" << std::endl;
//     LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
//     mainNode = m_partonMatcher->systematicDefine(mainNode,
//                                           "Electron_pass_triggers_NOSYS",
//                                           TtbarSingleLeptonVars::ElectronSelection::triggersnew,
//                                           {"trigPassed_HLT_e140_lhloose_L1eEM26M",
//     "trigPassed_HLT_e26_lhtight_ivarloose_L1eEM",
//     "trigPassed_HLT_e60_lhmedium_L1eEM"});


//     LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
//     mainNode = m_partonMatcher->systematicDefine(mainNode,
//                                           "Muon_pass_triggers_NOSYS",
//                                           TtbarSingleLeptonVars::MuonSelection::triggersnew,
//                                           {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH",
//   "trigPassed_HLT_mu50_L1MU14FCH"});
//   } else {
//     LOG(INFO) << "Sample is neither mc23d nor mc23e, using mc23e trigger definitions by default" << std::endl;
//     LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
//     mainNode = m_partonMatcher->systematicDefine(mainNode,
//                                           "Electron_pass_triggers_NOSYS",
//                                           TtbarSingleLeptonVars::ElectronSelection::triggersnew,
//                                           {"trigPassed_HLT_e140_lhloose_L1eEM26M",
//     "trigPassed_HLT_e26_lhtight_ivarloose_L1eEM",
//     "trigPassed_HLT_e60_lhmedium_L1eEM"});


//     LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
//     mainNode = m_partonMatcher->systematicDefine(mainNode,
//                                           "Muon_pass_triggers_NOSYS",
//                                           TtbarSingleLeptonVars::MuonSelection::triggersnew,
//                                           {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH",
//   "trigPassed_HLT_mu50_L1MU14FCH"});                                        
//                                         }
  
//   // create combined OR of electron OR muon trigger flags, signal region requires one lepton to have passed trigger
// LOG(INFO) << "Adding variable: SingleLepton_pass_triggers_NOSYS\n";

// auto orTriggers = [](bool e_pass, bool mu_pass) {
//     return (e_pass || mu_pass) ? 1 : 0;
// };


// mainNode = m_partonMatcher->systematicDefine(mainNode,
//                           "SingleLepton_pass_triggers_NOSYS",
//                           orTriggers,
//                           {"Electron_pass_triggers_NOSYS", "Muon_pass_triggers_NOSYS"});

                          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Define the Fully matched selection cuts ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fully matched definition                                       
  // LOG(INFO) << "Adding variable: is_fully_matched_NOSYS" << std::endl;     
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "is_fully_matched_NOSYS",
  //                                        TtbarDileptonVars::FullyMatched,
  //                                        {"event_jet_truth_idx",
  //                                         "event_jet_truth_candidates"
  //                                       });  
  
//   // First define the column with the correct value
//  mainNode = m_partonMatcher->systematicDefine(
//     mainNode,
//     "mcChannelNumber_NOSYS",
//     [](UInt_t mcID){ return (mcID == 601230) ? 1.0f : 0.0f; },
//     {"mcChannelNumber"} // column your lambda reads
// );

//   // Then register it in the systematic setup
//   mainNode = m_partonMatcher->systematicDefine(
//       mainNode,
//       "isTtbarDilepton_NOSYS",
//       [](const float& v){ return v; },
//       {"mcChannelNumber_NOSYS"} // use the column name, not a number
//   );
                  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////// btag //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // LOG(INFO) << "Adding variable: jet_GN2v01_Continuous_quantile_NOSYS\n";   
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "jet_GN2v01_Continuous_quantile_NOSYS",
  //                                        Variables::getPCBT,
  //                                        {"jet_GN2v01_FixedCutBEff_65_select", "jet_GN2v01_FixedCutBEff_70_select", "jet_GN2v01_FixedCutBEff_77_select", "jet_GN2v01_FixedCutBEff_85_select", "jet_GN2v01_FixedCutBEff_90_select"});
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////// Truth ////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
//   mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_t_afterFSR", "T_afterFSR");
//   mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_t_beforeFSR", "T_beforeFSR");
//   mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_tbar_afterFSR", "Tbar_afterFSR");
//   mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_tbar_beforeFSR", "Tbar_beforeFSR");
  
//   mainNode = m_partonMatcher->systematicDefine(
//     mainNode,
//     "truth_t_mass_afterFSR_NOSYS",
//     [](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& tlv) {
//         return tlv.M();
//     },
//     {"T_afterFSR_TLV_NOSYS"}
// );
//   mainNode = m_partonMatcher->systematicDefine(
//     mainNode,
//     "truth_tbar_mass_afterFSR_NOSYS",
//     [](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& tlv) {
//         return tlv.M();
//     },
//     {"Tbar_afterFSR_TLV_NOSYS"}
// );

//   LOG(INFO) << "Adding variable: truth_ttbar_mass_afterFSR_NOSYS\n";
//   mainNode = m_partonMatcher->systematicDefine(mainNode,
//                                          "truth_ttbar_mass_afterFSR_NOSYS",
//                                          Variables::invariantMassTwo,
//                                          {"T_afterFSR_TLV_NOSYS","Tbar_afterFSR_TLV_NOSYS"});  

//   LOG(INFO) << "Adding variable: truth_ttbar_mass_beforeFSR_NOSYS\n";
//   mainNode = m_partonMatcher->systematicDefine(mainNode,
//                                          "truth_ttbar_mass_beforeFSR_NOSYS",
//                                          Variables::invariantMassTwo,
//                                          {"T_beforeFSR_TLV_NOSYS","Tbar_beforeFSR_TLV_NOSYS"});    

  // mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_ttbar_afterFSR", "Ttbar_afterFSR");
  // mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_ttbar_beforeFSR", "Ttbar_beforeFSR");

  // LOG(INFO) << "Adding variable: truth_ttbar_mass_afterFSR_NOSYS\n";
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "truth_ttbar_mass_afterFSR_NOSYS",
  //                                        Variables::invariantMassTwo,
  //                                        {"T_afterFSR_TLV_NOSYS","Tbar_afterFSR_TLV_NOSYS"});  

  // LOG(INFO) << "Adding variable: truth_ttbar_mass_beforeFSR_NOSYS\n";
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "truth_ttbar_mass_beforeFSR_NOSYS",
  //                                        Variables::invariantMassTwo,
  //                                        {"T_beforeFSR_TLV_NOSYS","Tbar_beforeFSR_TLV_NOSYS"});                                          
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

ROOT::RDF::RNode ZllCalibrationFrame::defineVariablesNtuple(ROOT::RDF::RNode mainNode,
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
  // LOG(INFO) << "Adding variable: Truth_down_index\n";
  // mainNode = mainNode.Define("Truth_down_index",
  //                             Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
  //                             {"parton_truth_down_index_NOSYS","jets_ptOrderIndex_NOSYS"});
  // LOG(INFO) << "Adding variable: Truth_up_index\n";
  // mainNode = mainNode.Define("Truth_up_index",
  //                             Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
  //                             {"parton_truth_up_index_NOSYS","jets_ptOrderIndex_NOSYS"});
  LOG(INFO) << "Adding variable: Truth_had_b_index\n";
  mainNode = mainNode.Define("Truth_had_b_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_had_b_index_NOSYS","jets_ptOrderIndex_NOSYS"});
 
  LOG(INFO) << "Adding variable: Truth_had_bbar_index\n";
  mainNode = mainNode.Define("Truth_had_bbar_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_had_bbar_index_NOSYS","jets_ptOrderIndex_NOSYS"});
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

ROOT::RDF::RNode ZllCalibrationFrame::defineVariablesTruth(ROOT::RDF::RNode node,
                                                   const std::string& treeName,
                                                   const std::shared_ptr<Sample>& /*sample*/,
                                                   const UniqueSampleID& /*sampleID*/) {

  if (treeName != "parton") return node; // do not do anything for other trees ("truth" = parton truth)

  return node;
}

ROOT::RDF::RNode ZllCalibrationFrame::defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                           const std::string& /*treeName*/,
                                                           const std::shared_ptr<Sample>& /*sample*/,
                                                           const UniqueSampleID& /*sampleID*/) {
  return node;
}