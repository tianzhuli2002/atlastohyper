#include "PartonMatcher/PartonMatcher.h"
#include "FastFrames/DefineHelpers.h"
#include "FastFrames/UniqueSampleID.h"
#include "PartonMatcher/Variables.h"

ROOT::RDF::RNode PartonMatcher::defineVariables(ROOT::RDF::RNode mainNode,
                                                const std::shared_ptr<Sample>& /*sample*/,
                                                const UniqueSampleID& /*id*/) {

  // You can also use the UniqueSampleID object to apply a custom defione
  // based on the sample and the subsample
  //   sample->name(): is the name of the sample defined in the config
  //   id.dsid() returns sample DSID
  //   id.campaign() returns sample campaign
  //   id.simulation() return simulation flavour
  // You can use it in your functions to apply only per sample define
  // Add all the variables needed for the analysis
  LOG(INFO) << "You are running custom class: MyFrame"<<std::endl;    

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Object multiplicity ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: nJets_good_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nJets_good_NOSYS",
                                         Variables::nGoodJets,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS"});   

  LOG(INFO) << "Adding variable: jet_Ordered_TLV_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "jet_Ordered_TLV_NOSYS",
                                         static_cast<std::vector<TLV>(*)(const std::vector<TLV>&,
                                      const std::vector<char>&)>(DefineHelpers::sortedPassedVector),
                                         {"jet_TLV_NOSYS", "jet_select_baselineJvt_NOSYS"});                                           

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Pt indices ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: jets_ptOrderIndex_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "jets_ptOrderIndex_NOSYS",
                                         Variables::indexPtOrderedGoodOneSelec,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Bjets multiplicity and properties ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: nBjets_GN2v01_100_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nBjets_GN2v01_100_NOSYS",
                                         Variables::nBJets(100),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile"});

  LOG(INFO) << "Adding variable: nBjets_GN2v01_90_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nBjets_GN2v01_90_NOSYS",
                                         Variables::nBJets(90),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile"});

  LOG(INFO) << "Adding variable: nBjets_GN2v01_85_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nBjets_GN2v01_85_NOSYS",
                                         Variables::nBJets(85),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile"});

  LOG(INFO) << "Adding variable: nBjets_GN2v01_77_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nBjets_GN2v01_77_NOSYS",
                                         Variables::nBJets(77),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile"});

  LOG(INFO) << "Adding variable: nBjets_GN2v01_70_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nBjets_GN2v01_70_NOSYS",
                                         Variables::nBJets(70),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile"});
  
  LOG(INFO) << "Adding variable: nBjets_GN2v01_65_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nBjets_GN2v01_65_NOSYS",
                                         Variables::nBJets(65),
                                         {"jets_ptOrderIndex_NOSYS", "jet_GN2v01_Continuous_quantile"});

  LOG(INFO) << "Adding variable: GN2v01_PCBT_quantile_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "GN2v01_PCBT_quantile_NOSYS",
                                         Variables::pTSortedContainer<int>,
                                         {"jet_GN2v01_Continuous_quantile", "jets_ptOrderIndex_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Parton TLVs ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // W bosons
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_W_beforeFSR_from_t", "W_from_t_beforeFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_W_afterFSR_from_t", "W_from_t_afterFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_W_beforeFSR_from_tbar", "W_from_tbar_beforeFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_W_afterFSR_from_tbar", "W_from_tbar_afterFSR");

  // Top quarks
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_t_afterFSR", "T_afterFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_t_beforeFSR", "T_beforeFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_tbar_afterFSR", "Tbar_afterFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_tbar_beforeFSR", "Tbar_beforeFSR");

  // B quarks
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_b_afterFSR_from_t", "b_afterFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_b_beforeFSR_from_t", "b_beforeFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_bbar_afterFSR_from_tbar", "bbar_afterFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_bbar_beforeFSR_from_tbar", "bbar_beforeFSR");

  // Light quarks
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_Wdecay1_afterFSR_from_t", "q1_from_t_afterFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_Wdecay1_beforeFSR_from_t", "q1_from_t_beforeFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_Wdecay2_afterFSR_from_t", "q2_from_t_afterFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_Wdecay2_beforeFSR_from_t", "q2_from_t_beforeFSR");

  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_Wdecay1_afterFSR_from_tbar", "q1_from_tbar_afterFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_Wdecay1_beforeFSR_from_tbar", "q1_from_tbar_beforeFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_Wdecay2_afterFSR_from_tbar", "q2_from_tbar_afterFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_Wdecay2_beforeFSR_from_tbar", "q2_from_tbar_beforeFSR");

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Parton Kinematics ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: truth_t_mass_afterFSR_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "truth_t_mass_afterFSR_NOSYS",
                                         Variables::invariantMassThree,
                                         {"q1_from_t_afterFSR_TLV_NOSYS", "q2_from_t_afterFSR_TLV_NOSYS","b_afterFSR_TLV_NOSYS"});  

  LOG(INFO) << "Adding variable: truth_t_mass_beforeFSR_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "truth_t_mass_beforeFSR_NOSYS",
                                         Variables::invariantMassThree,
                                         {"q1_from_t_beforeFSR_TLV_NOSYS", "q2_from_t_beforeFSR_TLV_NOSYS","b_beforeFSR_TLV_NOSYS"});    

  LOG(INFO) << "Adding variable: truth_tbar_mass_afterFSR_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "truth_tbar_mass_afterFSR_NOSYS",
                                         Variables::invariantMassThree,
                                         {"q1_from_tbar_afterFSR_TLV_NOSYS", "q2_from_tbar_afterFSR_TLV_NOSYS","bbar_afterFSR_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: truth_tbar_mass_beforeFSR_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "truth_tbar_mass_beforeFSR_NOSYS",
                                         Variables::invariantMassThree,
                                         {"q1_from_tbar_beforeFSR_TLV_NOSYS", "q2_from_tbar_beforeFSR_TLV_NOSYS","bbar_beforeFSR_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: truth_W_mass_from_t_afterFSR_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "truth_W_mass_from_t_afterFSR_NOSYS",
                                         Variables::invariantMassTwo,
                                         {"q1_from_t_afterFSR_TLV_NOSYS", "q2_from_t_afterFSR_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: truth_W_mass_from_t_beforeFSR_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "truth_W_mass_from_t_beforeFSR_NOSYS",
                                         Variables::invariantMassTwo,
                                         {"q1_from_t_beforeFSR_TLV_NOSYS", "q2_from_t_beforeFSR_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: truth_W_mass_from_tbar_afterFSR_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "truth_W_mass_from_tbar_afterFSR_NOSYS",
                                         Variables::invariantMassTwo,
                                         {"q1_from_tbar_afterFSR_TLV_NOSYS", "q2_from_tbar_afterFSR_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: truth_W_mass_from_tbar_beforeFSR_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "truth_W_mass_from_tbar_beforeFSR_NOSYS",
                                         Variables::invariantMassTwo,
                                         {"q1_from_tbar_beforeFSR_TLV_NOSYS", "q2_from_tbar_beforeFSR_TLV_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Parton-reco matching ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: q1_index_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "q1_index_NOSYS",
                                         Variables::matchPartonRecoParticles,
                                         {"q1_from_t_beforeFSR_TLV_NOSYS", "jet_Ordered_TLV_NOSYS"});
  
  LOG(INFO) << "Adding variable: q2_index_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "q2_index_NOSYS",
                                         Variables::matchPartonRecoParticles,
                                         {"q2_from_t_beforeFSR_TLV_NOSYS", "jet_Ordered_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: b_index_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "b_index_NOSYS",
                                         Variables::matchPartonRecoParticles,
                                         {"b_beforeFSR_TLV_NOSYS", "jet_Ordered_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: q1bar_index_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "q1bar_index_NOSYS",
                                         Variables::matchPartonRecoParticles,
                                         {"q1_from_tbar_beforeFSR_TLV_NOSYS", "jet_Ordered_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: q2bar_index_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "q2bar_index_NOSYS",
                                         Variables::matchPartonRecoParticles,
                                         {"q2_from_tbar_beforeFSR_TLV_NOSYS", "jet_Ordered_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: bbar_index_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "bbar_index_NOSYS",
                                         Variables::matchPartonRecoParticles,
                                         {"bbar_beforeFSR_TLV_NOSYS", "jet_Ordered_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: W1_reco_mass_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "W1_reco_mass_NOSYS",
                                         Variables::recoInvariantMassTwo,
                                         {"q1_index_NOSYS", "q2_index_NOSYS", "jet_Ordered_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: W2_reco_mass_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "W2_reco_mass_NOSYS",
                                         Variables::recoInvariantMassTwo,
                                         {"q1bar_index_NOSYS", "q2bar_index_NOSYS", "jet_Ordered_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: Top1_reco_mass_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "Top1_reco_mass_NOSYS",
                                         Variables::recoInvariantMassThree,
                                         {"q1_index_NOSYS", "q2_index_NOSYS", "b_index_NOSYS", "jet_Ordered_TLV_NOSYS"});

  LOG(INFO) << "Adding variable: Top2_reco_mass_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "Top2_reco_mass_NOSYS",
                                         Variables::recoInvariantMassThree,
                                         {"q1bar_index_NOSYS", "q2bar_index_NOSYS", "bbar_index_NOSYS", "jet_Ordered_TLV_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Event classfication ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: W1_Matched_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "W1_Matched_NOSYS",
                                         Variables::isWMatchedCorrectly,
                                         {"q1_index_NOSYS", "q2_index_NOSYS"});
  
  LOG(INFO) << "Adding variable: W2_Matched_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "W2_Matched_NOSYS",
                                         Variables::isWMatchedCorrectly,
                                         {"q1bar_index_NOSYS", "q2bar_index_NOSYS"});

  LOG(INFO) << "Adding variable: W_Not_Overlap_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "W_Not_Overlap_NOSYS",
                                         Variables::isNotWOverlap,
                                         {"q1_index_NOSYS", "q2_index_NOSYS", "q1bar_index_NOSYS", "q2bar_index_NOSYS"});

  LOG(INFO) << "Adding variable: Top1_Matched_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "Top1_Matched_NOSYS",
                                         Variables::isTopMatchedCorrectly,
                                         {"q1_index_NOSYS", "q2_index_NOSYS", "b_index_NOSYS"});

  LOG(INFO) << "Adding variable: Top2_Matched_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "Top2_Matched_NOSYS",
                                         Variables::isTopMatchedCorrectly,
                                         {"q1bar_index_NOSYS", "q2bar_index_NOSYS", "bbar_index_NOSYS"});

  LOG(INFO) << "Adding variable: b_Not_Overlap_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "b_Not_Overlap_NOSYS",
                                         Variables::isNotBOverlap,
                                         {"b_index_NOSYS", "bbar_index_NOSYS"});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Trigger big OR ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Get branch names
  auto columns = mainNode.GetColumnNames();
  // Store the trigger branches
  std::vector<std::string> triggerBranches;
  for (const auto& column : columns) {
    if (column.find("trigPassed") != std::string::npos) {
      triggerBranches.push_back(column);
    }
  }

  auto firstTrigger = [](bool triggerDecision){return triggerDecision;};
  LOG(INFO) << "Adding variable: Passed_Any_Trigger_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "Passed_Any_Trigger_NOSYS",
                                         firstTrigger,
                                         {triggerBranches.at(0)});


  auto triggerOR = [](bool triggerDecision1, bool triggerDecision2){return triggerDecision1 || triggerDecision2;};
  LOG(INFO) << "Adding variable: Passed_Any_Trigger_NOSYS\n";
  for (std::size_t iTriggerChain{1}; iTriggerChain < triggerBranches.size(); iTriggerChain++) {
    
    mainNode = mainNode.Redefine("Passed_Any_Trigger_NOSYS",
                                 triggerOR,
                                 {"Passed_Any_Trigger_NOSYS", triggerBranches.at(iTriggerChain)});

  }

  // Uncomment to print variables!
  //auto display = mainNode.Display({"b_index_NOSYS", "q1_index_NOSYS", "q2_index_NOSYS", "Top1_Matched_NOSYS", "W1_Matched_NOSYS"}, 200);
  //display->Print();

  return mainNode;
}

ROOT::RDF::RNode PartonMatcher::defineVariablesNtuple(ROOT::RDF::RNode mainNode,
                                                      const std::shared_ptr<Sample>& sample,
                                                      const UniqueSampleID& id) {
  
  mainNode = this->defineVariables(mainNode, sample, id);
  return mainNode;
}

ROOT::RDF::RNode PartonMatcher::defineVariablesTruth(ROOT::RDF::RNode node,
                                                     const std::string& /*sample*/,
                                                     const std::shared_ptr<Sample>& /*sample*/,
                                                     const UniqueSampleID& /*sampleID*/) {
  return node;
}

ROOT::RDF::RNode PartonMatcher::defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                           const std::string& /*treeName*/,
                                                           const std::shared_ptr<Sample>& /*sample*/,
                                                           const UniqueSampleID& /*sampleID*/) {
  return node;
}
  
