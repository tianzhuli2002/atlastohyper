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

  std::string jet_OR_flag = "jet_select_passesOR_NOSYS"; // "jet_select_or_NOSYS"

  LOG(INFO) << "Adding variable: nJets_good_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "nJets_good_NOSYS",
                                         Variables::nGoodJets,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS",jet_OR_flag});                                                  

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Pt indices ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  LOG(INFO) << "Adding variable: jets_ptOrderIndex_NOSYS\n";
  mainNode = MainFrame::systematicDefine(mainNode,
                                         "jets_ptOrderIndex_NOSYS",
                                         Variables::indexPtOrderedGoodTwoSelec,
                                         {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS",jet_OR_flag});

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Parton TLVs ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_W_beforeFSR_from_t", "W_from_t_beforeFSR");
  mainNode = Variables::createPartonTLV(this, mainNode, "Ttbar_History_MC_W_afterFSR_from_t", "W_from_t_afterFSR");
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
  
