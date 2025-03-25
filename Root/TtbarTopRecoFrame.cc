#include "PartonMatcher/TtbarTopRecoFrame.h"
#include "FastFrames/DefineHelpers.h"
#include "PartonMatcher/Variables.h"

ROOT::RDF::RNode TtbarTopRecoFrame::defineVariables(ROOT::RDF::RNode mainNode,
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
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////// Trigger configuration ////////////////////////////////////////////////////
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

    // Define jet, electron and muon triggers
    std::vector<std::string> jetTrigggerBranches;
    std::vector<std::string> electronTriggerBranches;
    std::vector<std::string> muonTriggerBranches;
    // jet
    for (const auto& trigger : triggerBranches) {
        if (trigger.find("HLT_e") == std::string::npos && trigger.find("HLT_mu") == std::string::npos) {
            jetTrigggerBranches.push_back(trigger);
        }
    }
    // electron
    for (const auto& trigger : triggerBranches) {
        if (trigger.find("HLT_e") != std::string::npos) {
            electronTriggerBranches.push_back(trigger);
        }
    }
    // muon
    for (const auto& trigger : triggerBranches) {
        if (trigger.find("HLT_mu") != std::string::npos) {
            muonTriggerBranches.push_back(trigger);
        }
    }

    if (jetTrigggerBranches.size() == 0 && electronTriggerBranches.size() == 0 && muonTriggerBranches.size() == 0) {
        LOG(ERROR) << "No trigger branches found";
        throw std::invalid_argument("");
    } else if (jetTrigggerBranches.size() != 0 && (electronTriggerBranches.size() != 0 || muonTriggerBranches.size() != 0)) {
        LOG(ERROR) << "Ntuple with bad configuration! Jet triggers found, but also electron or muon triggers found";
        throw std::invalid_argument("");
    } else if (jetTrigggerBranches.size() == 0 && electronTriggerBranches.size() != 0 && muonTriggerBranches.size() == 0) {
        LOG(ERROR) << "Ntuple with bad configuration! Electron triggers found, but non muon triggers found";
        throw std::invalid_argument("");
    } else if (jetTrigggerBranches.size() == 0 && electronTriggerBranches.size() == 0 && muonTriggerBranches.size() != 0) {
        LOG(ERROR) << "Ntuple with bad configuration! Muon triggers found, but non electron triggers found";
        throw std::invalid_argument("");
    } else if (jetTrigggerBranches.size() != 0 && electronTriggerBranches.size() == 0 && muonTriggerBranches.size() == 0) {
        LOG(INFO) << "Good trigger config! Only jet triggers found.";
    } else if (jetTrigggerBranches.size() == 0 && electronTriggerBranches.size() != 0 && muonTriggerBranches.size() != 0) {
        LOG(INFO) << "Good trigger config! Only light-lepton triggers found.";
    } else {
        LOG(ERROR) << "Ntuple with bad configuration!";
        throw std::invalid_argument("");
    }

    auto firstTrigger = [](bool triggerDecision){return triggerDecision;};
    auto triggerOR = [](bool triggerDecision1, bool triggerDecision2){return triggerDecision1 || triggerDecision2;};
    if (jetTrigggerBranches.size() != 0){
        
        LOG(INFO) << "Adding variable: Passed_Any_Jet_Trigger_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
                                                "Passed_Any_Jet_Trigger_NOSYS",
                                                firstTrigger,
                                                {jetTrigggerBranches.at(0)});
        
        LOG(INFO) << "Adding variable: Passed_Any_Jet_Trigger_NOSYS\n";
        for (std::size_t iTriggerChain{1}; iTriggerChain < jetTrigggerBranches.size(); iTriggerChain++) {
    
        mainNode = mainNode.Redefine("Passed_Any_Jet_Trigger_NOSYS",
                                        triggerOR,
                                        {"Passed_Any_Jet_Trigger_NOSYS", jetTrigggerBranches.at(iTriggerChain)});
    
        }
    } else {
        LOG(INFO) << "Adding variable: Passed_Any_Muon_Trigger_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
                                                "Passed_Any_Muon_Trigger_NOSYS",
                                                firstTrigger,
                                                {muonTriggerBranches.at(0)});

        LOG(INFO) << "Adding variable: Passed_Any_Electron_Trigger_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
                                                "Passed_Any_Electron_Trigger_NOSYS",
                                                firstTrigger,
                                                {electronTriggerBranches.at(0)});

        for (std::size_t iTriggerChain{1}; iTriggerChain < muonTriggerBranches.size(); iTriggerChain++) {
    
        mainNode = mainNode.Redefine("Passed_Any_Muon_Trigger_NOSYS",
                                        triggerOR,
                                        {"Passed_Any_Muon_Trigger_NOSYS", muonTriggerBranches.at(iTriggerChain)});
    
        }

        for (std::size_t iTriggerChain{1}; iTriggerChain < electronTriggerBranches.size(); iTriggerChain++) {
        
        mainNode = mainNode.Redefine("Passed_Any_Electron_Trigger_NOSYS",
                                        triggerOR,
                                        {"Passed_Any_Electron_Trigger_NOSYS", electronTriggerBranches.at(iTriggerChain)});
    
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////// Object multiplicity ////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    LOG(INFO) << "Adding variable: nJets_good_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "nJets_good_NOSYS",
                                            Variables::nGoodJets,
                                            {"jet_pt_NOSYS", "jet_select_baselineJvt_NOSYS"});   

    LOG(INFO) << "Adding variable: jet_Ordered_TLV_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "jet_Ordered_TLV_NOSYS",
                                            static_cast<std::vector<TLV>(*)(const std::vector<TLV>&,
                                        const std::vector<char>&)>(DefineHelpers::sortedPassedVector),
                                            {"jet_TLV_NOSYS", "jet_select_baselineJvt_NOSYS"}); 


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////// Parton TLVs ////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // W bosons
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_W_beforeFSR_from_t", "W_from_t_beforeFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_W_afterFSR_from_t", "W_from_t_afterFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_W_beforeFSR_from_tbar", "W_from_tbar_beforeFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_W_afterFSR_from_tbar", "W_from_tbar_afterFSR");

    // Top quarks
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_t_afterFSR", "T_afterFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_t_beforeFSR", "T_beforeFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_tbar_afterFSR", "Tbar_afterFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_tbar_beforeFSR", "Tbar_beforeFSR");

    // B quarks
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_b_afterFSR_from_t", "b_afterFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_b_beforeFSR_from_t", "b_beforeFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_bbar_afterFSR_from_tbar", "bbar_afterFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_bbar_beforeFSR_from_tbar", "bbar_beforeFSR");

    // Light quarks
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_Wdecay1_afterFSR_from_t", "q1_from_t_afterFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_Wdecay1_beforeFSR_from_t", "q1_from_t_beforeFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_Wdecay2_afterFSR_from_t", "q2_from_t_afterFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_Wdecay2_beforeFSR_from_t", "q2_from_t_beforeFSR");

    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_Wdecay1_afterFSR_from_tbar", "q1_from_tbar_afterFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_Wdecay1_beforeFSR_from_tbar", "q1_from_tbar_beforeFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_Wdecay2_afterFSR_from_tbar", "q2_from_tbar_afterFSR");
    mainNode = Variables::createPartonTLV(m_partonMatcher, mainNode, "Ttbar_History_MC_Wdecay2_beforeFSR_from_tbar", "q2_from_tbar_beforeFSR");

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////// Parton Kinematics ////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    LOG(INFO) << "Adding variable: truth_t_mass_afterFSR_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "truth_t_mass_afterFSR_NOSYS",
                                            Variables::invariantMassThree,
                                            {"q1_from_t_afterFSR_TLV_NOSYS", "q2_from_t_afterFSR_TLV_NOSYS","b_afterFSR_TLV_NOSYS"});  

    LOG(INFO) << "Adding variable: truth_t_mass_beforeFSR_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "truth_t_mass_beforeFSR_NOSYS",
                                            Variables::invariantMassThree,
                                            {"q1_from_t_beforeFSR_TLV_NOSYS", "q2_from_t_beforeFSR_TLV_NOSYS","b_beforeFSR_TLV_NOSYS"});    

    LOG(INFO) << "Adding variable: truth_tbar_mass_afterFSR_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "truth_tbar_mass_afterFSR_NOSYS",
                                            Variables::invariantMassThree,
                                            {"q1_from_tbar_afterFSR_TLV_NOSYS", "q2_from_tbar_afterFSR_TLV_NOSYS","bbar_afterFSR_TLV_NOSYS"});

    LOG(INFO) << "Adding variable: truth_tbar_mass_beforeFSR_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "truth_tbar_mass_beforeFSR_NOSYS",
                                            Variables::invariantMassThree,
                                            {"q1_from_tbar_beforeFSR_TLV_NOSYS", "q2_from_tbar_beforeFSR_TLV_NOSYS","bbar_beforeFSR_TLV_NOSYS"});

    LOG(INFO) << "Adding variable: truth_W_mass_from_t_afterFSR_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "truth_W_mass_from_t_afterFSR_NOSYS",
                                            Variables::invariantMassTwo,
                                            {"q1_from_t_afterFSR_TLV_NOSYS", "q2_from_t_afterFSR_TLV_NOSYS"});

    LOG(INFO) << "Adding variable: truth_W_mass_from_t_beforeFSR_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "truth_W_mass_from_t_beforeFSR_NOSYS",
                                            Variables::invariantMassTwo,
                                            {"q1_from_t_beforeFSR_TLV_NOSYS", "q2_from_t_beforeFSR_TLV_NOSYS"});

    LOG(INFO) << "Adding variable: truth_W_mass_from_tbar_afterFSR_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "truth_W_mass_from_tbar_afterFSR_NOSYS",
                                            Variables::invariantMassTwo,
                                            {"q1_from_tbar_afterFSR_TLV_NOSYS", "q2_from_tbar_afterFSR_TLV_NOSYS"});

    LOG(INFO) << "Adding variable: truth_W_mass_from_tbar_beforeFSR_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                            "truth_W_mass_from_tbar_beforeFSR_NOSYS",
                                            Variables::invariantMassTwo,
                                            {"q1_from_tbar_beforeFSR_TLV_NOSYS", "q2_from_tbar_beforeFSR_TLV_NOSYS"});

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////// HyPER performance ////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // For AllHad 
    if (jetTrigggerBranches.size() != 0){
        LOG(INFO) << "Adding variable: Number_W_correctly_reconstructed_HyPER_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Number_W_correctly_reconstructed_HyPER_NOSYS",
            HyPERPerformance::countCorrectlyReconstructedW,
            {"hyper_TtbarAllHadronic_W1_Indices_NOSYS", "hyper_TtbarAllHadronic_W2_Indices_NOSYS", "event_jet_truth_idx", "event_jet_truth_candidates"});

        LOG(INFO) << "Adding variable: Number_W_correctly_reconstructed_KLF_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Number_W_correctly_reconstructed_KLF_NOSYS",
            HyPERPerformance::countCorrectlyReconstructedW_KLF,
            {"klfitter_lj1_from_top1_jetIndex_NOSYS", "klfitter_lj2_from_top1_jetIndex_NOSYS",
            "klfitter_lj1_from_top2_jetIndex_NOSYS", "klfitter_lj2_from_top2_jetIndex_NOSYS",
            "event_jet_truth_idx", "event_jet_truth_candidates"});

        LOG(INFO) << "Adding variable: Number_Top_correctly_reconstructed_HyPER_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Number_Top_correctly_reconstructed_HyPER_NOSYS",
            HyPERPerformance::countCorrectlyReconstructedTop,
            {"hyper_TtbarAllHadronic_Top1_Indices_NOSYS","hyper_TtbarAllHadronic_Top2_Indices_NOSYS","hyper_TtbarAllHadronic_W1_Indices_NOSYS", "hyper_TtbarAllHadronic_W2_Indices_NOSYS", "event_jet_truth_idx", "event_jet_truth_candidates"});
    
        LOG(INFO) << "Adding variable: Number_Top_correctly_reconstructed_KLF_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Number_Top_correctly_reconstructed_KLF_NOSYS",
            HyPERPerformance::countCorrectlyReconstructedTop_KLF,
            {"klfitter_lj1_from_top1_jetIndex_NOSYS","klfitter_lj2_from_top1_jetIndex_NOSYS", "klfitter_b_from_top1_jetIndex_NOSYS",
            "klfitter_lj1_from_top2_jetIndex_NOSYS","klfitter_lj2_from_top2_jetIndex_NOSYS", "klfitter_b_from_top2_jetIndex_NOSYS",
            "event_jet_truth_idx", "event_jet_truth_candidates"});

        // Check if the algorithms converged for two tops.
        LOG(INFO) << "Adding variable: Two_Tops_Reco_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Two_Tops_Reco_NOSYS",
            HyPERPerformance::HyPERConvergedTwoTops,
            {"hyper_TtbarAllHadronic_Top1_Indices_NOSYS","hyper_TtbarAllHadronic_Top2_Indices_NOSYS","hyper_TtbarAllHadronic_W1_Indices_NOSYS", "hyper_TtbarAllHadronic_W2_Indices_NOSYS"});

        LOG(INFO) << "Adding variable: Two_Tops_Reco_KLF_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Two_Tops_Reco_KLF_NOSYS",
            HyPERPerformance::KLFConvergedTwoTops,
            {"klfitter_lj1_from_top1_jetIndex_NOSYS","klfitter_lj2_from_top1_jetIndex_NOSYS", "klfitter_b_from_top1_jetIndex_NOSYS",
            "klfitter_lj1_from_top2_jetIndex_NOSYS","klfitter_lj2_from_top2_jetIndex_NOSYS", "klfitter_b_from_top2_jetIndex_NOSYS"});

        // Mass variables
        LOG(INFO) << "Adding variable: Top1_mass_KLF_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Top1_mass_KLF_NOSYS",
            Variables::recoInvariantMassThree_KLF<unsigned int>,
            {"klfitter_lj1_from_top1_jetIndex_NOSYS","klfitter_lj2_from_top1_jetIndex_NOSYS", "klfitter_b_from_top1_jetIndex_NOSYS","jet_TLV_NOSYS"});
        
        LOG(INFO) << "Adding variable: Top2_mass_KLF_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Top2_mass_KLF_NOSYS",
            Variables::recoInvariantMassThree_KLF<unsigned int>,
            {"klfitter_lj1_from_top2_jetIndex_NOSYS","klfitter_lj2_from_top2_jetIndex_NOSYS", "klfitter_b_from_top2_jetIndex_NOSYS","jet_TLV_NOSYS"});
        
        LOG(INFO) << "Adding variable: W1_mass_KLF_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "W1_mass_KLF_NOSYS",
            Variables::recoInvariantMassTwo_KLF<unsigned int>,
            {"klfitter_lj1_from_top1_jetIndex_NOSYS","klfitter_lj2_from_top1_jetIndex_NOSYS","jet_TLV_NOSYS"});
        
        LOG(INFO) << "Adding variable: W2_mass_KLF_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "W2_mass_KLF_NOSYS",
            Variables::recoInvariantMassTwo_KLF<unsigned int>,
            {"klfitter_lj1_from_top2_jetIndex_NOSYS","klfitter_lj2_from_top2_jetIndex_NOSYS","jet_TLV_NOSYS"});

        LOG(INFO) << "Adding variable: Top1_mass_HyPER_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Top1_mass_HyPER_NOSYS",
            Variables::recoInvariantMassThree_HyPER<long>,
            {"hyper_TtbarAllHadronic_Top1_Indices_NOSYS","jet_TLV_NOSYS"});
        
        LOG(INFO) << "Adding variable: Top2_mass_HyPER_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "Top2_mass_HyPER_NOSYS",
            Variables::recoInvariantMassThree_HyPER<long>,
            {"hyper_TtbarAllHadronic_Top2_Indices_NOSYS","jet_TLV_NOSYS"});

        LOG(INFO) << "Adding variable: W1_mass_HyPER_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "W1_mass_HyPER_NOSYS",
            Variables::recoInvariantMassTwo_HyPER<long>,
            {"hyper_TtbarAllHadronic_W1_Indices_NOSYS","jet_TLV_NOSYS"});
        
        LOG(INFO) << "Adding variable: W2_mass_HyPER_NOSYS\n";
        mainNode = m_partonMatcher->systematicDefine(mainNode,
            "W2_mass_HyPER_NOSYS",
            Variables::recoInvariantMassTwo_HyPER<long>,
            {"hyper_TtbarAllHadronic_W2_Indices_NOSYS","jet_TLV_NOSYS"});

    }

    return mainNode;
}

ROOT::RDF::RNode TtbarTopRecoFrame::defineVariablesNtuple(ROOT::RDF::RNode mainNode,
          const std::shared_ptr<Sample>& sample,
          const UniqueSampleID& id) {

mainNode = this->defineVariables(mainNode, sample, id);
return mainNode;
}

ROOT::RDF::RNode TtbarTopRecoFrame::defineVariablesTruth(ROOT::RDF::RNode node,
         const std::string& /*sample*/,
         const std::shared_ptr<Sample>& /*sample*/,
         const UniqueSampleID& /*sampleID*/) {
return node;
}

ROOT::RDF::RNode TtbarTopRecoFrame::defineVariablesNtupleTruth(ROOT::RDF::RNode node,
               const std::string& /*treeName*/,
               const std::shared_ptr<Sample>& /*sample*/,
               const UniqueSampleID& /*sampleID*/) {
return node;
}