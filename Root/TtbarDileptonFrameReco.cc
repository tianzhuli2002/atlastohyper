#include "PartonMatcher/TtbarDileptonFrameReco.h"
#include "FastFrames/UniqueSampleID.h"
#include "PartonMatcher/Variables.h"
#include "FastFrames/DefineHelpers.h"
#include <algorithm>
#include <vector>

ROOT::RDF::RNode TtbarDileptonFrameReco::defineVariables(ROOT::RDF::RNode mainNode,
                                              const std::shared_ptr<Sample>& sample,
                                              const UniqueSampleID& id) {

  // Add all the variables needed for the analysisx
  LOG(INFO) << "You are running custom class: TtbarDileptonFrameReco" << std::endl;

  // Provide local fallbacks for optional HyPER / nuflows branches so that
  // downstream Snapshot/ntuple-selection does not abort if those columns
  // are missing from the input files. Keep these fallbacks local to this
  // custom frame (avoid changing MainFrame.cc as requested).
  {
    auto cols = mainNode.GetColumnNames();
    auto has_col = [&](const std::string &n){ for (const auto &c : cols) if (c==n) return true; return false; };

    // Only provide fallbacks for optional INPUT branches (not for
    // ntuple output columns that will be defined later). This avoids
    // creating placeholder columns that later cause duplicate-Define
    // errors when the frame attempts to Define the real outputs.
    const std::vector<std::pair<std::string,std::string>> fallbacks = {
      {"hyper_TtbarDiLepton_Top1_Indices_NOSYS","vec_int"},
      {"hyper_TtbarDiLepton_Top1_IDs_NOSYS","vec_int"},
      {"hyper_TtbarDiLepton_Top2_Indices_NOSYS","vec_int"},
      {"hyper_TtbarDiLepton_Top2_IDs_NOSYS","vec_int"},
      {"nuflows_nu_out_NOSYS","vec_float"},
      {"nuflows_loglik_NOSYS","float"}
    };

    for (const auto &fb : fallbacks) {
      const auto &name = fb.first;
      const auto &kind = fb.second;
      if (has_col(name)) continue;
      LOG(INFO) << "Defining fallback " << name << " (empty) in reco frame" << std::endl;
      if (kind=="vec_int") {
        mainNode = mainNode.Define(name, [](){ ROOT::VecOps::RVec<int> out; return out; });
      } else if (kind=="vec_float") {
        mainNode = mainNode.Define(name, [](){ ROOT::VecOps::RVec<float> out; return out; });
      } else if (kind=="int") {
        mainNode = mainNode.Define(name, [](){ return 0; });
      } else { // float
        mainNode = mainNode.Define(name, [](){ return 0.0f; });
      }
      // refresh cols so subsequent checks see this column
      cols = mainNode.GetColumnNames();
    }
  }

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
  using sorted_particle_2sel = ROOT::VecOps::RVec<TLV>(*)(const ROOT::VecOps::RVec<TLV>&,const ROOT::VecOps::RVec<char>&, const ROOT::VecOps::RVec<char>&);
  LOG(INFO) << "Adding variable: sorted_bjet_TLV_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "sorted_bjet_TLV_NOSYS",
                                          static_cast<sorted_particle_2sel>(DefineHelpers::sortedPassedVector),
                                          {"jet_TLV_NOSYS", "jet_select_baselineJvt_NOSYS","jet_GN2v01_FixedCutBEff_77_select"});

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
                                          TtbarDileptonVars::regionName,
                                          {"n_electrons_NOSYS", "n_muons_NOSYS", "sorted_el_TLV_NOSYS", "sorted_mu_TLV_NOSYS"});
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Analysis specific variables ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  LOG(INFO) << "Adding variable: n_bjets_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "n_bjets_NOSYS",
                                          numberOfBJets30,
                                          {"sorted_bjet_TLV_NOSYS"});
  //LOG(INFO) << "Adding variable: PT diff"                                        

  // ttbar lepton pair invariant mass
  LOG(INFO) << "Adding variable: ttbar_mass_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "ttbar_mass_NOSYS",
                                          TtbarDileptonVars::ttbarLeptonPairInvariantMass,
                                          {"region_name_NOSYS", "sorted_mu_TLV_NOSYS", "sorted_el_TLV_NOSYS"});
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// MATCHING STUFF ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // minimum lepton-jet pair invariant mass                                          
  LOG(INFO) << "Adding variable: min_mlb_pair_masses_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "min_mlb_pair_masses_NOSYS",
                                         TtbarDileptonVars::minimumLeptonJetPairInvariantMass,
                                         {"sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS","sorted_bjet_TLV_NOSYS","jet_TLV_NOSYS","jet_GN2v01_FixedCutBEff_77_select","n_electrons_NOSYS","n_muons_NOSYS"});
  LOG(INFO) << "Adding variable: pass_min_mlb_pair_masses_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "pass_min_mlb_pair_masses_NOSYS",
                                         TtbarDileptonVars::pass_min_lb_invar_mass,
                                         {"sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS","sorted_bjet_TLV_NOSYS","jet_TLV_NOSYS","jet_GN2v01_FixedCutBEff_77_select","n_electrons_NOSYS","n_muons_NOSYS"});
  
  LOG(INFO) << "Fully matched mlb pairing masses added" << std::endl;
  // Add pairing-by-min-sum-of-squared-mlj outputs (includes correctness flag using truth arrays)
  LOG(INFO) << "Adding variable: pairing_minDeltaR_out_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_minDeltaR_out_NOSYS",
                                          TtbarDileptonVars::PairLeptonsToJets_MinSumDeltaR,
                                          {"el_TLV_NOSYS",
                                           "mu_TLV_NOSYS",
                                           "jet_TLV_NOSYS",
                                           "event_jet_truth_idx",
                                           "event_jet_truth_candidates",
                                           "el_charge",
                                            "mu_charge",
                                           "n_electrons_NOSYS",
                                           "n_muons_NOSYS"});

  // expose scalar extractors for easy histogramming
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_deltaR_lep0_jet_index_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>0? v[0] : -1.0f; },
                                          {"pairing_minDeltaR_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_deltaR_lep1_jet_index_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>1? v[1] : -1.0f; },
                                          {"pairing_minDeltaR_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_deltaR_mjl0_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>2? v[2] : -999999.0f; },
                                          {"pairing_minDeltaR_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_deltaR_mjl1_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>3? v[3] : -999999.0f; },
                                          {"pairing_minDeltaR_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_deltaR_region_code_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>4? v[4] : -9.0f; },
                                          {"pairing_minDeltaR_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_deltaR_correct_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>5? v[5] : -1.0f; },
                                          {"pairing_minDeltaR_out_NOSYS"});

  LOG(INFO) << "Adding variable: pairing_minSumSq_out_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_minSumSq_out_NOSYS",
                                          TtbarDileptonVars::PairLeptonsToTopJets_MinSumSq,
                                          {"el_TLV_NOSYS",
                                           "mu_TLV_NOSYS",
                                           "jet_TLV_NOSYS",
                                           
                                           "event_jet_truth_idx",
                                           "event_jet_truth_candidates",
                                           "el_charge",
                                            "mu_charge",
                                           "n_electrons_NOSYS",
                                           "n_muons_NOSYS"});

  // expose scalar extractors for easy histogramming
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_lep0_jet_index_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>0? v[0] : -1.0f; },
                                          {"pairing_minSumSq_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_lep1_jet_index_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>1? v[1] : -1.0f; },
                                          {"pairing_minSumSq_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_mjl0_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>2? v[2] : -999999.0f; },
                                          {"pairing_minSumSq_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_mjl1_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>3? v[3] : -999999.0f; },
                                          {"pairing_minSumSq_out_NOSYS"});                                        
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_region_code_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>4? v[4] : -9.0f; },
                                          {"pairing_minSumSq_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_correct_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>5? v[5] : -1.0f; },
                                          {"pairing_minSumSq_out_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_mlb_sum_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>6? v[6] : -9999999.0f; },
                                          {"pairing_minSumSq_out_NOSYS"});                                       
LOG(INFO) << "Added variable: mT_ttbar" << std::endl;
mainNode = m_partonMatcher->systematicDefine(mainNode,
                                       "mT_ttbar_NOSYS",
                                       TtbarDileptonVars::mT_ttbar_truth_from_truthIndices,
                                       {"jet_TLV_NOSYS",
                                        "el_TLV_NOSYS",
                                        "mu_TLV_NOSYS",
                                        "event_jet_truth_idx",
                                        "event_jet_truth_candidates",
                                        "el_charge",
                                        "mu_charge",
                                        "met_met_NOSYS",
                                        "met_phi_NOSYS"});
                                        
  // Chi2-based pairing using truth-driven mlb templates (defaults inside function)
  LOG(INFO) << "Adding variable: pairing_chi2_out_all_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_chi2_out_all_NOSYS",
                                           TtbarDileptonVars::PairLeptonsToJets_Chi2TruthAll,
                                          {"el_TLV_NOSYS",
                                           "mu_TLV_NOSYS",
                                           "jet_TLV_NOSYS",
                                           "met_met_NOSYS",
                                            "met_phi_NOSYS",
                                           "event_jet_truth_idx",
                                           "event_jet_truth_candidates",
                                            "el_charge",
                                            "mu_charge",
                                           "n_electrons_NOSYS",
                                           "n_muons_NOSYS"
                                           }
                                          );

  // expose scalar extractors for the chi2 pairing
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_chi2_lep0_jet_index_all_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>0? v[0] : -1.0f; },
                                          {"pairing_chi2_out_all_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_chi2_lep1_jet_index_all_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>1? v[1] : -1.0f; },
                                          {"pairing_chi2_out_all_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_chi2_mjl0_all_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>2? v[2] : -999999.0f; },
                                          {"pairing_chi2_out_all_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_chi2_mjl1_all_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>3? v[3] : -999999.0f; },
                                          {"pairing_chi2_out_all_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_chi2_chi2_all_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>4? v[4] : -1.0f; },
                                          {"pairing_chi2_out_all_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_chi2_correct_all_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>5? v[5] : -1.0f; },
                                          {"pairing_chi2_out_all_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_chi2_truthvalue_all_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>6? v[6] : -1.0f; },
                                          {"pairing_chi2_out_all_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "pairing_chi2_foundvaluediff_all_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>7? v[7] : -1.0f; },
                                          {"pairing_chi2_out_all_NOSYS"});

    // HyPER/nuflows presence check: if inputs don't provide these branches
    // skip all HyPER-dependent defines to keep processing robust.
    {
      auto cols = mainNode.GetColumnNames();
      auto has_col = [&](const std::string &n){ for (const auto &c : cols) if (c == n) return true; return false; };
      const bool has_hyper_inputs = has_col("hyper_TtbarDiLepton_Top1_Indices_NOSYS") || has_col("hyper_TtbarDiLepton_Top1_IDs_NOSYS") || has_col("nuflows_nu_out_NOSYS") || has_col("nuflows_loglik_NOSYS");
      if (has_hyper_inputs) {
        // Unified mttbar per-assignment methods
        auto cols2 = mainNode.GetColumnNames();
        auto has_col2 = [&](const std::string &n){ for (const auto &c : cols2) if (c == n) return true; return false; };

        if (!has_col2("mttbar_HYPER_NOSYS")) {
          mainNode = m_partonMatcher->systematicDefine(mainNode, "mttbar_HYPER_NOSYS",
      [](const std::vector<long>& top1_idx, const std::vector<long>& top1_IDs,
       const std::vector<long>& top2_idx, const std::vector<long>& top2_IDs,
       const ROOT::VecOps::RVec<TLV>& jet_TLV,
       const ROOT::VecOps::RVec<TLV>& el_TLV,
       const ROOT::VecOps::RVec<TLV>& mu_TLV,
       const ROOT::VecOps::RVec<float>& nuflows)->float {
        auto nus = TtbarDileptonVars::buildNuFromNuFlows(nuflows);
        std::vector<long> tops; tops.insert(tops.end(), top1_idx.begin(), top1_idx.end()); tops.insert(tops.end(), top2_idx.begin(), top2_idx.end());
        std::vector<long> ids; ids.insert(ids.end(), top1_IDs.begin(), top1_IDs.end()); ids.insert(ids.end(), top2_IDs.begin(), top2_IDs.end());
        auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::HYPER, tops, ids, jet_TLV, el_TLV, mu_TLV, nus.first, nus.second);
        if (!asg.valid) return -1.0f;
        std::vector<TLV> leps;
        if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
        else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
        else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
        if (leps.size() != 2) return -1.0f;
        auto b1 = jet_TLV.size() > static_cast<size_t>(asg.bjet1_idx) && asg.bjet1_idx>=0 ? jet_TLV[asg.bjet1_idx] : TLV(0,0,0,0);
        auto b2 = jet_TLV.size() > static_cast<size_t>(asg.bjet2_idx) && asg.bjet2_idx>=0 ? jet_TLV[asg.bjet2_idx] : TLV(0,0,0,0);
        auto l1 = leps.size() > static_cast<size_t>(asg.lep1_idx) && asg.lep1_idx>=0 ? leps[asg.lep1_idx] : TLV(0,0,0,0);
        auto l2 = leps.size() > static_cast<size_t>(asg.lep2_idx) && asg.lep2_idx>=0 ? leps[asg.lep2_idx] : TLV(0,0,0,0);
        auto l1m = ROOT::Math::PtEtaPhiMVector(l1.Pt(), l1.Eta(), l1.Phi(), l1.M());
        auto l2m = ROOT::Math::PtEtaPhiMVector(l2.Pt(), l2.Eta(), l2.Phi(), l2.M());
        auto b1m = ROOT::Math::PtEtaPhiMVector(b1.Pt(), b1.Eta(), b1.Phi(), b1.M());
        auto b2m = ROOT::Math::PtEtaPhiMVector(b2.Pt(), b2.Eta(), b2.Phi(), b2.M());
        return TtbarDileptonVars::compute_mttbar(l1m, l2m, b1m, b2m, nus.first, nus.second);
      },
      {"hyper_TtbarDiLepton_Top1_Indices_NOSYS","hyper_TtbarDiLepton_Top1_IDs_NOSYS","hyper_TtbarDiLepton_Top2_Indices_NOSYS","hyper_TtbarDiLepton_Top2_IDs_NOSYS","jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS","nuflows_nu_out_NOSYS"});
        }

        if (!has_col2("mttbar_CHI2_NOSYS")) {
          mainNode = m_partonMatcher->systematicDefine(mainNode, "mttbar_CHI2_NOSYS",
      [](const ROOT::VecOps::RVec<TLV>& jet_TLV, const ROOT::VecOps::RVec<TLV>& el_TLV, const ROOT::VecOps::RVec<TLV>& mu_TLV, const ROOT::VecOps::RVec<float>& nuflows)->float{
        auto nus = TtbarDileptonVars::buildNuFromNuFlows(nuflows);
        auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::CHI2, {}, {}, jet_TLV, el_TLV, mu_TLV, nus.first, nus.second);
        if (!asg.valid) return -1.0f;
        std::vector<TLV> leps;
        if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
        else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
        else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
        if (leps.size() != 2) return -1.0f;
        auto b1 = jet_TLV.size() > static_cast<size_t>(asg.bjet1_idx) && asg.bjet1_idx>=0 ? jet_TLV[asg.bjet1_idx] : TLV(0,0,0,0);
        auto b2 = jet_TLV.size() > static_cast<size_t>(asg.bjet2_idx) && asg.bjet2_idx>=0 ? jet_TLV[asg.bjet2_idx] : TLV(0,0,0,0);
        auto l1 = leps.size() > static_cast<size_t>(asg.lep1_idx) && asg.lep1_idx>=0 ? leps[asg.lep1_idx] : TLV(0,0,0,0);
        auto l2 = leps.size() > static_cast<size_t>(asg.lep2_idx) && asg.lep2_idx>=0 ? leps[asg.lep2_idx] : TLV(0,0,0,0);
        auto l1m = ROOT::Math::PtEtaPhiMVector(l1.Pt(), l1.Eta(), l1.Phi(), l1.M());
        auto l2m = ROOT::Math::PtEtaPhiMVector(l2.Pt(), l2.Eta(), l2.Phi(), l2.M());
        auto b1m = ROOT::Math::PtEtaPhiMVector(b1.Pt(), b1.Eta(), b1.Phi(), b1.M());
        auto b2m = ROOT::Math::PtEtaPhiMVector(b2.Pt(), b2.Eta(), b2.Phi(), b2.M());
        return TtbarDileptonVars::compute_mttbar(l1m, l2m, b1m, b2m, nus.first, nus.second);
      },
      {"jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS","nuflows_nu_out_NOSYS"});
        }

        if (!has_col2("mttbar_PT_NOSYS")) {
          mainNode = m_partonMatcher->systematicDefine(mainNode, "mttbar_PT_NOSYS",
      [](const ROOT::VecOps::RVec<TLV>& jet_TLV, const ROOT::VecOps::RVec<TLV>& el_TLV, const ROOT::VecOps::RVec<TLV>& mu_TLV)->float{
        auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::PT, {}, {}, jet_TLV, el_TLV, mu_TLV, ROOT::Math::PtEtaPhiMVector(0,0,0,0), ROOT::Math::PtEtaPhiMVector(0,0,0,0));
        if (!asg.valid) return -1.0f;
        std::vector<TLV> leps;
        if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
        else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
        else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
        if (leps.size() != 2) return -1.0f;
        auto b1 = jet_TLV.size() > static_cast<size_t>(asg.bjet1_idx) && asg.bjet1_idx>=0 ? jet_TLV[asg.bjet1_idx] : TLV(0,0,0,0);
        auto b2 = jet_TLV.size() > static_cast<size_t>(asg.bjet2_idx) && asg.bjet2_idx>=0 ? jet_TLV[asg.bjet2_idx] : TLV(0,0,0,0);
        auto l1 = leps.size() > static_cast<size_t>(asg.lep1_idx) && asg.lep1_idx>=0 ? leps[asg.lep1_idx] : TLV(0,0,0,0);
        auto l2 = leps.size() > static_cast<size_t>(asg.lep2_idx) && asg.lep2_idx>=0 ? leps[asg.lep2_idx] : TLV(0,0,0,0);
        auto l1m = ROOT::Math::PtEtaPhiMVector(l1.Pt(), l1.Eta(), l1.Phi(), l1.M());
        auto l2m = ROOT::Math::PtEtaPhiMVector(l2.Pt(), l2.Eta(), l2.Phi(), l2.M());
        auto b1m = ROOT::Math::PtEtaPhiMVector(b1.Pt(), b1.Eta(), b1.Phi(), b1.M());
        auto b2m = ROOT::Math::PtEtaPhiMVector(b2.Pt(), b2.Eta(), b2.Phi(), b2.M());
        return TtbarDileptonVars::compute_mttbar(l1m, l2m, b1m, b2m, ROOT::Math::PtEtaPhiMVector(0,0,0,0), ROOT::Math::PtEtaPhiMVector(0,0,0,0));
      },
      {"jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS"});
        }

        if (!has_col2("mttbar_DELTAR_NOSYS")) {
          mainNode = m_partonMatcher->systematicDefine(mainNode, "mttbar_DELTAR_NOSYS",
      [](const ROOT::VecOps::RVec<TLV>& jet_TLV, const ROOT::VecOps::RVec<TLV>& el_TLV, const ROOT::VecOps::RVec<TLV>& mu_TLV)->float{
        auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::DELTAR, {}, {}, jet_TLV, el_TLV, mu_TLV, ROOT::Math::PtEtaPhiMVector(0,0,0,0), ROOT::Math::PtEtaPhiMVector(0,0,0,0));
        if (!asg.valid) return -1.0f;
        std::vector<TLV> leps;
        if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
        else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
        else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
        if (leps.size() != 2) return -1.0f;
        auto b1 = jet_TLV.size() > static_cast<size_t>(asg.bjet1_idx) && asg.bjet1_idx>=0 ? jet_TLV[asg.bjet1_idx] : TLV(0,0,0,0);
        auto b2 = jet_TLV.size() > static_cast<size_t>(asg.bjet2_idx) && asg.bjet2_idx>=0 ? jet_TLV[asg.bjet2_idx] : TLV(0,0,0,0);
        auto l1 = leps.size() > static_cast<size_t>(asg.lep1_idx) && asg.lep1_idx>=0 ? leps[asg.lep1_idx] : TLV(0,0,0,0);
        auto l2 = leps.size() > static_cast<size_t>(asg.lep2_idx) && asg.lep2_idx>=0 ? leps[asg.lep2_idx] : TLV(0,0,0,0);
        auto l1m = ROOT::Math::PtEtaPhiMVector(l1.Pt(), l1.Eta(), l1.Phi(), l1.M());
        auto l2m = ROOT::Math::PtEtaPhiMVector(l2.Pt(), l2.Eta(), l2.Phi(), l2.M());
        auto b1m = ROOT::Math::PtEtaPhiMVector(b1.Pt(), b1.Eta(), b1.Phi(), b1.M());
        auto b2m = ROOT::Math::PtEtaPhiMVector(b2.Pt(), b2.Eta(), b2.Phi(), b2.M());
        return TtbarDileptonVars::compute_mttbar(l1m, l2m, b1m, b2m, ROOT::Math::PtEtaPhiMVector(0,0,0,0), ROOT::Math::PtEtaPhiMVector(0,0,0,0));
      },
      {"jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS"});
        }

        // Add NUWEIGHT-based mttbar variants when scalar neutrino branches exist
        const bool has_nu_scalar_inputs = has_col2("nu_pt_NOSYS") && has_col2("nu_eta_NOSYS") && has_col2("nu_phi_NOSYS") && has_col2("nu_e_NOSYS") && has_col2("nubar_pt_NOSYS") && has_col2("nubar_eta_NOSYS") && has_col2("nubar_phi_NOSYS") && has_col2("nubar_e_NOSYS");
        if (has_nu_scalar_inputs) {
          if (!has_col2("mttbar_HYPER_NUWEIGHT_NOSYS")) {
            mainNode = m_partonMatcher->systematicDefine(mainNode, "mttbar_HYPER_NUWEIGHT_NOSYS",
        [](const std::vector<long>& top1_idx, const std::vector<long>& top1_IDs,
           const std::vector<long>& top2_idx, const std::vector<long>& top2_IDs,
           const ROOT::VecOps::RVec<TLV>& jet_TLV,
           const ROOT::VecOps::RVec<TLV>& el_TLV,
           const ROOT::VecOps::RVec<TLV>& mu_TLV,
           float nu_pt, float nu_eta, float nu_phi, float nu_e,
           float nubar_pt, float nubar_eta, float nubar_phi, float nubar_e)->float {
              auto nu1 = TtbarDileptonVars::buildNuFromBranches(nu_pt, nu_eta, nu_phi, nu_e);
              auto nu2 = TtbarDileptonVars::buildNuFromBranches(nubar_pt, nubar_eta, nubar_phi, nubar_e);
              std::vector<long> tops; tops.insert(tops.end(), top1_idx.begin(), top1_idx.end()); tops.insert(tops.end(), top2_idx.begin(), top2_idx.end());
              std::vector<long> ids; ids.insert(ids.end(), top1_IDs.begin(), top1_IDs.end()); ids.insert(ids.end(), top2_IDs.begin(), top2_IDs.end());
              auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::HYPER, tops, ids, jet_TLV, el_TLV, mu_TLV, nu1, nu2);
              if (!asg.valid) return -1.0f;
              std::vector<TLV> leps;
              if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
              else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
              else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
              if (leps.size() != 2) return -1.0f;
              auto b1 = jet_TLV.size() > static_cast<size_t>(asg.bjet1_idx) && asg.bjet1_idx>=0 ? jet_TLV[asg.bjet1_idx] : TLV(0,0,0,0);
              auto b2 = jet_TLV.size() > static_cast<size_t>(asg.bjet2_idx) && asg.bjet2_idx>=0 ? jet_TLV[asg.bjet2_idx] : TLV(0,0,0,0);
              auto l1 = leps.size() > static_cast<size_t>(asg.lep1_idx) && asg.lep1_idx>=0 ? leps[asg.lep1_idx] : TLV(0,0,0,0);
              auto l2 = leps.size() > static_cast<size_t>(asg.lep2_idx) && asg.lep2_idx>=0 ? leps[asg.lep2_idx] : TLV(0,0,0,0);
              auto l1m = ROOT::Math::PtEtaPhiMVector(l1.Pt(), l1.Eta(), l1.Phi(), l1.M());
              auto l2m = ROOT::Math::PtEtaPhiMVector(l2.Pt(), l2.Eta(), l2.Phi(), l2.M());
              auto b1m = ROOT::Math::PtEtaPhiMVector(b1.Pt(), b1.Eta(), b1.Phi(), b1.M());
              auto b2m = ROOT::Math::PtEtaPhiMVector(b2.Pt(), b2.Eta(), b2.Phi(), b2.M());
              return TtbarDileptonVars::compute_mttbar(l1m, l2m, b1m, b2m, nu1, nu2);
        },
        {"hyper_TtbarDiLepton_Top1_Indices_NOSYS","hyper_TtbarDiLepton_Top1_IDs_NOSYS","hyper_TtbarDiLepton_Top2_Indices_NOSYS","hyper_TtbarDiLepton_Top2_IDs_NOSYS","jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS",
         "nu_pt_NOSYS","nu_eta_NOSYS","nu_phi_NOSYS","nu_e_NOSYS","nubar_pt_NOSYS","nubar_eta_NOSYS","nubar_phi_NOSYS","nubar_e_NOSYS"});
          }

          if (!has_col2("mttbar_CHI2_NUWEIGHT_NOSYS")) {
            mainNode = m_partonMatcher->systematicDefine(mainNode, "mttbar_CHI2_NUWEIGHT_NOSYS",
        [](const ROOT::VecOps::RVec<TLV>& jet_TLV, const ROOT::VecOps::RVec<TLV>& el_TLV, const ROOT::VecOps::RVec<TLV>& mu_TLV,
           float nu_pt, float nu_eta, float nu_phi, float nu_e,
           float nubar_pt, float nubar_eta, float nubar_phi, float nubar_e)->float{
              auto nu1 = TtbarDileptonVars::buildNuFromBranches(nu_pt, nu_eta, nu_phi, nu_e);
              auto nu2 = TtbarDileptonVars::buildNuFromBranches(nubar_pt, nubar_eta, nubar_phi, nubar_e);
              auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::CHI2, {}, {}, jet_TLV, el_TLV, mu_TLV, nu1, nu2);
              if (!asg.valid) return -1.0f;
              std::vector<TLV> leps;
              if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
              else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
              else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
              if (leps.size() != 2) return -1.0f;
              auto b1 = jet_TLV.size() > static_cast<size_t>(asg.bjet1_idx) && asg.bjet1_idx>=0 ? jet_TLV[asg.bjet1_idx] : TLV(0,0,0,0);
              auto b2 = jet_TLV.size() > static_cast<size_t>(asg.bjet2_idx) && asg.bjet2_idx>=0 ? jet_TLV[asg.bjet2_idx] : TLV(0,0,0,0);
              auto l1 = leps.size() > static_cast<size_t>(asg.lep1_idx) && asg.lep1_idx>=0 ? leps[asg.lep1_idx] : TLV(0,0,0,0);
              auto l2 = leps.size() > static_cast<size_t>(asg.lep2_idx) && asg.lep2_idx>=0 ? leps[asg.lep2_idx] : TLV(0,0,0,0);
              auto l1m = ROOT::Math::PtEtaPhiMVector(l1.Pt(), l1.Eta(), l1.Phi(), l1.M());
              auto l2m = ROOT::Math::PtEtaPhiMVector(l2.Pt(), l2.Eta(), l2.Phi(), l2.M());
              auto b1m = ROOT::Math::PtEtaPhiMVector(b1.Pt(), b1.Eta(), b1.Phi(), b1.M());
              auto b2m = ROOT::Math::PtEtaPhiMVector(b2.Pt(), b2.Eta(), b2.Phi(), b2.M());
              return TtbarDileptonVars::compute_mttbar(l1m, l2m, b1m, b2m, nu1, nu2);
        },
        {"jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS",
         "nu_pt_NOSYS","nu_eta_NOSYS","nu_phi_NOSYS","nu_e_NOSYS","nubar_pt_NOSYS","nubar_eta_NOSYS","nubar_phi_NOSYS","nubar_e_NOSYS"});
          }

          if (!has_col2("mttbar_PT_NUWEIGHT_NOSYS")) {
            mainNode = m_partonMatcher->systematicDefine(mainNode, "mttbar_PT_NUWEIGHT_NOSYS",
        [](const ROOT::VecOps::RVec<TLV>& jet_TLV, const ROOT::VecOps::RVec<TLV>& el_TLV, const ROOT::VecOps::RVec<TLV>& mu_TLV,
           float nu_pt, float nu_eta, float nu_phi, float nu_e,
           float nubar_pt, float nubar_eta, float nubar_phi, float nubar_e)->float{
              auto nu1 = TtbarDileptonVars::buildNuFromBranches(nu_pt, nu_eta, nu_phi, nu_e);
              auto nu2 = TtbarDileptonVars::buildNuFromBranches(nubar_pt, nubar_eta, nubar_phi, nubar_e);
              auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::PT, {}, {}, jet_TLV, el_TLV, mu_TLV, nu1, nu2);
              if (!asg.valid) return -1.0f;
              std::vector<TLV> leps;
              if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
              else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
              else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
              if (leps.size() != 2) return -1.0f;
              auto b1 = jet_TLV.size() > static_cast<size_t>(asg.bjet1_idx) && asg.bjet1_idx>=0 ? jet_TLV[asg.bjet1_idx] : TLV(0,0,0,0);
              auto b2 = jet_TLV.size() > static_cast<size_t>(asg.bjet2_idx) && asg.bjet2_idx>=0 ? jet_TLV[asg.bjet2_idx] : TLV(0,0,0,0);
              auto l1 = leps.size() > static_cast<size_t>(asg.lep1_idx) && asg.lep1_idx>=0 ? leps[asg.lep1_idx] : TLV(0,0,0,0);
              auto l2 = leps.size() > static_cast<size_t>(asg.lep2_idx) && asg.lep2_idx>=0 ? leps[asg.lep2_idx] : TLV(0,0,0,0);
              auto l1m = ROOT::Math::PtEtaPhiMVector(l1.Pt(), l1.Eta(), l1.Phi(), l1.M());
              auto l2m = ROOT::Math::PtEtaPhiMVector(l2.Pt(), l2.Eta(), l2.Phi(), l2.M());
              auto b1m = ROOT::Math::PtEtaPhiMVector(b1.Pt(), b1.Eta(), b1.Phi(), b1.M());
              auto b2m = ROOT::Math::PtEtaPhiMVector(b2.Pt(), b2.Eta(), b2.Phi(), b2.M());
              return TtbarDileptonVars::compute_mttbar(l1m, l2m, b1m, b2m, nu1, nu2);
        },
        {"jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS",
         "nu_pt_NOSYS","nu_eta_NOSYS","nu_phi_NOSYS","nu_e_NOSYS","nubar_pt_NOSYS","nubar_eta_NOSYS","nubar_phi_NOSYS","nubar_e_NOSYS"});
          }

          if (!has_col2("mttbar_DELTAR_NUWEIGHT_NOSYS")) {
            mainNode = m_partonMatcher->systematicDefine(mainNode, "mttbar_DELTAR_NUWEIGHT_NOSYS",
        [](const ROOT::VecOps::RVec<TLV>& jet_TLV, const ROOT::VecOps::RVec<TLV>& el_TLV, const ROOT::VecOps::RVec<TLV>& mu_TLV,
           float nu_pt, float nu_eta, float nu_phi, float nu_e,
           float nubar_pt, float nubar_eta, float nubar_phi, float nubar_e)->float{
              auto nu1 = TtbarDileptonVars::buildNuFromBranches(nu_pt, nu_eta, nu_phi, nu_e);
              auto nu2 = TtbarDileptonVars::buildNuFromBranches(nubar_pt, nubar_eta, nubar_phi, nubar_e);
              auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::DELTAR, {}, {}, jet_TLV, el_TLV, mu_TLV, nu1, nu2);
              if (!asg.valid) return -1.0f;
              std::vector<TLV> leps;
              if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
              else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
              else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
              if (leps.size() != 2) return -1.0f;
              auto b1 = jet_TLV.size() > static_cast<size_t>(asg.bjet1_idx) && asg.bjet1_idx>=0 ? jet_TLV[asg.bjet1_idx] : TLV(0,0,0,0);
              auto b2 = jet_TLV.size() > static_cast<size_t>(asg.bjet2_idx) && asg.bjet2_idx>=0 ? jet_TLV[asg.bjet2_idx] : TLV(0,0,0,0);
              auto l1 = leps.size() > static_cast<size_t>(asg.lep1_idx) && asg.lep1_idx>=0 ? leps[asg.lep1_idx] : TLV(0,0,0,0);
              auto l2 = leps.size() > static_cast<size_t>(asg.lep2_idx) && asg.lep2_idx>=0 ? leps[asg.lep2_idx] : TLV(0,0,0,0);
              auto l1m = ROOT::Math::PtEtaPhiMVector(l1.Pt(), l1.Eta(), l1.Phi(), l1.M());
              auto l2m = ROOT::Math::PtEtaPhiMVector(l2.Pt(), l2.Eta(), l2.Phi(), l2.M());
              auto b1m = ROOT::Math::PtEtaPhiMVector(b1.Pt(), b1.Eta(), b1.Phi(), b1.M());
              auto b2m = ROOT::Math::PtEtaPhiMVector(b2.Pt(), b2.Eta(), b2.Phi(), b2.M());
              return TtbarDileptonVars::compute_mttbar(l1m, l2m, b1m, b2m, nu1, nu2);
        },
        {"jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS",
         "nu_pt_NOSYS","nu_eta_NOSYS","nu_phi_NOSYS","nu_e_NOSYS","nubar_pt_NOSYS","nubar_eta_NOSYS","nubar_phi_NOSYS","nubar_e_NOSYS"});
          }
        }
      } else {
        LOG(INFO) << "Skipping HyPER-dependent defines: HyPER/nuflows branches not present in input" << std::endl;
      }
    }

  LOG(INFO) << "Adding variable: parton_truth_had_b_index_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode, "parton_truth_had_b_index_NOSYS",
                                          TtbarDileptonVars::getTruthIndexBHad,
                                          {"event_jet_truth_idx",
                                           "event_jet_truth_candidates"
                                           });
  LOG(INFO) << "Adding variable: parton_truth_had_bbar_index_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode, "parton_truth_had_bbar_index_NOSYS",
                                          TtbarDileptonVars::getTruthIndexBbarHad,
                                          {"event_jet_truth_idx",
                                           "event_jet_truth_candidates"
                                           });                                         
  // Fully matched mlb pairing masses added
  mainNode = m_partonMatcher->systematicDefine(mainNode, "FullyMatchedMLB_NOSYS",
                                          TtbarDileptonVars::FullyMatched_MLB,
                                          {"jet_TLV_NOSYS",                                    
                                          "el_TLV_NOSYS",
                                          "mu_TLV_NOSYS",
                                          "el_charge", "mu_charge",
                                          "event_jet_truth_idx",
                                          "event_jet_truth_candidates"
                                        });

  LOG(INFO) << "fully matched mlb truth assigned" << std::endl;
  mainNode= m_partonMatcher->systematicDefine(mainNode, "mlb_truth_b_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>0?v[0]:-1.0f;  },
                                          {"FullyMatchedMLB_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "mlb_truth_bbar_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>1?v[1]:-1.0f;  },
                                          {"FullyMatchedMLB_NOSYS"});        
  mainNode = m_partonMatcher->systematicDefine(mainNode, "mlb_alt_b_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>2?v[2]:-1.0f;  },
                                          {"FullyMatchedMLB_NOSYS"});        
  mainNode = m_partonMatcher->systematicDefine(mainNode, "mlb_alt_bbar_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>3?v[3]:-1.0f;  },
                                          {"FullyMatchedMLB_NOSYS"});
  mainNode = m_partonMatcher->systematicDefine(mainNode, "mlb_pTdiff_true_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>4?v[4]:-1.0f;  },
                                          {"FullyMatchedMLB_NOSYS"});                                        
  mainNode = m_partonMatcher->systematicDefine(mainNode, "mlb_sum_deltaR_true_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>5?v[5]:-1.0f;  },
                                          {"FullyMatchedMLB_NOSYS"});           
  mainNode = m_partonMatcher->systematicDefine(mainNode, "true_b_pT_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>5?v[6]:-1.0f;  },
                                          {"FullyMatchedMLB_NOSYS"});           
  mainNode = m_partonMatcher->systematicDefine(mainNode, "true_bbar_pT_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>5?v[7]:-1.0f;  },
                                          {"FullyMatchedMLB_NOSYS"});           
  mainNode = m_partonMatcher->systematicDefine(mainNode, "true_jet_deltaphi_NOSYS",
                                          [](const ROOT::VecOps::RVec<float>& v){ return v.size()>5?v[8]:-1.0f;  },
                                          {"FullyMatchedMLB_NOSYS"});

  // Invariant mass of the llbb system
  LOG(INFO) << "Adding variable: invariant_mass_llbb_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "invariant_mass_llbb_NOSYS",
                                         TtbarDileptonVars::invariantMassllbb,
                                         {"sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS","sorted_bjet_TLV_NOSYS","jet_TLV_NOSYS","jet_GN2v01_FixedCutBEff_77_select","n_electrons_NOSYS","n_muons_NOSYS"});
// Signal region definition
  LOG(INFO) << "Adding variable: signal_NOSYS" << std::endl;
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "signal_NOSYS",
                                         TtbarDileptonVars::Signal,
                                         {"sorted_el_TLV_NOSYS",
                                          "sorted_mu_TLV_NOSYS",
                                          "sorted_el_charge_NOSYS",
                                          "sorted_mu_charge_NOSYS",
                                          "n_jets_NOSYS",
                                          "n_bjets_NOSYS",
                                          "ttbar_mass_NOSYS",
                                          "min_mlb_pair_masses_NOSYS",
                                          "met_met_NOSYS"
                                        });                    
  LOG(INFO) << "Adding variable: selected_jet_indices_NOSYS" << std::endl;        
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "selected_jet_indices_NOSYS",
  //                                        selectTwoJetIndices_BtagPriority,
  //                                        {"jet_TLV_NOSYS",
  //                                         "jet_GN2v01_FixedCutBEff_77_select"
  //                                       });                                    
  // LOG(INFO) << "Adding variable: selected_match_truth_bpair_NOSYS" << std::endl;        
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                        "selected_match_truth_bpair_NOSYS",
  //                                        TtbarDileptonVars::SelectedJetsMatchTruthBPair,
  //                                        {"jets_ptOrderIndex_NOSYS",
  //                                         "event_jet_truth_idx",
  //                                         "event_jet_truth_candidates"
  //                                       });                                
  
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
  //std::string campaign = id.campaign();
  //std::cout << "Campaign is: " << campaign << std::endl;
  if (id.campaign() == "mc23a") {
    LOG(INFO) << "Sample is mc23a, using mc23a trigger definitions" << std::endl;
    LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "Electron_pass_triggers_NOSYS",
                                          TtbarSingleLeptonVars::ElectronSelection::triggers,
                                          {"trigPassed_HLT_e60_lhmedium_L1EM22VHI", "trigPassed_HLT_e140_lhloose_L1EM22VHI","trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI",
                                          "trigPassed_HLT_e26_lhtight_nod0_ivarloose","trigPassed_HLT_e140_lhloose_nod0","trigPassed_HLT_e60_lhmedium_nod0"});


    LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "Muon_pass_triggers_NOSYS",
                                          TtbarSingleLeptonVars::MuonSelection::triggers,
                                          {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH", "trigPassed_HLT_mu50_L1MU14FCH",
                                          "trigPassed_HLT_mu50","trigPassed_HLT_mu26_ivarmedium"});                                    
    } else if (id.campaign() == "mc23d") {
    LOG(INFO) << "Sample is mc23d, using mc23d trigger definitions" << std::endl;
    LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "Electron_pass_triggers_NOSYS",
                                          TtbarSingleLeptonVars::ElectronSelection::triggers,
                                          {"trigPassed_HLT_e60_lhmedium_L1EM22VHI", "trigPassed_HLT_e140_lhloose_L1EM22VHI","trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI",
                                          "trigPassed_HLT_e26_lhtight_nod0_ivarloose","trigPassed_HLT_e140_lhloose_nod0","trigPassed_HLT_e60_lhmedium_nod0"});


    LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "Muon_pass_triggers_NOSYS",
                                          TtbarSingleLeptonVars::MuonSelection::triggers,
                                          {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH", "trigPassed_HLT_mu50_L1MU14FCH",
                                          "trigPassed_HLT_mu50","trigPassed_HLT_mu26_ivarmedium"});
                                        
  } else if (id.campaign() == "mc23e") {
    LOG(INFO) << "Sample is mc23e, using mc23e trigger definitions" << std::endl;
    LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "Electron_pass_triggers_NOSYS",
                                          TtbarSingleLeptonVars::ElectronSelection::triggersnew,
                                          {"trigPassed_HLT_e140_lhloose_L1eEM26M",
    "trigPassed_HLT_e26_lhtight_ivarloose_L1eEM26M",
    "trigPassed_HLT_e60_lhmedium_L1eEM26M"});
    
  

    LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "Muon_pass_triggers_NOSYS",
                                          TtbarSingleLeptonVars::MuonSelection::triggersnew,
                                          {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH",
  "trigPassed_HLT_mu50_L1MU14FCH"});
  } 
  else if (id.campaign() == "mc20a" || id.campaign() == "mc20d" || id.campaign() == "mc20e") {
      LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "Electron_pass_triggers_NOSYS",
                                         TtbarSingleLeptonVars::ElectronSelection::triggers,
                                         {"trigPassed_HLT_e60_lhmedium_L1EM22VHI", "trigPassed_HLT_e140_lhloose_L1EM22VHI","trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI",
                                         "trigPassed_HLT_e26_lhtight_nod0_ivarloose","trigPassed_HLT_e140_lhloose_nod0","trigPassed_HLT_e60_lhmedium_nod0"});

    LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "Muon_pass_triggers_NOSYS",
                                         TtbarSingleLeptonVars::MuonSelection::triggers,
                                         {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH", "trigPassed_HLT_mu50_L1MU14FCH",
                                         "trigPassed_HLT_mu50","trigPassed_HLT_mu26_ivarmedium"});

  }
  else {
    LOG(INFO) << "Sample is neither mc23d nor mc23e, using mc23e trigger definitions by default" << std::endl;
    
    LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                        "Electron_pass_triggers_NOSYS",
                        TtbarSingleLeptonVars::ElectronSelection::triggersnew,
                        {"trigPassed_HLT_e140_lhloose_L1eEM26M",
    "trigPassed_HLT_e26_lhtight_ivarloose_L1eEM26M",
    "trigPassed_HLT_e60_lhmedium_L1eEM26M"});


    LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
    mainNode = m_partonMatcher->systematicDefine(mainNode,
                                          "Muon_pass_triggers_NOSYS",
                                          TtbarSingleLeptonVars::MuonSelection::triggersnew,
                                          {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH",
  "trigPassed_HLT_mu50_L1MU14FCH"});}
  //   LOG(INFO) << "Adding variable: Electron_pass_triggers_NOSYS\n";
  //   mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                         "Electron_pass_triggers_NOSYS",
  //                                         TtbarSingleLeptonVars::ElectronSelection::triggersnew,
  //                                         {"trigPassed_HLT_e140_lhloose_L1eEM26M",
  //   "trigPassed_HLT_e26_lhtight_ivarloose_L1eEM",
  //   "trigPassed_HLT_e60_lhmedium_L1eEM"});


  //   LOG(INFO) << "Adding variable: Muon_pass_triggers_NOSYS\n";
  //   mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                                         "Muon_pass_triggers_NOSYS",
  //                                         TtbarSingleLeptonVars::MuonSelection::triggersnew,
  //                                         {"trigPassed_HLT_mu24_ivarmedium_L1MU14FCH",
  // "trigPassed_HLT_mu50_L1MU14FCH"});                                        
  //                                       }
  
  // create combined OR of electron OR muon trigger flags, signal region requires one lepton to have passed trigger
LOG(INFO) << "Adding variable: SingleLepton_pass_triggers_NOSYS\n";

auto orTriggers = [](bool e_pass, bool mu_pass) {
    return (e_pass || mu_pass) ? 1 : 0;
};


mainNode = m_partonMatcher->systematicDefine(mainNode,
                          "SingleLepton_pass_triggers_NOSYS",
                          orTriggers,
                          {"Electron_pass_triggers_NOSYS", "Muon_pass_triggers_NOSYS"});

                          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////// Define the Fully matched selection cuts ////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fully matched definition                                       
  LOG(INFO) << "Adding variable: is_fully_matched_NOSYS" << std::endl;     
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "is_fully_matched_NOSYS",
                                         TtbarDileptonVars::FullyMatched,
                                         {"event_jet_truth_idx",
                                          "event_jet_truth_candidates"
                                        });  
  
  // First define the column with the correct value
 mainNode = m_partonMatcher->systematicDefine(
    mainNode,
    "mcChannelNumber_NOSYS",
    [](UInt_t mcID){ return (mcID == 601230) ? 1.0f : 0.0f; },
    {"mcChannelNumber"} // column your lambda reads
);

  // Then register it in the systematic setup
  mainNode = m_partonMatcher->systematicDefine(
      mainNode,
      "isTtbarDilepton_NOSYS",
      [](const float& v){ return v; },
      {"mcChannelNumber_NOSYS"} // use the column name, not a number
  );

  // LOG(INFO) << "id.dsid\n";   
  // std::cout << id.dsid() << std::endl; //returns sample DSID
                  
  // LOG(INFO) << "mcChannelNumber\n";
  // std::cout << id.mcChannelNumber() << std::endl; //returns sample MC channel number
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////// btag //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  LOG(INFO) << "Adding variable: jet_GN2v01_Continuous_quantile_NOSYS\n";   
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                                         "jet_GN2v01_Continuous_quantile_NOSYS",
                                         Variables::getPCBT,
                                         {"jet_GN2v01_FixedCutBEff_65_select", "jet_GN2v01_FixedCutBEff_70_select", "jet_GN2v01_FixedCutBEff_77_select", "jet_GN2v01_FixedCutBEff_85_select", "jet_GN2v01_FixedCutBEff_90_select"});
  
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
// Get branch names
  auto columns = mainNode.GetColumnNames();
  // Store the trigger branches
  bool isHyPERAvailable = false;
  for (const auto& column : columns) {
    
    if (column.find("hyper") != std::string::npos) {
      isHyPERAvailable = true;
    }
  }

  if (isHyPERAvailable) {
  LOG(INFO) << "Ensuring hyper_TtbarDiLepton top indices/IDs available\n";
  // If hyper_* columns are present they will be used directly; otherwise
  // create empty fallbacks earlier. We still expose standardized names
  // matching input: `hyper_TtbarDiLepton_Top*_Indices_NOSYS` and `*_IDs_NOSYS`.
  // No need to create HyPER_* aliases here.
  }
                                //muon as placeholder for taus

  // LOG(INFO) << "Adding variable: HyPER_WHad_Indices_NOSYS\n";
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                             "HyPER_WHad_Indices_NOSYS",
  //                             Variables::getSelecHyPERIndices,
  //                             { "hyper_TtbarLJets_WHad_Indices_NOSYS",
  //                               "hyper_TtbarLJets_TopHad_IDs_NOSYS",
  //                               "jets_ptOrderIndex_NOSYS",
  //                               "el_ptOrderIndex_NOSYS",
  //                               "mu_ptOrderIndex_NOSYS",
  //                               "mu_ptOrderIndex_NOSYS"});
  // LOG(INFO) << "Adding variable: HyPER_WLep_Indices_NOSYS\n";
  // mainNode = m_partonMatcher->systematicDefine(mainNode,
  //                             "HyPER_WLep_Indices_NOSYS",
  //                             Variables::getSelecHyPERIndices,
  //                             { "hyper_TtbarLJets_WLep_Indices_NOSYS",
  //                               "hyper_TtbarLJets_TopLep_IDs_NOSYS",
  //                               "jets_ptOrderIndex_NOSYS",
  //                               "el_ptOrderIndex_NOSYS",
  //                               "mu_ptOrderIndex_NOSYS",
  //                               "mu_ptOrderIndex_NOSYS"});

  return mainNode;
}

ROOT::RDF::RNode TtbarDileptonFrameReco::defineVariablesNtuple(ROOT::RDF::RNode mainNode,
                                                    const std::shared_ptr<Sample>& sample,
                                                    const UniqueSampleID& id) {
  // Provide a dummy jet_SV1_masssvx column if the input doesn't have it,
  // to avoid downstream MissingColumn errors in some configs. Do this before
  // calling defineVariables so any subsequent systematicDefine can rely on it.
  {
    auto cols = mainNode.GetColumnNames();
    bool hasSVcol = false;
    for (const auto &c : cols) if (c == "jet_SV1_masssvx") { hasSVcol = true; break; }
    if (!hasSVcol) {
      LOG(INFO) << "Defining fallback jet_SV1_masssvx (zeros) to avoid missing-column errors" << std::endl;
      mainNode = mainNode.Define("jet_SV1_masssvx",
                [](){ ROOT::VecOps::RVec<float> out; return out; });
    }
  }

  mainNode = this->defineVariables(mainNode, sample, id);

  // Define the collections to store.
  // Electrons
  LOG(INFO) << "Adding variable: Electron_pt\n";
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("el_pt_NOSYS")) {
      if (has_col_nt("Electron_pt")) mainNode = mainNode.Redefine("Electron_pt",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"el_pt_NOSYS","el_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Electron_pt",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"el_pt_NOSYS","el_ptOrderIndex_NOSYS"});
    }
  }
  LOG(INFO) << "Adding variable: Electron_eta\n";
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("el_eta")) {
      if (has_col_nt("Electron_eta")) mainNode = mainNode.Redefine("Electron_eta",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"el_eta","el_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Electron_eta",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"el_eta","el_ptOrderIndex_NOSYS"});
    }
  }
  LOG(INFO) << "Adding variable: Electron_phi\n";
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("el_phi")) {
      if (has_col_nt("Electron_phi")) mainNode = mainNode.Redefine("Electron_phi",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"el_phi","el_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Electron_phi",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"el_phi","el_ptOrderIndex_NOSYS"});
    }
  }
  LOG(INFO) << "Adding variable: Electron_e\n";
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("el_e_NOSYS")) {
      if (has_col_nt("Electron_e")) mainNode = mainNode.Redefine("Electron_e",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"el_e_NOSYS","el_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Electron_e",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"el_e_NOSYS","el_ptOrderIndex_NOSYS"});
    }
  }
  LOG(INFO) << "Adding variable: Electron_q\n";
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("el_charge")) {
      if (has_col_nt("Electron_q")) mainNode = mainNode.Redefine("Electron_q",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"el_charge","el_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Electron_q",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"el_charge","el_ptOrderIndex_NOSYS"});
    }
  }
  
  // Muons
  LOG(INFO) << "Adding variable: Muon_pt\n";
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("mu_pt_NOSYS")) {
      if (has_col_nt("Muon_pt")) mainNode = mainNode.Redefine("Muon_pt",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"mu_pt_NOSYS","mu_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Muon_pt",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"mu_pt_NOSYS","mu_ptOrderIndex_NOSYS"});
    }
    if (has_col_nt("mu_eta")) {
      if (has_col_nt("Muon_eta")) mainNode = mainNode.Redefine("Muon_eta",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"mu_eta","mu_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Muon_eta",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"mu_eta","mu_ptOrderIndex_NOSYS"});
    }
    if (has_col_nt("mu_phi")) {
      if (has_col_nt("Muon_phi")) mainNode = mainNode.Redefine("Muon_phi",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"mu_phi","mu_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Muon_phi",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"mu_phi","mu_ptOrderIndex_NOSYS"});
    }
    if (has_col_nt("mu_e_NOSYS")) {
      if (has_col_nt("Muon_e")) mainNode = mainNode.Redefine("Muon_e",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"mu_e_NOSYS","mu_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Muon_e",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"mu_e_NOSYS","mu_ptOrderIndex_NOSYS"});
    }
    if (has_col_nt("mu_charge")) {
      if (has_col_nt("Muon_q")) mainNode = mainNode.Redefine("Muon_q",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"mu_charge","mu_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Muon_q",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"mu_charge","mu_ptOrderIndex_NOSYS"});
    }
  }

  // Jets
  LOG(INFO) << "Adding variable: Jet_pt\n";
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("jet_pt_NOSYS")) {
      if (has_col_nt("Jet_pt")) mainNode = mainNode.Redefine("Jet_pt",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"jet_pt_NOSYS","jets_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Jet_pt",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"jet_pt_NOSYS","jets_ptOrderIndex_NOSYS"});
    }
    if (has_col_nt("jet_eta")) {
      if (has_col_nt("Jet_eta")) mainNode = mainNode.Redefine("Jet_eta",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_eta","jets_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Jet_eta",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_eta","jets_ptOrderIndex_NOSYS"});
    }
    if (has_col_nt("jet_phi")) {
      if (has_col_nt("Jet_phi")) mainNode = mainNode.Redefine("Jet_phi",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_phi","jets_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Jet_phi",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_phi","jets_ptOrderIndex_NOSYS"});
    }
    if (has_col_nt("jet_e_NOSYS")) {
      if (has_col_nt("Jet_e")) mainNode = mainNode.Redefine("Jet_e",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"jet_e_NOSYS","jets_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Jet_e",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(1000.0f),
                              {"jet_e_NOSYS","jets_ptOrderIndex_NOSYS"});
    }
  }
  LOG(INFO) << "Adding variable: Jet_SVmass\n";
  LOG(INFO) << "Jet_SVmass ntuple branch disabled for this config" << std::endl;

  //btag  
  LOG(INFO) << "Adding variable: Jet_GN2v01_score\n";
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("jet_GN2v01")) {
      if (has_col_nt("Jet_GN2v01_score")) mainNode = mainNode.Redefine("Jet_GN2v01_score",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_GN2v01","jets_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Jet_GN2v01_score",
                              TtbarSingleLeptonVars::CreatePtOrderedContainer<float>(),
                              {"jet_GN2v01","jets_ptOrderIndex_NOSYS"});
    }
  }

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
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("met_met_NOSYS")) {
      if (has_col_nt("Met_pt")) mainNode = mainNode.Redefine("Met_pt",
                              [](float met_pt) { return met_pt/1000.f; },
                              {"met_met_NOSYS"});
      else mainNode = mainNode.Define("Met_pt",
                              [](float met_pt) { return met_pt/1000.f; },
                              {"met_met_NOSYS"});
    }
    if (has_col_nt("met_phi_NOSYS")) {
      if (has_col_nt("Met_phi")) mainNode = mainNode.Redefine("Met_phi",
                              [](float met_phi) { return met_phi/1.f; },
                              {"met_phi_NOSYS"});
      else mainNode = mainNode.Define("Met_phi",
                              [](float met_phi) { return met_phi/1.f; },
                              {"met_phi_NOSYS"});
    }
  }

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
  {
    auto cols_nt = mainNode.GetColumnNames();
    auto has_col_nt = [&](const std::string &n){ for (const auto &c : cols_nt) if (c==n) return true; return false; };
    if (has_col_nt("parton_truth_had_b_index_NOSYS")) {
      if (has_col_nt("Truth_had_b_index")) mainNode = mainNode.Redefine("Truth_had_b_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_had_b_index_NOSYS","jets_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Truth_had_b_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_had_b_index_NOSYS","jets_ptOrderIndex_NOSYS"});
    }
    if (has_col_nt("parton_truth_had_bbar_index_NOSYS")) {
      if (has_col_nt("Truth_had_bbar_index")) mainNode = mainNode.Redefine("Truth_had_bbar_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_had_bbar_index_NOSYS","jets_ptOrderIndex_NOSYS"});
      else mainNode = mainNode.Define("Truth_had_bbar_index",
                              Variables::ConverPartonMatchtNtupleIndexToPtIndex(),
                              {"parton_truth_had_bbar_index_NOSYS","jets_ptOrderIndex_NOSYS"});
    }
  }
  
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

  // DiLepton HyPER ntuple variables (truth match helpers, HyPER b-index, convergence, correctness)
  LOG(INFO) << "Adding variable: Top1_Truth_Reco_Indices_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "Top1_Truth_Reco_Indices_NOSYS",
                              HyPERPerformanceDiLepton::getRecoIndicesTop1,
                              {"event_jet_truth_idx", "event_electron_truth_idx", "event_muon_truth_idx",
                               "event_jet_truth_candidates", "event_electron_truth_candidates", "event_muon_truth_candidates",
                               "event_nLeptons"});

  LOG(INFO) << "Adding variable: Top2_Truth_Reco_Indices_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "Top2_Truth_Reco_Indices_NOSYS",
                              HyPERPerformanceDiLepton::getRecoIndicesTop2,
                              {"event_jet_truth_idx", "event_electron_truth_idx", "event_muon_truth_idx",
                               "event_jet_truth_candidates", "event_electron_truth_candidates", "event_muon_truth_candidates",
                               "event_nLeptons"});

  // HyPER b-jet index for each top (uses HyPER IDs to locate the b)
  LOG(INFO) << "Adding variable: HyPER_Top1_bIndex_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_Top1_bIndex_NOSYS",
                              HyPERPerformanceSingleLepton::getHyPERBLepIndex,
                              {"HyPER_Top1_Indices_NOSYS", "hyper_TtbarDiLepton_Top1_IDs_NOSYS"});

  LOG(INFO) << "Adding variable: HyPER_Top2_bIndex_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_Top2_bIndex_NOSYS",
                              HyPERPerformanceSingleLepton::getHyPERBLepIndex,
                              {"HyPER_Top2_Indices_NOSYS", "hyper_TtbarDiLepton_Top2_IDs_NOSYS"});

  // HyPER convergence flags
  LOG(INFO) << "Adding variable: HyPER_Top1_converged_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_Top1_converged_NOSYS",
                              HyPERPerformanceDiLepton::HyPERTopConverged,
                              {"HyPER_Top1_Indices_NOSYS", "hyper_TtbarDiLepton_Top1_IDs_NOSYS"});

  LOG(INFO) << "Adding variable: HyPER_Top2_converged_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_Top2_converged_NOSYS",
                              HyPERPerformanceDiLepton::HyPERTopConverged,
                              {"HyPER_Top2_Indices_NOSYS", "hyper_TtbarDiLepton_Top2_IDs_NOSYS"});

  LOG(INFO) << "Adding variable: HyPER_two_tops_converged_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_two_tops_converged_NOSYS",
                              [](const std::vector<long>& a, const std::vector<long>& b, const std::vector<long>& c, const std::vector<long>& d){
                                  return HyPERPerformanceDiLepton::HyPERConvergedTwoLeptonicTops(a,b,c,d);
                              },
                              {"HyPER_Top1_Indices_NOSYS","hyper_TtbarDiLepton_Top1_IDs_NOSYS","HyPER_Top2_Indices_NOSYS","hyper_TtbarDiLepton_Top2_IDs_NOSYS"});

  // Correctness checks: compare HyPER b-index with truth-derived top indices
  LOG(INFO) << "Adding variable: HyPER_Top1_correct_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_Top1_correct_NOSYS",
                              HyPERPerformanceDiLepton::correctlyReconstructedTop,
                              {"HyPER_Top1_bIndex_NOSYS","Top1_Truth_Reco_Indices_NOSYS"});

  LOG(INFO) << "Adding variable: HyPER_Top2_correct_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_Top2_correct_NOSYS",
                              HyPERPerformanceDiLepton::correctlyReconstructedTop,
                              {"HyPER_Top2_bIndex_NOSYS","Top2_Truth_Reco_Indices_NOSYS"});

  LOG(INFO) << "Adding variable: HyPER_both_tops_correct_NOSYS\n";
  mainNode = m_partonMatcher->systematicDefine(mainNode,
                              "HyPER_both_tops_correct_NOSYS",
                              [](bool a, bool b){ return a && b; },
                              {"HyPER_Top1_correct_NOSYS","HyPER_Top2_correct_NOSYS"});

  // Note: avoid redundant non-systematic alias Defines here.
  // The systematic *_NOSYS columns are exposed directly by the ntupling
  // configuration. Keeping these alias Defines caused duplicate-Define
  // errors and bloated the frame. Upstream Snapshot/ntuple selection
  // should reference the *_NOSYS columns or define lightweight aliases
  // centrally if required.

  // Per-assignment-method correctness flags (NOSYS) - CHI2/PT/DELTAR use assignment comparison against truth
  {
    auto cols3 = mainNode.GetColumnNames();
    auto has_col3 = [&](const std::string &n){ for (const auto &c : cols3) if (c == n) return true; return false; };
    if (!has_col3("CHI2_correct_NOSYS")) {
      mainNode = m_partonMatcher->systematicDefine(mainNode, "CHI2_correct_NOSYS",
        [](const ROOT::VecOps::RVec<TLV>& jet_TLV, const ROOT::VecOps::RVec<TLV>& el_TLV, const ROOT::VecOps::RVec<TLV>& mu_TLV,
           const ROOT::VecOps::RVec<float>& nuflows, const std::vector<int>& top1_truth, const std::vector<int>& top2_truth)->bool{
          auto nus = TtbarDileptonVars::buildNuFromNuFlows(nuflows);
          auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::CHI2, {}, {}, jet_TLV, el_TLV, mu_TLV, nus.first, nus.second);
          if (!asg.valid) return false;
          int b1 = asg.bjet1_idx; int b2 = asg.bjet2_idx;
          if (top1_truth.size()==0 || top2_truth.size()==0) return false;
          int t1b = top1_truth.size()>0? top1_truth[0] : -1;
          int t2b = top2_truth.size()>0? top2_truth[0] : -1;
          if (t1b<0 || t2b<0) return false;
          return (b1==t1b && b2==t2b) || (b1==t2b && b2==t1b);
        },
        {"jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS","nuflows_nu_out_NOSYS","Top1_Truth_Reco_Indices_NOSYS","Top2_Truth_Reco_Indices_NOSYS"});
    }
  }

  {
    auto cols4 = mainNode.GetColumnNames();
    auto has_col4 = [&](const std::string &n){ for (const auto &c : cols4) if (c == n) return true; return false; };
    if (!has_col4("PT_correct_NOSYS")) {
      mainNode = m_partonMatcher->systematicDefine(mainNode, "PT_correct_NOSYS",
        [](const ROOT::VecOps::RVec<TLV>& jet_TLV, const ROOT::VecOps::RVec<TLV>& el_TLV, const ROOT::VecOps::RVec<TLV>& mu_TLV,
           const std::vector<int>& top1_truth, const std::vector<int>& top2_truth)->bool{
          auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::PT, {}, {}, jet_TLV, el_TLV, mu_TLV, ROOT::Math::PtEtaPhiMVector(0,0,0,0), ROOT::Math::PtEtaPhiMVector(0,0,0,0));
          if (!asg.valid) return false;
          int b1 = asg.bjet1_idx; int b2 = asg.bjet2_idx;
          if (top1_truth.size()==0 || top2_truth.size()==0) return false;
          int t1b = top1_truth.size()>0? top1_truth[0] : -1;
          int t2b = top2_truth.size()>0? top2_truth[0] : -1;
          if (t1b<0 || t2b<0) return false;
          return (b1==t1b && b2==t2b) || (b1==t2b && b2==t1b);
        },
        {"jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS","Top1_Truth_Reco_Indices_NOSYS","Top2_Truth_Reco_Indices_NOSYS"});
    }
  }

  {
    auto cols5 = mainNode.GetColumnNames();
    auto has_col5 = [&](const std::string &n){ for (const auto &c : cols5) if (c == n) return true; return false; };
    if (!has_col5("DELTAR_correct_NOSYS")) {
      mainNode = m_partonMatcher->systematicDefine(mainNode, "DELTAR_correct_NOSYS",
        [](const ROOT::VecOps::RVec<TLV>& jet_TLV, const ROOT::VecOps::RVec<TLV>& el_TLV, const ROOT::VecOps::RVec<TLV>& mu_TLV,
           const std::vector<int>& top1_truth, const std::vector<int>& top2_truth)->bool{
          auto asg = TtbarDileptonVars::getAssignment(TtbarDileptonVars::Method::DELTAR, {}, {}, jet_TLV, el_TLV, mu_TLV, ROOT::Math::PtEtaPhiMVector(0,0,0,0), ROOT::Math::PtEtaPhiMVector(0,0,0,0));
          if (!asg.valid) return false;
          int b1 = asg.bjet1_idx; int b2 = asg.bjet2_idx;
          if (top1_truth.size()==0 || top2_truth.size()==0) return false;
          int t1b = top1_truth.size()>0? top1_truth[0] : -1;
          int t2b = top2_truth.size()>0? top2_truth[0] : -1;
          if (t1b<0 || t2b<0) return false;
          return (b1==t1b && b2==t2b) || (b1==t2b && b2==t1b);
        },
        {"jet_TLV_NOSYS","sorted_el_TLV_NOSYS","sorted_mu_TLV_NOSYS","Top1_Truth_Reco_Indices_NOSYS","Top2_Truth_Reco_Indices_NOSYS"});
    }
  }

  // HYPER already has correctness helpers; alias to both-tops-correct
  mainNode = m_partonMatcher->systematicDefine(mainNode, "HYPER_correct_NOSYS",
    [](bool both_correct){ return both_correct; },
    {"HyPER_both_tops_correct_NOSYS"});

  return mainNode;
}

ROOT::RDF::RNode TtbarDileptonFrameReco::defineVariablesTruth(ROOT::RDF::RNode node,
                                                   const std::string& treeName,
                                                   const std::shared_ptr<Sample>& /*sample*/,
                                                   const UniqueSampleID& /*sampleID*/) {

  if (treeName != "parton") return node; // do not do anything for other trees ("truth" = parton truth)

  return node;
}

ROOT::RDF::RNode TtbarDileptonFrameReco::defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                           const std::string& /*treeName*/,
                                                           const std::shared_ptr<Sample>& /*sample*/,
                                                           const UniqueSampleID& /*sampleID*/) {
  return node;
}