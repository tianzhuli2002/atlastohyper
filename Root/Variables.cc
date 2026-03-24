#include<vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "Math/Vector4D.h"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "FastFrames/MainFrame.h"
#include "FastFrames/DefineHelpers.h"
#include "PartonMatcher/Variables.h"
#include "TRandom.h"

// Variables used in both channels
namespace Variables{
    /**
     * @brief Calculate the pt ordered indices of the objects that pass one selecttion.
     * @return std::vector<std::size_t> The indices of the objects that pass the selection.
     * @param ptVec The pT vector of the objects.
     * @param selection The first selection.
     */
    std::vector<std::size_t> indexPtOrderedGoodOneSelec(const std::vector<float>& ptVec,
                                        const std::vector<char>& selection){
        return DefineHelpers::sortedPassedIndices(ptVec,selection);
    }

    /**
     * @brief Calculate the pt ordered indices of the objects that pass two selections.
     * @return std::vector<std::size_t> The indices of the objects that pass the two selections.
     * @param ptVec The pT vector of the objects.
     * @param selection1 The first selection.
     * @param selection2 The second selection.
     */
    std::vector<std::size_t> indexPtOrderedGoodTwoSelec(const std::vector<float>& ptVec,
                                        const std::vector<char>& selection1,
                                        const std::vector<char>& selection2){
        return DefineHelpers::sortedPassedIndices(ptVec,selection1,selection2);
    }

    /**
     * @brief Calculate the number of jets that pass one selection and have a pT greater than 20 GeV.
     * @return int The number of jets that pass the selections.
     * @param ptVec The pT vector of the jets.
     * @param selected The first selection.
     */
    int nGoodJets(const std::vector<float>& ptVec,
                        const std::vector<char>& selected){
        return DefineHelpers::numberOfObjects(ptVec,20000,selected);
    }

    /**
     * @brief Create a TLV for a parton from the four momentum branches.
     * @return ROOT::RDF::RNode The node with the new TLV branch.
     * @param mainFrame The main frame.
     * @param mainNode The current main node.
     * @param partonBaseName The base name of the parton branches.
     * @param newBranchName The name of the new TLV branch.
     */
    ROOT::RDF::RNode createPartonTLV(MainFrame* mainFrame ,ROOT::RDF::RNode mainNode, const std::string& partonBaseName, const std::string& newBranchName){
        std::vector<std::string> kinematics{"_pt","_eta","_phi","_m"};
        
        auto buildTLV = [] (float pt, float eta, float phi, float m) {
            return ROOT::Math::PtEtaPhiMVector(pt,eta,phi,m);
        };
        
        std::string variableName = newBranchName + "_TLV_NOSYS";
        std::string ptName = partonBaseName + kinematics[0];
        std::string etaName = partonBaseName + kinematics[1];
        std::string phiName = partonBaseName + kinematics[2];
        std::string mName = partonBaseName + kinematics[3];

        LOG(INFO) << "Adding variable: " << variableName << "\n";
        mainNode = mainFrame->systematicDefine(mainNode,
                                         variableName,
                                         buildTLV,
                                         std::vector<std::string>{ptName,etaName,phiName,mName});

        return mainNode;
    }

    /**
     * @brief Calculate the invariant mass of two TLVs.
     * @param tlv1 The first TLV.
     * @param tlv2 The second TLV.
     * @return float The invariant mass of the two TLVs.
     */
    float invariantMassTwo(const ROOT::Math::PtEtaPhiMVector& tlv1, const ROOT::Math::PtEtaPhiMVector& tlv2){
        return (tlv1+tlv2).M()/1000.0f;
    }

    /**
     * @brief Calculate the invariant mass of three TLVs.
     * @param tlv1 The first TLV.
     * @param tlv2 The second TLV.
     * @param tlv3 The third TLV.
     * @return float The invariant mass of the three TLVs.
     */
    float invariantMassThree(const ROOT::Math::PtEtaPhiMVector& tlv1, const ROOT::Math::PtEtaPhiMVector& tlv2, const ROOT::Math::PtEtaPhiMVector& tlv3){
        return (tlv1+tlv2+tlv3).M()/1000.0f;
    }

    /**
     * @brief Calculate the invariant mass of two objects identified from indices, the objects are picked up from a vector of TLVs.
     * @param tlvVec The TLV vector.
     * @param index1 The index of the first object.
     * @param index2 The index of the second object.
     * @return float The invariant mass of the three TLVs.
     */
    float recoInvariantMassTwo(int index1, int index2, const std::vector<ROOT::Math::PtEtaPhiEVector>& tlvVector){
        ROOT::Math::PtEtaPhiEVector tlv1 = index1 >= 0 ? tlvVector.at(index1) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        ROOT::Math::PtEtaPhiEVector tlv2 = index2 >= 0 ? tlvVector.at(index2) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        
        return (tlv1+tlv2).M()/1000.0f;
    }

    /**
     * @brief Calculate the invariant mass of three objects identified from indices, the objects are picked up from a vector of TLVs.
     * @param tlvVec The TLV vector.
     * @param index1 The index of the first object.
     * @param index2 The index of the second object.
     * @param index3 The index of the third object.
     * @return float The invariant mass of the three TLVs.
     */
    float recoInvariantMassThree(int index1, int index2, int index3, const std::vector<ROOT::Math::PtEtaPhiEVector>& tlvVector){
        ROOT::Math::PtEtaPhiEVector tlv1 = index1 >= 0 ? tlvVector.at(index1) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        ROOT::Math::PtEtaPhiEVector tlv2 = index2 >= 0 ? tlvVector.at(index2) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        ROOT::Math::PtEtaPhiEVector tlv3 = index3 >= 0 ? tlvVector.at(index3) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);

        return (tlv1+tlv2+tlv3).M()/1000.0f;
    }

    /**
     * @brief Match a parton to a list of reco particles.
     * @param parton The parton TLV.
     * @param recoParticlesOrderedByPt The reco particles ordered by pT.
     * @return int The index of the matched reco particle.
     */
    int matchPartonRecoParticles(const ROOT::Math::PtEtaPhiMVector& parton, // Notice difference between E and M !!!
                                const std::vector<ROOT::Math::PtEtaPhiEVector>& recoParticlesOrderedByPt){
        float minDeltaR = 10000.0f;
        int matchedIndex = -1;
        for(std::size_t i = 0; i < recoParticlesOrderedByPt.size(); i++){
            float deltaR = ROOT::Math::VectorUtil::DeltaRapidityPhi(parton, recoParticlesOrderedByPt.at(i));
            if(deltaR < minDeltaR){
                minDeltaR = deltaR;
                matchedIndex = i;
            }
        }
        return minDeltaR <= 0.4 ? matchedIndex : -1;
    }

    /**
     * @brief Check if a W is matched correctly.
     * @param q1Index The index of the first quark.
     * @param q2Index The index of the second quark.
     */
    bool isWMatchedCorrectly(int q1Index, int q2Index){
        // TODO: An invariant mass check can be added here.
        if (q1Index == -1 || q2Index == -1){
            return false;
        }
        return q1Index != q2Index;
    }

    /**
     * @brief Check if two W overlap.
     * @param q1Index The index of the first quark from W1.
     * @param q2Index The index of the second quark from W1.
     * @param q1barIndex The index of the first quark from W2.
     * @param q2barIndex The index of the second quark from W2.
     */
    bool isNotWOverlap(int q1Index, int q2Index, int q1barIndex, int q2barIndex){
        bool q1NotOverlap = q1Index != q1barIndex && q1Index != q2barIndex;
        bool q2ONotverlap = q2Index != q1barIndex && q2Index != q2barIndex;
        return q1NotOverlap && q2ONotverlap;
    }

    /**
     * @brief Check if a Top is matched correctly.
     * @param q1Index The index of the first quark.
     * @param q2Index The index of the second quark.
     * @param bIndex The index of the b quark.
     */
    bool isTopMatchedCorrectly(int q1Index, int q2Index, int bIndex){
        // TODO: An invariant mass check can be added here.
        if (q1Index == -1 || q2Index == -1 || bIndex == -1){
            return false;
        }
        return q1Index != q2Index && q1Index != bIndex && q2Index != bIndex;
    }

    /**
     * @brief Check if two b quarks overlap.
     * @param bIndex The index of the b quark from top1.
     * @param bbarIndex The index of the b quark from top2.
     */
    bool isNotBOverlap(int bIndex, int bbarIndex){
        return bIndex != bbarIndex;
    }

    /**
     * @brief Create PCBT vector from the WPs information.
     * @param FCBEff_65 Is the jet there at the 65% WP (boolean).
     * @param FCBEff_70 Is the jet there at the 70% WP (boolean).
     * @param FCBEff_77 Is the jet there at the 77% WP (boolean).
     * @param FCBEff_85 Is the jet there at the 85% WP (boolean).
     * @param FCBEff_90 Is the jet there at the 90% WP (boolean).
     * @return std::vector<int> The output container.
     */
    std::vector<int> getPCBT(const std::vector<char>& FCBEff_65, const std::vector<char>& FCBEff_70, const std::vector<char>& FCBEff_77, const std::vector<char>& FCBEff_85, const std::vector<char>& FCBEff_90){
        std::vector<int> PCBT_output;

        size_t n = FCBEff_65.size();

        for (size_t i = 0; i < n; ++i) {
            if (FCBEff_65.at(i)) {
                PCBT_output.push_back(6);
            } else if (FCBEff_70.at(i)) {
                PCBT_output.push_back(5);
            } else if (FCBEff_77.at(i)) {
                PCBT_output.push_back(4);
            } else if (FCBEff_85.at(i)) {
                PCBT_output.push_back(3);
            } else if (FCBEff_90.at(i)) {
                PCBT_output.push_back(2);
            } else {
                PCBT_output.push_back(1);
            }
        }

        return PCBT_output;
    }

    /**
     * @brief Adjust the HyPER indices in function of the selections applied on the object containers
     * @return std::vector<long> The updated version of the HyPER indices.
     * @param edgeOrder 2 for edges, 3+ for hyperedges (number of final states associated with the (hyper)edge)
     * @param HyPER_indices Initial HyPER indices vector.
     * @param HyPER_IDs Corresponding HyPER IDs vector.
     * @param jet_idx Post selections jet indices.
     * @param el_idx Post selections electron indices.
     * @param mu_idx Post selections muon indices.
     */
    std::vector<long> getSelecHyPERIndices(const std::vector<long>& HyPER_indices,
                                        const std::vector<long>& HyPER_IDs,
                                        const std::vector<size_t>& jet_idx,
                                        const std::vector<size_t>& el_idx,
                                        const std::vector<size_t>& mu_idx,
                                        const std::vector<size_t>& tau_idx){
        
        std::vector<long> new_HyPER_indices = {};

        int edgeOrder = HyPER_indices.size();

        for(long i = 0; i < edgeOrder; ++i){
            
            long type;

            if(edgeOrder == 2){
                type = HyPER_IDs[i+1];
            }else{
                type = HyPER_IDs[i];
            }
            
            switch(type){
                case 1:
                    //Jet
                    for(size_t j = 0; j < jet_idx.size(); ++j){
                        if(static_cast<long>(jet_idx[j]) == HyPER_indices[i]){
                            new_HyPER_indices.push_back(j);
                        }
                    }

                    if(static_cast<long>(new_HyPER_indices.size()) < i){
                        new_HyPER_indices.push_back(-1);
                    }

                    break;
                case 2:
                    //Electron
                    for(size_t j = 0; j < el_idx.size(); ++j){
                        if(static_cast<long>(el_idx[j]) == HyPER_indices[i]){
                            new_HyPER_indices.push_back(j);
                        }
                    }

                    if(static_cast<long>(new_HyPER_indices.size()) < i){
                        new_HyPER_indices.push_back(-1);
                    }

                    break;
                case 3:
                    //Muon
                    for(size_t j = 0; j < mu_idx.size(); ++j){
                        if(static_cast<long>(mu_idx[j]) == HyPER_indices[i]){
                            new_HyPER_indices.push_back(j);
                        }
                    }

                    if(static_cast<long>(new_HyPER_indices.size()) < i){
                        new_HyPER_indices.push_back(-1);
                    }

                    break;
                case 4:
                    //MET or Neutrino
                    new_HyPER_indices.push_back(HyPER_indices[i]); //Not expecting any "fake" neutrino or neutrino removal via selection
                    break;
                case 5:
                    //Tau
                    for(size_t j = 0; j < tau_idx.size(); ++j){
                        if(static_cast<long>(tau_idx[j]) == HyPER_indices[i]){
                            new_HyPER_indices.push_back(j);
                        }
                    }

                    if(static_cast<long>(new_HyPER_indices.size()) < i){
                        new_HyPER_indices.push_back(-1);
                    }

                    break;
                default:
                    new_HyPER_indices.push_back(-1);
                    break;
            }
        }

        return new_HyPER_indices;
    }

    // Overload to accept integer-typed HyPER branches (convert to long and forward)
    std::vector<long> getSelecHyPERIndices(const std::vector<int>& HyPER_indices,
                                        const std::vector<int>& HyPER_IDs,
                                        const std::vector<size_t>& jet_idx,
                                        const std::vector<size_t>& el_idx,
                                        const std::vector<size_t>& mu_idx,
                                        const std::vector<size_t>& tau_idx){
        std::vector<long> HyPER_indices_long(HyPER_indices.begin(), HyPER_indices.end());
        std::vector<long> HyPER_IDs_long(HyPER_IDs.begin(), HyPER_IDs.end());
        return getSelecHyPERIndices(HyPER_indices_long, HyPER_IDs_long, jet_idx, el_idx, mu_idx, tau_idx);
    }

    /**
     * @brief Turn Char_t variable to bool.
     * @param myVariable variable that should be used as boolean
     */
    bool charToBool(Char_t myVariable){
        return static_cast<int>(myVariable) == 1;
    }
}
// ----------------- mttbar reconstruction helpers -----------------
namespace TtbarDileptonVars {

    float compute_mttbar(const ROOT::Math::PtEtaPhiMVector& lep1,
                         const ROOT::Math::PtEtaPhiMVector& lep2,
                         const ROOT::Math::PtEtaPhiMVector& bjet1,
                         const ROOT::Math::PtEtaPhiMVector& bjet2,
                         const ROOT::Math::PtEtaPhiMVector& nu1,
                         const ROOT::Math::PtEtaPhiMVector& nu2){
        auto top1 = lep1 + bjet1 + nu1;
        auto top2 = lep2 + bjet2 + nu2;
        auto ttbar = top1 + top2;
        return ttbar.M();
    }

    ROOT::Math::PtEtaPhiMVector buildNuFromBranches(float pt, float eta, float phi, float e){
        return ROOT::Math::PtEtaPhiMVector(pt, eta, phi, e);
    }

    std::pair<ROOT::Math::PtEtaPhiMVector, ROOT::Math::PtEtaPhiMVector>
    buildNuFromNuFlows(const ROOT::VecOps::RVec<float>& nuflows_nu_out){
        // Expect ordering: [nu1_px, nu1_py, nu1_pz, nu1_E, nu2_px, ...]
        if (nuflows_nu_out.size() < 8) return {ROOT::Math::PtEtaPhiMVector(0,0,0,0), ROOT::Math::PtEtaPhiMVector(0,0,0,0)};
        auto buildFromPxPyPzE = [](const float px, const float py, const float pz, const float E){
            float pt = std::hypot(px, py);
            float p = std::hypot(pt, pz);
            float eta = 0.0f;
            if (p > 0) eta = 0.5f * std::log((p + pz) / std::max(1e-9f, p - pz));
            float phi = std::atan2(py, px);
            return ROOT::Math::PtEtaPhiMVector(pt, eta, phi, std::sqrt(std::max(0.0f, E*E - p*p)));
        };
        auto nu1 = buildFromPxPyPzE(nuflows_nu_out[0], nuflows_nu_out[1], nuflows_nu_out[2], nuflows_nu_out[3]);
        auto nu2 = buildFromPxPyPzE(nuflows_nu_out[4], nuflows_nu_out[5], nuflows_nu_out[6], nuflows_nu_out[7]);
        return {nu1, nu2};
    }

    Assignment assignHyPER(const std::vector<long>& top_indices,
                           const std::vector<long>& top_IDs,
                           const ROOT::VecOps::RVec<TLV>& jet_TLV,
                           const ROOT::VecOps::RVec<TLV>& el_TLV,
                           const ROOT::VecOps::RVec<TLV>& mu_TLV){
        Assignment out;
        // top_indices expected as flat vector for one top: [i0,i1,i2] etc.; top_IDs holds corresponding type IDs
        if (top_indices.size() < 3 || top_IDs.size() < 3) return out;
        // Simple parser: find for top1 first 3 entries, then top2 next 3
        auto parseTriple = [&](size_t base)->Assignment{
            Assignment a; a.valid = false;
            int jet_idx = -1, lep_idx = -1;
            for(size_t k=0;k<3 && base+k<top_IDs.size() && base+k<top_indices.size();++k){
                long id = top_IDs[base+k];
                long idx = top_indices[base+k];
                if (id == 1) jet_idx = static_cast<int>(idx);
                if (id == 2 || id == 3) lep_idx = static_cast<int>(idx);
            }
            a.bjet1_idx = jet_idx; a.lep1_idx = lep_idx; a.valid = (jet_idx>=0 && lep_idx>=0);
            return a;
        };
        // parse first and second triples
        Assignment a1 = parseTriple(0);
        Assignment a2 = parseTriple(3);
        if (!a1.valid || !a2.valid) return out;
        out.bjet1_idx = a1.bjet1_idx;
        out.lep1_idx  = a1.lep1_idx;
        out.bjet2_idx = a2.bjet1_idx;
        out.lep2_idx  = a2.lep1_idx;
        out.valid = (out.bjet1_idx>=0 && out.bjet2_idx>=0 && out.lep1_idx>=0 && out.lep2_idx>=0 && out.bjet1_idx != out.bjet2_idx);
        return out;
    }

    Assignment assignPT(const ROOT::VecOps::RVec<TLV>& jet_TLV,
                        const ROOT::VecOps::RVec<TLV>& el_TLV,
                        const ROOT::VecOps::RVec<TLV>& mu_TLV){
        Assignment out;
        if (jet_TLV.size() < 2) return out;
        out.bjet1_idx = 0; out.bjet2_idx = 1;
        out.lep1_idx = 0; out.lep2_idx = 1;
        out.valid = true;
        return out;
    }

    Assignment assignDeltaR(const ROOT::VecOps::RVec<TLV>& jet_TLV,
                            const ROOT::VecOps::RVec<TLV>& el_TLV,
                            const ROOT::VecOps::RVec<TLV>& mu_TLV){
        Assignment out;
        // build lepton vector (e/mu mix): use first two available
        std::vector<TLV> leps;
        if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
        else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
        else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
        if (leps.size() != 2 || jet_TLV.empty()) return out;

        auto findClosest = [&](const TLV& lep)->int{
            float best = std::numeric_limits<float>::infinity(); int besti = -1;
            for (size_t j=0;j<jet_TLV.size();++j){
                float dR = ROOT::Math::VectorUtil::DeltaR(lep, jet_TLV[j]);
                if (dR < best){ best = dR; besti = static_cast<int>(j); }
            }
            return besti;
        };

        int j0 = findClosest(leps[0]);
        int j1 = findClosest(leps[1]);
        if (j0 == j1){ // resolve conflict: pick second-best for second lepton
            int second_best = -1; float second_best_val = std::numeric_limits<float>::infinity();
            for (size_t j=0;j<jet_TLV.size();++j){ if (static_cast<int>(j)==j0) continue; float dR = ROOT::Math::VectorUtil::DeltaR(leps[1], jet_TLV[j]); if (dR < second_best_val){ second_best_val = dR; second_best = static_cast<int>(j); } }
            j1 = second_best;
        }
        if (j0 < 0 || j1 < 0) return out;
        out.bjet1_idx = j0; out.bjet2_idx = j1; out.lep1_idx = 0; out.lep2_idx = 1; out.valid = true;
        return out;
    }

    Assignment assignChi2(const ROOT::VecOps::RVec<TLV>& jet_TLV,
                          const ROOT::VecOps::RVec<TLV>& el_TLV,
                          const ROOT::VecOps::RVec<TLV>& mu_TLV,
                          const ROOT::Math::PtEtaPhiMVector& nu1,
                          const ROOT::Math::PtEtaPhiMVector& nu2){
        Assignment out;
        if (jet_TLV.size() < 2) return out;
        // lepton ordering
        std::vector<TLV> leps;
        if (!el_TLV.empty() && !mu_TLV.empty()){ leps.push_back(el_TLV[0]); leps.push_back(mu_TLV[0]); }
        else if (el_TLV.size()>=2){ leps.push_back(el_TLV[0]); leps.push_back(el_TLV[1]); }
        else if (mu_TLV.size()>=2){ leps.push_back(mu_TLV[0]); leps.push_back(mu_TLV[1]); }
        if (leps.size() != 2) return out;

        const float m_top = 172500.0f; // MeV
        double bestScore = std::numeric_limits<double>::infinity();
        int best_i=-1, best_j=-1;
        for (size_t i=0;i<jet_TLV.size();++i){
            for (size_t j=0;j<jet_TLV.size();++j){
                if (i==j) continue;
                // Ensure all inputs to compute_mttbar are of type PtEtaPhiMVector by converting
                auto lep0_m = ROOT::Math::PtEtaPhiMVector(leps[0].Pt(), leps[0].Eta(), leps[0].Phi(), leps[0].M());
                auto lep1_m = ROOT::Math::PtEtaPhiMVector(leps[1].Pt(), leps[1].Eta(), leps[1].Phi(), leps[1].M());
                auto jet_i_m = ROOT::Math::PtEtaPhiMVector(jet_TLV[i].Pt(), jet_TLV[i].Eta(), jet_TLV[i].Phi(), jet_TLV[i].M());
                auto jet_j_m = ROOT::Math::PtEtaPhiMVector(jet_TLV[j].Pt(), jet_TLV[j].Eta(), jet_TLV[j].Phi(), jet_TLV[j].M());
                float m1 = compute_mttbar(lep0_m, lep1_m, jet_i_m, jet_j_m, nu1, nu2);
                // swap for second pairing considered below by using same neutrinos
                float m2 = compute_mttbar(lep1_m, lep0_m, jet_j_m, jet_i_m, nu2, nu1);
                double chi2 = (m1 - m_top)*(m1 - m_top) + (m2 - m_top)*(m2 - m_top);
                if (chi2 < bestScore){ bestScore = chi2; best_i = static_cast<int>(i); best_j = static_cast<int>(j); }
            }
        }
        if (best_i<0 || best_j<0) return out;
        out.bjet1_idx = best_i; out.bjet2_idx = best_j; out.lep1_idx = 0; out.lep2_idx = 1; out.valid = true;
        return out;
    }

    Assignment getAssignment(Method method,
                             const std::vector<long>& top_indices,
                             const std::vector<long>& top_IDs,
                             const ROOT::VecOps::RVec<TLV>& jet_TLV,
                             const ROOT::VecOps::RVec<TLV>& el_TLV,
                             const ROOT::VecOps::RVec<TLV>& mu_TLV,
                             const ROOT::Math::PtEtaPhiMVector& nu1,
                             const ROOT::Math::PtEtaPhiMVector& nu2){
        switch(method){
            case Method::HYPER: return assignHyPER(top_indices, top_IDs, jet_TLV, el_TLV, mu_TLV);
            case Method::CHI2:  return assignChi2(jet_TLV, el_TLV, mu_TLV, nu1, nu2);
            case Method::PT:    return assignPT(jet_TLV, el_TLV, mu_TLV);
            case Method::DELTAR: return assignDeltaR(jet_TLV, el_TLV, mu_TLV);
            default: return Assignment();
        }
    }

}


// Ensure a definition of ZllCalibrationVars::regionName exists in the global scope
// (some translation-units expect this exact symbol). Provide a minimal deterministic
// implementation so the library can be loaded at runtime.
std::string ZllCalibrationVars::regionName(
    const ROOT::VecOps::RVec<TLV>& sorted_elTLV,
    const ROOT::VecOps::RVec<TLV>& sorted_muTLV,
    const ROOT::VecOps::RVec<float>& sorted_elCharge,
    const ROOT::VecOps::RVec<float>& sorted_muCharge) {
    size_t nEl = sorted_elTLV.size();
    size_t nMu = sorted_muTLV.size();
    if (nEl == 2 && nMu == 0) return "ee";
    if (nMu == 2 && nEl == 0) return "mumu";
    if (nEl == 1 && nMu == 1) return "emu";
    return "other";
}

namespace TtbarSingleLeptonVars{
    /**
     * @brief Calculate the number of good leptons.
     * @param ptVec The pT vector of the leptons.
     * @param selected1 The first selection.
     * @param selected2 The second selection.
     * @return int The number of good leptons.
     */
    int nGoodLeptons(const std::vector<float>& ptVec,
                        const std::vector<char>& selected1,
                        const std::vector<char>& selected2){
        return DefineHelpers::numberOfObjects(ptVec,25000,selected1,selected2);
    }

    /**
     * @brief Calculate the number of good jets.
     * @param ptVec The pT vector of the jets.
     * @param selected1 The first selection.
     * @param selected2 The second selection.
     * @return int The number of good jets.
     */
    int nGoodJets(const std::vector<float>& ptVec,
                        const std::vector<char>& selected1,
                        const std::vector<char>& selected2){
        return DefineHelpers::numberOfObjects(ptVec,20000,selected1,selected2);
    }

    /**
     * @brief Calculate the indices of pT ordered objects passing two selections.
     * @param ptVec The pT vector of the objects.
     * @param selection1 The first selection.
     * @param selection2 The second selection.
     * @return std::vector<std::size_t> The pT ordered indices of the objects passing the selections.
     */
    std::vector<std::size_t> indexPtOrderedGoodTwoSelec(const std::vector<float>& ptVec,
                                        const std::vector<char>& selection1,
                                        const std::vector<char>& selection2){
        return DefineHelpers::sortedPassedIndices(ptVec,selection1,selection2);
    }

    /**
     * @brief Return the unique element in the KLF index vector.
     * @param klfIndex The KLF index vector.
     * @return int The unique element in the KLF index vector.
     */
    int fromKLFToNtupleIndex(const std::vector<unsigned int>& klfIndex){
        if (klfIndex.size() == 0) return -1;
        return klfIndex[0];
    }

    // Selection cuts for the muon channel
    namespace MuonSelection {
        // Events firing any of the triggers
        bool triggers(const bool& T1_2022,const bool& T2_2022,
                    const bool& T1_2018,const bool& T2_2018){

            bool pass_2018 = T1_2018 || T2_2018;
            bool pass_2022 = T1_2022 || T2_2022;
            return pass_2018 || pass_2022;
        }

        // Basic selection with number of good objects, triggers and correct ntuple branch.
        bool basic(const char& nTupleFlag, int numberOfMuons, int numberOfJets,bool triggered, bool triggerMatched){
            return nTupleFlag && numberOfMuons == 1 && triggered && triggerMatched && numberOfJets >= 4;//numberOfJets == 4;
        }
        // Events firing any of the available muon triggers
    bool triggersnew(const bool& mu24, const bool& mu50){
        // Just OR the available triggers
        return mu24 || mu50;
    }

    // Basic selection with number of good objects, triggers and correct ntuple branch.
    bool basicnew(const char& nTupleFlag, int numberOfMuons, int numberOfJets,
               bool triggered, bool triggerMatched){
        return nTupleFlag && numberOfMuons == 1 && triggered && triggerMatched && numberOfJets >= 4;
    }

    }

    // Selection cuts for the electron channel
    namespace ElectronSelection {
        // Events firing any of the triggers
        bool triggers(const bool& T1_2022,const bool& T2_2022, const bool& T3_2022,
                    const bool& T1_2018,const bool& T2_2018, const bool& T3_2018){

            bool pass_2018 = T1_2018 || T2_2018 || T3_2018; 
            bool pass_2022 = T1_2022 || T2_2022 || T3_2022;
            return pass_2018 || pass_2022;
        }

        // Basic selection with number of good objects, triggers and correct ntuple branch.
        bool basic(const char& nTupleFlag, int numberOfElectrons, int numberOfJets,bool triggered, bool triggerMatched){
            return nTupleFlag && numberOfElectrons == 1 && triggered && triggerMatched && numberOfJets >= 4;//numberOfJets == 4;
        }
        // Events firing any of the available electron triggers
    bool triggersnew(const bool& e26, const bool& e60, const bool& e140){
        // OR all available triggers
        return e26 || e60 || e140;
    }

    // Basic selection with number of good objects, triggers and correct ntuple branch.
    bool basicnew(const char& nTupleFlag, int numberOfElectrons, int numberOfJets,
               bool triggered, bool triggerMatched){
        return nTupleFlag && numberOfElectrons == 1 && triggered && triggerMatched && numberOfJets >= 4;
    }

   }

   //Selection cuts on jets
   namespace JetsSelection {
        bool fullyMatched(Int_t truth_down_index, Int_t truth_up_index, Int_t truth_bhad_index, Int_t truth_blep_index){
            bool have_down = (truth_down_index != -1);
            bool have_up   = (truth_up_index != -1);
            bool have_bhad = (truth_bhad_index != -1);
            bool have_blep = (truth_blep_index != -1);

            return (have_down && have_up && have_bhad && have_blep);
        } 
   }

    /**
     * @brief Get index of true down jet
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @return Int_t index of true down jet
     */
    Int_t getPartonTruthDown(const std::vector<int>& event_jet_truth_idx){

        Int_t myIndex = -1;

        if((event_jet_truth_idx[2] != -1) && (event_jet_truth_idx[5] == -1)){

            myIndex = event_jet_truth_idx[2];
        
        }else if((event_jet_truth_idx[2] == -1) && (event_jet_truth_idx[5] != -1)){
            
            myIndex = event_jet_truth_idx[5];
        }

        return myIndex;
    }

    /**
     * @brief Get index of true up jet
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @return Int_t index of true up jet
     */
    Int_t getPartonTruthUp(const std::vector<int>& event_jet_truth_idx){

        Int_t myIndex = -1;

        if((event_jet_truth_idx[1] != -1) && (event_jet_truth_idx[4] == -1)){

            myIndex = event_jet_truth_idx[1];
        
        }else if((event_jet_truth_idx[1] == -1) && (event_jet_truth_idx[4] != -1)){
            
            myIndex = event_jet_truth_idx[4];
        }
        
        return myIndex;
    }

    /**
     * @brief Get index of true hadronic b jet
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @return Int_t index of true hadronic b jet
     */
    Int_t getPartonTruthBHad(const std::vector<int>& event_jet_truth_idx){

        Int_t myIndex = -1;

        if((event_jet_truth_idx[0] != -1) && (event_jet_truth_idx[1] != -1) && (event_jet_truth_idx[2] != -1)){

            myIndex = event_jet_truth_idx[0];
        
        }else if((event_jet_truth_idx[3] != -1) && (event_jet_truth_idx[4] != -1) && (event_jet_truth_idx[5] != -1)){
            
            myIndex = event_jet_truth_idx[3];
        }
        
        return myIndex;
    }

    /**
     * @brief Get index of true leptonic b jet
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @return Int_t index of true leptonic b jet
     */
    Int_t getPartonTruthBLep(const std::vector<int>& event_jet_truth_idx){

        Int_t myIndex = -1;

        if((event_jet_truth_idx[0] != -1) && (event_jet_truth_idx[1] == -1) && (event_jet_truth_idx[2] == -1)){

            myIndex = event_jet_truth_idx[0];
        
        }else if((event_jet_truth_idx[3] == -1) && (event_jet_truth_idx[4] == -1) && (event_jet_truth_idx[5] == -1)){
            
            myIndex = event_jet_truth_idx[3];
        }
        
        return myIndex;
    }
}


namespace TtbarDileptonVars{
    float sumOfCharges( const ROOT::VecOps::RVec<float>& chargeV) { return ROOT::VecOps::Sum(chargeV); }
  
    ROOT::VecOps::RVec<float> sortedPassedChargeVector25(
        const ROOT::VecOps::RVec<float>& chargeV,
        const ROOT::VecOps::RVec<float>& ptV,
        const ROOT::VecOps::RVec<char>& selection){

        // Create a vector of decisions to test pT > 25 GeV.
        ROOT::VecOps::RVec<char> gt25 = ptV > 25000;

        // Get the pt-sorted indices of the leptons that pass the selection
        auto passedIndices = DefineHelpers::sortedPassedIndices(ptV, gt25 ,selection);
        // Create a vector to hold the sorted charges
        ROOT::VecOps::RVec<float> sortedCharges;

        // Loop over the passed indices and fill the sorted charges vector
        for (const auto& index : passedIndices) {
            sortedCharges.push_back(chargeV[index]);
        }

        return sortedCharges;
    }

    std::string regionName(
        std::size_t nElectrons,
        std::size_t nMuons,
        const ROOT::VecOps::RVec<TLV>& sorted_elTLV,
        const ROOT::VecOps::RVec<TLV>& sorted_muTLV) {

            // Get the leading lepton pT. Note: params are (nElectrons, nMuons, elTLV, muTLV)
            float leadingMuonPt = nMuons > 0 ? sorted_muTLV.at(0).Pt() : 0.0f;
            float leadingElectronPt = nElectrons > 0 ? sorted_elTLV.at(0).Pt() : 0.0f;
            float leadingLeptonPt = std::max(leadingMuonPt, leadingElectronPt);
            if (leadingLeptonPt < 28000 ) return "other";
            if (nMuons + nElectrons == 2) return "2l";
            return "other";
    }

    

    int recoIndexTruthBJet(const ROOT::VecOps::RVec<TLV>& recoBJets,
        const ROOT::Math::PtEtaPhiMVector& truthBJet) {

        // Check the inputs are not empty
        if (recoBJets.size() == 0 || truthBJet.Pt() == 0) {
        return -1;
        }

        // Initialize the index to -1 (no match)
        int index = -1;
        double maxDeltaR = 0.4;

        // Loop over the reco b-jets
        for (std::size_t i = 0; i < recoBJets.size(); ++i) {
        // Calculate the deltaR between the reco b-jet and the truth b-jet
        double deltaR = ROOT::Math::VectorUtil::DeltaR(recoBJets[i], truthBJet);
        // Check if the deltaR is less than the maximum allowed
        if (deltaR < maxDeltaR) {
        index = i;
        maxDeltaR = deltaR;
        }
        }

        return index;
    }

    float ttbarLeptonPairInvariantMass(
        const std::string& regionName,
        const ROOT::VecOps::RVec<TLV>& muonTLV,
        const ROOT::VecOps::RVec<TLV>& electronTLV){

        // If region is other, return 0.0f
        if (regionName == "other") return 0.0f;
        
        if (muonTLV.size() + electronTLV.size() == 2) {
            // For the 2l region, just return the invariant mass of the two leptons.
            if (muonTLV.size() == 2 && electronTLV.size() == 0) {
                return ROOT::Math::VectorUtil::InvariantMass(muonTLV[0], muonTLV[1]);
            }
            if (electronTLV.size() == 2 && muonTLV.size() == 0) {
                return ROOT::Math::VectorUtil::InvariantMass(electronTLV[0], electronTLV[1]);
            }
            if (muonTLV.size() == 1 && electronTLV.size() == 1) {
                return ROOT::Math::VectorUtil::InvariantMass(muonTLV[0], electronTLV[0]);
            }
        }
        return 0.0f;
        }



    float minimumLeptonJetPairInvariantMass(const ROOT::VecOps::RVec<TLV>& sorted_el_TLV,
                                            const ROOT::VecOps::RVec<TLV>& sorted_mu_TLV,
                                            const ROOT::VecOps::RVec<TLV>& sorted_bjet_TLV,
                                            const ROOT::VecOps::RVec<TLV>& jet_TLV,
                                            const ROOT::VecOps::RVec<char>& jet_GN2v01_FixedCutBEff_77_select,
                                            unsigned long n_electrons,
                                            unsigned long n_muons) {
                                           //ROOT::VecOps::RVec<float> out(4, 999999.0f);
                                            float out = -1.0f;//999999.0f;
                                           // Build the lepton pair in a consistent order.
                                           ROOT::VecOps::RVec<TLV> leps;
                                           if (n_electrons == 2 && n_muons == 0) {
                                             // two electrons
                                             if (sorted_el_TLV.size() < 2) return out;
                                             leps = {sorted_el_TLV[0], sorted_el_TLV[1]};
                                           } else if (n_muons == 2 && n_electrons == 0) {
                                             // two muons
                                             if (sorted_mu_TLV.size() < 2) return out;
                                             leps = {sorted_mu_TLV[0], sorted_mu_TLV[1]};
                                           } else if (n_muons == 1 && n_electrons == 1) {
                                             // e+mu: keep electron first then muon to be consistent with previous usages
                                             if (sorted_el_TLV.size() < 1 || sorted_mu_TLV.size() < 1) return out;
                                             leps = {sorted_el_TLV[0], sorted_mu_TLV[0]};
                                           } else {
                                             // Not a dilepton event we care about
                                             return out;
                                           }

                                           // Choose the two b-jet candidates per spec:
                                           // - if >=2 b-tagged jets, take the two highest-pT b-tagged jets (sorted_bjet_TLV assumed pT-sorted)
                                           // - if only one b-tagged jet, take that b-tagged jet and the highest-pT non-b-tagged jet
                                           ROOT::VecOps::RVec<TLV> bjets;
                                           if (sorted_bjet_TLV.size() >= 2) {
                                             bjets = {sorted_bjet_TLV[0], sorted_bjet_TLV[1]};
                                           } else if (sorted_bjet_TLV.size() == 1) {
                                             // find highest-pt non-b-tagged jet from jet_TLV using the provided btag selector
                                             TLV highestNonB; bool found=false;
                                             for (size_t i=0;i<jet_TLV.size();++i) {
                                               if (i >= jet_GN2v01_FixedCutBEff_77_select.size()) continue;
                                               if (jet_GN2v01_FixedCutBEff_77_select.at(i)) continue; // skip b-tagged
                                               if (!found || jet_TLV[i].Pt() > highestNonB.Pt()) { highestNonB = jet_TLV[i]; found = true; }
                                             }
                                             if (!found) return out; // no non-b jet found
                                             bjets = {sorted_bjet_TLV[0], highestNonB};
                                           } else {
                                             // fewer than 1 b-jet -> cannot form mlb pairs
                                             return out;
                                           }

                                           // Compute the four mlb values for the two possible assignments.
                                           float mlb_A_0 = ROOT::Math::VectorUtil::InvariantMass(leps[0], bjets[0]);
                                           float mlb_A_1 = ROOT::Math::VectorUtil::InvariantMass(leps[1], bjets[1]);
                                           float mlb_B_0 = ROOT::Math::VectorUtil::InvariantMass(leps[0], bjets[1]);
                                           float mlb_B_1 = ROOT::Math::VectorUtil::InvariantMass(leps[1], bjets[0]);

                                           float Amax = std::max(mlb_A_0, mlb_A_1);
                                           float Bmax = std::max(mlb_B_0, mlb_B_1);
                                           return std::min(Amax, Bmax);
                                           //return out;
                                         }
ROOT::VecOps::RVec<int> selectTwoJetIndices_BtagPriority(const ROOT::VecOps::RVec<TLV>& jet_TLV,
                                 const ROOT::VecOps::RVec<char>& jet_GN2v01_FixedCutBEff_77_select)
{
  ROOT::VecOps::RVec<int> bjet_indices;
  bjet_indices.reserve(jet_TLV.size());

  // Collect b-tagged jet indices
  for (size_t i = 0; i < jet_TLV.size(); ++i) {
    if (i >= jet_GN2v01_FixedCutBEff_77_select.size()) continue;
    if (jet_GN2v01_FixedCutBEff_77_select[i])
      bjet_indices.push_back(i);
  }

  // Sort b-tagged jets by pT (descending)
  std::sort(bjet_indices.begin(), bjet_indices.end(),
            [&](int a, int b) {
              return jet_TLV[a].Pt() > jet_TLV[b].Pt();
            });

  ROOT::VecOps::RVec<int> out;

  // Case 1: ≥2 b-tagged jets
  if (bjet_indices.size() >= 2) {
    out.push_back(bjet_indices[0]);
    out.push_back(bjet_indices[1]);
    return out;
  }

  // Case 2: exactly 1 b-tagged jet
  if (bjet_indices.size() == 1) {
    int bidx = bjet_indices[0];

    int best_nonb = -1;
    float best_pt = -1.0f;

    for (size_t i = 0; i < jet_TLV.size(); ++i) {
      if (i >= jet_GN2v01_FixedCutBEff_77_select.size()) continue;
      if (jet_GN2v01_FixedCutBEff_77_select[i]) continue; // skip b-tagged

      float pt = jet_TLV[i].Pt();
      if (pt > best_pt) {
        best_pt = pt;
        best_nonb = static_cast<int>(i);
      }
    }

    if (best_nonb < 0) return ROOT::VecOps::RVec<int>{};

    out.push_back(bidx);
    out.push_back(best_nonb);
    return out;
  }

  // Case 3: no b-tagged jets
  return ROOT::VecOps::RVec<int>{};
}

int SelectedJetsMatchTruthBPair(
    const ROOT::VecOps::RVec<int>& selected_jet_indices,
    const std::vector<int>& event_jet_truth_idx,
    const std::vector<int>& event_jet_truth_candidates)
{
  if (selected_jet_indices.size() != 2)
    return 0;

  // Inline truth-match validity check
  auto jet_ok = [&](size_t idx) -> bool {
    return event_jet_truth_idx.size() > idx &&
           event_jet_truth_candidates.size() > idx &&
           event_jet_truth_idx[idx] != -1 &&
           event_jet_truth_candidates[idx] == 1;
  };

  // Truth b and b̄ availability
  bool b_ok    = jet_ok(0);
  bool bbar_ok = jet_ok(3);

  // No valid truth b jets
  if (!b_ok && !bbar_ok)
    return 0;

  int sel0 = selected_jet_indices[0];
  int sel1 = selected_jet_indices[1];

  int t0 = event_jet_truth_idx[sel0];
  int t1 = event_jet_truth_idx[sel1];

  bool sel0_is_b    = b_ok    && (t0 == event_jet_truth_idx[0]);
  bool sel0_is_bbar = bbar_ok && (t0 == event_jet_truth_idx[3]);
  bool sel1_is_b    = b_ok    && (t1 == event_jet_truth_idx[0]);
  bool sel1_is_bbar = bbar_ok && (t1 == event_jet_truth_idx[3]);

  int n_matches =
      sel0_is_b + sel0_is_bbar +
      sel1_is_b + sel1_is_bbar;

  // Correct b–b̄ pair (any order)
  if (n_matches == 2)
    return 1;

  // No match at all
  if (n_matches == 0)
    return 0;

  // Partial or ambiguous match
  return -1;
}

                                     
//     float compute_mllbb(const ROOT::VecOps::RVec<TLV>& sorted_el_TLV,
//                     const ROOT::VecOps::RVec<TLV>& sorted_mu_TLV,
//                     const ROOT::VecOps::RVec<TLV>& sorted_bjet_TLV,
//                     const ROOT::VecOps::RVec<TLV>& jet_TLV,
//                     const ROOT::VecOps::RVec<char>& jet_GN2v01_FixedCutBEff_77_select,
//                     unsigned long n_electrons,
//                     unsigned long n_muons)
// {
//   // Build dilepton system
//   ROOT::VecOps::RVec<TLV> leps;

//   if (n_electrons == 2 && n_muons == 0) {
//     if (sorted_el_TLV.size() < 2) return -1.0f;
//     leps = {sorted_el_TLV[0], sorted_el_TLV[1]};
//   } else if (n_muons == 2 && n_electrons == 0) {
//     if (sorted_mu_TLV.size() < 2) return -1.0f;
//     leps = {sorted_mu_TLV[0], sorted_mu_TLV[1]};
//   } else if (n_muons == 1 && n_electrons == 1) {
//     if (sorted_el_TLV.size() < 1 || sorted_mu_TLV.size() < 1) return -1.0f;
//     leps = {sorted_el_TLV[0], sorted_mu_TLV[0]};
//   } else {
//     return -1.0f;
//   }

// //   // Select the two jets (b-tag priority)
// //   ROOT::VecOps::RVec<TLV> jets =
// //     selectTwoJets_BtagPriority(sorted_bjet_TLV,
// //                                jet_TLV,
// //                                jet_GN2v01_FixedCutBEff_85_select);

// //   if (jets.size() != 2) return -1.0f;
//   auto jet_idx = selectTwoJetIndices_BtagPriority(jet_TLV, jet_GN2v01_FixedCutBEff_77_select);
// if (jet_idx.size() == 2) {
//   TLV llbb = leps[0] + leps[1] + jet_TLV[jet_idx[0]] + jet_TLV[jet_idx[1]];
//   // Build llbb four-vector
//   //TLV llbb = leps[0] + leps[1] + jets[0] + jets[1];
// } else {
//   return -1.0f;
// }

//   return llbb.M();
// }                                     
        

    bool pass_min_lb_invar_mass(const ROOT::VecOps::RVec<TLV>& sorted_el_TLV,
                                            const ROOT::VecOps::RVec<TLV>& sorted_mu_TLV,
                                            const ROOT::VecOps::RVec<TLV>& sorted_bjet_TLV,
                                            const ROOT::VecOps::RVec<TLV>& jet_TLV,
                                            const ROOT::VecOps::RVec<char>& jet_GN2v01_FixedCutBEff_85_select,
                                            unsigned long n_electrons,
                                            unsigned long n_muons) {
                                           //ROOT::VecOps::RVec<float> out(4, 999999.0f);
                                            float out = -1.0f;//999999.0f;
                                           // Build the lepton pair in a consistent order.
                                           ROOT::VecOps::RVec<TLV> leps;
                                           if (n_electrons == 2 && n_muons == 0) {
                                             // two electrons
                                             if (sorted_el_TLV.size() < 2) return out;
                                             leps = {sorted_el_TLV[0], sorted_el_TLV[1]};
                                           } else if (n_muons == 2 && n_electrons == 0) {
                                             // two muons
                                             if (sorted_mu_TLV.size() < 2) return out;
                                             leps = {sorted_mu_TLV[0], sorted_mu_TLV[1]};
                                           } else if (n_muons == 1 && n_electrons == 1) {
                                             // e+mu: keep electron first then muon to be consistent with previous usages
                                             if (sorted_el_TLV.size() < 1 || sorted_mu_TLV.size() < 1) return out;
                                             leps = {sorted_el_TLV[0], sorted_mu_TLV[0]};
                                           } else {
                                             // Not a dilepton event we care about
                                             return out;
                                           }

                                           // Choose the two b-jet candidates per spec:
                                           // - if >=2 b-tagged jets, take the two highest-pT b-tagged jets (sorted_bjet_TLV assumed pT-sorted)
                                           // - if only one b-tagged jet, take that b-tagged jet and the highest-pT non-b-tagged jet
                                           ROOT::VecOps::RVec<TLV> bjets;
                                           if (sorted_bjet_TLV.size() >= 2) {
                                             bjets = {sorted_bjet_TLV[0], sorted_bjet_TLV[1]};
                                           } else if (sorted_bjet_TLV.size() == 1) {
                                             // find highest-pt non-b-tagged jet from jet_TLV using the provided btag selector
                                             TLV highestNonB; bool found=false;
                                             for (size_t i=0;i<jet_TLV.size();++i) {
                                               if (i >= jet_GN2v01_FixedCutBEff_85_select.size()) continue;
                                               if (jet_GN2v01_FixedCutBEff_85_select.at(i)) continue; // skip b-tagged
                                               if (!found || jet_TLV[i].Pt() > highestNonB.Pt()) { highestNonB = jet_TLV[i]; found = true; }
                                             }
                                             if (!found) return out; // no non-b jet found
                                             bjets = {sorted_bjet_TLV[0], highestNonB};
                                           } else {
                                             // fewer than 1 b-jet -> cannot form mlb pairs
                                             return out;
                                           }

                                           // Compute the four mlb values for the two possible assignments.
                                           float mlb_A_0 = ROOT::Math::VectorUtil::InvariantMass(leps[0], bjets[0]);
                                           float mlb_A_1 = ROOT::Math::VectorUtil::InvariantMass(leps[1], bjets[1]);
                                           float mlb_B_0 = ROOT::Math::VectorUtil::InvariantMass(leps[0], bjets[1]);
                                           float mlb_B_1 = ROOT::Math::VectorUtil::InvariantMass(leps[1], bjets[0]);

                                           float Amax = std::max(mlb_A_0, mlb_A_1);
                                           float Bmax = std::max(mlb_B_0, mlb_B_1);
                                           if (std::min(Amax, Bmax) < 150000) return true;
                                             return false;
                                           //return out;
                                         }
                                         
        ROOT::VecOps::RVec<float> FullyMatched_MLB(
        const ROOT::VecOps::RVec<TLV>& jet_TLVs,
        const ROOT::VecOps::RVec<TLV>& el_TLVs,
        const ROOT::VecOps::RVec<TLV>& mu_TLVs,
        const ROOT::VecOps::RVec<float>& el_charge,
        const ROOT::VecOps::RVec<float>& mu_charge,
        const ROOT::VecOps::RVec<int>& event_jet_truth_idx,
        const ROOT::VecOps::RVec<int>& event_jet_truth_candidates
        ){
            ROOT::VecOps::RVec<float> out(9, -1.0f); // [0]=true pairing sum, [1]=alt pairing sum

            int idx_b    = event_jet_truth_idx[0];
            int idx_bbar = event_jet_truth_idx[3];

            if (idx_b < 0 || idx_bbar < 0) return out;
            if (event_jet_truth_candidates[0] != 1 || event_jet_truth_candidates[3] != 1) return out;
            if (idx_b >= (int)jet_TLVs.size() || idx_bbar >= (int)jet_TLVs.size()) return out;
            if (idx_b == idx_bbar) return out; // must be distinct reco jets

            TLV jet_b    = jet_TLVs[idx_b];
            TLV jet_bbar = jet_TLVs[idx_bbar];

            // require opposite-sign leptons
            float total_charge = 0;
            if (!el_charge.empty()) total_charge += el_charge[0];
            if (!mu_charge.empty()) total_charge += mu_charge[0];
            if (total_charge != 0.0f) return out;

            // assign lepton1 and lepton2 consistently
            TLV leptonplus, leptonminus;
            if (!el_charge.empty() && el_charge[0] > 0) {
                leptonplus = el_TLVs[0];
                //std::cout << "el_TLVs Pt: " << el_TLVs << "\n";
                leptonminus = mu_TLVs[0];
            } else if (!mu_charge.empty() && mu_charge[0] > 0) {
                leptonplus = mu_TLVs[0];
                leptonminus = el_TLVs[0];
            } else {
                // fallback: just assign in order if positive/negative not known
                if (!el_TLVs.empty() && !mu_TLVs.empty()) { return out; }
                //     lepton1 = el_TLVs[0];
                //     lepton2 = mu_TLVs[0];
                // } else {
                //     return out; // sanity check
                // }
            }

            // Compute sums of mlb for the true pairing (lepton1->b, lepton2->bbar)
            // and the alternative pairing (swap b/bbar). Units follow VectorUtil::InvariantMass.
            
            float mlb_true_b = ROOT::Math::VectorUtil::InvariantMass(leptonplus, jet_b);
            float mlb_true_bbar = ROOT::Math::VectorUtil::InvariantMass(leptonminus, jet_bbar);

            float pTvis0 = (leptonplus + jet_b).Pt();
            float pTvis1 = (leptonminus + jet_bbar).Pt();
            float pTdiff_true = pTvis0 - pTvis1; // MeV
            float sum_deltaR_true = ROOT::Math::VectorUtil::DeltaR(leptonplus, jet_b)
                              + ROOT::Math::VectorUtil::DeltaR(leptonminus, jet_bbar);
            
            float mlb_alt_b = ROOT::Math::VectorUtil::InvariantMass(leptonplus, jet_bbar);
            float mlb_alt_bbar = ROOT::Math::VectorUtil::InvariantMass(leptonminus, jet_b);
            
            out[0] = mlb_true_b;
            out[1] = mlb_true_bbar;
            out[2] = mlb_alt_b;
            out[3] = mlb_alt_bbar;
            out[4] = pTdiff_true;
            out[5] = sum_deltaR_true;
            out[6] = jet_b.Pt();
            out[7] = jet_bbar.Pt();
            out[8] = deltaPhi(jet_b, jet_bbar);
            return out;
    }          

    float invariantMassllbb(
        const ROOT::VecOps::RVec<TLV>& sorted_el_TLV,
        const ROOT::VecOps::RVec<TLV>& sorted_mu_TLV,
        const ROOT::VecOps::RVec<TLV>& sorted_bjet_TLV,
        const ROOT::VecOps::RVec<TLV>& jet_TLV,
        const ROOT::VecOps::RVec<char>& jet_GN2v01_FixedCutBEff_77_select,
        unsigned long n_electrons,
        unsigned long n_muons)
    {
        float mllbb = -1.0f;

        // --- build the lepton pair ---
        ROOT::VecOps::RVec<TLV> leps;
        if (n_electrons == 2 && n_muons == 0) {
            if (sorted_el_TLV.size() < 2) return mllbb;
            leps = {sorted_el_TLV[0], sorted_el_TLV[1]};
        } else if (n_muons == 2 && n_electrons == 0) {
            if (sorted_mu_TLV.size() < 2) return mllbb;
            leps = {sorted_mu_TLV[0], sorted_mu_TLV[1]};
        } else if (n_muons == 1 && n_electrons == 1) {
            if (sorted_el_TLV.size() < 1 || sorted_mu_TLV.size() < 1) return mllbb;
            leps = {sorted_el_TLV[0], sorted_mu_TLV[0]};
        } else {
            return mllbb;
        }

        // --- select two b-jet candidates ---
        ROOT::VecOps::RVec<TLV> bjets;
        if (sorted_bjet_TLV.size() >= 2) {
            bjets = {sorted_bjet_TLV[0], sorted_bjet_TLV[1]};
        } else if (sorted_bjet_TLV.size() == 1) {
            TLV highestNonB; bool found = false;
            for (size_t i = 0; i < jet_TLV.size(); ++i) {
                if (i >= jet_GN2v01_FixedCutBEff_77_select.size()) continue;
                // skip b-tagged jets
                if (jet_GN2v01_FixedCutBEff_77_select.at(i)) continue;
                // also skip the matched b-jet index if it was found (defensive)
                //if (matchedBJetIndex != -1 && static_cast<int>(i) == matchedBJetIndex) continue;
                if (!found || jet_TLV[i].Pt() > highestNonB.Pt()) {
                    highestNonB = jet_TLV[i];
                    found = true;
                }
            }
            if (!found) return mllbb;
            bjets = {sorted_bjet_TLV[0], highestNonB};
        } else {
            return mllbb;
        }

        // --- compute invariant mass of system ---
        TLV total = leps[0] + leps[1] + bjets[0] + bjets[1];
        mllbb = total.M();
        return mllbb;
    }
float Signal(
    const ROOT::VecOps::RVec<TLV>& elTLV,
    const ROOT::VecOps::RVec<TLV>& muTLV,
    const ROOT::VecOps::RVec<float>& elChargeV,
    const ROOT::VecOps::RVec<float>& muChargeV,
    unsigned long nJets,
    unsigned long nBJets,
    float dilepton_mass_NOSYS,
    float lepton_bjet_invariant_mass,
    float met_met_NOSYS
) {
    const unsigned long nElectrons = elTLV.size();
    const unsigned long nMuons     = muTLV.size();

    // --------------------------------------------------
    // Exactly two leptons
    // --------------------------------------------------
    if (nElectrons + nMuons != 2) return 0;

    // --------------------------------------------------
    // Opposite-sign leptons
    // --------------------------------------------------
    float sumCharge = 0.0f;
    for (float q : elChargeV) sumCharge += q;
    for (float q : muChargeV) sumCharge += q;
    if (sumCharge != 0.0f) return 0;

    // --------------------------------------------------
    // Lepton pT cuts (TLVs are pT-sorted)
    // --------------------------------------------------
    float pt0 = 0.0f;
    float pt1 = 0.0f;

    if (nElectrons == 2) {
        pt0 = elTLV[0].Pt();
        pt1 = elTLV[1].Pt();
    }
    else if (nMuons == 2) {
        pt0 = muTLV[0].Pt();
        pt1 = muTLV[1].Pt();
    }
    else {
        // eμ channel
        pt0 = std::max(elTLV[0].Pt(), muTLV[0].Pt());
        pt1 = std::min(elTLV[0].Pt(), muTLV[0].Pt());
    }

    // both leptons pT > 25 GeV
    if (pt0 < 25000 || pt1 < 25000) return 0;

    // leading lepton pT > 28 GeV
    if (pt0 < 28000) return 0;

    // --------------------------------------------------
    // Jet and b-tag requirements
    // --------------------------------------------------
    if (nJets < 2) return 0;
    if (nBJets < 1) return 0;

    // --------------------------------------------------
    // Dilepton mass cuts
    // --------------------------------------------------
    if (dilepton_mass_NOSYS < 15000) return 0;

    bool isSF = (nElectrons == 2 || nMuons == 2);
    if (isSF) {
        // MET > 45 GeV
        if (met_met_NOSYS < 45000) return 0;

        // Z veto
        if (dilepton_mass_NOSYS > 81000 &&
            dilepton_mass_NOSYS < 101000)
            return 0;
    }

    // --------------------------------------------------
    // m_lb cut
    // --------------------------------------------------
    if (lepton_bjet_invariant_mass > 150000) return 0;

    return 1;
}



    // matched under the strict criteria, 0 otherwise.
    int FullyMatched(const std::vector<int>& event_jet_truth_idx,
                 const std::vector<int>& event_jet_truth_candidates){
    // Helper for safe access
    auto jet_ok = [&](size_t idx)->bool{
        return event_jet_truth_idx.size() > idx &&
               event_jet_truth_candidates.size() > idx &&
               event_jet_truth_idx[idx] != -1 &&
               event_jet_truth_candidates[idx] == 1;
    };

    // Require both b and b̄ jets to be well-matched
    bool b_matched = jet_ok(0);
    bool bbar_matched = jet_ok(3);
    //if (!(b_matched && bbar_matched)) return 0;
     if (!b_matched && !bbar_matched)
        return 0;

    // Case 2: Exactly one matched → -1
    if (b_matched ^ bbar_matched)
        return -1;

    // Case 3: Both matched → check if they match different truth jets
    if (event_jet_truth_idx[0] == event_jet_truth_idx[3])
        return 0;  // both match the *same* truth jet → invalid

    return 1;  // fully matched
}

// Compute transverse mass for vis-system + MET:
// mT = sqrt( (Et_vis + MET)^2 - |pT_vis + MET_vec|^2 )
float transverseMassSystem_safe(const TLV& vis, float met_pt, float met_phi) {
    float m_vis = vis.M();
    float pT_vis = vis.Pt();
    float Et_vis = std::sqrt(std::max(0.0f, m_vis*m_vis + pT_vis*pT_vis));

    float met_px = met_pt * std::cos(met_phi);
    float met_py = met_pt * std::sin(met_phi);

    float pvis_px = vis.Px();
    float pvis_py = vis.Py();

    float px_sum = pvis_px + met_px;
    float py_sum = pvis_py + met_py;
    float sumEt = Et_vis + met_pt;

    float arg = sumEt*sumEt - (px_sum*px_sum + py_sum*py_sum);
    if (arg <= 0.0f) return 0.0f; // guard numerical noise -> return 0 rather than NaN
    return std::sqrt(arg);
}


Int_t getTruthIndexBHad(const std::vector<int>& event_jet_truth_idx,
                        const std::vector<int>& event_jet_truth_candidates){

        Int_t myIndex = -1;

        if (event_jet_truth_idx.size() > 3 && event_jet_truth_candidates.size() > 3) {
            if (event_jet_truth_idx[0] != -1  && event_jet_truth_candidates[0] == 1) {
                myIndex = event_jet_truth_idx[0];
            }
        }

        return myIndex;
    }

Int_t getTruthIndexBbarHad(const std::vector<int>& event_jet_truth_idx,
const std::vector<int>& event_jet_truth_candidates){

        Int_t myIndex = -1;

        if (event_jet_truth_idx.size() > 3 && event_jet_truth_candidates.size() > 3) {
            if (event_jet_truth_idx[3] != -1  && event_jet_truth_candidates[3] == 1) {
                myIndex = event_jet_truth_idx[3];
            }
        }
        return myIndex;
    }

// Compute mT_ttbar using truth-matched b and bbar (no ΔR). Returns -1 on failure.
float mT_ttbar_truth_from_truthIndices(
    const ROOT::VecOps::RVec<TLV>& jet_TLVs,
    const ROOT::VecOps::RVec<TLV>& el_TLVs,
    const ROOT::VecOps::RVec<TLV>& mu_TLVs,
    const ROOT::VecOps::RVec<int>& event_jet_truth_idx,
    const ROOT::VecOps::RVec<int>& event_jet_truth_candidates,
    const ROOT::VecOps::RVec<float>& el_charge,
    const ROOT::VecOps::RVec<float>& mu_charge,
    float met_pt,
    float met_phi)
{
    // Validate truth arrays: expecting b at slot 0 and bbar at slot 3, candidate==1
    if (event_jet_truth_idx.size() <= 3 || event_jet_truth_candidates.size() <= 3) return -1.0f;
    if (event_jet_truth_idx[0] == -1 || event_jet_truth_idx[3] == -1) return -1.0f;
    if (event_jet_truth_candidates[0] != 1 || event_jet_truth_candidates[3] != 1) return -1.0f;

    int idx_b    = event_jet_truth_idx[0];
    int idx_bbar = event_jet_truth_idx[3];

    if (idx_b < 0 || idx_bbar < 0) return -1.0f;
    if (idx_b >= static_cast<int>(jet_TLVs.size()) || idx_bbar >= static_cast<int>(jet_TLVs.size())) return -1.0f;
    if (idx_b == idx_bbar) return -1.0f;

    // require we have opposite-sign leptons (same check used in FullyMatched_MLB)
    float total_charge = 0.0f;
    if (!el_charge.empty()) total_charge += el_charge[0];
    if (!mu_charge.empty()) total_charge += mu_charge[0];
    if (total_charge != 0.0f) return -1.0f;

    // Build lepton TLVs in the same consistent ordering used elsewhere (electron first if positive)
    TLV lep0, lep1;
    if (!el_TLVs.empty() && !mu_TLVs.empty()) {
        if (el_charge[0] > 0) { lep0 = el_TLVs[0]; lep1 = mu_TLVs[0]; }
        else { lep0 = mu_TLVs[0]; lep1 = el_TLVs[0]; }
    } else {
        return -1.0f;
    }

    TLV jet_b    = jet_TLVs[idx_b];
    TLV jet_bbar = jet_TLVs[idx_bbar];

    TLV vis = lep0 + lep1 + jet_b + jet_bbar;

    // Important: keep units consistent. If TLVs are MeV, pass met_pt in MeV too.
    float mT = transverseMassSystem_safe(vis, met_pt, met_phi);
    return mT;
}

#include <limits>

// Helper: compute deltaPhi in range (-pi, +pi]
inline float deltaPhi(const TLV& a,
                           const TLV& b) {
    return std::fabs(ROOT::Math::VectorUtil::DeltaPhi(a, b)); // [-pi, pi]
}

inline float deltaR(const TLV &a, const TLV &b) {
    return ROOT::Math::VectorUtil::DeltaR(a, b);
}

// Returns: [ jet_for_lep0, jet_for_lep1, dR_lep0, dR_lep1, sum_dR, correctness_flag ]
ROOT::VecOps::RVec<float> PairLeptonsToJets_MinSumDeltaR(
    const ROOT::VecOps::RVec<TLV>& el_TLV,
    const ROOT::VecOps::RVec<TLV>& mu_TLV,
    const ROOT::VecOps::RVec<TLV>& jet_TLV,
    const ROOT::VecOps::RVec<int>& event_jet_truth_idx,
    const ROOT::VecOps::RVec<int>& event_jet_truth_candidates,
    const ROOT::VecOps::RVec<float>& el_charge,
    const ROOT::VecOps::RVec<float>& mu_charge,
    unsigned long n_electrons,
    unsigned long n_muons
) {
    const float SENTINEL_NEG = -1.0f;
    ROOT::VecOps::RVec<float> out(6, SENTINEL_NEG);
    // out[0] = jet index for lep0
    // out[1] = jet index for lep1
    // out[2] = dR(lepton0, jet_for_lep0)
    // out[3] = dR(lepton1, jet_for_lep1)
    // out[4] = sum dR
    // out[5] = correctness: -1 no truth, 0 wrong, 1 correct

    // require 1 electron + 1 muon (as in your previous code)
    if (!(n_muons == 1 && n_electrons == 1)) return out;
    if (el_TLV.empty() || mu_TLV.empty()) return out;
    if (jet_TLV.size() < 2) return out;

    // order leptons so leps[0] = positive, leps[1] = negative (same convention as before)
    ROOT::VecOps::RVec<TLV> leps;
    if (el_charge[0] > 0) leps = {el_TLV[0], mu_TLV[0]};
    else                  leps = {mu_TLV[0], el_TLV[0]};

    // If exactly two jets, just compare the two possible assignments
    float best_sum = std::numeric_limits<float>::infinity();
    int best_i = -1, best_j = -1;
    float best_dR0 = -1.0f, best_dR1 = -1.0f;

    const size_t nj = jet_TLV.size();
    // loop over ordered pairs i != j (we preserve ordering because lepton0/lepton1 are distinct)
    for (size_t i = 0; i < nj; ++i) {
        for (size_t j = 0; j < nj; ++j) {
            if (i == j) continue;
            float dR0 = deltaR(leps[0], jet_TLV[i]);
            float dR1 = deltaR(leps[1], jet_TLV[j]);
            float sum = dR0 + dR1;
            if (sum < best_sum) {
                best_sum = sum;
                best_i = static_cast<int>(i);
                best_j = static_cast<int>(j);
                best_dR0 = dR0;
                best_dR1 = dR1;
            }
        }
    }

    if (best_i < 0 || best_j < 0) return out;

    out[0] = static_cast<float>(best_i);
    out[1] = static_cast<float>(best_j);
    out[2] = best_dR0;
    out[3] = best_dR1;
    out[4] = best_sum;

    // --- truth matching: accept swapped order as correct by default ---
    int truth_b = -1, truth_bbar = -1;
    bool haveTruth = false;
    if (event_jet_truth_idx.size() > 3 && event_jet_truth_candidates.size() > 3) {
        if (event_jet_truth_idx[0] != -1 && event_jet_truth_idx[3] != -1
            && event_jet_truth_candidates[0] == 1 && event_jet_truth_candidates[3] == 1) {
            truth_b = event_jet_truth_idx[0];
            truth_bbar = event_jet_truth_idx[3];
            if (truth_b >= 0 && truth_b < static_cast<int>(nj) &&
                truth_bbar >= 0 && truth_bbar < static_cast<int>(nj) &&
                truth_b != truth_bbar) {
                haveTruth = true;
            }
        }
    }

    if (!haveTruth) {
        out[5] = -1.0f; // no truth info
        return out;
    }

    int reco_j0 = static_cast<int>(out[0]);
    int reco_j1 = static_cast<int>(out[1]);
    bool correct = (reco_j0 == truth_b && reco_j1 == truth_bbar); //||(reco_j0 == truth_bbar && reco_j1 == truth_b));
    out[5] = correct ? 1.0f : 0.0f;

    return out;
}


 ROOT::VecOps::RVec<float> PairLeptonsToTopJets_MinSumSq(
    const ROOT::VecOps::RVec<TLV>& el_TLV,
    const ROOT::VecOps::RVec<TLV>& mu_TLV,
    const ROOT::VecOps::RVec<TLV>& jet_TLV,
    //const ROOT::VecOps::RVec<char>& jet_btag_select,
    const ROOT::VecOps::RVec<int>& event_jet_truth_idx,
    const ROOT::VecOps::RVec<int>& event_jet_truth_candidates,
    const ROOT::VecOps::RVec<float>& el_charge,
    const ROOT::VecOps::RVec<float>& mu_charge,
    unsigned long n_electrons,
    unsigned long n_muons)
{
    const float SENTINEL_NEG = -1.0f;
    ROOT::VecOps::RVec<float> out(7, SENTINEL_NEG);

    // Build lepton pair in consistent order (same convention used elsewhere)
    ROOT::VecOps::RVec<TLV> leps;
    if (!(n_muons == 1 && n_electrons == 1)) return out;
    if (el_TLV.empty() || mu_TLV.empty()) return out;
    if (el_charge[0] > 0) {
    // electron is positive → assign as lepton 0
    leps = {el_TLV[0], mu_TLV[0]};
} else {
    // muon must then be positive → assign as lepton 0
    leps = {mu_TLV[0], el_TLV[0]};
}    
    // Need at least two jets
    if (jet_TLV.size() < 2) return out;

    // Find indices of the two highest-pT jets
    size_t i0 = 0, i1 = 1;
    if (jet_TLV.size() > 2) {
        // simple top-2 selection
        if (jet_TLV[1].Pt() > jet_TLV[0].Pt()) { i0 = 1; i1 = 0; }
        for (size_t i = 2; i < jet_TLV.size(); ++i) {
            if (jet_TLV[i].Pt() > jet_TLV[i0].Pt()) { i1 = i0; i0 = i; }
            else if (jet_TLV[i].Pt() > jet_TLV[i1].Pt()) { i1 = i; }
        }
    } else {
        // exactly two jets
        if (jet_TLV[1].Pt() > jet_TLV[0].Pt()) { i0 = 1; i1 = 0; }
        else { i0 = 0; i1 = 1; }
    }

    // compute the two pairing objectives (sum of squared masses)
    float mA0 = ROOT::Math::VectorUtil::InvariantMass(leps[0], jet_TLV[i0]) ; // A pairing: (lep0->i0)
    float mA1 = ROOT::Math::VectorUtil::InvariantMass(leps[1], jet_TLV[i1]) ; // (lep1->i1)
    float objA = mA0*mA0 + mA1*mA1;

    float mB0 = ROOT::Math::VectorUtil::InvariantMass(leps[0], jet_TLV[i1]) ; // B pairing: swap
    float mB1 = ROOT::Math::VectorUtil::InvariantMass(leps[1], jet_TLV[i0]) ;
    float objB = mB0*mB0 + mB1*mB1;

    // choose minimal-sum-of-squares pairing
    size_t jet_for_lep0 = (objA <= objB) ? i0 : i1;
    size_t jet_for_lep1 = (objA <= objB) ? i1 : i0;
    float   m_jl0 = (objA <= objB) ? mA0 : mB0;
    float   m_jl1 = (objA <= objB) ? mA1 : mB1;

    // Fill output: cast jet indices to float for compatibility with RVec<float>
    out[0] = static_cast<float>(jet_for_lep0);
    out[1] = static_cast<float>(jet_for_lep1);
    out[2] = m_jl0;
    out[3] = m_jl1;
    out[6] = m_jl0 + m_jl1;
    // out[4] = regionCode;

    // If truth vectors are not present or malformed, set -1.
    // event_jet_truth_idx layout expected to have b at slot 0 and bbar at slot 3.
    int truth_b = -1;
    int truth_bbar = -1;
    bool haveTruth = false;
    if (event_jet_truth_idx.size() > 3 && event_jet_truth_candidates.size() > 3) {
        if (event_jet_truth_idx[0] != -1 && event_jet_truth_idx[3] != -1
            && event_jet_truth_candidates[0] == 1 && event_jet_truth_candidates[3] == 1) {
            truth_b = event_jet_truth_idx[0];
            truth_bbar = event_jet_truth_idx[3];
            // sanity: indices must be within reco jets
            if (truth_b >= 0 && truth_b < static_cast<int>(jet_TLV.size()) &&
                truth_bbar >= 0 && truth_bbar < static_cast<int>(jet_TLV.size()) &&
                truth_b != truth_bbar) {
                haveTruth = true;
            }
        }
    }
    if (!haveTruth) {
        //out[5] = -1.0f; // no truth info / not evaluable
        return out;
    }
    
    int reco_j0 = static_cast<int>(out[0]);
    int reco_j1 = static_cast<int>(out[1]);
    bool correct = ((reco_j0 == truth_b && reco_j1 == truth_bbar )); // (reco_j0 == truth_bbar && reco_j1 == truth_b));
    out[5] = correct ? 1.0f : 0.0f;
    return out;
} 

ROOT::VecOps::RVec<float> PairLeptonsToJets_Chi2TruthAll(
    const ROOT::VecOps::RVec<TLV>& el_TLV,
    const ROOT::VecOps::RVec<TLV>& mu_TLV,
    const ROOT::VecOps::RVec<TLV>& jet_TLV,
    float met_met,
    float met_phi,
    const ROOT::VecOps::RVec<int>& event_jet_truth_idx,
    const ROOT::VecOps::RVec<int>& event_jet_truth_candidates,
    const ROOT::VecOps::RVec<float>& el_charge,
    const ROOT::VecOps::RVec<float>& mu_charge,
    unsigned long n_electrons,
    unsigned long n_muons
) {
    
    float mean_mlb  = 92000.0f;  // Default
    float sigma_mlb = 21000.0f;

    if (jet_TLV.size() == 2) {
        mean_mlb  = 96374.0f;
        sigma_mlb = 27445.0f;
    }
    else if (jet_TLV.size() == 3) {
        mean_mlb  = 93600.0f;
        sigma_mlb = 21131.0f;
    }
    else if (jet_TLV.size() == 4) {
        mean_mlb  = 93000.0f;
        sigma_mlb = 26653.0f;
    }
    else if (jet_TLV.size() >= 5) {
        mean_mlb  = 90600.0f;
        sigma_mlb = 24341.0f;
    }

    const float mean_mllbb  = 357000.0f;  // MeV
    const float sigma_mllbb = 140000.0f;  // MeV

    const float mean_pTdiff = -363.5f;       // MeV
    const float sigma_pTdiff = 65000.0f;  // MeV

    const float mean_mT     = 389760.0f;  // MeV
    const float sigma_mT    = 88627.1f;  // MeV

    const float mean_deltaR = 3.567f;
    const float sigma_deltaR = 1.211f;

    const float w_mlb   = 1.0f;
    const float w_mllbb = 0.0f;
    const float w_pT    = 1.0f;
    const float w_mT    = 1.0f;
    const float w_deltaR = 1.0f;

    // Output layout:
    // 0: jet index for lepton 0 (positive lepton)
    // 1: jet index for lepton 1 (negative lepton)
    // 2: m(l0, j0)
    // 3: m(l1, j1)
    // 4: chi2 of chosen assignment
    // 5: correctness flag  (-1 no truth, 0 wrong, 1 correct)
    // 6: chi2 of truth (b, bbar) assignment (if available)
    // 7: delta chi2 (chosen - truth)  (if truth available)
    // 8: chi2 of swapped truth (bbar, b) assignment (if available)
    const float SENTINEL_NEG = -1.0f;
    ROOT::VecOps::RVec<float> out(9, SENTINEL_NEG);

    // Require exactly one e and one mu.
    if (!(n_electrons == 1 && n_muons == 1)) return out;
    if (el_TLV.empty() || mu_TLV.empty()) return out;
    if (jet_TLV.size() < 2) return out;

    // Charge vectors must be present; require opposite sign.
    if (el_charge.size() < 1 || mu_charge.size() < 1) return out;
    if (el_charge[0] * mu_charge[0] >= 0.0f) {
        // Same-sign or zero charge; leave sentinels.
        return out;
    }

    // Order leptons so leps[0] is positive (defines our convention).
    ROOT::VecOps::RVec<TLV> leps(2);
    if (el_charge[0] > 0.0f) {
        leps[0] = el_TLV[0];
        leps[1] = mu_TLV[0];
    } else {
        leps[0] = mu_TLV[0];
        leps[1] = el_TLV[0];
    }

    // Precompute MET components (assumed MeV).
    const float met_px = met_met * std::cos(met_phi);
    const float met_py = met_met * std::sin(met_phi);

    struct Chi2Terms {
        float chi2;
        float mlb0;
        float mlb1;
        float mllbb;
        float pTdiff;
        float mT_ttbar;
        float sum_deltaR;
    };

    // Helper to compute chi2 (and keep terms) for assignment (ii, jj).
    auto evaluate_pair = [&](size_t ii, size_t jj) -> Chi2Terms {
        Chi2Terms info{};
        info.mlb0 = ROOT::Math::VectorUtil::InvariantMass(leps[0], jet_TLV[ii]);
        info.mlb1 = ROOT::Math::VectorUtil::InvariantMass(leps[1], jet_TLV[jj]);

        const TLV vis0 = leps[0] + jet_TLV[ii];
        const TLV vis1 = leps[1] + jet_TLV[jj];
        info.pTdiff = vis0.Pt() - vis1.Pt();
        info.sum_deltaR = deltaR(leps[0], jet_TLV[ii]) + deltaR(leps[1], jet_TLV[jj]);

        const TLV vis = leps[0] + leps[1] + jet_TLV[ii] + jet_TLV[jj];
        info.mllbb = vis.M();
        const float m_vis = vis.M();
        const float pT_vis = vis.Pt();
        const float Et_vis = std::sqrt(std::max(0.0f, m_vis * m_vis + pT_vis * pT_vis));
        const float px_sum = vis.Px() + met_px;
        const float py_sum = vis.Py() + met_py;
        const float arg = (Et_vis + met_met) * (Et_vis + met_met) - (px_sum * px_sum + py_sum * py_sum);
        info.mT_ttbar = (arg > 0.0f) ? std::sqrt(arg) : 0.0f;

        const float term_mlb0 = (info.mlb0 - mean_mlb) / sigma_mlb;
        const float term_mlb1 = (info.mlb1 - mean_mlb) / sigma_mlb;
        const float term_mllbb = (info.mllbb - mean_mllbb) / sigma_mllbb;
        const float term_pT = (info.pTdiff - mean_pTdiff) / sigma_pTdiff;
        const float term_mT = (info.mT_ttbar - mean_mT) / sigma_mT;
        const float term_deltaR = (info.sum_deltaR - mean_deltaR) / sigma_deltaR;

        info.chi2 = w_mlb   * (term_mlb0 * term_mlb0 + term_mlb1 * term_mlb1)
                  + w_mllbb * (term_mllbb * term_mllbb)
                  + w_pT    * (term_pT * term_pT)
                  + w_mT    * (term_mT * term_mT)
                  + w_deltaR * (term_deltaR * term_deltaR);
        return info;
    };

    // Scan all ordered jet pairs.
    float best_chi2 = std::numeric_limits<float>::infinity();
    int best_i = -1;
    int best_j = -1;
    Chi2Terms best_terms{};

    const size_t nj = jet_TLV.size();
    for (size_t i = 0; i < nj; ++i) {
        for (size_t j = 0; j < nj; ++j) {
            if (i == j) continue;
            Chi2Terms info = evaluate_pair(i, j);
            if (info.chi2 < best_chi2) {
                best_chi2 = info.chi2;
                best_i = static_cast<int>(i);
                best_j = static_cast<int>(j);
                best_terms = info;
            }
        }
    }

    if (best_i < 0 || best_j < 0) {
        // No valid pair found.
        return out;
    }

    out[0] = static_cast<float>(best_i);
    out[1] = static_cast<float>(best_j);
    out[2] = best_terms.mlb0;
    out[3] = best_terms.mlb1;
    out[4] = best_terms.chi2;

    // --- Truth handling (assumed: event_jet_truth_idx already in the same ordering as jet_TLV) ---
    int truth_b = -1;
    int truth_bbar = -1;
    bool haveTruth = false;

    if (event_jet_truth_idx.size() >= 4 && event_jet_truth_candidates.size() >= 4) {
        const bool validB    = (event_jet_truth_idx[0] != -1 && event_jet_truth_candidates[0] == 1);
        const bool validBbar = (event_jet_truth_idx[3] != -1 && event_jet_truth_candidates[3] == 1);
        if (validB && validBbar) {
            truth_b    = event_jet_truth_idx[0];
            truth_bbar = event_jet_truth_idx[3];
            // NOTE: this assumes truth indices already correspond to positions in jet_TLV.
            if (truth_b >= 0 && truth_bbar >= 0 &&
                truth_b < static_cast<int>(nj) &&
                truth_bbar < static_cast<int>(nj) &&
                truth_b != truth_bbar) {
                haveTruth = true;
            }
        }
    }

    if (!haveTruth) {
        // out[5], out[6], out[7], out[8] remain sentinel -1.
        return out;
    }

    // Compute chi2 for truth orientations.
    const Chi2Terms truth_terms = evaluate_pair(static_cast<size_t>(truth_b),
                                                static_cast<size_t>(truth_bbar));
    const Chi2Terms truth_terms_swapped = evaluate_pair(static_cast<size_t>(truth_bbar),
                                                        static_cast<size_t>(truth_b));

    out[6] = truth_terms.chi2;
    out[7] = best_terms.chi2 - truth_terms.chi2;
    out[8] = truth_terms_swapped.chi2;

    const bool correct =
        (best_i == truth_b && best_j == truth_bbar);
        //(best_i == truth_bbar && best_j == truth_b);
    out[5] = correct ? 1.0f : 0.0f;

    return out;
}
/**
     * @brief Create a TLV for a parton from the four momentum branches.
     * @return ROOT::RDF::RNode The node with the new TLV branch.
     * @param mainFrame The main frame.
     * @param mainNode The current main node.
     * @param partonBaseName The base name of the parton branches.
     * @param newBranchName The name of the new TLV branch.
     */
    ROOT::RDF::RNode createNeutrinoTLV(MainFrame* mainFrame ,ROOT::RDF::RNode mainNode, const std::string& partonBaseName, const std::string& newBranchName){
        std::vector<std::string> kinematics{"_pt","_eta","_phi"};
        
        auto buildTLV = [] (float pt, float eta, float phi) {
            return ROOT::Math::PtEtaPhiMVector(pt,eta,phi,0);
        };
        
        std::string variableName = newBranchName + "_TLV";
        std::string ptName = partonBaseName + kinematics[0];
        std::string etaName = partonBaseName + kinematics[1];
        std::string phiName = partonBaseName + kinematics[2];
        //std::string mName = partonBaseName + kinematics[3];

        LOG(INFO) << "Adding variable: " << variableName << "\n";
        mainNode = mainNode.Define(variableName,
                                         buildTLV,
                                         std::vector<std::string>{ptName,etaName,phiName});
        return mainNode;
    }

//     struct TruthNuInfo {
//     ROOT::Math::PtEtaPhiMVector nu_from_t;
//     ROOT::Math::PtEtaPhiMVector nu_from_tbar;
//     int pdg_from_t;
//     int pdg_from_tbar;TtbarDileptonVars::
// };TtbarDileptonVars::
TruthNuInfo ExtractTruthNeutrinos(
    float t_w1_pt,    float t_w1_eta,    float t_w1_phi,    float t_w1_m,    int t_w1_pdg,
    float t_w2_pt,    float t_w2_eta,    float t_w2_phi,    float t_w2_m,    int t_w2_pdg,
    float tb_w1_pt,   float tb_w1_eta,   float tb_w1_phi,   float tb_w1_m,   int tb_w1_pdg,
    float tb_w2_pt,   float tb_w2_eta,   float tb_w2_phi,   float tb_w2_m,   int tb_w2_pdg
) {
    TruthNuInfo out;

    auto makeTLV = [](float pt, float eta, float phi, float m) {
        return ROOT::Math::PtEtaPhiMVector(pt, eta, phi, m);
    };

    auto isNu = [](int pdg) {
        return std::abs(pdg) == 12 || std::abs(pdg) == 14 || std::abs(pdg) == 16;
    };

    // ---- From top ----
    if (isNu(t_w1_pdg)) {
        out.nu_from_t  = makeTLV(t_w1_pt, t_w1_eta, t_w1_phi, t_w1_m);
        out.pdg_from_t = t_w1_pdg;
    } else if (isNu(t_w2_pdg)) {
        out.nu_from_t  = makeTLV(t_w2_pt, t_w2_eta, t_w2_phi, t_w2_m);
        out.pdg_from_t = t_w2_pdg;
    } else {
        out.nu_from_t  = makeTLV(0,0,0,0);
        out.pdg_from_t = 0;
    }

    // ---- From antitop ----
    if (isNu(tb_w1_pdg)) {
        out.nu_from_tbar  = makeTLV(tb_w1_pt, tb_w1_eta, tb_w1_phi, tb_w1_m);
        out.pdg_from_tbar = tb_w1_pdg;
    } else if (isNu(tb_w2_pdg)) {
        out.nu_from_tbar  = makeTLV(tb_w2_pt, tb_w2_eta, tb_w2_phi, tb_w2_m);
        out.pdg_from_tbar = tb_w2_pdg;
    } else {
        out.nu_from_tbar  = makeTLV(0,0,0,0);
        out.pdg_from_tbar = 0;
    }

    return out;
}


NuFlavourInfo ExtractNeutrinoFlavours(
    float t_w1_pt,    float t_w1_eta,    float t_w1_phi,    float t_w1_m,    int t_w1_pdg,
    float t_w2_pt,    float t_w2_eta,    float t_w2_phi,    float t_w2_m,    int t_w2_pdg,
    float tb_w1_pt,   float tb_w1_eta,   float tb_w1_phi,   float tb_w1_m,   int tb_w1_pdg,
    float tb_w2_pt,   float tb_w2_eta,   float tb_w2_phi,   float tb_w2_m,   int tb_w2_pdg
) {
    NuFlavourInfo out;

    auto makeTLV = [](float pt, float eta, float phi, float m) {
        return ROOT::Math::PtEtaPhiMVector(pt, eta, phi, m);
    };

    // Helper to process any daughter
    auto process = [&](float pt, float eta, float phi, float m, int pdg) {
        int apdg = std::abs(pdg);
        if (apdg == 12 ) { //&& out.pdg_e == 0
            out.nu_e  = makeTLV(pt, eta, phi, m);
            out.pdg_e = pdg;
        } else if (apdg == 14 ) { //&& out.pdg_mu == 0
            out.nu_mu  = makeTLV(pt, eta, phi, m);
            out.pdg_mu = pdg;
        }
    };

    // Process all W daughters (top and antitop)
    process(t_w1_pt,  t_w1_eta,  t_w1_phi,  t_w1_m,  t_w1_pdg);
    process(t_w2_pt,  t_w2_eta,  t_w2_phi,  t_w2_m,  t_w2_pdg);
    process(tb_w1_pt, tb_w1_eta, tb_w1_phi, tb_w1_m, tb_w1_pdg);
    process(tb_w2_pt, tb_w2_eta, tb_w2_phi, tb_w2_m, tb_w2_pdg);

    // Fill zeros if absent
    if (out.pdg_e == 0) {
        out.nu_e = makeTLV(0,0,0,0);
    }
    if (out.pdg_mu == 0) {
        out.nu_mu = makeTLV(0,0,0,0);
    }

    return out;
}





}

namespace HyPERPerformance{
  
    int countCorrectlyReconstructedW(const std::vector<long> HyPER_W1_Indices,
        const std::vector<long> HyPER_W2_Indices,
        const std::vector<int> jetTruthIdx,
        const std::vector<int> nMatchedObjects){

            // In TCT v2.17.0 the truth Ws are [1,2] and [4,5]. The assumption is that truth jets are uniquely matched.
            // Let's mark events which are not matcheable with -2.
            // And the events where HyPER did not converge properly with -1.

            // First check all sizes are appropriate.
            if (jetTruthIdx.size() != 6 || nMatchedObjects.size() != 6) return -2;

            // Check both truth Ws are uniquely matched.
            bool truthW1Unique = nMatchedObjects.at(1) == 1 && nMatchedObjects.at(2) == 1;
            bool truthW2Unique = nMatchedObjects.at(4) == 1 && nMatchedObjects.at(5) == 1;
            if (!truthW1Unique || !truthW2Unique) return -2;

            // Define the labels
            std::vector<int> truthW1{jetTruthIdx.at(1),jetTruthIdx.at(2)};
            std::vector<int> truthW2{jetTruthIdx.at(4),jetTruthIdx.at(5)};

            // Check we have a truth matcheable event
            bool truthW1Matcheable = truthW1.at(0) != -1 && truthW1.at(1) != -1;
            bool truthW2Matcheable = truthW2.at(0) != -1 && truthW2.at(1) != -1;
            if (!truthW1Matcheable || !truthW2Matcheable) return -2;
            
            // Now check the HyPER predictions.
            if (HyPER_W1_Indices.size() != 2 || HyPER_W2_Indices.size() != 2) return -1;

            // Check HyPER at least has one valid prediction.
            bool validW1HyPER = allValidIndices<long>(HyPER_W1_Indices);
            bool validW2HyPER = allValidIndices<long>(HyPER_W2_Indices);
            if (!validW1HyPER && !validW2HyPER) return -1;

            // Check HyPER has unique predictions.
            bool WPredictionsOverlap = twoVectorsOverlap<long,long>(HyPER_W1_Indices,HyPER_W2_Indices);
            if (WPredictionsOverlap) return -1;

            int nCorrectW = 0;
            if (twoVectorsEqual<long,int>(HyPER_W1_Indices,truthW1)) nCorrectW++;
            if (twoVectorsEqual<long,int>(HyPER_W1_Indices,truthW2)) nCorrectW++;
            if (twoVectorsEqual<long,int>(HyPER_W2_Indices,truthW1)) nCorrectW++;
            if (twoVectorsEqual<long,int>(HyPER_W2_Indices,truthW2)) nCorrectW++;

            if (nCorrectW > 2) {
                throw std::invalid_argument("The number of correctly reconstructed Ws is greater than 2.");
            }

            return nCorrectW;
        }

    int countCorrectlyReconstructedW_KLF(const std::vector<unsigned int> KLF_W11_Index,
        const std::vector<unsigned int> KLF_W12_Index,
        const std::vector<unsigned int> KLF_W21_Index,
        const std::vector<unsigned int> KLF_W22_Index,
        const std::vector<int> jetTruthIdx,
        const std::vector<int> nMatchedObjects){

            // Form the indices pairs and use the HyPER function.
            long KLF_W11 = KLF_W11_Index.size() == 0 ? -1 : KLF_W11_Index.at(0);
            long KLF_W12 = KLF_W12_Index.size() == 0 ? -1 : KLF_W12_Index.at(0);
            long KLF_W21 = KLF_W21_Index.size() == 0 ? -1 : KLF_W21_Index.at(0);
            long KLF_W22 = KLF_W22_Index.size() == 0 ? -1 : KLF_W22_Index.at(0);
            std::vector<long> KLF_W1_Indices{KLF_W11,KLF_W12};
            std::vector<long> KLF_W2_Indices{KLF_W21,KLF_W22};

            return countCorrectlyReconstructedW(KLF_W1_Indices,KLF_W2_Indices,jetTruthIdx,nMatchedObjects);
        }

    int countCorrectlyReconstructedTop(const std::vector<long> HyPER_T1_Indices,
        const std::vector<long> HyPER_T2_Indices,
        const std::vector<long> HyPER_W1_Indices,
        const std::vector<long> HyPER_W2_Indices,
        const std::vector<int> jetTruthIdx,
        const std::vector<int> nMatchedObjects){

            // In TCT v2.17.0 the truth Tops are [0,1,2] and [3,4,5]. 0 and 3 are the b-jets. The assumption is that truth jets are uniquely matched.
            // Let's mark events which are not matcheable with -2.
            // And the events where HyPER did not converge properly with -1.

            // First check all sizes are appropriate.
            if (jetTruthIdx.size() != 6 || nMatchedObjects.size() != 6) return -2;

            // Check both truth Tops are uniquely matched.
            bool truthT1Unique = nMatchedObjects.at(0) == 1 && nMatchedObjects.at(1) == 1 && nMatchedObjects.at(2) == 1;
            bool truthT2Unique = nMatchedObjects.at(3) == 1 && nMatchedObjects.at(4) == 1 && nMatchedObjects.at(5) == 1;
            if (!truthT1Unique || !truthT2Unique) return -2;

            // Define the labels
            std::vector<int> truthT1{jetTruthIdx.at(0),jetTruthIdx.at(1),jetTruthIdx.at(2)};
            std::vector<int> truthT2{jetTruthIdx.at(3),jetTruthIdx.at(4),jetTruthIdx.at(5)};
            std::vector<int> truthW1{jetTruthIdx.at(1),jetTruthIdx.at(2)};
            std::vector<int> truthW2{jetTruthIdx.at(4),jetTruthIdx.at(5)};

            // Check we have a truth matcheable event
            bool truthW1Matcheable = truthT1.at(0) != -1 && truthT1.at(1) != -1 && truthT1.at(2) != -1;
            bool truthW2Matcheable = truthT2.at(0) != -1 && truthT2.at(1) != -1 && truthT2.at(2) != -1;
            if (!truthW1Matcheable || !truthW2Matcheable) return -2;

            // Now check the HyPER predictions.
            if (HyPER_T1_Indices.size() != 3 || HyPER_T2_Indices.size() != 3) return -1;
            if (HyPER_W1_Indices.size() != 2 || HyPER_W2_Indices.size() != 2) return -1;

            // Check HyPER at least has one valid Top prediction.
            bool validT1HyPER = allValidIndices<long>(HyPER_T1_Indices);
            bool validT2HyPER = allValidIndices<long>(HyPER_T2_Indices);
            if (!validT1HyPER && !validT2HyPER) return -1;
            // Check HyPER at least has one valid W prediction.
            bool validW1HyPER = allValidIndices<long>(HyPER_W1_Indices);
            bool validW2HyPER = allValidIndices<long>(HyPER_W2_Indices);
            if (!validW1HyPER && !validW2HyPER) return -1;

            // Check HyPER has unique predictions.
            bool TPredictionsOverlap = twoVectorsOverlap<long,long>(HyPER_T1_Indices,HyPER_T2_Indices);
            if (TPredictionsOverlap) return -1;
            bool WPredictionsOverlap = twoVectorsOverlap<long,long>(HyPER_W1_Indices,HyPER_W2_Indices);
            if (WPredictionsOverlap) return -1;

            int nCorrectTop = 0;
            if (twoVectorsEqual<long,int>(HyPER_T1_Indices,truthT1) && twoVectorsEqual<long,int>(HyPER_W1_Indices,truthW1)) nCorrectTop++;
            if (twoVectorsEqual<long,int>(HyPER_T1_Indices,truthT2) && twoVectorsEqual<long,int>(HyPER_W1_Indices,truthW2)) nCorrectTop++;
            if (twoVectorsEqual<long,int>(HyPER_T2_Indices,truthT1) && twoVectorsEqual<long,int>(HyPER_W2_Indices,truthW1)) nCorrectTop++;
            if (twoVectorsEqual<long,int>(HyPER_T2_Indices,truthT2) && twoVectorsEqual<long,int>(HyPER_W2_Indices,truthW2)) nCorrectTop++;

            if (nCorrectTop > 2) {
                throw std::invalid_argument("The number of correctly reconstructed Ws is greater than 2.");
            }

            return nCorrectTop;
        }

        int countCorrectlyReconstructedTop_KLF(const std::vector<unsigned int> KLF_W11_Index,
            const std::vector<unsigned int> KLF_W12_Index,
            const std::vector<unsigned int> KLF_b1_Index,
            const std::vector<unsigned int> KLF_W21_Index,
            const std::vector<unsigned int> KLF_W22_Index,
            const std::vector<unsigned int> KLF_b2_Index,
            const std::vector<int> jetTruthIdx,
            const std::vector<int> nMatchedObjects){

            // Form the indices pairs and use the HyPER function.
            long KLF_W11 = KLF_W11_Index.size() == 0 ? -1 : KLF_W11_Index.at(0);
            long KLF_W12 = KLF_W12_Index.size() == 0 ? -1 : KLF_W12_Index.at(0);
            long KLF_W21 = KLF_W21_Index.size() == 0 ? -1 : KLF_W21_Index.at(0);
            long KLF_W22 = KLF_W22_Index.size() == 0 ? -1 : KLF_W22_Index.at(0);
            long KLF_b1 = KLF_b1_Index.size() == 0 ? -1 : KLF_b1_Index.at(0);
            long KLF_b2 = KLF_b2_Index.size() == 0 ? -1 : KLF_b2_Index.at(0);
            std::vector<long> KLF_T1_Indices{KLF_b1,KLF_W11,KLF_W12};
            std::vector<long> KLF_T2_Indices{KLF_b2,KLF_W21,KLF_W22};
            std::vector<long> KLF_W1_Indices{KLF_W11,KLF_W12};
            std::vector<long> KLF_W2_Indices{KLF_W21,KLF_W22};

            return countCorrectlyReconstructedTop(KLF_T1_Indices,KLF_T2_Indices,KLF_W1_Indices,KLF_W2_Indices,jetTruthIdx,nMatchedObjects);
        }

        bool HyPERConvergedTwoTops(const std::vector<long> HyPER_T1_Indices,
            const std::vector<long> HyPER_T2_Indices,
            const std::vector<long> HyPER_W1_Indices,
            const std::vector<long> HyPER_W2_Indices){

            // Now check the HyPER predictions.
            if (HyPER_T1_Indices.size() != 3 || HyPER_T2_Indices.size() != 3) return false;
            if (HyPER_W1_Indices.size() != 2 || HyPER_W2_Indices.size() != 2) return false;

            // Check HyPER has two valid Top predictions.
            bool validT1HyPER = allValidIndices<long>(HyPER_T1_Indices);
            bool validT2HyPER = allValidIndices<long>(HyPER_T2_Indices);
            if (!validT1HyPER || !validT2HyPER) return false;
            // Check HyPER has two valid W predictions.
            bool validW1HyPER = allValidIndices<long>(HyPER_W1_Indices);
            bool validW2HyPER = allValidIndices<long>(HyPER_W2_Indices);
            if (!validW1HyPER || !validW2HyPER) return false;

            // Check HyPER has unique predictions.
            bool TPredictionsOverlap = twoVectorsOverlap<long,long>(HyPER_T1_Indices,HyPER_T2_Indices);
            if (TPredictionsOverlap) return false;
            bool WPredictionsOverlap = twoVectorsOverlap<long,long>(HyPER_W1_Indices,HyPER_W2_Indices);
            if (WPredictionsOverlap) return false;

            return true;
        }

        bool KLFConvergedTwoTops(const std::vector<unsigned int> KLF_W11_Index,
            const std::vector<unsigned int> KLF_W12_Index,
            const std::vector<unsigned int> KLF_b1_Index,
            const std::vector<unsigned int> KLF_W21_Index,
            const std::vector<unsigned int> KLF_W22_Index,
            const std::vector<unsigned int> KLF_b2_Index){

            // Form the indices pairs and use the HyPER function.
            long KLF_W11 = KLF_W11_Index.size() == 0 ? -1 : KLF_W11_Index.at(0);
            long KLF_W12 = KLF_W12_Index.size() == 0 ? -1 : KLF_W12_Index.at(0);
            long KLF_W21 = KLF_W21_Index.size() == 0 ? -1 : KLF_W21_Index.at(0);
            long KLF_W22 = KLF_W22_Index.size() == 0 ? -1 : KLF_W22_Index.at(0);
            long KLF_b1 = KLF_b1_Index.size() == 0 ? -1 : KLF_b1_Index.at(0);
            long KLF_b2 = KLF_b2_Index.size() == 0 ? -1 : KLF_b2_Index.at(0);
            std::vector<long> KLF_T1_Indices{KLF_b1,KLF_W11,KLF_W12};
            std::vector<long> KLF_T2_Indices{KLF_b2,KLF_W21,KLF_W22};
            std::vector<long> KLF_W1_Indices{KLF_W11,KLF_W12};
            std::vector<long> KLF_W2_Indices{KLF_W21,KLF_W22};

            return HyPERConvergedTwoTops(KLF_T1_Indices,KLF_T2_Indices,KLF_W1_Indices,KLF_W2_Indices);
        }

        int FindBIndexInHyPER(const std::vector<long>& HyPER_W_indices, const std::vector<long>& HyPER_Top_indices){
            // Sanity checks
            if (HyPER_W_indices.size() != 2 || HyPER_Top_indices.size() != 3) return -1;
            if (!allValidIndices<long>(HyPER_W_indices)) return -1; // Bad W
            if (!allValidIndices<long>(HyPER_Top_indices)) return -1; // Bad Top

            // Check W indices are contained in Top
            for (const auto& W_index : HyPER_W_indices){
                if (std::find(HyPER_Top_indices.begin(), HyPER_Top_indices.end(), W_index) == HyPER_Top_indices.end()){
                    throw std::invalid_argument("W index not found in Top indices.");
                }
            }

            // Find the b index: It is the one that is not in the W indices
            for (const auto& top_index : HyPER_Top_indices){
                if (std::find(HyPER_W_indices.begin(), HyPER_W_indices.end(), top_index) == HyPER_W_indices.end()){
                    return static_cast<int>(top_index);
                }
            }
            
            throw std::invalid_argument("The code should never get here!.");
        }
}

namespace HyPERPerformanceSingleLepton {

    std::vector<int> getRecoIndicesTopLep(const std::vector<int>& jetTruthIndices,
        const std::vector<int>& eTruthIndices,
        const std::vector<int>& muTruthIndices,
        const std::vector<int>& nMatchedJets,
        const std::vector<int>& nMatchedElectrons,
        const std::vector<int>& nMatchedMuons,
        int nLeptons){
        
        std::vector<int> defaultIndices{-1,-1,-1};
        // Check the sizes of the input vectors
        if (jetTruthIndices.size() != 6 || nMatchedJets.size() != 6) return defaultIndices;
        if (eTruthIndices.size() != 2 || nMatchedElectrons.size() != 2) return defaultIndices;
        if (muTruthIndices.size() != 2 || nMatchedMuons.size() != 2) return defaultIndices;
        if (nLeptons != 1) return defaultIndices;

        // Check that the lepton indices make sense.
        bool electronEvent = eTruthIndices.at(0) != -1 || eTruthIndices.at(1) != -1;
        bool muonEvent = muTruthIndices.at(0) != -1 || muTruthIndices.at(1) != -1;
        if (electronEvent && muonEvent) {
            LOG(DEBUG) << "Event contains both electron and muon matched, but at the same time it is a single lepton event.";
            return defaultIndices;
        }
        if (!electronEvent && !muonEvent) {
            LOG(DEBUG) << "Event does not contain any lepton matched, but at the same time it is a single lepton event.";
            return defaultIndices;
        }

        // Define the lepton flavour.
        std::vector<int> leptonTruthIndices = electronEvent ? eTruthIndices : muTruthIndices;
        std::vector<int> nMatchedLeptons = electronEvent ? nMatchedElectrons : nMatchedMuons;    

        bool tDecay = leptonTruthIndices.at(0) != -1;
        bool tBarDecay = leptonTruthIndices.at(1) != -1;
        if (tDecay && tBarDecay) {
            LOG(DEBUG) << "Event contains both t and tbar decay matched leptons. But at this point it should be Single-lepton event.";
            return defaultIndices;
        }

        // Define which top to look for. Check that only one reco lepton is matched.
        std::size_t decayIndex = tDecay ? 0 : 1;
        if (nMatchedLeptons.at(decayIndex) != 1) return defaultIndices;

        // At this point the W is matched correctly.
        defaultIndices = {-1, leptonTruthIndices.at(decayIndex), 0};

        // Now find the jet indices.
        std::size_t jetIndex = tDecay ? 0 : 3;
        if (jetTruthIndices.at(jetIndex) == -1) return defaultIndices;
        if (nMatchedJets.at(jetIndex) != 1) return defaultIndices;

        // Build the top and return it.
        defaultIndices.at(0) = jetTruthIndices.at(jetIndex);
        return defaultIndices;
    }

    std::vector<int> getRecoIndicesTopHad(const std::vector<int>& jetTruthIndices,
        const std::vector<int>& eTruthIndices,
        const std::vector<int>& muTruthIndices,
        const std::vector<int>& nMatchedJets,
        const std::vector<int>& nMatchedElectrons,
        const std::vector<int>& nMatchedMuons,
        int nLeptons){

        // First part, try to figure out which is the hadronic top.
        // Only return default indices if there is no way to find out.
        std::vector<int> defaultIndices{-1,-1,-1};
        // Check the sizes of the input vectors
        if (jetTruthIndices.size() != 6 || nMatchedJets.size() != 6) return defaultIndices;
        if (eTruthIndices.size() != 2 || nMatchedElectrons.size() != 2) return defaultIndices;
        if (muTruthIndices.size() != 2 || nMatchedMuons.size() != 2) return defaultIndices;
        if (nLeptons != 1) return defaultIndices;

        // Check that the lepton indices make sense.
        bool electronEvent = eTruthIndices.at(0) != -1 || eTruthIndices.at(1) != -1;
        bool muonEvent = muTruthIndices.at(0) != -1 || muTruthIndices.at(1) != -1;
        if (electronEvent && muonEvent) {
            LOG(DEBUG) << "Event contains both electron and muon matched, but at the same time it is a single lepton event.";
            return defaultIndices;
        }
        if (!electronEvent && !muonEvent) {
            LOG(DEBUG) << "Event does not contain any lepton matched, but at the same time it is a single lepton event.";
            return defaultIndices;
        }

        // Define the lepton flavour.
        std::vector<int> leptonTruthIndices = electronEvent ? eTruthIndices : muTruthIndices;
        std::vector<int> nMatchedLeptons = electronEvent ? nMatchedElectrons : nMatchedMuons;    

        bool tDecayLeptonic = leptonTruthIndices.at(0) != -1;
        bool tBarDecayLeptonic = leptonTruthIndices.at(1) != -1;
        if (tDecayLeptonic && tBarDecayLeptonic) {
            LOG(DEBUG) << "Event contains both t and tbar decay matched leptons. But at this point it should be Single-lepton event.";
            return defaultIndices;
        }

        // Now find the jet indices.
        std::size_t bJetIndex;
        std::size_t q1JetIndex;
        std::size_t q2JetIndex;
        if (tDecayLeptonic){
            bJetIndex = 3;
            q1JetIndex = 4;
            q2JetIndex = 5;
        } else {
            bJetIndex = 0;
            q1JetIndex = 1;
            q2JetIndex = 2;
        } 
        
        // Form the W
        if (jetTruthIndices.at(q1JetIndex) == -1 || jetTruthIndices.at(q2JetIndex) == -1) return defaultIndices;
        if (nMatchedJets.at(q1JetIndex) != 1 || nMatchedJets.at(q2JetIndex) != 1) return defaultIndices;
        defaultIndices = {-1, jetTruthIndices.at(q1JetIndex), jetTruthIndices.at(q2JetIndex)};

        // Find the b-jet and form the top if possible.
        if (jetTruthIndices.at(bJetIndex) == -1) return defaultIndices;
        if (nMatchedJets.at(bJetIndex) != 1) return defaultIndices;
        defaultIndices.at(0) = jetTruthIndices.at(bJetIndex);

        return defaultIndices;
    }

    bool isTopLepReconstructable(const std::vector<int>& TruthIndices){
        // Check the sizes of the input vectors
        if (TruthIndices.size() != 3) return false;

        // Check that the lepton indices are all valid.
        return HyPERPerformance::allValidIndices<int>(TruthIndices);
    }

    bool isTopHadReconstructable(const std::vector<int>& TruthIndices){
        //std::cout << "Hello World!" << std::endl;
        // Check the sizes of the input vectors
        if (TruthIndices.size() != 3) return false;

       //std::cout<< "TruthIndices: " << TruthIndices.at(0) << " " << TruthIndices.at(1) << " " << TruthIndices.at(2) << std::endl;

        // Check that the lepton indices are all valid.
        if (!HyPERPerformance::allValidIndices<int>(TruthIndices)) return false;

        //std::cout<< "All indices are valid." << std::endl;

        // Check that the indices are unique.
        return HyPERPerformance::allElementsUnique<int>(TruthIndices);
    }

    int getTopLepLeptonFlavour(const std::vector<int>& eTruthIndices,
        const std::vector<int>& muTruthIndices,
        const std::vector<int>& nMatchedElectrons,
        const std::vector<int>& nMatchedMuons,
        int nLeptons){

        int defaultFlavour{-1};
        // Check the sizes of the input vectors
        if (eTruthIndices.size() != 2 || nMatchedElectrons.size() != 2) return defaultFlavour;
        if (muTruthIndices.size() != 2 || nMatchedMuons.size() != 2) return defaultFlavour;
        if (nLeptons != 1) return defaultFlavour;

        // Check that the lepton indices make sense.
        bool electronEvent = eTruthIndices.at(0) != -1 || eTruthIndices.at(1) != -1;
        bool muonEvent = muTruthIndices.at(0) != -1 || muTruthIndices.at(1) != -1;
        if (electronEvent && muonEvent) {
            LOG(DEBUG) << "Event contains both electron and muon matched, but at the same time it is a single lepton event.";
            return defaultFlavour;
        }
        if (!electronEvent && !muonEvent) {
            LOG(DEBUG) << "Event does not contain any lepton matched, but at the same time it is a single lepton event.";
            return defaultFlavour;
        }

        // Define the lepton flavour.
        std::vector<int> leptonTruthIndices = electronEvent ? eTruthIndices : muTruthIndices;
        std::vector<int> nMatchedLeptons = electronEvent ? nMatchedElectrons : nMatchedMuons;    

        bool tDecay = leptonTruthIndices.at(0) != -1;
        bool tBarDecay = leptonTruthIndices.at(1) != -1;
        if (tDecay && tBarDecay) {
            LOG(DEBUG) << "Event contains both t and tbar decay matched leptons. But at this point it should be Single-lepton event.";
            return defaultFlavour;
        }

        // Define which top to look for. Check that only one reco lepton is matched.
        std::size_t decayIndex = tDecay ? 0 : 1;
        if (nMatchedLeptons.at(decayIndex) != 1) return defaultFlavour;

        defaultFlavour = electronEvent ? 2 : 3;
        return defaultFlavour;
    }

    bool HyPERTopLepConverged(const std::vector<long>& HyPERTopLepIndices,
        const std::vector<long>& HyPERTopLepIDs){
            // Check the sizes of the vectors are correct.
            if (HyPERTopLepIndices.size() != 3 || HyPERTopLepIDs.size() != 3) return false;

            // Check there is only one jet.
            std::size_t nJets = 0;
            for (const auto& ID : HyPERTopLepIDs){
                if (ID == 1) nJets++;
            }
            if (nJets != 1) return false;

            // Check there is only one MET.
            std::size_t nMET = 0;
            for (const auto& ID : HyPERTopLepIDs){
                if (ID == 4) nMET++;
            }
            if (nMET != 1) return false;

            // Check there is only one lepton.
            std::size_t nLepton = 0;
            for (const auto& ID : HyPERTopLepIDs){
                if (ID == 2 || ID == 3) nLepton++;
            }
            if (nLepton != 1) return false;

            // Check all indices are valid.
            if (!HyPERPerformance::allValidIndices<long>(HyPERTopLepIndices)) return false;

            return true;
    }

    bool HyPERTopHadConverged(const std::vector<long>& HyPERTopHadIndices,
        const std::vector<long>& HyPERTopHadIDs){

            // Check the sizes of the vectors are correct.
            if (HyPERTopHadIndices.size() != 3 || HyPERTopHadIDs.size() != 3) return false;

            // Check all are jets.
            for (const auto& ID : HyPERTopHadIDs){
                if (ID != 1) return false;
            }

            // Check all indices are valid.
            if (!HyPERPerformance::allValidIndices<long>(HyPERTopHadIndices)) return false;

            // Check all indices are unique.
            return HyPERPerformance::allElementsUnique<long>(HyPERTopHadIndices);
    }

    bool HyPERConvergedTwoTops(const bool topLepConverged,
        const bool topHadConverged,
        const std::vector<long>& HyPERTopLepIndices,
        const std::vector<long>& HyPERTopLepIDs,
        const std::vector<long>& HyPERTopHadIndices,
        const std::vector<long>& HyPERTopHadIDs,
        const std::vector<long>& HyPERWHadIndices){

        // If any of the two tops is not converged, return false.
        if (!topLepConverged || !topHadConverged) return false;

        // Check the sizes of the vectors are correct.
        if (HyPERTopLepIndices.size() != 3 || HyPERTopLepIDs.size() != 3) return false;
        if (HyPERTopHadIndices.size() != 3 || HyPERTopHadIDs.size() != 3) return false;
        if (HyPERWHadIndices.size() != 2) return false;

        // Look for the b-jet in the topLep indices.
        int topLepBJetIndex = -1;
        for (std::size_t i = 0; i < 3; i++){
            if (HyPERTopLepIDs.at(i) == 1){
                topLepBJetIndex = HyPERTopLepIndices.at(i);
                break;
            }
        }
        if (topLepBJetIndex == -1) {
            LOG(ERROR) << "The topLep indices do not contain a b-jet. At this point they should.";
            throw std::invalid_argument("");
        }

        // Look for the b-jet in the topHad indices.
        if (!HyPERPerformance::allValidIndices<long>(HyPERWHadIndices)){
            LOG(ERROR) << "The WHad indices are not valid. While the TopHad indices are.";
            throw std::invalid_argument("");
        }

        std::vector<int> topHadBJetCandidates {}; // This should only contain one element at the end.
        // But we allow for more possibilites to be safe and check.
        for (std::size_t tophad = 0; tophad < 3; tophad++){
            bool matchFound = false;
            int topHadBJetIndex = HyPERTopHadIndices.at(tophad);
            for (std::size_t whad = 0; whad < 2; whad++){
                if (HyPERTopHadIndices.at(tophad) == HyPERWHadIndices.at(whad)){
                    matchFound = true;
                    break;
                }
            }
            if (!matchFound) topHadBJetCandidates.push_back(topHadBJetIndex);
        }

        if (topHadBJetCandidates.size() != 1) {
            LOG(ERROR) << "There is more than one index in topHad that is not in WHad.";
            throw std::invalid_argument("");
        }
        if (topHadBJetCandidates.at(0) == -1) {
            LOG(ERROR) << "The topHad indices do not contain a b-jet. At this point they should.";
            throw std::invalid_argument("");
        }

        // Check that both b-jets are different, if they are reconstruction was successful.
        return (topLepBJetIndex != topHadBJetCandidates.at(0));

    }

    bool correctlyReconstructedWHad(std::vector<int> TopHadTruthRecoIndices,
                                    const std::vector<long>& HyPERWHadIndices){ 

        // Check the sizes of the vectors are correct.
        if (TopHadTruthRecoIndices.size() != 3 || HyPERWHadIndices.size() != 2) return false;

        // Check HyPER indices are valid.
        if (!HyPERPerformance::allValidIndices<long>(HyPERWHadIndices)) return false;
        if (!HyPERPerformance::allElementsUnique<long>(HyPERWHadIndices)) return false;

        // Check truth indices are valid.
        if (!HyPERPerformance::allValidIndices<int>(TopHadTruthRecoIndices)) return false;
        if (!HyPERPerformance::allElementsUnique<int>(TopHadTruthRecoIndices)) return false;
        
        // Select the truth W indices.
        std::vector<long> truthWHadIndices{TopHadTruthRecoIndices.at(1),TopHadTruthRecoIndices.at(2)};

        return HyPERPerformance::twoVectorsEqual<long,long>(truthWHadIndices,HyPERWHadIndices);
    }

    bool correctlyReconstructedWLep(std::vector<int> TopLepTruthRecoIndices,
                                    const std::vector<long>& HyPERWLepIndices,
                                    int truthLeptonFlavour,
                                    const std::vector<long>& HyPERTopLepIDs){ 

         // Check the sizes of the vectors are correct.
        if (TopLepTruthRecoIndices.size() != 3 || HyPERWLepIndices.size() != 2 || HyPERTopLepIDs.size() != 3) return false;
        // Check HyPER indices are valid.
        if (!HyPERPerformance::allValidIndices<long>(HyPERWLepIndices)) return false;

        // Check truth indices are valid.
        if (!HyPERPerformance::allValidIndices<int>(TopLepTruthRecoIndices)) return false;

        // Check the lepton flavour is correct.
        if (truthLeptonFlavour != 2 && truthLeptonFlavour != 3) return false;
        bool matchedFlavour = false;
        int nMatches = 0;
        for (const auto& ID : HyPERTopLepIDs){
            if (ID == truthLeptonFlavour) {
                matchedFlavour = true;
                nMatches++;
            }
        }

        if (!matchedFlavour) return false;
        if (nMatches != 1) {
            LOG(ERROR) << "The lepton flavour is not unique in the HyPERTopLepIDs.";
            throw std::invalid_argument("");
        }

        // Select the truth W indices.
        std::vector<long> truthWLepIndices{TopLepTruthRecoIndices.at(1),TopLepTruthRecoIndices.at(2)};

        return HyPERPerformance::twoVectorsEqual<long,long>(truthWLepIndices,HyPERWLepIndices);
    }

    int getHyPERBHadIndex(const std::vector<long>& HyPERTopHadIndices,
        const std::vector<long>& HyPERTopHadIDs,
        const std::vector<long>& HyPERWHadIndices){

        // Check the sizes of the vectors are correct.
        if (HyPERTopHadIndices.size() != 3 || HyPERTopHadIDs.size() != 3) return -1;
        if (HyPERWHadIndices.size() != 2) return -1;

        // Check all indices are valid.
        if (!HyPERPerformance::allValidIndices<long>(HyPERTopHadIndices)) return -1;
        if (!HyPERPerformance::allValidIndices<long>(HyPERTopHadIDs)) return -1;
        if (!HyPERPerformance::allValidIndices<long>(HyPERTopHadIndices)) return -1;

        // Check all IDs are equal to 1.
        for (const auto& ID : HyPERTopHadIDs){
            if (ID != 1) return -1;
        }

        std::vector<int> topHadBJetCandidates {}; // This should only contain one element at the end.
        // But we allow for more possibilites to be safe and check.
        for (std::size_t tophad = 0; tophad < 3; tophad++){
            bool matchFound = false;
            int topHadBJetIndex = HyPERTopHadIndices.at(tophad);
            for (std::size_t whad = 0; whad < 2; whad++){
                if (HyPERTopHadIndices.at(tophad) == HyPERWHadIndices.at(whad)){
                    matchFound = true;
                    break;
                }
            }
            if (!matchFound) topHadBJetCandidates.push_back(topHadBJetIndex);
        }

        if (topHadBJetCandidates.size() != 1) {
            LOG(ERROR) << "There is more than one index in topHad that is not in WHad.";
            throw std::invalid_argument("");
        }
        if (topHadBJetCandidates.at(0) == -1) {
            LOG(ERROR) << "The topHad indices do not contain a b-jet. At this point they should.";
            throw std::invalid_argument("");
        }
        
        return topHadBJetCandidates.at(0);
    }

    int getHyPERBLepIndex(const std::vector<long>& HyPERTopLepIndices,
            const std::vector<long>& HyPERTopLepIDs){

        // Check the sizes of the vectors are correct.
        if (HyPERTopLepIndices.size() != 3 || HyPERTopLepIDs.size() != 3) return -1;
        // Check all indices are valid.
        if (!HyPERPerformance::allValidIndices<long>(HyPERTopLepIndices)) return -1;
        if (!HyPERPerformance::allValidIndices<long>(HyPERTopLepIDs)) return -1;

        // Check the position of the jet in the IDs. Allow for more than one jet and throw error if found more than one.
        int topLepBJetIndex = -1;
        std::size_t nBJets = 0;
        for (std::size_t i = 0; i < 3; i++){
            if (HyPERTopLepIDs.at(i) == 1){
                topLepBJetIndex = HyPERTopLepIndices.at(i);
                nBJets++;
            }
        }

        if (topLepBJetIndex == -1) {
            LOG(ERROR) << "The topLep indices do not contain a b-jet. At this point they should.";
            throw std::invalid_argument("");
        }
        if (nBJets != 1) {
            LOG(ERROR) << "The topLep indices contain more than one b-jet. At this point they should.";
            throw std::invalid_argument("");
        }

        return topLepBJetIndex;
    }

    bool correctlyReconstructedTopLep(const bool correctWLep,
        int HyPERBLepIndex,
        const std::vector<int>& TopLepTruthRecoIndices){

            // Check the W is correctly reconstructed.
            if (!correctWLep) return false;

            // Check the size of the array are correct.
            if (TopLepTruthRecoIndices.size() != 3) return false;

            if (HyPERBLepIndex == -1) return false;

            // Check the b-jet is correctly matched.
            return HyPERBLepIndex == TopLepTruthRecoIndices.at(0);

    }

    bool correctlyReconstructedTopHad(const bool correctWHad,
        int HyPERBHadIndex,
        const std::vector<int>& TopHadTruthRecoIndices){

            // Check the W is correctly reconstructed.
            if (!correctWHad) return false;

            // Check the size of the array are correct.
            if (TopHadTruthRecoIndices.size() != 3) return false;

            if (HyPERBHadIndex == -1) return false;

            // Check the b-jet is correctly matched.
            return HyPERBHadIndex == TopHadTruthRecoIndices.at(0);
    }

    bool correctlyReconstructedTopLep_KLF(int KLFBLepIndex,
                                          const std::vector<int>& TopLepTruthRecoIndices){

        // Use HyPER function.
        return correctlyReconstructedTopLep(true,KLFBLepIndex,TopLepTruthRecoIndices);
    }

    bool correctlyReconstructedWHad_KLF(std::vector<int> TopHadTruthRecoIndices,
                                        int KLW1Index,
                                        int KLW2Index){

        // Form the indices.
        std::vector<long> KLFWHadIndices{KLW1Index,KLW2Index};
        
        // Use HyPER function.
        return correctlyReconstructedWHad(TopHadTruthRecoIndices, KLFWHadIndices);
    }

    bool correctlyReconstructedTopHad_KLF(const bool correctWHad,
    int KLFBHadIndex,
    const std::vector<int>& TopHadTruthRecoIndices){

        // Use HyPER function.
        return correctlyReconstructedTopHad(correctWHad,KLFBHadIndex,TopHadTruthRecoIndices);
    }
}

namespace HyPERPerformanceDiLepton {

    std::vector<int> getRecoIndicesTop1(const std::vector<int>& jetTruthIndices,
        const std::vector<int>& eTruthIndices,
        const std::vector<int>& muTruthIndices,
        const std::vector<int>& nMatchedJets,
        const std::vector<int>& nMatchedElectrons,
        const std::vector<int>& nMatchedMuons,
        int nLeptons){

        std::vector<int> defaultIndices{-1,-1,-1};
        // Basic checks
        if (jetTruthIndices.size() != 6 || nMatchedJets.size() != 6) return defaultIndices;
        if (eTruthIndices.size() != 2 || nMatchedElectrons.size() != 2) return defaultIndices;
        if (muTruthIndices.size() != 2 || nMatchedMuons.size() != 2) return defaultIndices;
        if (nLeptons != 2) return defaultIndices;

        // For dilepton both tops should have a matched lepton.
        if ((eTruthIndices.at(0) == -1 && muTruthIndices.at(0) == -1) ||
            (eTruthIndices.at(1) == -1 && muTruthIndices.at(1) == -1)) return defaultIndices;

        // Check that each top has exactly one matched lepton (electron or muon).
        if ((nMatchedElectrons.at(0) + nMatchedMuons.at(0)) != 1) return defaultIndices;
        if ((nMatchedElectrons.at(1) + nMatchedMuons.at(1)) != 1) return defaultIndices;

        // Determine lepton truth index for top1 (use electron if present, otherwise muon)
        int lepIndexTop1 = eTruthIndices.at(0) != -1 ? eTruthIndices.at(0) : muTruthIndices.at(0);

        // b-jet for top1 is stored in jetTruthIndices position 0
        std::size_t bJetPos = 0;
        if (jetTruthIndices.at(bJetPos) == -1) return defaultIndices;
        if (nMatchedJets.at(bJetPos) != 1) return defaultIndices;

        defaultIndices.at(0) = jetTruthIndices.at(bJetPos);
        defaultIndices.at(1) = lepIndexTop1;
        defaultIndices.at(2) = 0; // neutrino placeholder

        return defaultIndices;
    }

    std::vector<int> getRecoIndicesTop2(const std::vector<int>& jetTruthIndices,
        const std::vector<int>& eTruthIndices,
        const std::vector<int>& muTruthIndices,
        const std::vector<int>& nMatchedJets,
        const std::vector<int>& nMatchedElectrons,
        const std::vector<int>& nMatchedMuons,
        int nLeptons){

        std::vector<int> defaultIndices{-1,-1,-1};
        // Basic checks
        if (jetTruthIndices.size() != 6 || nMatchedJets.size() != 6) return defaultIndices;
        if (eTruthIndices.size() != 2 || nMatchedElectrons.size() != 2) return defaultIndices;
        if (muTruthIndices.size() != 2 || nMatchedMuons.size() != 2) return defaultIndices;
        if (nLeptons != 2) return defaultIndices;

        // For dilepton both tops should have a matched lepton.
        if ((eTruthIndices.at(0) == -1 && muTruthIndices.at(0) == -1) ||
            (eTruthIndices.at(1) == -1 && muTruthIndices.at(1) == -1)) return defaultIndices;

        // Check that each top has exactly one matched lepton (electron or muon).
        if ((nMatchedElectrons.at(0) + nMatchedMuons.at(0)) != 1) return defaultIndices;
        if ((nMatchedElectrons.at(1) + nMatchedMuons.at(1)) != 1) return defaultIndices;

        // Determine lepton truth index for top2 (use electron if present, otherwise muon)
        int lepIndexTop2 = eTruthIndices.at(1) != -1 ? eTruthIndices.at(1) : muTruthIndices.at(1);

        // b-jet for top2 is stored in jetTruthIndices position 3
        std::size_t bJetPos = 3;
        if (jetTruthIndices.at(bJetPos) == -1) return defaultIndices;
        if (nMatchedJets.at(bJetPos) != 1) return defaultIndices;

        defaultIndices.at(0) = jetTruthIndices.at(bJetPos);
        defaultIndices.at(1) = lepIndexTop2;
        defaultIndices.at(2) = 0; // neutrino placeholder

        return defaultIndices;
    }

    bool HyPERTopConverged(const std::vector<long>& HyPERTopIndices,
        const std::vector<long>& HyPERTopIDs){

        // Same logic as HyPERTopLepConverged but generalised name for dilepton use.
        if (HyPERTopIndices.size() != 3 || HyPERTopIDs.size() != 3) return false;

        // Check there is only one jet.
        std::size_t nJets = 0;
        for (const auto& ID : HyPERTopIDs){
            if (ID == 1) nJets++;
        }
        if (nJets != 1) return false;

        // Check there is only one MET.
        std::size_t nMET = 0;
        for (const auto& ID : HyPERTopIDs){
            if (ID == 4) nMET++;
        }
        if (nMET != 1) return false;

        // Check there is only one lepton.
        std::size_t nLepton = 0;
        for (const auto& ID : HyPERTopIDs){
            if (ID == 2 || ID == 3) nLepton++;
        }
        if (nLepton != 1) return false;

        // Check all indices are valid.
        if (!HyPERPerformance::allValidIndices<long>(HyPERTopIndices)) return false;

        return true;
    }

    bool HyPERConvergedTwoLeptonicTops(const std::vector<long>& HyPERTop1Indices,
        const std::vector<long>& HyPERTop1IDs,
        const std::vector<long>& HyPERTop2Indices,
        const std::vector<long>& HyPERTop2IDs){

        // Individual convergence checks
        if (!HyPERTopConverged(HyPERTop1Indices,HyPERTop1IDs)) return false;
        if (!HyPERTopConverged(HyPERTop2Indices,HyPERTop2IDs)) return false;

        // Extract b-jet indices
        int b1 = HyPERPerformanceSingleLepton::getHyPERBLepIndex(HyPERTop1Indices,HyPERTop1IDs);
        int b2 = HyPERPerformanceSingleLepton::getHyPERBLepIndex(HyPERTop2Indices,HyPERTop2IDs);

        if (b1 == -1 || b2 == -1) return false;

        // Require distinct b-jets
        return (b1 != b2);
    }

    bool correctlyReconstructedTop(int HyPERBIndex,
        const std::vector<int>& TopTruthRecoIndices){

        if (TopTruthRecoIndices.size() != 3) return false;
        if (HyPERBIndex == -1) return false;
        return HyPERBIndex == TopTruthRecoIndices.at(0);
    }

}