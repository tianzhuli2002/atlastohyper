#include<vector>
#include <string>
#include <algorithm>
#include "Math/Vector4D.h"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "FastFrames/MainFrame.h"
#include "FastFrames/DefineHelpers.h"
#include "PartonMatcher/Variables.h"

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
            return nTupleFlag && numberOfMuons == 1 && triggered && triggerMatched && numberOfJets == 4;
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
            return nTupleFlag && numberOfElectrons == 1 && triggered && triggerMatched && numberOfJets == 4;
        }
    }
}