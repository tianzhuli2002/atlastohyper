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

    /**
     * @brief Turn Char_t variable to bool.
     * @param myVariable variable that should be used as boolean
     */
    bool charToBool(Char_t myVariable){
        return static_cast<int>(myVariable) == 1;
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
            return nTupleFlag && numberOfMuons == 1 && triggered && triggerMatched && numberOfJets >= 4;//numberOfJets == 4;
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

}