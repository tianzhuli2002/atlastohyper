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
            // First check all sizes are appropriate.
            if (HyPER_W1_Indices.size() != 2 || HyPER_W2_Indices.size() != 2 || jetTruthIdx.size() != 6 || nMatchedObjects.size() != 6) return -1;

            // Check HyPER at least has one valid prediction.
            bool validW1HyPER = allValidIndices<long>(HyPER_W1_Indices);
            bool validW2HyPER = allValidIndices<long>(HyPER_W2_Indices);
            if (!validW1HyPER && !validW2HyPER) return -1;

            // Check both truth Ws are uniquely matched.
            bool truthW1Unique = nMatchedObjects.at(1) == 1 && nMatchedObjects.at(2) == 1;
            bool truthW2Unique = nMatchedObjects.at(4) == 1 && nMatchedObjects.at(5) == 1;
            if (!truthW1Unique || !truthW2Unique) return -1;

            // Check HyPER has unique predictions.
            bool WPredictionsOverlap = twoVectorsOverlap<long,long>(HyPER_W1_Indices,HyPER_W2_Indices);
            if (WPredictionsOverlap) return -1;

            // Define the labels
            std::vector<int> truthW1{jetTruthIdx.at(1),jetTruthIdx.at(2)};
            std::vector<int> truthW2{jetTruthIdx.at(4),jetTruthIdx.at(5)};

            // Check we have a truth matcheable event
            bool truthW1Matcheable = truthW1.at(0) != -1 && truthW1.at(1) != -1;
            bool truthW2Matcheable = truthW2.at(0) != -1 && truthW2.at(1) != -1;
            if (!truthW1Matcheable || !truthW2Matcheable) return -1;

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

    int countCorrectlyReconstructedTop(const std::vector<long> HyPER_T1_Indices,
        const std::vector<long> HyPER_T2_Indices,
        const std::vector<long> HyPER_W1_Indices,
        const std::vector<long> HyPER_W2_Indices,
        const std::vector<int> jetTruthIdx,
        const std::vector<int> nMatchedObjects){

            // In TCT v2.17.0 the truth Tops are [0,1,2] and [3,4,5]. 0 and 3 are the b-jets. The assumption is that truth jets are uniquely matched.
            // First check all sizes are appropriate.
            if (HyPER_T1_Indices.size() != 3 || HyPER_T2_Indices.size() != 3 || jetTruthIdx.size() != 6 || nMatchedObjects.size() != 6) return -1;
            if (HyPER_W1_Indices.size() != 2 || HyPER_W2_Indices.size() != 2) return -1;

            // Check HyPER at least has one valid Top prediction.
            bool validT1HyPER = allValidIndices<long>(HyPER_T1_Indices);
            bool validT2HyPER = allValidIndices<long>(HyPER_T2_Indices);
            if (!validT1HyPER && !validT2HyPER) return -1;
            // Check HyPER at least has one valid W prediction.
            bool validW1HyPER = allValidIndices<long>(HyPER_W1_Indices);
            bool validW2HyPER = allValidIndices<long>(HyPER_W2_Indices);
            if (!validW1HyPER && !validW2HyPER) return -1;

            // Check both truth Ws are uniquely matched.
            bool truthT1Unique = nMatchedObjects.at(0) == 1 && nMatchedObjects.at(1) == 1 && nMatchedObjects.at(2) == 1;
            bool truthT2Unique = nMatchedObjects.at(3) == 1 && nMatchedObjects.at(4) == 1 && nMatchedObjects.at(5) == 1;
            if (!truthT1Unique || !truthT2Unique) return -1;

            // Check HyPER has unique predictions.
            bool TPredictionsOverlap = twoVectorsOverlap<long,long>(HyPER_T1_Indices,HyPER_T2_Indices);
            if (TPredictionsOverlap) return -1;
            bool WPredictionsOverlap = twoVectorsOverlap<long,long>(HyPER_W1_Indices,HyPER_W2_Indices);
            if (WPredictionsOverlap) return -1;

            // Define the labels
            std::vector<int> truthT1{jetTruthIdx.at(0),jetTruthIdx.at(1),jetTruthIdx.at(2)};
            std::vector<int> truthT2{jetTruthIdx.at(3),jetTruthIdx.at(4),jetTruthIdx.at(5)};
            std::vector<int> truthW1{jetTruthIdx.at(1),jetTruthIdx.at(2)};
            std::vector<int> truthW2{jetTruthIdx.at(4),jetTruthIdx.at(5)};

            // Check we have a truth matcheable event
            bool truthW1Matcheable = truthT1.at(0) != -1 && truthT1.at(1) != -1 && truthT1.at(2) != -1;
            bool truthW2Matcheable = truthT2.at(0) != -1 && truthT2.at(1) != -1 && truthT2.at(2) != -1;
            if (!truthW1Matcheable || !truthW2Matcheable) return -1;

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
}