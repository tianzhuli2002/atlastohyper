#pragma once
// Inlcude the math library
#include<vector>
#include <string>
#include <algorithm>
#include "Math/Vector4D.h"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

// Variables used in both channels
namespace Variables{
    /**
     * @brief Calculate the pt ordered indices of the objects that pass one selecttion.
     * @return std::vector<std::size_t> The indices of the objects that pass the selection.
     * @param ptVec The pT vector of the objects.
     * @param selection The first selection.
     */
    std::vector<std::size_t> indexPtOrderedGoodOneSelec(const std::vector<float>& ptVec,
                                        const std::vector<char>& selection);

    /**
     * @brief Calculate the pt ordered indices of the objects that pass two selections.
     * @return std::vector<std::size_t> The indices of the objects that pass the two selections.
     * @param ptVec The pT vector of the objects.
     * @param selection1 The first selection.
     * @param selection2 The second selection.
     */
    std::vector<std::size_t> indexPtOrderedGoodTwoSelec(const std::vector<float>& ptVec,
                                        const std::vector<char>& selection1,
                                        const std::vector<char>& selection2);
    /**
     * @brief Calculate the number of jets that pass one selection and have a pT greater than 20 GeV.
     * @return int The number of jets that pass the selections.
     * @param ptVec The pT vector of the jets.
     * @param selected The first selection.
     */
    int nGoodJets(const std::vector<float>& ptVec,
                        const std::vector<char>& selected);

    /**
     * @brief Create a TLV for a parton from the four momentum branches.
     * @return ROOT::RDF::RNode The node with the new TLV branch.
     * @param mainFrame The main frame.
     * @param mainNode The current main node.
     * @param partonBaseName The base name of the parton branches.
     * @param newBranchName The name of the new TLV branch.
     */
    ROOT::RDF::RNode createPartonTLV(MainFrame* mainFrame ,ROOT::RDF::RNode mainNode,
                                     const std::string& partonBaseName, const std::string& newBranchName);

    /**
     * @brief Create a container from the indices of pT sorted objects and another container.
     * @param T The type of the output container (should match the input).
     * @param pTOrderedJetIndices The indices of the jets ordered by pT.
     * @param inputContainer The input container.
     * @return std::vector<T> The output container.
     */
    template<typename T>
    std::vector<T> pTSortedContainer(const std::vector<T>& inputContainer, const std::vector<std::size_t>& pTOrderedJetIndices){
        std::vector<T> outputContainer;
        for(const auto& i : pTOrderedJetIndices){
            outputContainer.push_back(inputContainer.at(i));
        }
        return outputContainer;
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
    std::vector<int> getPCBT(const std::vector<char>& FCBEff_65, const std::vector<char>& FCBEff_70, const std::vector<char>& FCBEff_77, const std::vector<char>& FCBEff_85, const std::vector<char>& FCBEff_90);

    /** 
     * @brief For the parton match branches convert from the ntuple index to the pT index.
     * @param defaultValue The value to use if the index is not found. Unmatched particle.
     * @return The converted index.
     */
    class ConverPartonMatchtNtupleIndexToPtIndex{
        public:
        ConverPartonMatchtNtupleIndexToPtIndex(int defaultValue = -2) : m_defaultValue(defaultValue) {}
        // Overload operator ()
        int operator()(int ntupleIndex, const std::vector<std::size_t>& ptIndices){
            if (ntupleIndex == -1) return ntupleIndex; // This means that at the TCT level the event has been found to be di-lepton or ambigous matching.
            // Try to find the index in the ptIndices vector. If found return the position. If not return the default value.
            auto it = std::find(ptIndices.begin(),ptIndices.end(),ntupleIndex);
            if (it != ptIndices.end()){
                return std::distance(ptIndices.begin(), it);
            } else {
                return m_defaultValue;
            }
        }

        private:
        int m_defaultValue = -2;
    };

    /**
     * @brief For the parton match branches convert from the ntuple index to the pT index.
     * @param defaultValue The value to use if the index is not found. Unmatched particle.
     * @return The converted index.
     */
    template <typename T>
    class ConvertNtupleIndexToPtIndex{
        public:
        ConvertNtupleIndexToPtIndex(T defaultValue) : m_defaultValue(defaultValue) {}
        // Overload operator ()
        T operator()(T ntupleIndex, const std::vector<std::size_t>& ptIndices){
            // Try to find the index in the ptIndices vector. If found return the position. If not return the default value.
            auto it = std::find(ptIndices.begin(),ptIndices.end(),ntupleIndex);
            if (it != ptIndices.end()){
                return std::distance(ptIndices.begin(), it);
            } else {
                return m_defaultValue;
            }
        }

        private:
        T m_defaultValue;
    };

     /**
     * @brief For the parton match branches convert from the ntuple index to the pT index.
     * @param defaultValue The value to use if the index is not found. Unmatched particle.
     * @return The converted index.
     */
    template <typename T>
    class ConvertOneElementVectorNtupleIndexToPtIndex{
        public:
        ConvertOneElementVectorNtupleIndexToPtIndex(T defaultValue) : m_defaultValue(defaultValue) {}
        // Overload operator ()
        T operator()(std::vector<T> ntupleIndex, const std::vector<std::size_t>& ptIndices){
            if (ntupleIndex.size() == 0) return m_defaultValue;
            // Try to find the index in the ptIndices vector. If found return the position. If not return the default value.
            auto it = std::find(ptIndices.begin(),ptIndices.end(),ntupleIndex[0]);
            if (it != ptIndices.end()){
                return std::distance(ptIndices.begin(), it);
            } else {
                return m_defaultValue;
            }
        }

        private:
        T m_defaultValue;
    };

    /**
     * @brief Count number of b-jets above certain GN2v01 working point
     * @param pTOrderedJetIndices The indices of the jets ordered by pT.
     * @param jetBTagVec The vector storing the pseudocontinuos b-tagging scores.
     * @param jetBTagWP The working point for b-tagging.
     * @return int The number of b-jets above the working point.
     */
    class nBJets{
        public:
        nBJets(const std::size_t jetBTagWP) : m_jetBTagWP(jetBTagWP){}
        int getCuantile(){
            switch (m_jetBTagWP){
            case 100:
                return 1;
            case 90:
                return 2;
            case 85:
                return 3;
            case 77:
                return 4;
            case 70:
                return 5;
            case 65:
                return 6;
            default:
                LOG(ERROR) << "Invalid b-tagging working point for this algorithm.";
                throw std::invalid_argument("");
            }
        }

        int operator()(const std::vector<std::size_t>& pTOrderedJetIndices, const std::vector<int>& jetBTagVec){
            int nBJets = 0;
            for(const auto& i : pTOrderedJetIndices){
                if(jetBTagVec.at(i) >= this->getCuantile()){
                    nBJets++;
                }
            }
            return nBJets;
        }

        private:
        const std::size_t m_jetBTagWP;
    };

    /**
     * @brief Calculate the invariant mass of two TLVs.
     * @param tlv1 The first TLV.
     * @param tlv2 The second TLV.
     * @return float The invariant mass of the two TLVs.
     */
    float invariantMassTwo(const ROOT::Math::PtEtaPhiMVector& tlv1,
                           const ROOT::Math::PtEtaPhiMVector& tlv2);

    /**
     * @brief Calculate the invariant mass of three TLVs.
     * @param tlv1 The first TLV.
     * @param tlv2 The second TLV.
     * @param tlv3 The third TLV.
     * @return float The invariant mass of the three TLVs.
     */
    float invariantMassThree(const ROOT::Math::PtEtaPhiMVector& tlv1,
                            const ROOT::Math::PtEtaPhiMVector& tlv2,
                            const ROOT::Math::PtEtaPhiMVector& tlv3);

    /**
     * @brief Calculate the invariant mass of two objects identified from indices, the objects are picked up from a vector of TLVs.
     * @param tlvVec The TLV vector.
     * @param index1 The index of the first object.
     * @param index2 The index of the second object.
     * @return float The invariant mass of the three TLVs.
     */
    float recoInvariantMassTwo(int index1, int index2, const std::vector<ROOT::Math::PtEtaPhiEVector>& tlvVector);

    /**
     * @brief Calculate the invariant mass of three objects identified from indices, the objects are picked up from a vector of TLVs.
     * @param tlvVec The TLV vector.
     * @param index1 The index of the first object.
     * @param index2 The index of the second object.
     * @param index3 The index of the third object.
     * @return float The invariant mass of the three TLVs.
     */
    float recoInvariantMassThree(int index1, int index2, int index3, const std::vector<ROOT::Math::PtEtaPhiEVector>& tlvVector);
    
    /**
     * @brief Calculate the invariant mass of three objects identified from indices, the objects are picked up from a vector of TLVs.
     * The indices are stored in vectors. This is a KLF version.
     * @param tlvVec The TLV vector.
     * @param index1 The index of the first object.
     * @param index2 The index of the second object.
     * @param index3 The index of the third object.
     * @return float The invariant mass of the three TLVs.
     */
    template <typename T>
    float recoInvariantMassThree_KLF(const std::vector<T>& index1, const std::vector<T>& index2, const std::vector<T>& index3, const std::vector<ROOT::Math::PtEtaPhiEVector>& tlvVector){
        ROOT::Math::PtEtaPhiEVector tlv1 = index1.size() > 0 ? tlvVector.at(index1[0]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        ROOT::Math::PtEtaPhiEVector tlv2 = index2.size() > 0 ? tlvVector.at(index2[0]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        ROOT::Math::PtEtaPhiEVector tlv3 = index3.size() > 0 ? tlvVector.at(index3[0]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);

        return (tlv1+tlv2+tlv3).M()/1000.0f;
    }

    /**
     * @brief Calculate the invariant mass of two objects identified from indices, the objects are picked up from a vector of TLVs.
     * The indices are stored in vectors. This is a KLF version.
     * @param tlvVec The TLV vector.
     * @param index1 The index of the first object.
     * @param index2 The index of the second object.
     * @param index3 The index of the third object.
     * @return float The invariant mass of the two TLVs.
     */
    template <typename T>
    float recoInvariantMassTwo_KLF(const std::vector<T>& index1, const std::vector<T>& index2, const std::vector<ROOT::Math::PtEtaPhiEVector>& tlvVector){
        ROOT::Math::PtEtaPhiEVector tlv1 = index1.size() > 0 ? tlvVector.at(index1[0]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        ROOT::Math::PtEtaPhiEVector tlv2 = index2.size() > 0 ? tlvVector.at(index2[0]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);

        return (tlv1+tlv2).M()/1000.0f;
    }

    /**
     * @brief Calculate the invariant mass of three objects identified from indices, the objects are picked up from a vector of TLVs.
     * The indices are stored in vectors. This is a HyPER version.
     * @param tlvVec The TLV vector.
     * @param indices The vector holding the indices.
     * @return float The invariant mass of the three TLVs.
     */
    template <typename T>
    float recoInvariantMassThree_HyPER(const std::vector<T>& indices, const std::vector<ROOT::Math::PtEtaPhiEVector>& tlvVector){
        ROOT::Math::PtEtaPhiEVector tlv1 = indices.at(0) != -1 ? tlvVector.at(indices[0]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        ROOT::Math::PtEtaPhiEVector tlv2 = indices.at(1) != -1 ? tlvVector.at(indices[1]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        ROOT::Math::PtEtaPhiEVector tlv3 = indices.at(2) != -1 ? tlvVector.at(indices[2]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);

        return (tlv1+tlv2+tlv3).M()/1000.0f;
    }

    /**
     * @brief Calculate the invariant mass of two objects identified from indices, the objects are picked up from a vector of TLVs.
     * The indices are stored in vectors. This is a HyPER version.
     * @param tlvVec The TLV vector.
      * @param indices The vector holding the indices.
     * @return float The invariant mass of the two TLVs.
     */
    template <typename T>
    float recoInvariantMassTwo_HyPER(const std::vector<T>& indices, const std::vector<ROOT::Math::PtEtaPhiEVector>& tlvVector){
        ROOT::Math::PtEtaPhiEVector tlv1 = indices.at(0) != -1 ? tlvVector.at(indices[0]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);
        ROOT::Math::PtEtaPhiEVector tlv2 = indices.at(1) != -1 ? tlvVector.at(indices[1]) : ROOT::Math::PtEtaPhiEVector(0,0,0,0);

        return (tlv1+tlv2).M()/1000.0f;
    }

    /**
     * @brief Match a parton to a list of reco particles.
     * @param parton The parton TLV.
     * @param recoParticlesOrderedByPt The reco particles ordered by pT.
     * @return int The index of the matched reco particle.
     */
    int matchPartonRecoParticles(const ROOT::Math::PtEtaPhiMVector& parton, // Notice difference between E and M !!!
                                const std::vector<ROOT::Math::PtEtaPhiEVector>& recoParticlesOrderedByPt);

    /**
     * @brief Check if a W is matched correctly.
     * @param q1Index The index of the first quark.
     * @param q2Index The index of the second quark.
     */
    bool isWMatchedCorrectly(int q1Index, int q2Index);

    /**
     * @brief Check if two W overlap.
     * @param q1Index The index of the first quark from W1.
     * @param q2Index The index of the second quark from W1.
     * @param q1barIndex The index of the first quark from W2.
     * @param q2barIndex The index of the second quark from W2.
     */
    bool isNotWOverlap(int q1Index, int q2Index, int q1barIndex, int q2barIndex);

    /**
     * @brief Check if a Top is matched correctly.
     * @param q1Index The index of the first quark.
     * @param q2Index The index of the second quark.
     * @param bIndex The index of the b quark.
     */
    bool isTopMatchedCorrectly(int q1Index, int q2Index, int bIndex);

    /**
     * @brief Check if two b quarks overlap.
     * @param bIndex The index of the b quark from top1.
     * @param bbarIndex The index of the b quark from top2.
     */
    bool isNotBOverlap(int bIndex, int bbarIndex);

    /**
     * @brief Adjust the HyPER indices in function of the selections applied on the object containers
     * @return std::vector<long> The updated version of the HyPER indices.
     * @param HyPER_indices Initial HyPER indices vector.
     * @param HyPER_IDs Corresponding HyPER IDs vector.
     * @param jet_idx Post selections jet indices.
     * @param el_idx Post selections electron indices.
     * @param mu_idx Post selections muon indices.
     * @param tau_idx Post selections tau indices.
     */
    std::vector<long> getSelecHyPERIndices(const std::vector<long>& HyPER_indices,
                                        const std::vector<long>& HyPER_IDs,
                                        const std::vector<size_t>& jet_idx,
                                        const std::vector<size_t>& el_idx,
                                        const std::vector<size_t>& mu_idx,
                                        const std::vector<size_t>& tau_idx);

    /**
     * @brief Turn Char_t variable to bool.
     * @param myVariable variable that should be used as boolean
     */
    bool charToBool(Char_t myVariable);

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
                        const std::vector<char>& selected2);

    /**
     * @brief Calculate the number of good jets.
     * @param ptVec The pT vector of the jets.
     * @param selected1 The first selection.
     * @param selected2 The second selection.
     * @return int The number of good jets.
     */
    int nGoodJets(const std::vector<float>& ptVec,
                        const std::vector<char>& selected1,
                        const std::vector<char>& selected2);

    /**
     * @brief Calculate the indices of pT ordered objects passing two selections.
     * @param ptVec The pT vector of the objects.
     * @param selection1 The first selection.
     * @param selection2 The second selection.
     * @return std::vector<std::size_t> The pT ordered indices of the objects passing the selections.
     */
    std::vector<std::size_t> indexPtOrderedGoodTwoSelec(const std::vector<float>& ptVec,
                                        const std::vector<char>& selection1,
                                        const std::vector<char>& selection2);
    /**
     * @brief Return the unique element in the KLF index vector.
     * @param klfIndex The KLF index vector.
     * @return int The unique element in the KLF index vector.
     */
    int fromKLFToNtupleIndex(const std::vector<unsigned int>& klfIndex);

    /**
     * @brief Get index of true down jet
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @return Int_t index of true down jet
     */
    Int_t getPartonTruthDown(const std::vector<int>& event_jet_truth_idx);

    /**
     * @brief Get index of true up jet
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @return Int_t index of true up jet
     */
    Int_t getPartonTruthUp(const std::vector<int>& event_jet_truth_idx);

    /**
     * @brief Get index of true hadronic b jet
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @return Int_t index of true hadronic b jet
     */
    Int_t getPartonTruthBHad(const std::vector<int>& event_jet_truth_idx);

    /**
     * @brief Get index of true leptonic b jet
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @return Int_t index of true leptonic b jet
     */
    Int_t getPartonTruthBLep(const std::vector<int>& event_jet_truth_idx);

    /**
     * @brief Create a container from a variable stored in a column. The container is created with the indices of the pT ordered objects.
     * @param scale The unit scale to apply to the variable.
     * @return T The type of the variables inside the container to register.
     */
    template <typename T>
    class CreatePtOrderedContainer{
        public:
        CreatePtOrderedContainer(float scale = 1.0f) : m_scale(scale) {}
        // Overload operator ()
        std::vector<T> operator()(const std::vector<T>& inputCollection, const std::vector<std::size_t>& ptIndices){
            std::vector<T> container{};
            for (std::size_t index : ptIndices){
                if (inputCollection.size() <= index) {
                    LOG(ERROR) << "An index in the ptIndices vector is greater than the size of the passed collection." << std::endl;
                    throw std::runtime_error("");
                }
                else container.push_back(inputCollection[index]/m_scale);
            }
            return container;
        }

        private:
        float m_scale;
    };

    // Selection cuts for the muon channel
    namespace MuonSelection {
        // Events firing any of the triggers
        bool triggers(const bool& T1_2022,const bool& T2_2022,
                    const bool& T1_2018,const bool& T2_2018);

        // Basic selection with number of good objects, triggers and correct ntuple branch.
        bool basic(const char& nTupleFlag, int numberOfMuons,
                   int numberOfJets,bool triggered, bool triggerMatched);
    }

    // Selection cuts for the electron channel
    namespace ElectronSelection {
        // Events firing any of the triggers
        bool triggers(const bool& T1_2022,const bool& T2_2022, const bool& T3_2022,
                    const bool& T1_2018,const bool& T2_2018, const bool& T3_2018);

        // Basic selection with number of good objects, triggers and correct ntuple branch.
        bool basic(const char& nTupleFlag, int numberOfElectrons,
                   int numberOfJets,bool triggered, bool triggerMatched);
    }

    //Selection cuts on jets
    namespace JetsSelection {
        bool fullyMatched(Int_t truth_down_index, Int_t truth_up_index, Int_t truth_bhad_index, Int_t truth_blep_index);
    }
}


namespace HyPERPerformance{
    template <typename T>
    bool allValidIndices(const std::vector<T>& indices){
        for (const auto& index : indices){
            if (index == -1) return false;
        }
        return true;
    }

    template <typename T, typename U>
    bool twoVectorsEqual(std::vector<T> vec1, std::vector<U> vec2){
        if (vec1.size() != vec2.size()) return false;
        // Order the vectors
        std::sort(vec1.begin(), vec1.end());
        std::sort(vec2.begin(), vec2.end());
        //Compare
        for (std::size_t i = 0; i < vec1.size(); i++){
            if (vec1[i] != vec2[i]) return false;
        }
        return true;
    }

    template <typename T, typename U>
    bool twoVectorsOverlap(const std::vector<T>& vec1, const std::vector<U>& vec2){
        if (vec1.size() != vec2.size()) throw std::invalid_argument("The vectors must have the same size.");
        // Compare
        for (const auto& index1 : vec1){
            for (const auto& index2 : vec2){
                if (index1 == index2) return true;
            }
        }
        return false;
    }

    /**
     * @brief Count the number of correctly reconstructed W bosons.
     * @param HyPER_W1_Indices The HyPER indices of the W1 jets.
     * @param HyPER_W2_Indices The HyPER indices of the W2 jets.
     * @param jetTruthIdx The indices of the truth matched jets.
     * @param nMatchedObjects The number of matched objects per truth object.
     */
    int countCorrectlyReconstructedW(const std::vector<long> HyPER_W1_Indices,
                                     const std::vector<long> HyPER_W2_Indices,
                                     const std::vector<int> jetTruthIdx,
                                     const std::vector<int> nMatchedObjects);

    /**
     * @brief Count the number of correctly reconstructed W bosons.
     * @param KLF_W11_Index The KLF indices of the W1 jet1.
     * @param KLF_W12_Index The KLF indices of the W1 jet2.
     * @param KLF_W21_Index The KLF indices of the W2 jet1.
     * @param KLF_W22_Index The KLF indices of the W2 jet2.
     * @param jetTruthIdx The indices of the truth matched jets.
     * @param nMatchedObjects The number of matched objects per truth object.
     */
    int countCorrectlyReconstructedW_KLF(const std::vector<unsigned int> KLF_W11_Index,
        const std::vector<unsigned int> KLF_W12_Index,
        const std::vector<unsigned int> KLF_W21_Index,
        const std::vector<unsigned int> KLF_W22_Index,
        const std::vector<int> jetTruthIdx,
        const std::vector<int> nMatchedObjects);

    /**
     * @brief Count the number of correctly reconstructed Top quarks.
     * @param HyPER_T1_Indices The HyPER indices of the Top1 jets.
     * @param HyPER_T2_Indices The HyPER indices of the Top2 jets.
     * @param HyPER_W1_Indices The HyPER indices of the W1 jets.
     * @param HyPER_W2_Indices The HyPER indices of the W2 jets.
     * @param jetTruthIdx The indices of the truth matched jets.
     * @param nMatchedObjects The number of matched objects per truth object.
     */
    int countCorrectlyReconstructedTop(const std::vector<long> HyPER_T1_Indices,
        const std::vector<long> HyPER_T2_Indices,
        const std::vector<long> HyPER_W1_Indices,
        const std::vector<long> HyPER_W2_Indices,
        const std::vector<int> jetTruthIdx,
        const std::vector<int> nMatchedObjects);

    /**
     * @brief Count the number of correctly reconstructed Top quarks.
     * @param KLF_W11_Index The KLF indices of the W1 jet1.
     * @param KLF_W12_Index The KLF indices of the W1 jet2.
     * @param KLF_b1_Index The KLF indices of the b1 jet.
     * @param KLF_W21_Index The KLF indices of the W2 jet1.
     * @param KLF_W22_Index The KLF indices of the W2 jet2.
     * @param KLF_b2_Index The KLF indices of the b2 jet.
     * @param jetTruthIdx The indices of the truth matched jets.
     * @param nMatchedObjects The number of matched objects per truth object.
     */
    int countCorrectlyReconstructedTop_KLF(const std::vector<unsigned int> KLF_W11_Index,
        const std::vector<unsigned int> KLF_W12_Index,
        const std::vector<unsigned int> KLF_b1_Index,
        const std::vector<unsigned int> KLF_W21_Index,
        const std::vector<unsigned int> KLF_W22_Index,
        const std::vector<unsigned int> KLF_b2_Index,
        const std::vector<int> jetTruthIdx,
        const std::vector<int> nMatchedObjects);

    /**
     * @brief Check if the HyPER algorithm converged for Two Tops.
     * @param HyPER_T1_Indices The HyPER indices of the Top1 jets.
     * @param HyPER_T2_Indices The HyPER indices of the Top2 jets.
     * @param HyPER_W1_Indices The HyPER indices of the W1 jets.
     * @param HyPER_W2_Indices The HyPER indices of the W2 jets.
     * @return bool True if the algorithm converged and gives two valid tops.
     */
    bool HyPERConvergedTwoTops(const std::vector<long> HyPER_T1_Indices,
        const std::vector<long> HyPER_T2_Indices,
        const std::vector<long> HyPER_W1_Indices,
        const std::vector<long> HyPER_W2_Indices);
    
    /**
     * @brief Check if the KLF algorithm converged for Two Tops.
     */
    bool KLFConvergedTwoTops(const std::vector<unsigned int> KLF_W11_Index,
        const std::vector<unsigned int> KLF_W12_Index,
        const std::vector<unsigned int> KLF_b1_Index,
        const std::vector<unsigned int> KLF_W21_Index,
        const std::vector<unsigned int> KLF_W22_Index,
        const std::vector<unsigned int> KLF_b2_Index);

    class ConvertHyPERIndexToPtIndex{
        public:
        ConvertHyPERIndexToPtIndex(std::size_t index) : m_index(index) {}
        // Overload operator ()
        int operator()(const std::vector<long>& HyPER_indices, const std::vector<std::size_t>& ptIndices){
            // Sanity checks
            if (m_index == -1) throw std::invalid_argument("The index is not set.");
            if (HyPER_indices.size() <= static_cast<std::size_t>(m_index)) return -1;
            long nTupleIndex = HyPER_indices.at(m_index);
            if (nTupleIndex == -1) return -1;

            Variables::ConvertNtupleIndexToPtIndex<int> converter(-1);
            return converter(nTupleIndex,ptIndices);   
        }

        private:
        int m_index = -1;
    };
    
    int FindBIndexInHyPER(const std::vector<long>& HyPER_W_indices, const std::vector<long>& HyPER_Top_indices);
}