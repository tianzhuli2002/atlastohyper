#pragma once
// Inlcude the math library
#include<vector>
#include <string>
#include <algorithm>
#include "Math/Vector4D.h"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "ROOT/RVec.hxx"

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

    /** 
     * @brief Convert a one-element RVector of type T to a simple type of type V.
     * @param InputVector The input vector.
     */
    template <typename T, typename V>
    V convertOneElementRVectorToIntegerType(const ROOT::VecOps::RVec<T>& InputVector){
        if (InputVector.size() == 0) return -1;
        else if (InputVector.size() == 1) return static_cast<V>(InputVector[0]);
        else {
            LOG(ERROR) << "The input vector has more than one element." << std::endl;
            throw std::runtime_error("");
        }
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
        // Events firing any of the available muon triggers
    bool triggersnew(const bool& mu24, const bool& mu50);

    // Basic selection with number of good objects, triggers and correct ntuple branch.
    bool basicnew(const char& nTupleFlag, int numberOfMuons, int numberOfJets,
               bool triggered, bool triggerMatched);
           
    }

    // Selection cuts for the electron channel
    namespace ElectronSelection {
        // Events firing any of the triggers
        bool triggers(const bool& T1_2022,const bool& T2_2022, const bool& T3_2022,
                    const bool& T1_2018,const bool& T2_2018, const bool& T3_2018);

        // Basic selection with number of good objects, triggers and correct ntuple branch.
        bool basic(const char& nTupleFlag, int numberOfElectrons,
                   int numberOfJets,bool triggered, bool triggerMatched);
        // Events firing any of the available electron triggers

    bool triggersnew(const bool& e26, const bool& e60, const bool& e140);

    // Basic selection with number of good objects, triggers and correct ntuple branch.
    bool basicnew(const char& nTupleFlag, int numberOfElectrons, int numberOfJets,
               bool triggered, bool triggerMatched);
           
    }

    //Selection cuts on jets
    namespace JetsSelection {
        bool fullyMatched(Int_t truth_down_index, Int_t truth_up_index, Int_t truth_bhad_index, Int_t truth_blep_index);
    }
}

namespace TtbarDileptonVars{
    /**
     * @brief Function to get the sum of the charges of the leptons.
     * @param chargeV Vector of charges of the leptons.
     * @return float Sum of the charges.
     */
    float sumOfCharges( const ROOT::VecOps::RVec<float>& chargeV);

    /**
     * @brief Function to get the pt-sorted vector of charges of the leptons that pass a selection.
     * And have a pT > 7 GeV.
     * @param chargeV Vector of charges of the leptons.
     * @param ptV Vector of pT of the leptons.
     * @param selection Vector of selection flags for the leptons.
     * @return ROOT::VecOps::RVec<float> Vector of charges of the leptons that pass the selection, sorted by pt.
     */
    ROOT::VecOps::RVec<float> sortedPassedChargeVector25(
        const ROOT::VecOps::RVec<float>& chargeV,
        const ROOT::VecOps::RVec<float>& ptV,
        const ROOT::VecOps::RVec<char>& selection);

    /**
     * @brief Function to determine the region based on the pT of the leading lepton
     * and the number of leptons of each type.
     * @param nElectrons Number of electrons.
     * @param nMuons Number of muons.
     * @param elTLV Vector of TLorentzVectors of the electrons (pt-sorted).
     * @param muTLV Vector of TLorentzVectors of the muons (pt-sorted).
     * @return std::string Region name.
     */
    std::string regionName(
        std::size_t nElectrons,
        std::size_t nMuons,
        const ROOT::VecOps::RVec<TLV>& elTLV,
        const ROOT::VecOps::RVec<TLV>& muTLV);

    /**
     * @brief Functor class to create a TLorentzVector for a truth particle that must have a given particle ID.
     * If the truth particle does not have the given ID, it returns a TLorentzVector with zero values.
     * @param particleID Particle ID of the truth particle.
     */
    class makeTruthTLV {
    using TLV = ROOT::Math::PtEtaPhiMVector;
    public:
        makeTruthTLV(int particleID) : m_particleID(particleID) {}

        TLV operator() (float pt,
                    float eta,
                    float phi,
                    float m,
                    int pdgId) const {
            // Create a vector to hold the TLorentzVectors
            TLV tlv(0,0,0,0);

            // Check if the particle ID matches the given ID
            if (pdgId != m_particleID) return tlv;

            // If it matches, fill the TLorentzVector with the given values.
            tlv.SetCoordinates(pt, eta, phi, m);
            return tlv;
        }

    private:
        int m_particleID;
    };
    
    /**
     * @brief Function to get the index of the reco tagged b-jet that matches the truth b-jet.
     * Returns -1 if no match is found.
     * @param recoBJets Vector of TLorentzVectors of the reco b-jets.
     * @param truthBJet TLorentzVector of the truth b-jet.
     * @return int Index of the reco b-jet that matches the truth b-jet.
     */
    int recoIndexTruthBJet(const ROOT::VecOps::RVec<TLV>& recoBJets,
                        const ROOT::Math::PtEtaPhiMVector& truthBJet);

    /**
     * @brief Invariant mass of the ttbar lepton pair.
     * -1 if no pair is found.
     * @param regionName Name of the region.
     * @param muonTLV Vector of TLorentzVectors of the muons.
     * @param electronTLV Vector of TLorentzVectors of the electrons.
     * @return float Invariant mass of the ttbar lepton pair.
     */
    float ttbarLeptonPairInvariantMass(
        const std::string& regionName,
        const ROOT::VecOps::RVec<TLV>& muonTLV,
        const ROOT::VecOps::RVec<TLV>& electronTLV
        );
    /**
     * @brief Minimum invariant mass of lepton-jet pairs in the event.
     * @param sorted_el_TLV Vector of TLorentzVectors of the electrons (pt
     * sorted).
     * @param sorted_mu_TLV Vector of TLorentzVectors of the muons (pt
     * sorted).
     * @param sorted_bjet_TLV Vector of TLorentzVectors of the b-jets (pt
     * sorted).
     * @param jet_TLV Vector of TLorentzVectors of the jets.
     * @param jet_GN2v01_FixedCutBEff_85_select Vector of selection flags for the jets.
     * @param n_electrons Number of electrons.
     * @param n_muons Number of muons.
     * @return float Minimum invariant mass of lepton-jet pairs in the event.
     */
    /**
     * @brief Minimum invariant mass of lepton-jet pairs in the event.
     * @param sorted_el_TLV Vector of TLorentzVectors of the electrons (pt
     * sorted).
     * @param sorted_mu_TLV Vector of TLorentzVectors of the muons (
     *  pt sorted).
     * @param sorted_bjet_TLV Vector of TLorentzVectors of the b-jets (pt
     * sorted).
     * @param jet_TLV Vector of TLorentzVectors of the jets.        
     * @param jet_GN2v01_FixedCutBEff_85_select Vector of selection flags for the jets.
     * @param n_electrons Number of electrons.
     * @param n_muons Number of muons.
     * @return float Minimum invariant mass of lepton-jet pairs in the event.
     */
    float minimumLeptonJetPairInvariantMass(
        const ROOT::VecOps::RVec<TLV>& sorted_el_TLV,
        const ROOT::VecOps::RVec<TLV>& sorted_mu_TLV,
        const ROOT::VecOps::RVec<TLV>& sorted_bjet_TLV,
        const ROOT::VecOps::RVec<TLV>& jet_TLV,
        const ROOT::VecOps::RVec<char>& jet_GN2v01_FixedCutBEff_85_select,
        unsigned long n_electrons,
        unsigned long n_muons);
    
        /**
     * @brief Minimum invariant mass of lepton-jet pairs in the event.
     * @param sorted_el_TLV Vector of TLorentzVectors of the electrons (pt
     * sorted).
     * @param sorted_mu_TLV Vector of TLorentzVectors of the muons (
     *  pt sorted).
     * @param sorted_bjet_TLV Vector of TLorentzVectors of the b-jets (pt
     * sorted).
     * @param jet_TLV Vector of TLorentzVectors of the jets.        
     * @param jet_GN2v01_FixedCutBEff_85_select Vector of selection flags for the jets.
     * @param n_electrons Number of electrons.
     * @param n_muons Number of muons.
     * @return float Minimum invariant mass of lepton-jet pairs in the event.
     */
ROOT::VecOps::RVec<int> selectTwoJetIndices_BtagPriority(const ROOT::VecOps::RVec<TLV>& jet_TLV,
                                 const ROOT::VecOps::RVec<char>& jet_GN2v01_FixedCutBEff_77_select);

                                 
     int SelectedJetsMatchTruthBPair(
    const ROOT::VecOps::RVec<int>& selected_jet_indices,
    const std::vector<int>& event_jet_truth_idx,
    const std::vector<int>& event_jet_truth_candidates);


    bool pass_min_lb_invar_mass(
        const ROOT::VecOps::RVec<TLV>& sorted_el_TLV,
        const ROOT::VecOps::RVec<TLV>& sorted_mu_TLV,
        const ROOT::VecOps::RVec<TLV>& sorted_bjet_TLV,
        const ROOT::VecOps::RVec<TLV>& jet_TLV,
        const ROOT::VecOps::RVec<char>& jet_GN2v01_FixedCutBEff_85_select,
        unsigned long n_electrons,
        unsigned long n_muons);


    /**
     * @brief Determine if the event is fully matched.
     * @param jet_TLVs Vector of TLorentzVectors of the jets.
     * @param el_TLVs Vector of TLorentzVectors of the electrons.
     * @param mu_TLVs Vector of TLorentzVectors of the muons.
     * @param el_ChargeV Vector of charges of the electrons.
     * @param mu_ChargeV Vector of charges of the muons.
     * @param event_jet_truth_idx Vector of truth indices of the jets.
     * @param event_jet_truth_candidates Vector of truth candidates of the jets.
     * @return ROOT::VecOps::RVec<float> Vector with the invariant mass of the truth matched lepton+b-jet pairs, 
     * as well as many other variables witht the truth matching information.
     */
    ROOT::VecOps::RVec<float> FullyMatched_MLB(
        const ROOT::VecOps::RVec<TLV>& jet_TLVs,
        const ROOT::VecOps::RVec<TLV>& el_TLVs,
        const ROOT::VecOps::RVec<TLV>& mu_TLVs,
        const ROOT::VecOps::RVec<float>& el_ChargeV,
        const ROOT::VecOps::RVec<float>& mu_ChargeV,
        const ROOT::VecOps::RVec<int>& event_jet_truth_idx,
        const ROOT::VecOps::RVec<int>& event_jet_truth_candidates
    );

    /** 
     * @brief Invariant mass of the llbb system.
     * @param sorted_el_TLV Vector of TLorentzVectors of the electrons (pt
     * sorted).
     * @param sorted_mu_TLV Vector of TLorentzVectors of the muons (
     *  pt sorted).
     * @param sorted_bjet_TLV Vector of TLorentzVectors of the b-jets (pt
     * sorted).
     * @param jet_TLV Vector of TLorentzVectors of the jets.    
     * @param jet_GN2v01_FixedCutBEff_85_select Vector of selection flags for the jets.
     * @param n_electrons Number of electrons.
     * @param n_muons Number of muons.
     * @return float Invariant mass of the llbb system. 
    */
    float invariantMassllbb(
        const ROOT::VecOps::RVec<TLV>& sorted_el_TLV,
        const ROOT::VecOps::RVec<TLV>& sorted_mu_TLV,
        const ROOT::VecOps::RVec<TLV>& sorted_bjet_TLV,
        const ROOT::VecOps::RVec<TLV>& jet_TLV,
        const ROOT::VecOps::RVec<char>& jet_GN2v01_FixedCutBEff_85_select,
        unsigned long n_electrons,
        unsigned long n_muons);    
    /**
     *@brief Determines whether region is Signal or not
     * @param elChargeV Vector of electron charges
     * @param muChargeV Vector of muon charges
     * @param ttbar_mass_NOSYS Invariant mass of the ttbar system without systematic variations
     * @return float 1 if the event is in the signal region, 0 otherwise
     */
    // float Signal(
    //     const ROOT::VecOps::RVec<float>& elChargeV,
    //     const ROOT::VecOps::RVec<float>& muChargeV,
    //     float ttbar_mass_NOSYS,
    //     float lepton_bjet_invariant_mass,
    //     unsigned long nElectrons,
    //     unsigned long nMuons,
    //     float met_met_NOSYS
    // );
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
);
    /**
     * @brief Check if the event is fully matched.
     * @param event_jet_truth_idx Vector of truth indices of the jets.
     * @param event_jet_truth_candidates Vector of truth candidates of the jets.
     * @return int 1 if the event is fully matched, 0 otherwise.
     */
    int FullyMatched(const std::vector<int>& event_jet_truth_idx,
                     const std::vector<int>& event_jet_truth_candidates
    );

    /**
     * @brief Get index of true b hadron
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @param event_jet_truth_candidates vector of int describing the topology-dependent parton-level object corresponding to a detector-level jet index
     * @return Int_t index of true b hadron
     */
    Int_t getTruthIndexBHad(const std::vector<int>& event_jet_truth_idx,
    const std::vector<int>& event_jet_truth_candidates);

    /**
     * @brief Get index of true bbar hadron
     * @param event_jet_truth_idx  vector of int describing the detector-level jet index corresponding to a topology-dependent parton-level object 
     * @param event_jet_truth_candidates vector of int describing the topology-dependent parton-level object corresponding to a detector-level jet index
     * @return Int_t index of true bbar hadron
     */
    Int_t getTruthIndexBbarHad(const std::vector<int>& event_jet_truth_idx,
    const std::vector<int>& event_jet_truth_candidates);

    /**
     * @brief Calculate the transverse mass of a system given its visible TLV and the MET.
     * @param vis The visible TLV of the system.
     * @param met_pt The MET pT.
     * @param met_phi The MET phi.
     * @return float The transverse mass of the system.
     */
    float transverseMassSystem_safe(const TLV& vis, float met_pt, float met_phi
    );
    /**
     * @brief Calculate the transverse mass of the ttbar system using truth matching information.
     * @param jet_TLVs Vector of TLorentzVectors of the jets.
     * @param el_TLVs Vector of TLorentzVectors of the electrons.
     * @param mu_TLVs Vector of TLorentzVectors of the muons.
     * @param event_jet_truth_idx Vector of truth indices of the jets.
     * @param event_jet_truth_candidates Vector of truth candidates of the jets.
     * @param el_charge Vector of charges of the electrons.
     * @param mu_charge Vector of charges of the muons.
     * @param met_pt The MET pT.
     * @param met_phi The MET phi.
     * @return float The transverse mass of the ttbar system.
     */
    float mT_ttbar_truth_from_truthIndices(
    const ROOT::VecOps::RVec<TLV>& jet_TLVs,
    const ROOT::VecOps::RVec<TLV>& el_TLVs,
    const ROOT::VecOps::RVec<TLV>& mu_TLVs,
    const ROOT::VecOps::RVec<int>& event_jet_truth_idx,
    const ROOT::VecOps::RVec<int>& event_jet_truth_candidates,
    const ROOT::VecOps::RVec<float>& el_charge,
    const ROOT::VecOps::RVec<float>& mu_charge,
    float met_pt,
    float met_phi
);
    
 /**
  * @brief Calculate deltaPhi between two TLVs.
  * @param a The first TLV.
  * @param b The second TLV.
  * @return float The deltaPhi between the two TLVs.
  */
inline float deltaPhi(const TLV& a,
                           const TLV& b
                        );

/**
 * @brief Calculate deltaR between two TLVs.
 * @param a The first TLV.
 * @param b The second TLV.
 * @return float The deltaR between the two TLVs.
 */
inline float deltaR(const TLV &a, const TLV &b);

/** 
 * @brief Pair leptons to jets by minimizing the sum of deltaR between lepton and jet.
 * @param el_TLV Vector of TLorentzVectors of the electrons.
 * @param mu_TLV Vector of TLorentzVectors of the muons.
 * @param jet_TLV Vector of TLorentzVectors of the jets.
 * @param event_jet_truth_idx Vector of truth indices of the jets.
 * @param event_jet_truth_candidates Vector of truth candidates of the jets.
 * @param el_charge Vector of charges of the electrons.
 * @param mu_charge Vector of charges of the muons.
 * @param n_electrons Number of electrons.
 * @param n_muons Number of muons.
 * @return ROOT::VecOps::RVec<float> Vector with the sum of deltaR between lepton and jet for the best pairing.
*/
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
); 
/** 
 * @brief Pair leptons to top jets by minimizing the sum of squared differences between reconstructed and true top masses.
 * @param el_TLV Vector of TLorentzVectors of the electrons.
 * @param mu_TLV Vector of TLorentzVectors of the muons.
 * @param jet_TLV Vector of TLorentzVectors of the jets.
 * @param event_jet_truth_idx Vector of truth indices of the jets.  
 * @param event_jet_truth_candidates Vector of truth candidates of the jets.
 * @param el_charge Vector of charges of the electrons.
 * @param mu_charge Vector of charges of the muons.
 * @param n_electrons Number of electrons.
 * @param n_muons Number of muons.
 * @return ROOT::VecOps::RVec<float> if whether pairing was successful (1.0) or not (0.0).
 */
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
    unsigned long n_muons);               

    /** 
     * @brief Pair leptons to jets by minimizing a Chi2 based on m_lb.
     * @param el_TLV Vector of TLorentzVectors of the electrons.
     * @param mu_TLV Vector of TLorentzVectors of the muons.
     * @param jet_TLV Vector of TLorentzVectors of the jets.
     * @param met_met MET magnitude.
     * @param met_phi MET phi.
     * @param event_jet_truth_idx Vector of truth indices of the jets.
     * @param event_jet_truth_candidates Vector of truth candidates of the jets.
     * @param el_charge Vector of charges of the electrons.
     * @param mu_charge Vector of charges of the muons.
     * @param n_electrons Number of electrons.
     * @param n_muons Number of muons.
     * @return ROOT::VecOps::RVec<float> Vector with the Chi2 value for the best pairing.
    */
    // Chi2-based pairing that searches jet->lepton assignments and minimizes
    // sum( (m_lb - mean)^2 / sigma^2 ). Defaults for mean/sigma are in MeV.
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
    );
        /**
     * @brief Create a TLV for a parton from the four momentum branches.
     * @return ROOT::RDF::RNode The node with the new TLV branch.
     * @param mainFrame The main frame.
     * @param mainNode The current main node.
     * @param partonBaseName The base name of the parton branches.
     * @param newBranchName The name of the new TLV branch.
     */
    ROOT::RDF::RNode createNeutrinoTLV(MainFrame* mainFrame ,ROOT::RDF::RNode mainNode, const std::string& partonBaseName, const std::string& newBranchName);


    struct TruthNuInfo {
    ROOT::Math::PtEtaPhiMVector nu_from_t;
    ROOT::Math::PtEtaPhiMVector nu_from_tbar;
    int pdg_from_t;
    int pdg_from_tbar;
};
TruthNuInfo ExtractTruthNeutrinos(
    float t_w1_pt,    float t_w1_eta,    float t_w1_phi,    float t_w1_m,    int t_w1_pdg,
    float t_w2_pt,    float t_w2_eta,    float t_w2_phi,    float t_w2_m,    int t_w2_pdg,
    float tb_w1_pt,   float tb_w1_eta,   float tb_w1_phi,   float tb_w1_m,   int tb_w1_pdg,
    float tb_w2_pt,   float tb_w2_eta,   float tb_w2_phi,   float tb_w2_m,   int tb_w2_pdg
);

struct NuFlavourInfo {
    ROOT::Math::PtEtaPhiMVector nu_e;
    ROOT::Math::PtEtaPhiMVector nu_mu;

    int pdg_e;   // ±12
    int pdg_mu;  // ±14
};

NuFlavourInfo ExtractNeutrinoFlavours(
    float t_w1_pt,    float t_w1_eta,    float t_w1_phi,    float t_w1_m,    int t_w1_pdg,
    float t_w2_pt,    float t_w2_eta,    float t_w2_phi,    float t_w2_m,    int t_w2_pdg,
    float tb_w1_pt,   float tb_w1_eta,   float tb_w1_phi,   float tb_w1_m,   int tb_w1_pdg,
    float tb_w2_pt,   float tb_w2_eta,   float tb_w2_phi,   float tb_w2_m,   int tb_w2_pdg
);

}
namespace ZllCalibrationVars{
       std::string regionName(
        const ROOT::VecOps::RVec<TLV>& sorted_elTLV,
        const ROOT::VecOps::RVec<TLV>& sorted_muTLV,
        const ROOT::VecOps::RVec<float>& sorted_elCharge,
        const ROOT::VecOps::RVec<float>& sorted_muCharge);
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

    template <typename T>
    bool allElementsUnique(const std::vector<T>& vec){
        std::set<T> uniqueElements(vec.begin(), vec.end());
        return uniqueElements.size() == vec.size();
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

namespace HyPERPerformanceSingleLepton{
    /**
     * @brief Form the Reco indices of TopLep and TopHad. 
     * The indices will be [blep, lepton , met = 0 always, bhad, q1, q2]
     * @param jetTruthIndices The vector with the truth indices of the jets.
     * @param eTruthIndices The vector with the truth indices of the electrons.
     * @param muTruthIndices The vector with the truth indices of the muons.
     * @param nMatchedJets Vector with number of matched jets per truth object.
     * @param nMatchedElectrons Vector with number of matched electrons per truth object.
     * @param nMatchedMuons Vector with number of matched muons per truth object.
     * @param nLeptons The number of truth leptons in the event. Always needs to be 1 for single lepton.
     * @return std::vector<int> The vector with the locations of the truth objects.
     */
    std::vector<int> getRecoIndicesTopLep(const std::vector<int>& jetTruthIndices,
                                        const std::vector<int>& eTruthIndices,
                                        const std::vector<int>& muTruthIndices,
                                        const std::vector<int>& nMatchedJets,
                                        const std::vector<int>& nMatchedElectrons,
                                        const std::vector<int>& nMatchedMuons,
                                        int nLeptons);

    /**
     * @brief Form the Reco indices of TopHad.
     * The indices will be [b, q1, q2]
     * @param jetTruthIndices The vector with the truth indices of the jets.
     * @param eTruthIndices The vector with the truth indices of the electrons.
     * @param muTruthIndices The vector with the truth indices of the muons.
     * @param nMatchedJets Vector with number of matched jets per truth object.
     * @param nMatchedElectrons Vector with number of matched electrons per truth object.
     * @param nMatchedMuons Vector with number of matched muons per truth object.
     * @param nLeptons The number of truth leptons in the event. Always needs to be 1 for single lepton.
     * @return std::vector<int> The vector with the locations of the truth objects.
     */
    std::vector<int> getRecoIndicesTopHad(const std::vector<int>& jetTruthIndices,
                                        const std::vector<int>& eTruthIndices,
                                        const std::vector<int>& muTruthIndices,
                                        const std::vector<int>& nMatchedJets,
                                        const std::vector<int>& nMatchedElectrons,
                                        const std::vector<int>& nMatchedMuons,
                                        int nLeptons);


    /**
     * @brief Functions to check whether the TopLep and TopHad are reconstrutable.
     * @param TruthIndices
     */
    bool isTopLepReconstructable(const std::vector<int>& TruthIndices);
    bool isTopHadReconstructable(const std::vector<int>& TruthIndices);
    
    /**
     * @brief Get the TopLep lepton flavour. 
     * The default/invalid value is -1. Electron = 2, Muon = 3.
     * @param eTruthIndices The vector with the truth indices of the electrons.
     * @param muTruthIndices The vector with the truth indices of the muons.
     * @param nMatchedElectrons Vector with number of matched electrons per truth object.
     * @param nMatchedMuons Vector with number of matched muons per truth object.
     * @param nLeptons The number of truth leptons in the event. Always needs to be 1 for single lepton.
    */
   int getTopLepLeptonFlavour(const std::vector<int>& eTruthIndices,
                                        const std::vector<int>& muTruthIndices,
                                        const std::vector<int>& nMatchedElectrons,
                                        const std::vector<int>& nMatchedMuons,
                                        int nLeptons);

    /**
     * @brief Check if the HyPER result for the TopLep makes sense.
     * @param HyPERTopLepIndices The HyPER indices of the TopLep.
     * @param HyPERTopLepIDs The HyPER IDs of the TopLep indices.
     */
    bool HyPERTopLepConverged(const std::vector<long>& HyPERTopLepIndices,
        const std::vector<long>& HyPERTopLepIDs);

    /**
     * @brief Check if the HyPER result for the TopHad makes sense.
     * @param HyPERTopHadIndices The HyPER indices of the TopHad.
     * @param HyPERTopHadIDs The HyPER IDs of the TopHad indices.
     */
    bool HyPERTopHadConverged(const std::vector<long>& HyPERTopHadIndices,
        const std::vector<long>& HyPERTopHadIDs);

    /**
     * @brief Check if we have two valid solutions for the tops.
     * @param topLepConverged
     * @param topHadConverged
     * @param HyPERTopLepIndices The HyPER indices of the TopLep.
     * @param HyPERTopLepIDs The HyPER IDs of the TopLep indices.
     * @param HyPERTopHadIndices The HyPER indices of the TopHad.
     * @param HyPERTopHadIDs The HyPER IDs of the TopHad indices.
     */
    bool HyPERConvergedTwoTops(const bool topLepConverged,
        const bool topHadConverged,
        const std::vector<long>& HyPERTopLepIndices,
        const std::vector<long>& HyPERTopLepIDs,
        const std::vector<long>& HyPERTopHadIndices,
        const std::vector<long>& HyPERTopHadIDs,
        const std::vector<long>& HyPERWHadIndices);
    
    /** 
     * @brief Check if the WHad is correctly reconstructed.
     * @param TopHadTruthRecoIndices The indices of the TopHad partciles from the truth record. Constructed via getRecoIndicesTopHad.
     * @param HyPERWHadIndices The HyPER indices of the WHad.
    */
    bool correctlyReconstructedWHad(std::vector<int> TopHadTruthRecoIndices,
                                     const std::vector<long>& HyPERWHadIndices);

    /** 
     * @brief Check if the WHad is correctly reconstructed.
     * @param TopLepTruthRecoIndices The indices of the TopLep particles from the truth record. Constructed via getRecoIndicesTopLep.
     * @param HyPERWLepIndices The HyPER indices of the WLep.
    */
    bool correctlyReconstructedWLep(std::vector<int> TopLepTruthRecoIndices,
                                     const std::vector<long>& HyPERWLepIndices,
                                     int truthLeptonFlavour,
                                     const std::vector<long>& HyPERTopLepIDs);

    /** 
     * @brief Check if the WHad is correctly reconstructed.
     * @param TopHadTruthRecoIndices The indices of the TopHad partciles from the truth record. Constructed via getRecoIndicesTopHad.
     * @param HyPERWHadIndices The HyPER indices of the WHad.
    */
   bool correctlyReconstructedWHad(std::vector<int> TopHadTruthRecoIndices,
    const std::vector<long>& HyPERWHadIndices);
    
    // These function return -1 if no valid solution is found.
    int getHyPERBHadIndex(const std::vector<long>& HyPERTopHadIndices,
                          const std::vector<long>& HyPERTopHadIDs,
                          const std::vector<long>& HyPERWHadIndices);

    int getHyPERBLepIndex(const std::vector<long>& HyPERTopLepIndices,
                          const std::vector<long>& HyPERTopLepIDs);

    /**
     * @brief Check if the TopLep was reconstructed correctly by HyPER.
     * @param correctWLep True if the WLep was reconstructed correctly.
     * @param HyPERBLepIndex The HyPER index of the BLep.
     * @param TopLepTruthRecoIndices The indices of the TopLep particles from the truth record. Constructed via getRecoIndicesTopLep.
     */
    bool correctlyReconstructedTopLep(const bool correctWLep,
        int HyPERBLepIndex,
        const std::vector<int>& TopLepTruthRecoIndices);

    /**
     * @brief Check if the TopHad was reconstructed correctly by HyPER.
     * @param correctWHad True if the WHad was reconstructed correctly.
     * @param HyPERBHadIndex The HyPER indicex of the BHad.
     * @param TopHadTruthRecoIndices The indices of the TopHad particles from the truth record. Constructed via getRecoIndicesTopHad.
     */
    bool correctlyReconstructedTopHad(const bool correctWHad,
        int HyPERBHadIndex,
        const std::vector<int>& TopHadTruthRecoIndices);
    
    /**
     * @brief These are the KLF versions of the functions above. 
     * The idea is to build the same inputs and use the HyPER logic.
     */
    bool correctlyReconstructedTopLep_KLF(int KLFBLepIndex,
                                          const std::vector<int>& TopLepTruthRecoIndices);

    bool correctlyReconstructedWHad_KLF(std::vector<int> TopHadTruthRecoIndices,
                                        int KLW1Index,
                                        int KLW2Index);

    bool correctlyReconstructedTopHad_KLF(const bool correctWHad,
        int KLFBHadIndex,
        const std::vector<int>& TopHadTruthRecoIndices);
}