#pragma once
// Inlcude the math library
#include<vector>
#include <string>
#include <algorithm>
#include "Math/Vector4D.h"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"

// Variables used in both channels
namespace Variables{
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
     * @brief Calculate the number of jets that pass two selections and have a pT greater than 20 GeV.
     * @return int The number of jets that pass the two selections.
     * @param ptVec The pT vector of the jets.
     * @param selected1 The first selection.
     * @param selected2 The second selection.
     */
    int nGoodJets(const std::vector<float>& ptVec,
                        const std::vector<char>& selected1,
                        const std::vector<char>& selected2){
        return DefineHelpers::numberOfObjects(ptVec,20000,selected1,selected2);
    }

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

}