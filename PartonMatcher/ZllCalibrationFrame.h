#pragma once
#include "PartonMatcher/PartonMatcher.h"

class ZllCalibrationFrame : public ChannelFrame {
public:

  explicit ZllCalibrationFrame(PartonMatcher* partonMatcherPtr) : ChannelFrame(partonMatcherPtr) {}

  virtual ~ZllCalibrationFrame() = default;

  virtual ROOT::RDF::RNode defineVariables(ROOT::RDF::RNode mainNode,
                                           const std::shared_ptr<Sample>& sample,
                                           const UniqueSampleID& id);
  
  virtual ROOT::RDF::RNode defineVariablesNtuple(ROOT::RDF::RNode mainNode,
                                                 const std::shared_ptr<Sample>& sample,
                                                 const UniqueSampleID& id);

  virtual ROOT::RDF::RNode defineVariablesTruth(ROOT::RDF::RNode node,
                                                const std::string& treeName,
                                                const std::shared_ptr<Sample>& sample,
                                                const UniqueSampleID& sampleID);

    
  virtual ROOT::RDF::RNode defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                      const std::string& treeName,
                                                      const std::shared_ptr<Sample>& sample,
                                                      const UniqueSampleID& sampleID);

};
