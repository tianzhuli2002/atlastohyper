#pragma once

#include "FastFrames/MainFrame.h"

#include "FastFrames/ConfigSetting.h"
#include "FastFrames/Sample.h"
#include "FastFrames/UniqueSampleID.h"

#include "Math/Vector4D.h"
#include "ROOT/RDataFrame.hxx"
#include "TClass.h"

#include <memory>

using TLV = ROOT::Math::PtEtaPhiEVector;

class UniqueSampleID;
class PartonMatcher;

class ChannelFrame {
  public:
  explicit ChannelFrame(PartonMatcher* partonMatcherPtr) : m_partonMatcher(partonMatcherPtr) {}

  virtual ~ChannelFrame() { m_partonMatcher = nullptr; }


  virtual ROOT::RDF::RNode defineVariables(ROOT::RDF::RNode mainNode,
                                           const std::shared_ptr<Sample>& sample,
                                           const UniqueSampleID& id) = 0;
  
  virtual ROOT::RDF::RNode defineVariablesNtuple(ROOT::RDF::RNode mainNode,
                                                 const std::shared_ptr<Sample>& sample,
                                                 const UniqueSampleID& id) = 0;

  virtual ROOT::RDF::RNode defineVariablesTruth(ROOT::RDF::RNode node,
                                                const std::string& truth,
                                                const std::shared_ptr<Sample>& sample,
                                                const UniqueSampleID& sampleID) = 0;
  
  virtual ROOT::RDF::RNode defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                      const std::string& treeName,
                                                      const std::shared_ptr<Sample>& sample,
                                                      const UniqueSampleID& sampleID) = 0;

  PartonMatcher* m_partonMatcher;
};

class PartonMatcher : public MainFrame {
public:

  explicit PartonMatcher() = default;

  virtual ~PartonMatcher() = default;

  virtual void init() override {
      MainFrame::init(); 
      classFactory(); 
  }

  virtual ROOT::RDF::RNode defineVariables(ROOT::RDF::RNode mainNode,
                                           const std::shared_ptr<Sample>& sample,
                                           const UniqueSampleID& id) override 
                                      {
                                        return m_channelFrame->defineVariables(mainNode, sample, id);
                                      }
  
  virtual ROOT::RDF::RNode defineVariablesNtuple(ROOT::RDF::RNode mainNode,
                                                 const std::shared_ptr<Sample>& sample,
                                                 const UniqueSampleID& id) override
                                                 {
                                                    return m_channelFrame->defineVariablesNtuple(mainNode, sample, id);
                                                 }

  virtual ROOT::RDF::RNode defineVariablesTruth(ROOT::RDF::RNode node,
                                                const std::string& truth,
                                                const std::shared_ptr<Sample>& sample,
                                                const UniqueSampleID& sampleID) override
                                                {
                                                  return m_channelFrame->defineVariablesTruth(node, truth, sample, sampleID);
                                                }
  
  virtual ROOT::RDF::RNode defineVariablesNtupleTruth(ROOT::RDF::RNode node,
                                                      const std::string& treeName,
                                                      const std::shared_ptr<Sample>& sample,
                                                      const UniqueSampleID& sampleID) override
                                                      {
                                                        return m_channelFrame->defineVariablesNtupleTruth(node, treeName, sample, sampleID);
                                                      }

  void classFactory();
    
private:
  std::unique_ptr<ChannelFrame> m_channelFrame;

  ClassDefOverride(PartonMatcher, 1);

};