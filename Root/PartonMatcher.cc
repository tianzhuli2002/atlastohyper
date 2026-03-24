#include "PartonMatcher/PartonMatcher.h"
#include "PartonMatcher/TtbarAllHadFrame.h"
#include "PartonMatcher/TtbarSingleLeptonFrame.h"
#include "PartonMatcher/TtbarTopRecoFrame.h"
#include "PartonMatcher/TtbarDileptonFrame.h"
#include "PartonMatcher/TtbarDileptonFramevyper.h"
#include "PartonMatcher/TtbarDileptonFrameReco.h"
#include "PartonMatcher/ZllCalibrationFrame.h"
void PartonMatcher::classFactory(){
    const std::string className = m_config->customOptions().getOption<std::string>("class_name","");
    LOG(INFO) << "Class name: " << className << std::endl;

    if (className == "") {
      LOG(ERROR) << "No class name provided";
      throw std::invalid_argument("");
    }

    if (className == "TtbarAllHadFrame") {
      m_channelFrame = std::make_unique<TtbarAllHadFrame>(this);
    } else if (className == "TtbarSingleLeptonFrame") {
      m_channelFrame = std::make_unique<TtbarSingleLeptonFrame>(this);
    } else if (className == "TtbarTopRecoFrame") {
      m_channelFrame = std::make_unique<TtbarTopRecoFrame>(this);
    } else if (className == "TtbarDileptonFrame") {
      m_channelFrame = std::make_unique<TtbarDileptonFrame>(this);
    } else if (className == "TtbarDileptonFrameReco") {
      m_channelFrame = std::make_unique<TtbarDileptonFrameReco>(this);
    } else if (className == "TtbarDileptonFramevyper") {
      m_channelFrame = std::make_unique<TtbarDileptonFramevyper>(this);
    } else if (className == "ZllCalibrationFrame") {
      m_channelFrame = std::make_unique<ZllCalibrationFrame>(this);
    } else {
      LOG(ERROR) << "Invalid class name provided";
      throw std::invalid_argument("");
    }
}