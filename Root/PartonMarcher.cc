#include "PartonMatcher/PartonMatcher.h"
#include "PartonMatcher/TtbarAllHadFrame.h"
#include "PartonMatcher/TtbarSingleLeptonFrame.h"
#include "PartonMatcher/TtbarTopRecoFrame.h"

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
    } else {
      LOG(ERROR) << "Invalid class name provided";
      throw std::invalid_argument("");
    }
}