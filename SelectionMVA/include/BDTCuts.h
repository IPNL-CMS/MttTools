#pragma once

#include <map>
#include <stdexcept>
#include <string>

enum class BDTType {
  BACKGROUND,
  SIGNAL
};

std::string BDTTypeToString(BDTType type);
BDTType stringToBDTType(const std::string& type);

class BDTCuts {
  public:
    BDTCuts(const std::string& cutsFile = "../SelectionMVA/bdt_cuts.yml");

    float getCut(BDTType bdtType, int btagCategory);

  private:
    void parse();

    std::string m_file;
    std::map<std::pair<BDTType, int>, float> m_cuts;
};

