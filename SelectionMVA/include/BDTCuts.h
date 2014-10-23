#pragma once

#include <map>
#include <stdexcept>
#include <string>

enum class BDTType {
  BACKGROUND,
  SIGNAL
};

struct BDTInfo {
  float cut;
  std::string weights;
};

std::string BDTTypeToString(BDTType type);
BDTType stringToBDTType(const std::string& type);

class BDTCuts {
  public:
    BDTCuts(const std::string& cutsFile = "../SelectionMVA/bdt_cuts.yml");

    float getCut(BDTType bdtType, int btagCategory);
    std::string getWeights(BDTType bdtType, int btagCategory);

  private:
    void parse();

    std::string m_file;
    std::map<std::pair<BDTType, int>, BDTInfo> m_cuts;
};

