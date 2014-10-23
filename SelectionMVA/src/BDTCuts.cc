#include <BDTCuts.h>

#include <yaml-cpp/yaml.h>

BDTCuts::BDTCuts(const std::string& cutsFile/* = "../SelectionMVA/bdt_cuts.yml"*/):
  m_file(cutsFile) {

    parse();
}

void BDTCuts::parse() {
  YAML::Node f = YAML::LoadFile(m_file);

  for (YAML::const_iterator it = f.begin(); it != f.end(); ++it) {
    std::string category = it->first.as<std::string>();

    YAML::Node cutNode = it->second;
    for (std::size_t i = 0; i < cutNode.size(); i++) {
      int btagCategory = cutNode[i]["b-tag"].as<int>();
      float cut = cutNode[i]["cut"].as<float>();
      std::string weights = cutNode[i]["weights"].as<std::string>();

      m_cuts[std::make_pair(stringToBDTType(category), btagCategory)].cut = cut;
      m_cuts[std::make_pair(stringToBDTType(category), btagCategory)].weights = weights;
    }
  }
}

float BDTCuts::getCut(BDTType bdtType, int btagCategory) {
  return m_cuts[std::make_pair(bdtType, btagCategory)].cut;
}

std::string BDTCuts::getWeights(BDTType bdtType, int btagCategory) {
  return m_cuts[std::make_pair(bdtType, btagCategory)].weights;
}

std::string BDTTypeToString(BDTType type) {

  switch (type) {
    case BDTType::BACKGROUND:
      return "background";

    case BDTType::SIGNAL:
      return "signal";
  }

  return "invalid_type";
}

BDTType stringToBDTType(const std::string& type) {

  if (type == "background")
    return BDTType::BACKGROUND;
  else if (type == "signal")
    return BDTType::SIGNAL;

  throw std::invalid_argument("type is not valid");
}

