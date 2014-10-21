#include <MVAReader.h>

MVAReader::MVAReader(const std::vector<std::string>& inputFiles):
  m_inputFiles(inputFiles) {
}

void MVAReader::initMVA(const std::string& weightFile) {
  initTrees();
  createReader(weightFile);
}

void MVAReader::createReader(const std::string& weightFile) {
  m_reader.reset(new TMVA::Reader());

  setupVariables();

  m_reader->BookMVA("mva", weightFile.c_str());
}

float MVAReader::evaluate(uint64_t entry) {

  for (auto& chain: m_chains) {
    chain->GetEntry(entry);
  }

  computeVariables();

  return m_reader->EvaluateMVA("mva");
}
