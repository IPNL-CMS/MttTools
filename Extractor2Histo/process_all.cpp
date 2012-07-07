#include <TSystem.h>
#include "Extractor2Histos.h"

#include <fstream>
#include <sstream>
#include <vector>

#include "tclap/CmdLine.h"

struct MCData {
  std::string inputFile;
  std::string outputFile;
  std::string name;
  bool semimu;
};

std::vector<MCData> parseMCInputs(const std::string& file, bool mc) {

  std::vector<MCData> mcDatas;

  fstream f(file.c_str());
  std::string line;
  MCData data;

  while (std::getline(f, line)) {
    if (!line.length() || line[0] == '#')
      continue;

    if (line == "end")
      break;

    std::stringstream ss(line);
    if (mc) {
      ss >> data.inputFile >> data.outputFile >> data.name;
    } else {
      ss >> data.inputFile >> data.outputFile >> data.semimu;
    }

    //FIXME?
    data.inputFile = "inputs/" + data.inputFile;

    mcDatas.push_back(data);
  }

  return mcDatas;
}

int main(int argc, char** argv)
{
  // MC
  // semimu

  try {

    TCLAP::CmdLine cmd("Convert Extractor output to histograms", ' ', "0.1");

    TCLAP::SwitchArg semimuArg("", "semimu", "Process semimu", cmd, false);
    TCLAP::SwitchArg semieArg("", "semie", "Process semie", cmd, false);
    TCLAP::SwitchArg dataArg("", "data", "Process data", cmd, false);

    cmd.parse(argc, argv);

    Extractor2Histos *a = NULL;

    std::vector<MCData> mcDatas;

    if (semimuArg.getValue()) {

      mcDatas = parseMCInputs("input_mu.list", true);

      for (std::vector<MCData>::iterator it = mcDatas.begin(); it != mcDatas.end(); ++it) {

        a = new Extractor2Histos(it->inputFile, it->outputFile, it->name, true, true);
        a->Loop();
        delete a;

      }
    }

    if (semieArg.getValue()) {

      mcDatas = parseMCInputs("input_e.list", true);

      for (std::vector<MCData>::iterator it = mcDatas.begin(); it != mcDatas.end(); ++it) {

        a = new Extractor2Histos(it->inputFile, it->outputFile, it->name, false, true);
        a->Loop();
        delete a;

      }
    }

    // Data
    if (dataArg.getValue()) {

      mcDatas = parseMCInputs("input_data.list", false);

      for (std::vector<MCData>::iterator it = mcDatas.begin(); it != mcDatas.end(); ++it) {

        a = new Extractor2Histos(it->inputFile, it->outputFile, "", it->semimu, false);
        a->Loop();
        delete a;

      }
    }

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}

