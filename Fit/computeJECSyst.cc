#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>
#include <unistd.h>

#include <tclap/CmdLine.h>
#include <json/json.h>

#include "Utils.h"

std::vector<std::string> JEC_PARAMS;
std::string base_path = "";

void process(std::vector<int>& masses, bool muonsOnly, int btag, const std::string& file, bool singleFile) {
  std::vector<pid_t> children; // For fork()

  for (std::vector<int>::iterator it = masses.begin(); it != masses.end(); ++it) {

    pid_t pid = fork();
    if (pid == 0) {
      for (std::vector<std::string>::const_iterator param = JEC_PARAMS.begin(); param != JEC_PARAMS.end(); ++param) {
        std::stringstream ss;
        ss << *it;

        pid_t child = fork();
        if (child == 0) {
          char** params = getSystCLParameters(ss.str(), file, singleFile, muonsOnly, btag, "--syst", (*param).c_str(), NULL);
          execv("./fitMtt", params);
          exit(0);
        } else {
          waitpid(child, NULL, 0);
        }
      }
      exit(0);
    } else {
      children.push_back(pid);
    }

  }

  // wait for all children
  for (std::vector<pid_t>::iterator it = children.begin(); it != children.end(); ++it) {
    waitpid(*it, NULL, 0);
  }
}

void fillParams() {
  JEC_PARAMS.push_back("JECup");
  JEC_PARAMS.push_back("JECdown");
}

void saveSystematic(int mass, int btag, double syst) {

  FILE* lock = fopen((base_path + "/systematics.lock").c_str(), "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  std::ifstream file((base_path + "/systematics.json").c_str());
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  root[getAnalysisUUID()][strMass][btagStr]["jec"] = syst;

  FILE* fd = fopen((base_path + "/systematics.json").c_str(), "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

double computeSystValue(double sigma, double up, double down) {
  double systup = fabs(up - sigma) / fabs(sigma);
  double systdown = fabs(down - sigma) / fabs(sigma);
  return (systup + systdown) / 2.;
}

void computeSyst(std::vector<int> masses, int btag) {
  std::cout << std::endl << std::endl;

  Json::Reader reader;
  Json::Value root;
  std::ifstream file((base_path + "/systematics_parameters.json").c_str());
  bool success = reader.parse(file, root);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "systematics_parameters.json" << "'. Exiting." << std::endl;
    exit(1);
  }

  root = root[getAnalysisUUID()];

  Json::Value refRoot;
  file.open((base_path + "/sigma_reference.json").c_str());
  success = reader.parse(file, refRoot);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "sigma_reference.json" << "'. Exiting." << std::endl;
    std::cerr << "You may need to run fitMtt first." << std::endl;
    exit(1);
  }

  refRoot = refRoot[getAnalysisUUID()];

  for (std::vector<int>::iterator it = masses.begin(); it != masses.end(); ++it) {
    std::stringstream ss;
    ss << *it;
    std::string strMass = ss.str();

    ss.clear(); ss.str(std::string());
    ss << btag;
    std::string btagStr = ss.str();

    if (! root.isMember(strMass) || ! refRoot.isMember(strMass)) {
      std::cerr << "ERROR: mass '" << *it << "' not found in JSON file. Exiting." << std::endl;
      exit(1);
    }

    std::cout << std::endl;
    std::cout << "#### M = " << *it << " GeV ###" << std::endl;
    std::cout << std::endl;

    double sigma_ref = refRoot[strMass][btagStr]["sigma"].asDouble();
    std::cout << "sigma_ref = " << sigma_ref << std::endl;

    double up = 0, down = 0;
    Json::Value jecNode = root[strMass][btagStr]["jec"];

    up = jecNode["JECup"]["sigma"].asDouble();
    down = jecNode["JECdown"]["sigma"].asDouble();
    double syst = computeSystValue(sigma_ref, up, down);
    saveSystematic(*it, btag, syst);
    std::cout << "JEC syst for MZ' = " << *it << " : " << syst << std::endl;
  }
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Compute JEC systematic", ' ', "0.1");
    
    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::SwitchArg muonsOnlyArg("", "muons-only", "Compute sigmaref using only semi-mu data", cmd);
    TCLAP::SwitchArg extractArg("", "dont-extract", "Don't run fitMtt. Only compute systematic with previous results", cmd);
    TCLAP::MultiArg<int> massArg("m", "mass", "Zprime mass", false, "integer", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "integer", cmd);

    cmd.parse(argc, argv);

    std::vector<int> masses = massArg.getValue();
    if (masses.size() == 0) {
      masses.push_back(750);
      masses.push_back(1000);
      masses.push_back(1250);
      masses.push_back(1500);
    }

    base_path = "analysis/" + getAnalysisUUID();

    fillParams();

    if (! extractArg.getValue())
      process(masses, muonsOnlyArg.getValue(), btagArg.getValue(), inputFileArg.isSet() ? inputFileArg.getValue() : inputListArg.getValue(), inputFileArg.isSet());

    computeSyst(masses, btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
