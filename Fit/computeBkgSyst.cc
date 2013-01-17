#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>
#include <unistd.h>

#include <tclap/CmdLine.h>
#include <json/json.h>

#include "Utils.h"

std::string buildConfigFileName(const std::string& fctName) {
  return "fit_pdf_" + fctName + ".json";
}

std::vector<std::string> BKG_FUNCTIONS;
std::string base_path = "";

void process(std::vector<int>& masses, bool onlyMuons, int btag, const std::string& file, bool singleFile, bool fixBackground) {
  std::vector<pid_t> children; // For fork()

  for (std::vector<int>::iterator it = masses.begin(); it != masses.end(); ++it) {

    pid_t pid = fork();
    if (pid == 0) {
      for (std::vector<std::string>::const_iterator param = BKG_FUNCTIONS.begin(); param != BKG_FUNCTIONS.end(); ++param) {
        std::stringstream ss;
        ss << *it;

        pid_t child = fork();
        if (child == 0) {
          std::string fileName = buildConfigFileName(*param);
          char** params = getSystCLParameters(ss.str(), file, singleFile, fixBackground, onlyMuons, btag, "--config-file", fileName.c_str(), NULL);
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

  std::cout << std::endl << std::endl;
}

void fillParams() {
  BKG_FUNCTIONS.push_back("falt");
  //BKG_FUNCTIONS.push_back("faltc");
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

  root[getAnalysisUUID()][strMass][btagStr]["background_pdf"] = syst;

  FILE* fd = fopen((base_path + "/systematics.json").c_str(), "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

double computeSystValue(double sigma_ref, double sigma) {
  return fabs(fabs(sigma) - fabs(sigma_ref));
}

void computeSyst(std::vector<int> masses, int btag) {

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

    double totalSyst = 0;
    int n = 0;
    Json::Value bkgNode = root[strMass][btagStr]["background"];

    Json::Value::Members members = bkgNode.getMemberNames();
    for (Json::Value::Members::iterator member = members.begin(); member != members.end(); ++member) {
      if (std::find_if(BKG_FUNCTIONS.begin(), BKG_FUNCTIONS.end(), [&member](std::string& s) {
            return buildConfigFileName(s) == *member;
          }) == BKG_FUNCTIONS.end())
          continue;

      double sigma = bkgNode[*member]["sigma"].asDouble();
      std::cout << "sigma for " << *member  << ": " << sigma << std::endl;
      totalSyst += computeSystValue(sigma_ref, sigma);
      n++;
    }

    totalSyst /= n;

    saveSystematic(*it, btag, totalSyst / fabs(sigma_ref));
    std::cout << "Background PDF syst for MZ' = " << *it << " : " << totalSyst << " pb; " << totalSyst / fabs(sigma_ref) * 100 << "%" << std::endl;
  }
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Compute bkg PDF systematic", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::SwitchArg muonsOnlyArg("", "muons-only", "Compute sigmaref using only semi-mu data", cmd);
    TCLAP::SwitchArg extractArg("", "dont-extract", "Don't run fitMtt. Only compute systematic with previous results", cmd);
    TCLAP::MultiArg<int> massArg("m", "mass", "Zprime mass", false, "integer", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "integer", cmd);

    TCLAP::SwitchArg fixBackgroundArg("", "fix-background", "Fix background when fitting", cmd);

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
      process(masses, muonsOnlyArg.getValue(), btagArg.getValue(), inputFileArg.isSet() ? inputFileArg.getValue() : inputListArg.getValue(), inputFileArg.isSet(), fixBackgroundArg.getValue());

    computeSyst(masses, btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
