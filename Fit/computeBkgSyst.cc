#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>

#include <tclap/CmdLine.h>
#include <json/json.h>

#include "Utils.h"

std::vector<std::string> BKG_FUNCTIONS;

void process(std::vector<int>& masses, bool onlyMuons, int btag) {
  std::vector<pid_t> children; // For fork()

  for (std::vector<int>::iterator it = masses.begin(); it != masses.end(); ++it) {

    pid_t pid = fork();
    if (pid == 0) {
      for (std::vector<std::string>::const_iterator param = BKG_FUNCTIONS.begin(); param != BKG_FUNCTIONS.end(); ++param) {
        std::stringstream ss;
        ss << *it;

        pid_t child = fork();
        if (child == 0) {
          char** params = getSystCLParameters(ss.str(), onlyMuons, btag, "--name", (*param).c_str(), NULL);
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
  //BKG_FUNCTIONS.push_back("faltC");
}

void saveSystematic(int mass, int btag, double syst) {

  FILE* lock = fopen("systematics.lock", "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  std::ifstream file("systematics.json");
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  root[strMass][btagStr]["background_pdf"] = syst;

  FILE* fd = fopen("systematics.json", "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

double computeSystValue(double sigma_ref, double sigma) {
  return fabs(sigma - sigma_ref) / fabs(sigma_ref);
}

void computeSyst(std::vector<int> masses, int btag) {

  Json::Reader reader;
  Json::Value root;
  std::ifstream file("systematics_parameters.json");
  bool success = reader.parse(file, root);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "systematics_parameters.json" << "'. Exiting." << std::endl;
    exit(1);
  }

  Json::Value refRoot;
  file.open("sigma_reference.json");
  success = reader.parse(file, refRoot);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "sigma_reference.json" << "'. Exiting." << std::endl;
    std::cerr << "You may need to run fitMtt first." << std::endl;
    exit(1);
  }

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
      double sigma = bkgNode[*member]["sigma"].asDouble();
      totalSyst += computeSystValue(sigma_ref, sigma);
      n++;
    }

    totalSyst /= n;

    saveSystematic(*it, btag, totalSyst);
    std::cout << "Background PDF syst for MZ' = " << *it << " : " << totalSyst << std::endl;
  }
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Compute bkg PDF systematic", ' ', "0.1");

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

    fillParams();

    if (! extractArg.getValue())
      process(masses, muonsOnlyArg.getValue(), btagArg.getValue());

    computeSyst(masses, btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
