#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>

#include <tclap/CmdLine.h>
#include <json/json.h>

#include "Utils.h"

std::vector<std::string> CB_PARAMS;

void process(std::vector<int>& masses, bool muonsOnly, int btag) {
  std::vector<pid_t> children; // For fork()

  for (std::vector<int>::iterator it = masses.begin(); it != masses.end(); ++it) {

    pid_t pid = fork();
    if (pid == 0) {
      for (std::vector<std::string>::const_iterator param = CB_PARAMS.begin(); param != CB_PARAMS.end(); ++param) {
        std::stringstream ss;
        ss << *it;

        const char* signs[] = {
          "up",
          "down"
        };
        for (int i = 0; i < 2; i++) {
          pid_t child = fork();
          if (child == 0) {
            char** params = getSystCLParameters(ss.str(), muonsOnly, btag, "--systCB", (*param).c_str(), "--syst-sign", signs[i], NULL);
            execv("./fitMtt", params);
            exit(0);
          } else {
            waitpid(child, NULL, 0);
          }
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

void fillParams(bool muonsOnly) {
  CB_PARAMS.push_back("muon_alpha");
  CB_PARAMS.push_back("muon_sigma");
  CB_PARAMS.push_back("muon_mean");
  if (! muonsOnly) {
    CB_PARAMS.push_back("electron_alpha");
    CB_PARAMS.push_back("electron_sigma");
    CB_PARAMS.push_back("electron_mean");
  }
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

  root[strMass][btagStr]["signal_pdf"] = syst;

  FILE* fd = fopen("systematics.json", "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

double computeSystValue(double sigma, double up, double down) {
  return (fabs(up - sigma) + fabs(down - sigma)) / (2. * fabs(sigma));
}

void computeSyst(std::vector<int> masses, int btag) {
  std::cout << std::endl << std::endl;

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

    double up = 0, down = 0, totalSyst = 0;
    Json::Value signalNode = root[strMass][btagStr]["signal"];
    Json::Value::Members members = signalNode.getMemberNames();
    for (Json::Value::Members::const_iterator member = members.begin(); member != members.end(); ++member) {

      // Check if the current node is in our array of variable to vary
      if (std::find(CB_PARAMS.begin(), CB_PARAMS.end(), *member) == CB_PARAMS.end()) {
        std::cout << "Node " << *member << " is in json file, but not in CB_PARAMS array. Removing it from syst. calculation" << std::endl;
        continue;
      }

      up = signalNode[*member]["sigma"]["up"].asDouble();
      down = signalNode[*member]["sigma"]["down"].asDouble();
      double syst = computeSystValue(sigma_ref, up, down);
      totalSyst += syst * syst;
      std::cout << *member << " : sigma up = " << up << "; sigma down = " << down << "; syst = " << syst << std::endl;
    }
    double syst = sqrt(totalSyst);
    saveSystematic(*it, btag, syst);
    std::cout << "Signal PDF syst for MZ' = " << *it << " : " << syst << std::endl;
  }
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Compute CrystalBall PDF systematic", ' ', "0.1");

    TCLAP::SwitchArg muonsOnlyArg("", "muons-only", "Compute sigmaref using only semi-mu data", cmd);
    TCLAP::SwitchArg extractArg("", "dont-extract", "Don't run fitMtt for each parameters. Only compute systematic with previous results", cmd);
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

    fillParams(muonsOnlyArg.getValue());

    if (! extractArg.getValue())
      process(masses, muonsOnlyArg.getValue(), btagArg.getValue());

    computeSyst(masses, btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
