#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>

#include <tclap/CmdLine.h>
#include <json/json.h>

#include "Utils.h"

void loadObservedLimits(const std::vector<int>& masses, int btag, std::vector<double>& obs) {
  Json::Reader reader;
  Json::Value root;
  std::ifstream file("likelihood_scan.json");
  bool success = reader.parse(file, root);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "likelihood_scan.json" << "'. Exiting." << std::endl;
    std::cerr << "Please run the likelihood scan using 'runLikelihoodScan'" << std::endl;
    exit(1);
  }

  root = root[getAnalysisUUID()];

  for (std::vector<int>::const_iterator it = masses.begin(); it != masses.end(); ++it) {
    std::stringstream ss;
    ss << (int) *it;
    std::string strMass = ss.str();

    ss.clear(); ss.str(std::string());
    ss << btag;
    std::string btagStr = ss.str();

    if (! root.isMember(strMass)) {
      std::cerr << "ERROR: mass '" << *it << "' not found in JSON file. Exiting." << std::endl;
      exit(1);
    }

    double observedLimit = root[strMass][btagStr]["scan_wsyst_cut_limit"].asDouble();
    obs.push_back(observedLimit);
  }
}

void loadExpectedLimits(const std::vector<int>& masses, int btag, std::vector<double>& exp, std::vector<double>& error_h_95, std::vector<double>& error_l_95, std::vector<double>& error_h_68, std::vector<double>& error_l_68) {
  Json::Reader reader;
  Json::Value root;
  std::ifstream file("expected_limits.json");
  bool success = reader.parse(file, root);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "expected_limits.json" << "'. Exiting." << std::endl;
    std::cerr << "Please run 'treatToyStuff'" << std::endl;
    exit(1);
  }

  root = root[getAnalysisUUID()];

  for (std::vector<int>::const_iterator it = masses.begin(); it != masses.end(); ++it) {
    std::stringstream ss;
    ss << (int) *it;
    std::string strMass = ss.str();

    ss.clear(); ss.str(std::string());
    ss << btag;
    std::string btagStr = ss.str();

    if (! root.isMember(strMass)) {
      std::cerr << "ERROR: mass '" << *it << "' not found in JSON file. Exiting." << std::endl;
      exit(1);
    }

    double expectedLimit = root[strMass][btagStr]["median"].asDouble();
    double eyl68 = root[strMass][btagStr]["widthM68"].asDouble();
    double eyh68 = root[strMass][btagStr]["widthP68"].asDouble();
    double eyl95 = root[strMass][btagStr]["widthM95"].asDouble();
    double eyh95 = root[strMass][btagStr]["widthP95"].asDouble();

    exp.push_back(expectedLimit);
    error_h_95.push_back(eyh95);
    error_l_95.push_back(eyl95);
    error_h_68.push_back(eyh68);
    error_l_68.push_back(eyl68);
  }
}

void extractLimits(std::vector<int> masses, bool showObserved, bool showExpected, int btag) {

  std::vector<double> observed; // Observed limits
  std::vector<double> expected; // Expected limits

  std::vector<double> error_high_95; // 95%
  std::vector<double> error_low_95; // 95%

  std::vector<double> error_high_68;
  std::vector<double> error_low_68;

  if (showObserved)
    loadObservedLimits(masses, btag, observed);

  if (showExpected)
    loadExpectedLimits(masses, btag, expected, error_high_95, error_low_95, error_high_68, error_low_68);

  int i = 0;
  for (std::vector<int>::iterator it = masses.begin(); it != masses.end(); ++it, i++) {
    std::cout << "#### M = " << *it << " GeV ####" << std::endl;

    if (showObserved) {
      std::cout << "Observed limit = " << observed[i] << " pb" << std::endl;
    }

    if (showExpected) {
      std::cout << "Expected limit = " << expected[i] << " pb" << std::endl;
      std::cout << "68% = -" << error_low_68[i] << "\t+" << error_high_68[i] << std::endl;
      std::cout << "95% = -" << error_low_95[i] << "\t+" << error_high_95[i] << std::endl;
    }

    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Extract observed and expected limits", ' ', "0.1");

    TCLAP::SwitchArg observedArg("", "observed", "Extract observed limit", cmd);
    TCLAP::SwitchArg expectedArg("", "expected", "Extract expected limit", cmd);

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

    extractLimits(masses, observedArg.getValue(), expectedArg.getValue(), btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
