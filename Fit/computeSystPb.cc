#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>

#include <tclap/CmdLine.h>
#include <json/json.h>

void computeSyst(std::vector<int> masses) {
  Json::Reader reader;
  Json::Value root;
  std::ifstream file("systematics.json");
  bool success = reader.parse(file, root);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "systematics.json" << "'. Exiting." << std::endl;
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

    if (! root.isMember(strMass) || ! refRoot.isMember(strMass)) {
      std::cerr << "ERROR: mass '" << *it << "' not found in JSON file. Exiting." << std::endl;
      exit(1);
    }

    std::cout << std::endl;
    std::cout << "#### M = " << *it << " GeV ###" << std::endl;
    std::cout << std::endl;

    double sigma_ref = refRoot[strMass]["sigma"].asDouble();
    std::cout << "sigma_ref = " << sigma_ref << std::endl;

    double bkg = 0, signal = 0, jec = 0;
    jec = root[strMass]["jec"].asDouble() * sigma_ref;
    signal = root[strMass]["signal_pdf"].asDouble() * sigma_ref;
    bkg = root[strMass]["background_pdf"].asDouble() * sigma_ref;

    double totalSyst = sqrt(jec * jec + signal * signal + bkg * bkg);

    std::cout << "syst_jec = " << jec << " | syst_bkg = " << bkg << " | syst_signal = " << signal << std::endl;
    std::cout << "Total syst for MZ' = " << *it << " : " << totalSyst << " pb" << std::endl;
  }
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Compute JEC systematic", ' ', "0.1");

    TCLAP::MultiArg<int> massArg("m", "mass", "Zprime mass", false, "integer", cmd);

    cmd.parse(argc, argv);

    std::vector<int> masses = massArg.getValue();
    if (masses.size() == 0) {
      masses.push_back(750);
      masses.push_back(1000);
      masses.push_back(1250);
      masses.push_back(1500);
    }

    computeSyst(masses);

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
