#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>

#include <tclap/CmdLine.h>
#include <json/json.h>

void process(std::vector<int>& masses, bool muonsOnly, bool onlyLumiSyst, int btag) {
  std::vector<pid_t> children; // For fork()

  for (std::vector<int>::iterator it = masses.begin(); it != masses.end(); ++it) {

    pid_t pid = fork();
    if (pid == 0) {
      std::stringstream ss;
      ss << *it;
      std::string mass = ss.str();

      ss.clear(); ss.str(std::string());
      ss << btag;
      std::string btagStr = ss.str();

      if (! muonsOnly) {
        if (!onlyLumiSyst)
          execl("./fitMtt", "fitMtt", "-m", mass.c_str(), "--scan", "--no-text-files", "--no-root-files", "--no-figs", "--b-tag", btagStr.c_str(), NULL);
        else
          execl("./fitMtt", "fitMtt", "-m", mass.c_str(), "--scan", "--no-text-files", "--no-root-files", "--no-figs", "--only-lumi-syst", "--b-tag", btagStr.c_str(), NULL);
      } else {
        if (!onlyLumiSyst)
          execl("./fitMtt", "fitMtt", "-m", mass.c_str(), "--scan", "--no-text-files", "--no-root-files", "--no-figs", "--muons-only", "--b-tag", btagStr.c_str(), NULL);
        else
          execl("./fitMtt", "fitMtt", "-m", mass.c_str(), "--scan", "--no-text-files", "--no-root-files", "--no-figs", "--muons-only", "--only-lumi-syst", "--b-tag", btagStr.c_str(), NULL);
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

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Run Likelood scan", ' ', "0.1");

    TCLAP::MultiArg<int> massArg("m", "mass", "Zprime mass", false, "integer", cmd);
    TCLAP::SwitchArg muonsOnlyArg("", "muons-only", "Compute sigmaref using only semi-mu data", cmd);
    TCLAP::SwitchArg onlyLumiSystArg("", "only-lumi-syst", "Only use luminosity error for systematics", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);

    cmd.parse(argc, argv);

    std::vector<int> masses = massArg.getValue();
    if (masses.size() == 0) {
      masses.push_back(750);
      masses.push_back(1000);
      masses.push_back(1250);
      masses.push_back(1500);
    }

    process(masses, muonsOnlyArg.getValue(), onlyLumiSystArg.getValue(), btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
