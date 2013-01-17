#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>
#include <unistd.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <fcntl.h>

#include <tclap/CmdLine.h>
#include <json/json.h>

#include <TRandom.h>
#include <TRandom2.h>

#include "Utils.h"

std::vector<std::string> JEC_PARAMS;
std::string base_path = "";

void process(int mass, bool muonsOnly, int btag, const std::string& file, bool singleFile, bool fixBackground, std::map<std::string, float>& crossSections) {
  std::vector<pid_t> children; // For fork()

  TRandom* random = new TRandom2(0);
  int shmid;

  // Create a shared memory area to store the cross sections
  key_t key = (key_t) random->Uniform(1024, 1000000);
  if ((shmid = shmget(key, sizeof(SHMFitResults) * 8, IPC_CREAT | 0666)) < 0) {
    perror("Can't create shared memory area");
    exit(1);
  }

  delete random;

  SHMFitResults* shm = NULL;
  if ((shm = static_cast<SHMFitResults*>(shmat(shmid, NULL, 0))) == (void *) -1) {
    perror("Can't map shared memory to local memory");
    exit(1);
  }

  std::stringstream strKey;
  strKey << key;

  for (std::vector<std::string>::const_iterator param = JEC_PARAMS.begin(); param != JEC_PARAMS.end(); ++param) {
    std::stringstream ss;
    ss << mass;

    std::cout << "Fitting for " << *param << std::endl;
    pid_t child = fork();
    if (child == 0) {
      
      int fd = open("/dev/null", O_WRONLY);
      dup2(fd, STDOUT_FILENO); 
      dup2(fd, STDIN_FILENO);
      char** params = getSystCLParameters(ss.str(), file, singleFile, fixBackground, muonsOnly, btag, "--shared-memory", "--shm-key", strKey.str().c_str(), "--syst", (*param).c_str(), NULL);
      if (execv("./fitMtt", params) < 0) {
        perror("Failed to launch fitMtt");
      }
      exit(0);
    } else {
      waitpid(child, NULL, 0);
    }
    std::cout << "Done. # events: " << shm->nSignalEvents << "; Cross-section: " << shm->sigma << " pb." << std::endl << std::endl;
    crossSections[*param] = shm->sigma;
  }
}

void fillParams() {
  JEC_PARAMS.push_back("nominal");
  JEC_PARAMS.push_back("JECup");
  JEC_PARAMS.push_back("JECdown");
}

void saveSystematic(double syst, const std::string& outputFile) {

  Json::Value root;
  root["result"]= syst;

  Json::StyledWriter writer;
  const std::string json = writer.write(root);

  std::ofstream file(outputFile.c_str());
  file << json;
  file.close();
}

double computeSystValue(double sigma, double up, double down) {
  double systup = fabs(up - sigma);
  double systdown = fabs(down - sigma);
  return (systup + systdown) / 2.;
}

void computeSyst(std::map<std::string, float>& crossSections, const std::string& outputFile) {

  double syst = computeSystValue(crossSections["nominal"], crossSections["JECup"], crossSections["JECdown"]);
  saveSystematic(syst, outputFile);
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Compute JEC systematic", ' ', "0.1");
    
    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::SwitchArg muonsOnlyArg("", "muons-only", "Compute sigmaref using only semi-mu data", cmd);
    TCLAP::SwitchArg extractArg("", "dont-extract", "Don't run fitMtt. Only compute systematic with previous results", cmd);
    //TCLAP::MultiArg<int> massArg("m", "mass", "Zprime mass", false, "integer", cmd);
    TCLAP::ValueArg<int> massArg("m", "mass", "Zprime mass", true, 0, "integer", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "integer", cmd);

    TCLAP::SwitchArg fixBackgroundArg("", "fix-background", "Fix background when fitting", cmd);

    TCLAP::ValueArg<std::string> outputArg("o", "output", "JSON output file", true, "", "string", cmd);

    cmd.parse(argc, argv);

    /*std::vector<int> masses = massArg.getValue();
    if (masses.size() == 0) {
      masses.push_back(750);
      masses.push_back(1000);
      masses.push_back(1250);
      masses.push_back(1500);
    }*/

    base_path = "analysis/" + getAnalysisUUID();

    fillParams();

    std::map<std::string, float> crossSections;
    process(massArg.getValue(), muonsOnlyArg.getValue(), btagArg.getValue(), inputFileArg.isSet() ? inputFileArg.getValue() : inputListArg.getValue(), inputFileArg.isSet(), fixBackgroundArg.getValue(), crossSections);

    computeSyst(crossSections, outputArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
