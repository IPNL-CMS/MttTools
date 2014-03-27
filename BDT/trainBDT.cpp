#include "TMVA/Factory.h"

#include <vector>
#include <fstream>

#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>

#include "tclap/CmdLine.h"

// Load libdpm at startup, on order to be sure that rfio files are working
#include <dlfcn.h>
struct Dummy
{
  Dummy()
  {
    dlopen("libdpm.so", RTLD_NOW|RTLD_GLOBAL);
  }
};
static Dummy foo;

void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

TChain* loadChain(const std::vector<std::string>& inputFiles, const std::string& treeName) {
  TChain* c = new TChain(treeName.c_str());
  for (const auto& inputFile: inputFiles)
    c->Add(inputFile.c_str());

  return c;
}

void process(const std::vector<std::string>& inputFiles, const std::string& name, const std::string& outputFile) {
  TChain* signal = loadChain(inputFiles, "signal");
  TChain* background = loadChain(inputFiles, "background");

  TFile* output = TFile::Open(outputFile.c_str(), "recreate");

  TMVA::Factory* factory = new TMVA::Factory(name.c_str(), output, "V");
  factory->AddSignalTree(signal, 1.);
  factory->AddBackgroundTree(background, 1.);

  factory->AddVariable("lightJet1p2_Pt");
  factory->AddVariable("leptonic_B_Pt");
  factory->AddVariable("leptonic_W_Pt");
  factory->AddVariable("leptonic_Top_Pt");
  factory->AddVariable("leptonic_Top_M");
  factory->AddVariable("hadronic_B_Pt");
  //factory->AddVariable("hadronic_W_Pt");
  factory->AddVariable("hadronic_W_M");
  factory->AddVariable("hadronic_Top_Pt");
  factory->AddVariable("hadronic_Top_M");

  factory->AddVariable("delta_phi_tops");
  factory->AddVariable("delta_phi_lightjets");
  factory->AddVariable("delta_phi_W");
  factory->AddVariable("delta_R_tops");
  factory->AddVariable("delta_R_lightjets");
  factory->AddVariable("delta_R_W");
  
  factory->AddVariable("leptonic_B_CSV");
  factory->AddVariable("hadronic_B_CSV");

  //factory->AddVariable("pt_tt_system");
  factory->AddVariable("ht_fraction");

  factory->SetWeightExpression("weight");

  factory->PrepareTrainingAndTestTree("", "", "V:VerboseLevel=Info");

  //factory->BookMethod( TMVA::Types::kBDT, "BDT", "V:BoostType=AdaBoost:nCuts=20:VarTransform=D");
  factory->BookMethod(TMVA::Types::kBDT, "BDT", "V:BoostType=AdaBoost:nCuts=20");
  factory->BookMethod(TMVA::Types::kMLP, "NN", "V:VarTransform=Norm");
  //factory->BookMethod(TMVA::Types::kPDERS, "PDERS", "V");

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  output->Close();
  delete output;

  delete signal;
  delete background;
}

int main(int argc, char** argv) {

  try {

    TCLAP::CmdLine cmd("Train BDT", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> nameArg("n", "name", "dataset name", true, "", "string", cmd);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", true, "", "string", cmd);

    cmd.parse(argc, argv);

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }
    
    process(inputFiles, nameArg.getValue(), outputFileArg.getValue());

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
