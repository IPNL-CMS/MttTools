#include "TMVA/Config.h"
#include "TMVA/Factory.h"

#include <vector>
#include <fstream>

#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>

#include "tclap/CmdLine.h"

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
int main(int argc, char** argv) {

  try {
    TCLAP::CmdLine cmd("Train the MVA algorithm", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", true, "", "string", cmd);

    TCLAP::SwitchArg chi2Arg("", "chi2", "Use chi2 sorting algorithm", false);
    TCLAP::SwitchArg mvaArg("", "mva", "Use MVA instead of chi2", false);
    TCLAP::SwitchArg kfArg("", "kf", "Use KF instead of chi2", false);
    TCLAP::SwitchArg hybridArg("", "hybrid", "Use hybrid method for sorting algorithm", false);
    std::vector<TCLAP::Arg*>  xorlist;
    xorlist.push_back(&chi2Arg);
    xorlist.push_back(&mvaArg);
    xorlist.push_back(&kfArg);
    xorlist.push_back(&hybridArg);
    cmd.xorAdd( xorlist );

    cmd.parse(argc, argv);

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    TChain* signal = loadChain(inputFiles, "signal");
    TChain* background = loadChain(inputFiles, "background");

    TFile* output = TFile::Open(outputFileArg.getValue().c_str(), "recreate");

    TMVA::Factory* factory = new TMVA::Factory("MVA_test", output, "V");
    factory->AddSignalTree(signal, 1.);
    factory->AddBackgroundTree(background, 1.);
    factory->SetWeightExpression("20000 * weight");

    factory->AddVariable("aplanarity");
    factory->AddVariable("circularity");
    factory->AddVariable("sphericity");
    factory->AddVariable("mean_csv");
    factory->AddVariable("neutrino_pt");
    //factory->AddVariable("neutrino_eta");
    factory->AddVariable("lepton_pt");
    factory->AddVariable("lepton_eta");
    factory->AddVariable("leptonic_B_pt");
    factory->AddVariable("leptonic_B_eta");
    factory->AddVariable("hadronic_B_pt");
    factory->AddVariable("hadronic_B_eta");
    factory->AddVariable("hadronic_first_jet_pt");
    factory->AddVariable("hadronic_first_jet_eta");
    factory->AddVariable("hadronic_second_jet_pt");
    factory->AddVariable("hadronic_second_jet_eta");
    factory->AddVariable("leptonic_W_pt");
    //factory->AddVariable("leptonic_W_eta");
    factory->AddVariable("leptonic_W_rapidity");
    //factory->AddVariable("leptonic_W_mass");
    factory->AddVariable("leptonic_W_transverse_mass");
    factory->AddVariable("hadronic_W_pt");
    //factory->AddVariable("hadronic_W_eta");
    factory->AddVariable("hadronic_W_rapidity");
    if (chi2Arg.isSet())
      factory->AddVariable("hadronic_W_mass");
    factory->AddVariable("hadronic_W_transverse_mass");
    factory->AddVariable("leptonic_T_pt");
    //factory->AddVariable("leptonic_T_eta");
    factory->AddVariable("leptonic_T_rapidity");
    if (chi2Arg.isSet())
      factory->AddVariable("leptonic_T_mass");
    factory->AddVariable("leptonic_T_transverse_mass");
    factory->AddVariable("hadronic_T_pt");
    //factory->AddVariable("hadronic_T_eta");
    factory->AddVariable("hadronic_T_rapidity");
    if (chi2Arg.isSet())
      factory->AddVariable("hadronic_T_mass");
    factory->AddVariable("hadronic_T_transverse_mass");
    factory->AddVariable("resonance_pt");
    factory->AddVariable("resonance_eta");
    //factory->AddVariable("resonance_rapidity");
    factory->AddVariable("neutrino_lepton_delta_R");
    //factory->AddVariable("neutrino_lepton_delta_eta");
    //factory->AddVariable("neutrino_lepton_delta_phi");
    //factory->AddVariable("neutrino_lepton_delta_rapidity");
    factory->AddVariable("neutrino_leptonic_B_delta_R");
    //factory->AddVariable("neutrino_leptonic_B_delta_eta");
    //factory->AddVariable("neutrino_leptonic_B_delta_phi");
    //factory->AddVariable("neutrino_leptonic_B_delta_rapidity");
    factory->AddVariable("neutrino_hadronic_B_delta_R");
    //factory->AddVariable("neutrino_hadronic_B_delta_eta");
    //factory->AddVariable("neutrino_hadronic_B_delta_phi");
    //factory->AddVariable("neutrino_hadronic_B_delta_rapidity");
    factory->AddVariable("neutrino_hadronic_first_jet_delta_R");
    //factory->AddVariable("neutrino_hadronic_first_jet_delta_eta");
    //factory->AddVariable("neutrino_hadronic_first_jet_delta_phi");
    //factory->AddVariable("neutrino_hadronic_first_jet_delta_rapidity");
    factory->AddVariable("neutrino_hadronic_second_jet_delta_R");
    //factory->AddVariable("neutrino_hadronic_second_jet_delta_eta");
    //factory->AddVariable("neutrino_hadronic_second_jet_delta_phi");
    //factory->AddVariable("neutrino_hadronic_second_jet_delta_rapidity");
    factory->AddVariable("lepton_leptonic_B_delta_R");
    //factory->AddVariable("lepton_leptonic_B_delta_eta");
    //factory->AddVariable("lepton_leptonic_B_delta_phi");
    //factory->AddVariable("lepton_leptonic_B_delta_rapidity");
    factory->AddVariable("lepton_hadronic_B_delta_R");
    //factory->AddVariable("lepton_hadronic_B_delta_eta");
    //factory->AddVariable("lepton_hadronic_B_delta_phi");
    //factory->AddVariable("lepton_hadronic_B_delta_rapidity");
    factory->AddVariable("lepton_hadronic_first_jet_delta_R");
    //factory->AddVariable("lepton_hadronic_first_jet_delta_eta");
    //factory->AddVariable("lepton_hadronic_first_jet_delta_phi");
    //factory->AddVariable("lepton_hadronic_first_jet_delta_rapidity");
    factory->AddVariable("lepton_hadronic_second_jet_delta_R");
    //factory->AddVariable("lepton_hadronic_second_jet_delta_eta");
    //factory->AddVariable("lepton_hadronic_second_jet_delta_phi");
    //factory->AddVariable("lepton_hadronic_second_jet_delta_rapidity");
    factory->AddVariable("leptonic_B_hadronic_B_delta_R");
    //factory->AddVariable("leptonic_B_hadronic_B_delta_eta");
    //factory->AddVariable("leptonic_B_hadronic_B_delta_phi");
    //factory->AddVariable("leptonic_B_hadronic_B_delta_rapidity");
    factory->AddVariable("leptonic_B_hadronic_first_jet_delta_R");
    //factory->AddVariable("leptonic_B_hadronic_first_jet_delta_eta");
    //factory->AddVariable("leptonic_B_hadronic_first_jet_delta_phi");
    //factory->AddVariable("leptonic_B_hadronic_first_jet_delta_rapidity");
    factory->AddVariable("leptonic_B_hadronic_second_jet_delta_R");
    //factory->AddVariable("leptonic_B_hadronic_second_jet_delta_eta");
    //factory->AddVariable("leptonic_B_hadronic_second_jet_delta_phi");
    //factory->AddVariable("leptonic_B_hadronic_second_jet_delta_rapidity");
    factory->AddVariable("hadronic_B_hadronic_first_jet_delta_R");
    //factory->AddVariable("hadronic_B_hadronic_first_jet_delta_eta");
    //factory->AddVariable("hadronic_B_hadronic_first_jet_delta_phi");
    //factory->AddVariable("hadronic_B_hadronic_first_jet_delta_rapidity");
    factory->AddVariable("hadronic_B_hadronic_second_jet_delta_R");
    //factory->AddVariable("hadronic_B_hadronic_second_jet_delta_eta");
    //factory->AddVariable("hadronic_B_hadronic_second_jet_delta_phi");
    //factory->AddVariable("hadronic_B_hadronic_second_jet_delta_rapidity");
    factory->AddVariable("hadronic_first_jet_hadronic_second_jet_delta_R");
    //factory->AddVariable("hadronic_first_jet_hadronic_second_jet_delta_eta");
    //factory->AddVariable("hadronic_first_jet_hadronic_second_jet_delta_phi");
    //factory->AddVariable("hadronic_first_jet_hadronic_second_jet_delta_rapidity");

    factory->PrepareTrainingAndTestTree("", "", "V:VerboseLevel=Info");

    factory->BookMethod(TMVA::Types::kBDT, "BDT", "V:NTrees=2000:BoostType=AdaBoost:AdaBoostBeta=0.6:UseNVars=6:nCuts=2000:MaxDepth=4");
    //factory->BookMethod(TMVA::Types::kMLP, "NN", "V:VarTransform=D");

    (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 1000;

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    output->Close();
    delete output;

    delete signal;
    delete background;

  } catch (TCLAP::ArgException& e) {
    std::cout << e.what() << std::endl;
  }

}
