//#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <algorithm>
#include <string>

#include <TChain.h>
#include <TFile.h>

#include <tclap/CmdLine.h>

#include "../PUReweighting/PUReweighter.h"

bool OVERRIDE_TYPE;
std::string OVERRIDED_TYPE;
PUProfile puProfile;

void reduce(TChain* mtt, TChain* event, const std::string& outputFile, bool isData, const std::string& type) {

  TString outputFileFormated = OVERRIDE_TYPE ? outputFile : TString::Format(outputFile.c_str(), (isData) ? "" : type.c_str());
  std::cout << outputFileFormated << std::endl;

  std::map<int, TTree*> outputTrees {
    {0, new TTree("dataset_0btag", "dataset for 0 b-tagged jet") },
      {1, new TTree("dataset_1btag", "dataset for 1 b-tagged jet") },
      {2, new TTree("dataset_2btag", "dataset for at least 2 b-tagged jets") }
  };

  float mtt_afterChi2, pt_1stJet, pt_2ndJet, bestSolChi2;
  int isSel, nBtaggedJets_CSVM;

  mtt->SetBranchAddress("mtt_AfterChi2", &mtt_afterChi2, NULL);
  mtt->SetBranchAddress("1stjetpt", &pt_1stJet, NULL);
  mtt->SetBranchAddress("2ndjetpt", &pt_2ndJet, NULL);
  mtt->SetBranchAddress("bestSolChi2", &bestSolChi2, NULL);
  mtt->SetBranchAddress("isSel", &isSel, NULL);
  if (puProfile == PUProfile::S6)
    mtt->SetBranchAddress("nBtaggedJets_TCHET", &nBtaggedJets_CSVM, NULL);
  else
    mtt->SetBranchAddress("nBtaggedJets_CSVM", &nBtaggedJets_CSVM, NULL);

  mtt->SetBranchStatus("*", 0);
  mtt->SetBranchStatus("mtt_AfterChi2", 1);
  mtt->SetBranchStatus("1stjetpt", 1);
  mtt->SetBranchStatus("2ndjetpt", 1);
  mtt->SetBranchStatus("bestSolChi2", 1);
  mtt->SetBranchStatus("isSel", 1);
  if (puProfile == PUProfile::S6)
    mtt->SetBranchStatus("nBtaggedJets_TCHET", 1);
  else
    mtt->SetBranchStatus("nBtaggedJets_CSVM", 1);

  float n_trueInteractions;

  event->SetBranchAddress("nTrueInteractions", &n_trueInteractions, NULL);

  event->SetBranchStatus("*", 0);
  event->SetBranchStatus("nTrueInteractions", 1);

  float weight;
  int lepton = (type == "semimu") ? 13 : 11;
  for (auto& outputTree: outputTrees) {
    outputTree.second->Branch("mtt", &mtt_afterChi2, "mtt/F");
    outputTree.second->Branch("weight", &weight, "weight/F");
    outputTree.second->Branch("lepton_type", &lepton, "lepton_type/I");
  }

  PUReweighter* puReweigher = NULL;
  if (! isData) {
    puReweigher = new PUReweighter(type == "semimu", puProfile);
  }

  int64_t entries = mtt->GetEntries();
  for (int64_t i = 0; i < entries; i++) {
    mtt->GetEntry(i);
    event->GetEntry(i);

    if (i % 1000000 == 0) {
      std::cout << "Processing event #" << i + 1 << " over " << entries << " (" << (float) i / entries * 100 << " %)" << std::endl;
    }

    // Selection
    if (isSel == 1 && pt_1stJet > 70. && pt_2ndJet > 50. && mtt_afterChi2 > 0. && bestSolChi2 < 500.) {
      // Good event

      weight = 1.;
      if (! isData) {
        weight = puReweigher->weight(n_trueInteractions);
      }

      int index;
      if (nBtaggedJets_CSVM == 0)
        index = 0;
      else if (nBtaggedJets_CSVM == 1)
        index = 1;
      else
        index = 2;

      outputTrees[index]->Fill();
    }
  }

  TFile* output = TFile::Open(outputFileFormated, "recreate");

  for (auto& tree: outputTrees) {
    tree.second->SetDirectory(output);
    tree.second->Write();
  }

  output->Close();
  delete output;

  delete puReweigher;
}

void loadChain(const std::vector<std::string>& inputFiles, bool isData, const std::string& suffix, TChain*& mtt, TChain*& event) {

  mtt = new TChain("Mtt");
  event = new TChain("event");

  if (isData) {
    for (const std::string& file: inputFiles) {
      mtt->Add(file.c_str());
      event->Add(file.c_str());
    }
  } else {
    for (const std::string& file: inputFiles) {
      TString f = OVERRIDE_TYPE ? file : TString::Format(file.c_str(), suffix.c_str());
      mtt->Add(f);
      event->Add(f);
    }
  }
}

void reduce(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isData) {

  TChain* mtt = NULL, *event = NULL;
  if (isData || OVERRIDE_TYPE) {
    
    std::string type;
    if (isData) {
      type = TString(outputFile).Contains("Mu") ? "semimu" : "semie";
    } else {
      type = OVERRIDED_TYPE;
    }

    loadChain(inputFiles, isData, type, mtt, event);
    reduce(mtt, event, outputFile, isData, type);
    delete mtt; delete event;

  } else {

      loadChain(inputFiles, false, "semimu", mtt, event);
      reduce(mtt, event, outputFile, false, "semimu");
      delete mtt; delete event;

      loadChain(inputFiles, false, "semie", mtt, event);
      reduce(mtt, event, outputFile, false, "semie");
      delete mtt; delete event;

  }
}

void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

int main(int argc, char** argv)
{
  try {
    TCLAP::CmdLine cmd("reduce Zprime dataset", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", true, "", "string", cmd);

    TCLAP::SwitchArg dataArg("", "data", "Is this data?", false);
    TCLAP::SwitchArg mcArg("", "mc", "Is this mc?", false);

    cmd.xorAdd(dataArg, mcArg);

    TCLAP::ValueArg<std::string> typeArg("", "type", "current inputfile type (semie or semimu)", false, "", "string", cmd);
    TCLAP::ValueArg<std::string> pileupArg("", "pileup", "PU profile used for MC production", false, "S7", "string", cmd);

    cmd.parse(argc, argv);

    if (typeArg.isSet()) {
      OVERRIDE_TYPE = true;
      OVERRIDED_TYPE = typeArg.getValue();
    }

    std::string p = pileupArg.getValue();
    std::transform(p.begin(), p.end(), p.begin(), ::tolower);
    if (p == "s6")
      puProfile = PUProfile::S6;
    else
      puProfile = PUProfile::S7;
    
    bool isData = dataArg.isSet();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    reduce(inputFiles, outputFileArg.getValue(), isData);    

  } catch (TCLAP::ArgException& e) {
    std::cout << e.what() << std::endl;
  }

}
