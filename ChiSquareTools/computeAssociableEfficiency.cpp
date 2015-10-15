//#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <algorithm>
#include <string>

#include <ios>

#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <TH1F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TCut.h>
#include <TProfile.h>

#include <tclap/CmdLine.h>

// Variables

float MC_mtt;

int genLeptonicBIndex, genHadronicBIndex, genFirstJetIndex, genSecondJetIndex;

bool CSV_MODE = true;

double binomialError(double efficiency, int n) {
  return sqrt((efficiency * (1 - efficiency)) / (n));
}

void printEff(double efficiency, int n) {
  if (CSV_MODE) {
    std::cout << std::left << std::setw(15) << efficiency << "\t" << binomialError(efficiency, n) << std::endl;
  } else {
    std::cout << efficiency * 100 << " +/- " << binomialError(efficiency, n) * 100 << std::endl;
  }
}

void loadChain(const std::vector<std::string>& inputFiles, TChain*& mtt, TChain*& event) {
  std::cout << "Opening files..." << std::endl;

  //mc = new TChain("MC");
  event = new TChain("event");
  //jets = new TChain("jet_PF");
  //MET = new TChain("MET_PF");
  //muons = new TChain("muon_PF");
  //electrons = new TChain("electron_PF");
  mtt = new TChain("Mtt");

  for (const std::string& file: inputFiles) {
    //mc->Add(file.c_str());
    event->Add(file.c_str());
    //jets->Add(file.c_str());
    //MET->Add(file.c_str());
    //muons->Add(file.c_str());
    //electrons->Add(file.c_str());
    mtt->Add(file.c_str());
  }

  //jets->SetCacheSize(30*1024*1024);
  event->SetCacheSize(30*1024*1024);
  mtt->SetCacheSize(30*1024*1024);

  std::cout << "... done." << std::endl;
}

void SetBranchAddress(TChain* chain, const char* branchName, void* address) {
  chain->SetBranchStatus(branchName, 1);
  chain->SetBranchAddress(branchName, address, NULL);
}

bool isMatched(int index, int jetsMcIndex[], TClonesArray* p4s, int nJets) {
  for (int i = 0; i < nJets; i++) {
    if (jetsMcIndex[i] == index) {
      TLorentzVector* p4 = static_cast<TLorentzVector*>((*p4s)[i]);
      return p4->Pt() > 30 && fabs(p4->Eta()) < 2.4;
    }
  }

  return false;
}


void process(const std::vector<std::string>& inputFiles) {

  //TH1::SetDefaultSumw2(true);
  gROOT->SetBatch(true);

  TChain* MC = NULL, *event = NULL, *jets = NULL, *MET = NULL, *muons = NULL, *electrons = NULL, *mtt = NULL;
  loadChain(inputFiles, mtt, event);

  //MC->SetBranchStatus("*", 0);
  event->SetBranchStatus("*", 0);
  //jets->SetBranchStatus("*", 0);
  //MET->SetBranchStatus("*", 0);
  //muons->SetBranchStatus("*", 0);
  //electrons->SetBranchStatus("*", 0);
  mtt->SetBranchStatus("*", 0);

  // Jets
/*  SetBranchAddress(jets, "n_jets", &n_jets);*/
  //SetBranchAddress(jets, "jet_4vector", &m_jet_lorentzvector);
  /*SetBranchAddress(jets, "jet_mcParticleIndex", &m_jet_MCIndex);*/

  // Event
  float generator_weight;
  SetBranchAddress(event, "generator_weight", &generator_weight);

  // MTT
  SetBranchAddress(mtt, "MC_leptonicBIndex", &genLeptonicBIndex);
  SetBranchAddress(mtt, "MC_hadronicBIndex", &genHadronicBIndex);
  SetBranchAddress(mtt, "MC_hadronicFirstJetIndex", &genFirstJetIndex);
  SetBranchAddress(mtt, "MC_hadronicSecondJetIndex", &genSecondJetIndex);

  int numComb, channel, isSel;
  bool eventIsAssociable;
  SetBranchAddress(mtt, "MC_channel", &channel);
  SetBranchAddress(mtt, "isSel", &isSel);
  SetBranchAddress(mtt, "numComb_chi2", &numComb);
  SetBranchAddress(mtt, "MC_mtt", &MC_mtt);
  SetBranchAddress(mtt, "eventIsAssociable", &eventIsAssociable);

  SetBranchAddress(mtt, "MC_mtt", &MC_mtt);

  uint64_t possibleAssociableEntries = 0;
  uint64_t associableEntries = 0;

  uint64_t possibleAssociableEntries_450_550 = 0;
  uint64_t associableEntries_450_550 = 0;

  uint64_t possibleAssociableEntries_550_650 = 0;
  uint64_t associableEntries_550_650 = 0;

  uint64_t possibleAssociableEntries_650_750 = 0;
  uint64_t associableEntries_650_750 = 0;

  uint64_t possibleAssociableEntries_750_850 = 0;
  uint64_t associableEntries_750_850 = 0;

  uint64_t entries = mtt->GetEntries();
  //uint64_t entries = 0;

  for (uint64_t entry = 0; entry < entries; entry++) {
    //MC->GetEntry(entry);
    event->GetEntry(entry);
    //jets->GetEntry(entry);
    //MET->GetEntry(entry);
    //muons->GetEntry(entry);
    //electrons->GetEntry(entry);
    mtt->GetEntry(entry);

    if (((entry) % 100000) == 0) {
      std::cout << "Processing entry " << entry + 1 << " out of " << entries << " (" << (entry + 1) / (float) entries * 100 << "%)" << std::endl;
    }

    if (genLeptonicBIndex == -1 || genHadronicBIndex == -1 || genFirstJetIndex == -1 || genHadronicBIndex == -1)
      continue;

    if (channel != 1 && channel != 2)
      continue;

    if (isSel != 1)
      continue;

    if (generator_weight > 0.) {
      possibleAssociableEntries++;
      if(eventIsAssociable) {
        associableEntries++;
      }    
    }

    if (MC_mtt >= 450 && MC_mtt < 550) {
      possibleAssociableEntries_450_550++;
      if(eventIsAssociable) {
        associableEntries_450_550++;
      }
    } else if (MC_mtt >= 550 && MC_mtt < 650) {
      possibleAssociableEntries_550_650++;
      if(eventIsAssociable) {
        associableEntries_550_650++;
      }
    } else if (MC_mtt >= 650 && MC_mtt < 750) {
      possibleAssociableEntries_650_750++;
      if(eventIsAssociable) {
        associableEntries_650_750++;
      }
    } else if (MC_mtt >= 750 && MC_mtt < 850) {
      possibleAssociableEntries_750_850++;
      if(eventIsAssociable) {
        associableEntries_750_850++;
      }
    }


  }

  std::cout << "Number of total entries: " << entries << std::endl;

  std::cout << "Associable efficiency" << std::endl;
  printEff((double) associableEntries / (possibleAssociableEntries), possibleAssociableEntries);

  std::cout << "Associable efficiency for Mtt Gen in [450 ; 550[ GeV" << std::endl;
  printEff((double) associableEntries_450_550 / (possibleAssociableEntries_450_550), possibleAssociableEntries_450_550);

  std::cout << "Associable efficiency for Mtt Gen in [550 ; 650[ GeV" << std::endl;
  printEff((double) associableEntries_550_650 / (possibleAssociableEntries_550_650), possibleAssociableEntries_550_650);

  std::cout << "Associable efficiency for Mtt Gen in [650 ; 750[ GeV" << std::endl;
  printEff((double) associableEntries_650_750 / (possibleAssociableEntries_650_750), possibleAssociableEntries_650_750);

  std::cout << "Associable efficiency for Mtt Gen in [750 ; 850[ GeV" << std::endl;
  printEff((double) associableEntries_750_850 / (possibleAssociableEntries_750_850), possibleAssociableEntries_750_850);

  delete MC; delete event; delete MET; delete jets; delete muons; delete electrons;

}

void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  std::ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

int main(int argc, char** argv)
{

  try {
    TCLAP::CmdLine cmd("compute various distribution to optimize chi square", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    //TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", false, "", "string", cmd);

    /*TCLAP::ValueArg<std::string> typeArg("", "type", "current inputfile type (semie or semimu)", false, "", "string", cmd);
      TCLAP::ValueArg<std::string> pileupArg("", "pileup", "PU profile used for MC production", false, "S7", "string", cmd);*/

    TCLAP::SwitchArg csvModeArg("", "csv", "csv mode", cmd);

    cmd.parse(argc, argv);

    CSV_MODE = csvModeArg.getValue();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    process(inputFiles);    

  } catch (TCLAP::ArgException& e) {
    std::cout << e.what() << std::endl;
  }

}

