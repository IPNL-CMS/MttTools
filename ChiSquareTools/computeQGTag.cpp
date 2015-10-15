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
static const int 	m_jets_MAX       = 200;

TClonesArray* m_jet_lorentzvector = new TClonesArray("TLorentzVector");
int  n_jets;
int  m_jet_MCIndex[m_jets_MAX];
float  m_jet_qgtag_likelihood[m_jets_MAX];


float MC_mtt;

int genLeptonicBIndex, genHadronicBIndex, genFirstJetIndex, genSecondJetIndex;

uint64_t chi2_matchedEntries1 = 0;
uint64_t chi2_matchedEntries2 = 0;
uint64_t chi2_matchedEntries3 = 0;
uint64_t chi2_matchedEntries4 = 0;

uint64_t chi2_matchedAndWellPlacedEntries1 = 0;
uint64_t chi2_matchedAndWellPlacedEntries2 = 0;
uint64_t chi2_matchedAndWellPlacedEntries3 = 0;
uint64_t chi2_matchedAndWellPlacedEntries4 = 0;

uint64_t jets_matchedEntries1 = 0;
uint64_t jets_matchedEntries2 = 0;
uint64_t jets_matchedEntries3 = 0;
uint64_t jets_matchedEntries4 = 0;

const int nBins = 15;
const double bins[] = {340, 360, 380, 400, 420, 460, 500, 550, 600, 650, 750, 850, 950, 1050, 1200, 1400, 1600};

//TProfile* chi2_matchedEntries_vs_mtt = new TProfile("matched_entries_vs_mtt", "Efficiency;Generated m_{tt};Number of matched jets", nBins, bins);
TH1F* chi2_selectedEntries_vs_mtt;
TH1F* chi2_matchedEntries_vs_mtt;
TH1F* chi2_efficiency_vs_mtt;

TH1F* jets_selectedEntries_vs_mtt;
TH1F* jets_matchedEntries_vs_mtt;
TH1F* jets_efficiency_vs_mtt;

TH1F* h_qgLikelihood_hadronicBJets;
TH1F* h_qgLikelihood_leptonicBJets;
TH1F* h_qgLikelihood_lightJets;
TH1F* h_qgLikelihood_allJets;
TH1F* h_qgLikelihood_allTTJets;
TH1F* h_qgLikelihood_allNoTTJets;

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

void loadChain(const std::vector<std::string>& inputFiles, TChain*& jets, TChain*& mtt) {
  std::cout << "Opening files..." << std::endl;

  //mc = new TChain("MC");
  //event = new TChain("event");
  jets = new TChain("jet_PF");
  //MET = new TChain("MET_PF");
  //muons = new TChain("muon_PF");
  //electrons = new TChain("electron_PF");
  mtt = new TChain("Mtt");

  for (const std::string& file: inputFiles) {
    //mc->Add(file.c_str());
    //event->Add(file.c_str());
    jets->Add(file.c_str());
    //MET->Add(file.c_str());
    //muons->Add(file.c_str());
    //electrons->Add(file.c_str());
    mtt->Add(file.c_str());
  }

  jets->SetCacheSize(30*1024*1024);
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


void computeEfficiencyChiSquare(int selectedLeptonicBIndex, int selectedHadronicBIndex, int selectedFirstJetIndex, int selectedSecondJetIndex) {
  int numberOfMatchedJets = 0;
  if (m_jet_MCIndex[selectedLeptonicBIndex] == genLeptonicBIndex || m_jet_MCIndex[selectedLeptonicBIndex] == genHadronicBIndex ||
      m_jet_MCIndex[selectedLeptonicBIndex] == genFirstJetIndex || m_jet_MCIndex[selectedLeptonicBIndex] == genSecondJetIndex)
    numberOfMatchedJets++;

  if (m_jet_MCIndex[selectedHadronicBIndex] == genLeptonicBIndex || m_jet_MCIndex[selectedHadronicBIndex] == genHadronicBIndex ||
      m_jet_MCIndex[selectedHadronicBIndex] == genFirstJetIndex || m_jet_MCIndex[selectedHadronicBIndex] == genSecondJetIndex)
    numberOfMatchedJets++;

  if (m_jet_MCIndex[selectedFirstJetIndex] == genLeptonicBIndex || m_jet_MCIndex[selectedFirstJetIndex] == genHadronicBIndex ||
      m_jet_MCIndex[selectedFirstJetIndex] == genFirstJetIndex || m_jet_MCIndex[selectedFirstJetIndex] == genSecondJetIndex)
    numberOfMatchedJets++;

  if (m_jet_MCIndex[selectedSecondJetIndex] == genLeptonicBIndex || m_jet_MCIndex[selectedSecondJetIndex] == genHadronicBIndex ||
      m_jet_MCIndex[selectedSecondJetIndex] == genFirstJetIndex || m_jet_MCIndex[selectedSecondJetIndex] == genSecondJetIndex)
    numberOfMatchedJets++;

  int numberOfMatchedAndWellPlacedJets = 0;
  if (m_jet_MCIndex[selectedLeptonicBIndex] == genLeptonicBIndex)
    numberOfMatchedAndWellPlacedJets++;

  if (m_jet_MCIndex[selectedHadronicBIndex] == genHadronicBIndex)
    numberOfMatchedAndWellPlacedJets++;

  if (m_jet_MCIndex[selectedFirstJetIndex] == genFirstJetIndex || m_jet_MCIndex[selectedFirstJetIndex] == genSecondJetIndex)
    numberOfMatchedAndWellPlacedJets++;

  if (m_jet_MCIndex[selectedSecondJetIndex] == genSecondJetIndex || m_jet_MCIndex[selectedSecondJetIndex] == genFirstJetIndex)
    numberOfMatchedAndWellPlacedJets++;

  switch (numberOfMatchedAndWellPlacedJets) {
    case 4:
      chi2_matchedAndWellPlacedEntries4++;

    case 3:
      chi2_matchedAndWellPlacedEntries3++;

    case 2:
      chi2_matchedAndWellPlacedEntries2++;

    case 1:
      chi2_matchedAndWellPlacedEntries1++;
  }

  switch (numberOfMatchedJets) {
    case 4:
      chi2_matchedEntries4++;
      chi2_matchedEntries_vs_mtt->Fill(MC_mtt);

    case 3:
      chi2_matchedEntries3++;

    case 2:
      chi2_matchedEntries2++;

    case 1:
      chi2_matchedEntries1++;
  }


}

void computeEfficiencyFourLeadingJets(int selectedLeptonicBIndex, int selectedHadronicBIndex, int selectedFirstJetIndex, int selectedSecondJetIndex) {
int numberOfMatchedJets = 0;
  if (m_jet_MCIndex[selectedLeptonicBIndex] == genLeptonicBIndex || m_jet_MCIndex[selectedLeptonicBIndex] == genHadronicBIndex ||
      m_jet_MCIndex[selectedLeptonicBIndex] == genFirstJetIndex || m_jet_MCIndex[selectedLeptonicBIndex] == genSecondJetIndex)
    numberOfMatchedJets++;

  if (m_jet_MCIndex[selectedHadronicBIndex] == genLeptonicBIndex || m_jet_MCIndex[selectedHadronicBIndex] == genHadronicBIndex ||
      m_jet_MCIndex[selectedHadronicBIndex] == genFirstJetIndex || m_jet_MCIndex[selectedHadronicBIndex] == genSecondJetIndex)
    numberOfMatchedJets++;

  if (m_jet_MCIndex[selectedFirstJetIndex] == genLeptonicBIndex || m_jet_MCIndex[selectedFirstJetIndex] == genHadronicBIndex ||
      m_jet_MCIndex[selectedFirstJetIndex] == genFirstJetIndex || m_jet_MCIndex[selectedFirstJetIndex] == genSecondJetIndex)
    numberOfMatchedJets++;

  if (m_jet_MCIndex[selectedSecondJetIndex] == genLeptonicBIndex || m_jet_MCIndex[selectedSecondJetIndex] == genHadronicBIndex ||
      m_jet_MCIndex[selectedSecondJetIndex] == genFirstJetIndex || m_jet_MCIndex[selectedSecondJetIndex] == genSecondJetIndex)
    numberOfMatchedJets++;

  switch (numberOfMatchedJets) {
    case 4:
      jets_matchedEntries4++;
      jets_matchedEntries_vs_mtt->Fill(MC_mtt);

    case 3:
      jets_matchedEntries3++;

    case 2:
      jets_matchedEntries2++;

    case 1:
      jets_matchedEntries1++;
  }
}

void process(const std::vector<std::string>& inputFiles, const std::string& outputFile) {

  //TH1::SetDefaultSumw2(true);
  gROOT->SetBatch(true);

  jets_selectedEntries_vs_mtt = new TH1F("jets_selected_entries_vs_mtt", "Efficiency;Generated m_{tt};Number of selected jets", nBins, bins);
  chi2_selectedEntries_vs_mtt = new TH1F("chi2_selected_entries_vs_mtt", "Efficiency;Generated m_{tt};Number of selected jets", nBins, bins);

  chi2_matchedEntries_vs_mtt = new TH1F("matched_entries_vs_mtt", "Efficiency;Generated m_{tt};Number of matched jets", nBins, bins);
  jets_matchedEntries_vs_mtt = new TH1F("jets_entries_vs_mtt", "Efficiency;Generated m_{tt};Number of matched jets", nBins, bins);

  chi2_efficiency_vs_mtt = new TH1F("chi2_efficiency_vs_mtt", "Efficiency;Generated m_{tt};Jet selection efficiency", nBins, bins);
  jets_efficiency_vs_mtt = new TH1F("jets_efficiency_vs_mtt", "Efficiency;Generated m_{tt};Jet selection efficiency", nBins, bins);

  h_qgLikelihood_hadronicBJets = new TH1F("h_qgLikelihood_hadronicBJets", "QG h_qgLikelihood_hadronicBJets", 50, 0, 1);
  h_qgLikelihood_leptonicBJets = new TH1F("h_qgLikelihood_leptonicBJets", "QG h_qgLikelihood_leptonicBJets", 50, 0, 1);
  h_qgLikelihood_lightJets = new TH1F("h_qgLikelihood_lightJets", "QG h_qgLikelihood_lightJets", 50, 0, 1);
  h_qgLikelihood_allJets = new TH1F("h_qgLikelihood_allJets", "QG h_qgLikelihood_allJets", 50, 0, 1);
  h_qgLikelihood_allTTJets = new TH1F("h_qgLikelihood_allTTJets", "QG h_qgLikelihood_allTTJets", 50, 0, 1);
  h_qgLikelihood_allNoTTJets = new TH1F("h_qgLikelihood_allNoTTJets", "QG h_qgLikelihood_allNoTTJets", 50, 0, 1);

  TChain* MC = NULL, *event = NULL, *jets = NULL, *MET = NULL, *muons = NULL, *electrons = NULL, *mtt = NULL;
  loadChain(inputFiles, jets, mtt);

  //MC->SetBranchStatus("*", 0);
  //event->SetBranchStatus("*", 0);
  jets->SetBranchStatus("*", 0);
  //MET->SetBranchStatus("*", 0);
  //muons->SetBranchStatus("*", 0);
  //electrons->SetBranchStatus("*", 0);
  mtt->SetBranchStatus("*", 0);

  // Jets
  SetBranchAddress(jets, "n_jets", &n_jets);
  SetBranchAddress(jets, "jet_4vector", &m_jet_lorentzvector);
  SetBranchAddress(jets, "jet_mcParticleIndex", &m_jet_MCIndex);
  SetBranchAddress(jets, "jet_qgtag_likelihood", &m_jet_qgtag_likelihood);

  // MTT
  SetBranchAddress(mtt, "MC_leptonicBIndex", &genLeptonicBIndex);
  SetBranchAddress(mtt, "MC_hadronicBIndex", &genHadronicBIndex);
  SetBranchAddress(mtt, "MC_hadronicFirstJetIndex", &genFirstJetIndex);
  SetBranchAddress(mtt, "MC_hadronicSecondJetIndex", &genSecondJetIndex);

  int selectedLeptonicBIndex, selectedHadronicBIndex, selectedFirstJetIndex, selectedSecondJetIndex;
  SetBranchAddress(mtt, "selectedLeptonicBIndex_AfterChi2", &selectedLeptonicBIndex);
  SetBranchAddress(mtt, "selectedHadronicBIndex_AfterChi2", &selectedHadronicBIndex);
  SetBranchAddress(mtt, "selectedHadronicFirstJetIndex_AfterChi2", &selectedFirstJetIndex);
  SetBranchAddress(mtt, "selectedHadronicSecondJetIndex_AfterChi2", &selectedSecondJetIndex);

  int numComb, channel, isSel;
  SetBranchAddress(mtt, "MC_channel", &channel);
  SetBranchAddress(mtt, "isSel", &isSel);
  SetBranchAddress(mtt, "numComb_chi2", &numComb);
  SetBranchAddress(mtt, "MC_mtt", &MC_mtt);

  SetBranchAddress(mtt, "MC_mtt", &MC_mtt);

  uint64_t jets_selectedEntries = 0;
  uint64_t chi2_selectedEntries = 0;

  uint64_t possibleMatchableEntries = 0;
  uint64_t matchableEntries = 0;

  uint64_t entries = mtt->GetEntries();
  //uint64_t entries = 0;

  for (uint64_t entry = 0; entry < entries; entry++) {
    //MC->GetEntry(entry);
    //event->GetEntry(entry);
    jets->GetEntry(entry);
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
    
    if (numComb) {
      possibleMatchableEntries++;
    }

    // Ensure gen particles are matched to one reco jets
    if (!isMatched(genLeptonicBIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets) ||
        !isMatched(genHadronicBIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets) ||
        !isMatched(genFirstJetIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets)  ||
        !isMatched(genSecondJetIndex, m_jet_MCIndex, m_jet_lorentzvector, n_jets))
      continue;

    h_qgLikelihood_hadronicBJets->Fill(m_jet_qgtag_likelihood[genHadronicBIndex]);
    h_qgLikelihood_leptonicBJets->Fill(m_jet_qgtag_likelihood[genLeptonicBIndex]);

    h_qgLikelihood_lightJets->Fill(m_jet_qgtag_likelihood[genFirstJetIndex]);
    h_qgLikelihood_lightJets->Fill(m_jet_qgtag_likelihood[genSecondJetIndex]);

    h_qgLikelihood_allTTJets->Fill(m_jet_qgtag_likelihood[genHadronicBIndex]);
    h_qgLikelihood_allTTJets->Fill(m_jet_qgtag_likelihood[genLeptonicBIndex]);
    h_qgLikelihood_allTTJets->Fill(m_jet_qgtag_likelihood[genFirstJetIndex]);
    h_qgLikelihood_allTTJets->Fill(m_jet_qgtag_likelihood[genSecondJetIndex]);

    for (int i = 0; i < n_jets; i++) {
      h_qgLikelihood_allJets->Fill(m_jet_qgtag_likelihood[i]);
      if (i != genLeptonicBIndex && i != genHadronicBIndex && i != genFirstJetIndex && i != genSecondJetIndex) {
        h_qgLikelihood_allNoTTJets->Fill(m_jet_qgtag_likelihood[i]);
      }
    } 

    jets_selectedEntries++;
    //jets_selectedEntries_vs_mtt->Fill(MC_mtt);

    //// Use the four leading jet as selected jets
    //computeEfficiencyFourLeadingJets(0, 1, 2, 3);

    //if (numComb < 1)
      //continue;

    matchableEntries++;

    chi2_selectedEntries++;
    //chi2_selectedEntries_vs_mtt->Fill(MC_mtt);

    //computeEfficiencyChiSquare(selectedLeptonicBIndex, selectedHadronicBIndex, selectedFirstJetIndex, selectedSecondJetIndex);
  }

/*  std::cout << "Number of total entries: " << entries << std::endl;*/
  //std::cout << "Number of selected entries: " << chi2_selectedEntries << std::endl << std::endl;

  //std::cout << "Matchable efficiency" << std::endl;
  //printEff((double) matchableEntries / (possibleMatchableEntries), possibleMatchableEntries);

  //std::cout << std::endl << "Chi square" << std::endl;
  //std::cout << "Matched and well placed: 1 / 2 / 3 / 4 jets" << std::endl;
  //printEff((double) chi2_matchedAndWellPlacedEntries1 / chi2_selectedEntries, chi2_selectedEntries);
  //printEff((double) chi2_matchedAndWellPlacedEntries2 / chi2_selectedEntries, chi2_selectedEntries);
  //printEff((double) chi2_matchedAndWellPlacedEntries3 / chi2_selectedEntries, chi2_selectedEntries);
  //printEff((double) chi2_matchedAndWellPlacedEntries4 / chi2_selectedEntries, chi2_selectedEntries);

  //std::cout << std::endl << "Matched: 1 / 2 / 3 / 4 jets" << std::endl;
  //printEff((double) chi2_matchedEntries1 / chi2_selectedEntries, chi2_selectedEntries);
  //printEff((double) chi2_matchedEntries2 / chi2_selectedEntries, chi2_selectedEntries);
  //printEff((double) chi2_matchedEntries3 / chi2_selectedEntries, chi2_selectedEntries);
  //printEff((double) chi2_matchedEntries4 / chi2_selectedEntries, chi2_selectedEntries);

  //std::cout << std::endl << "Four leading jets" << std::endl;
  //std::cout << std::endl << "Matched: 1 / 2 / 3 / 4 jets" << std::endl;
  //printEff((double) jets_matchedEntries1 / jets_selectedEntries, jets_selectedEntries);
  //printEff((double) jets_matchedEntries2 / jets_selectedEntries, jets_selectedEntries);
  //printEff((double) jets_matchedEntries3 / jets_selectedEntries, jets_selectedEntries);
  /*printEff((double) jets_matchedEntries4 / jets_selectedEntries, jets_selectedEntries);*/

  delete MC; delete event; delete MET; delete jets; delete muons; delete electrons;

  if (outputFile.length() > 0) {
    TFile* f = TFile::Open(outputFile.c_str(), "recreate");
/*    chi2_efficiency_vs_mtt->Divide(chi2_matchedEntries_vs_mtt, chi2_selectedEntries_vs_mtt);*/
    //jets_efficiency_vs_mtt->Divide(jets_matchedEntries_vs_mtt, jets_selectedEntries_vs_mtt);

    //// Compute binomial errors
    //for (int i = 1; i < nBins; i++) {
      //chi2_efficiency_vs_mtt->SetBinError(i, binomialError(chi2_efficiency_vs_mtt->GetBinContent(i), chi2_selectedEntries_vs_mtt->GetBinContent(i)));
      //jets_efficiency_vs_mtt->SetBinError(i, binomialError(jets_efficiency_vs_mtt->GetBinContent(i), jets_selectedEntries_vs_mtt->GetBinContent(i)));
    //}

    //chi2_efficiency_vs_mtt->Write();
    //jets_efficiency_vs_mtt->Write();

    /*h_mtt_res_chi2->Write();*/

    h_qgLikelihood_hadronicBJets->Write();
    h_qgLikelihood_leptonicBJets->Write();
    h_qgLikelihood_lightJets->Write();
    h_qgLikelihood_allJets->Write();
    h_qgLikelihood_allTTJets->Write();
    h_qgLikelihood_allNoTTJets->Write();


    f->Close();
    delete f;
  }
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

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", false, "", "string", cmd);

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

    process(inputFiles, outputFileArg.getValue());    

  } catch (TCLAP::ArgException& e) {
    std::cout << e.what() << std::endl;
  }

}

