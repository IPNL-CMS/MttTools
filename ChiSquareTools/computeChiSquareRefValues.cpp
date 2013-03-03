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
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <TH1F.h>
#include <TF1.h>
#include <TROOT.h>
#include <TCut.h>

#include "KinFit.h"
#include <tclap/CmdLine.h>

bool CSV_MODE = false;

//Breit-Wigner function
double bw(double* x, double* par)
{
  double arg1 = 14.0/22.0; // 2 over pi
  double arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  double arg3 = ((x[0]*x[0]) - (par[1]*par[1]))*((x[0]*x[0]) - (par[1]*par[1]));
  double arg4 = x[0]*x[0]*x[0]*x[0]*((par[2]*par[2])/(par[1]*par[1]));

  return par[0]*arg1*arg2/(arg3 + arg4);
}

void loadChain(const std::vector<std::string>& inputFiles, TChain*& mc, TChain*& event, TChain*& jets, TChain*& MET, TChain*& muons, TChain*& electrons, TChain*& mtt) {
  std::cout << "Opening files..." << std::endl;

  mc = new TChain("MC");
  event = new TChain("event");
  jets = new TChain("jet_PF");
  MET = new TChain("MET_PF");
  muons = new TChain("muon_PF");
  electrons = new TChain("electron_PF");
  //mtt = new TChain("Mtt");

  for (const std::string& file: inputFiles) {
    mc->Add(file.c_str());
    event->Add(file.c_str());
    jets->Add(file.c_str());
    MET->Add(file.c_str());
    muons->Add(file.c_str());
    electrons->Add(file.c_str());
    //mtt->Add(file.c_str());
  }

  std::cout << "... done." << std::endl;
}

void SetBranchAddress(TChain* chain, const char* branchName, void* address) {
  chain->SetBranchStatus(branchName, 1);
  chain->SetBranchAddress(branchName, address, NULL);
}

#define ID_B (5)
#define ID_T (6)

#define ID_E (11)
#define ID_NEUTRINO_E (12)
#define ID_MU (13)
#define ID_NEUTRINO_MU (14)
#define ID_TAU (15)
#define ID_NEUTRINO_TAU (16)

#define ID_W (24)

void process(const std::vector<std::string>& inputFiles, const std::string& outputFile) {

  gROOT->SetBatch(true);

  TChain* MC = NULL, *event = NULL, *jets = NULL, *MET = NULL, *muons = NULL, *electrons = NULL, *mtt = NULL;
  loadChain(inputFiles, MC, event, jets, MET, muons, electrons, mtt);

  MC->SetBranchStatus("*", 0);
  event->SetBranchStatus("*", 0);
  jets->SetBranchStatus("*", 0);
  MET->SetBranchStatus("*", 0);
  muons->SetBranchStatus("*", 0);
  electrons->SetBranchStatus("*", 0);

  static const int 	m_MCs_MAX        = 1000;

  TClonesArray *m_MC_lorentzvector;
  int   m_n_MCs;
  //int   m_MC_index[m_MCs_MAX];
  //int   m_MC_status[m_MCs_MAX];
  int   m_MC_type[m_MCs_MAX];
  int   m_MC_imot1[m_MCs_MAX];
 /* int   m_MC_imot2[m_MCs_MAX];
  int   m_MC_generation[m_MCs_MAX];
  float m_MC_E[m_MCs_MAX];
  float	m_MC_px[m_MCs_MAX];
  float	m_MC_py[m_MCs_MAX];
  float	m_MC_pz[m_MCs_MAX];
  float	m_MC_vx[m_MCs_MAX];
  float	m_MC_vy[m_MCs_MAX];
  float	m_MC_vz[m_MCs_MAX];
  float	m_MC_eta[m_MCs_MAX];
  float	m_MC_phi[m_MCs_MAX];*/

  m_MC_lorentzvector = new TClonesArray("TLorentzVector");

  SetBranchAddress(MC, "n_MCs",  &m_n_MCs);
  SetBranchAddress(MC, "MC_4vector", &m_MC_lorentzvector);
  //SetBranchAddress(MC, "MC_index",   &m_MC_index);
  SetBranchAddress(MC, "MC_type",    &m_MC_type);
  SetBranchAddress(MC, "MC_mot1",    &m_MC_imot1);
  //SetBranchAddress(MC, "MC_mot2",    &m_MC_imot2);
  //SetBranchAddress(MC, "MC_generation",   &m_MC_generation);
  /*SetBranchAddress(MC, "MC_e",   &m_MC_E);
  SetBranchAddress(MC, "MC_px",  &m_MC_px);
  SetBranchAddress(MC, "MC_py",  &m_MC_py);
  SetBranchAddress(MC, "MC_pz",  &m_MC_pz);
  SetBranchAddress(MC, "MC_vx",  &m_MC_vx);
  SetBranchAddress(MC, "MC_vy",  &m_MC_vy);
  SetBranchAddress(MC, "MC_vz",  &m_MC_vz);
  SetBranchAddress(MC, "MC_eta", &m_MC_eta);
  SetBranchAddress(MC, "MC_phi", &m_MC_phi);*/

  // Jets
  static const int 	m_jets_MAX       = 200;
 
  TClonesArray* m_jet_lorentzvector = new TClonesArray("TLorentzVector");
  int  n_jets;
  int  m_jet_MCIndex[m_jets_MAX];

  SetBranchAddress(jets, "n_jets", &n_jets);
  SetBranchAddress(jets, "jet_4vector", &m_jet_lorentzvector);
  SetBranchAddress(jets, "jet_mcParticleIndex", &m_jet_MCIndex);

  // Electrons
  static const int 	m_electrons_MAX  = 100;
  TClonesArray* m_ele_lorentzvector = new TClonesArray("TLorentzVector");
  int n_electrons;
  int m_ele_MCIndex[m_electrons_MAX];
  
  SetBranchAddress(electrons, "n_electrons", &n_electrons);
  SetBranchAddress(electrons, "electron_4vector", &m_ele_lorentzvector);
  SetBranchAddress(electrons, "electron_mcParticleIndex", &m_ele_MCIndex);

  // Muons
  static const int 	m_muons_MAX  = 100;
  TClonesArray* m_muo_lorentzvector = new TClonesArray("TLorentzVector");
  int n_muons;
  int m_muo_MCIndex[m_muons_MAX];

  SetBranchAddress(muons, "n_muons",  &n_muons);
  SetBranchAddress(muons, "muon_4vector", &m_muo_lorentzvector);
  SetBranchAddress(muons, "muon_mcParticleIndex", &m_muo_MCIndex);

  // MET
  TClonesArray* m_met_lorentzvector = new TClonesArray("TLorentzVector");
  SetBranchAddress(MET, "met_4vector", &m_met_lorentzvector);

  // Histograms
  TH1* h_deltaPtFirstJet  = new TH1F("deltapt_firstjet", "Delta pt first jet", 70, -1.5, 1.5);
  TH1* h_deltaPtSecondJet = new TH1F("deltapt_secondjet", "Delta pt second jet", 70, -1.5, 1.5);

  TH1* h_hadronicWMass = new TH1F("w_mass", "hadronic w mass", 70, 30, 140);
  TH1* h_hadronicWMass_large = new TH1F("w_mass_large", "hadronic w mass", 200, 0, 600);
  TH1* h_hadronicTopMass = new TH1F("hadronic_top_mass", "hadronic top mass", 100, 100, 300);
  TH1* h_hadronicTopMass_large = new TH1F("hadronic_top_mass_large", "hadronic top mass", 200, 0, 800);

  TH1* h_hadronicWMassMC = new TH1F("w_mass_mc", "hadronic w mass (gen)", 80, 50, 110);
  TH1* h_hadronicTopMassMC = new TH1F("hadronic_top_mass_mc", "hadronic top mass (gen)", 100, 150, 190);

  TH1* h_leptonicTopMassMu = new TH1F("leptonic_top_mass_mu", "leptonic top mass - mu channel", 60, 100, 500);
  TH1* h_leptonicTopMassMu_large = new TH1F("leptonic_top_mass_mu_large", "leptonic top mass - mu channel", 200, 0, 800);

  TH1* h_leptonicTopMassE = new TH1F("leptonic_top_mass_e", "leptonic top mass - e channel", 60, 100, 500);
  TH1* h_leptonicTopMassE_large = new TH1F("leptonic_top_mass_e_large", "leptonic top mass - e channel", 200, 0, 800);
  
  TH1* h_ptSystem  = new TH1F("pt_system", "tt system pt", 100, 0, 200);
  TH1* h_ptSystem_large  = new TH1F("pt_system_large", "tt system pt", 200, 0, 800);

  TH1* h_htFrac  = new TH1F("ht_frac", "HT frac", 50, 0, 1);

  // Wrong combinaison plot
  TH1* h_hadronicWMass_wrong = new TH1F("w_mass_wrong", "hadronic w mass (wrong selection)", 200, 0, 600);
  TH1* h_hadronicTopMass_wrong = new TH1F("hadronic_top_mass_wrong", "hadronic top mass (wrong selection)", 200, 0, 800);
  TH1* h_leptonicTopMassMu_wrong = new TH1F("leptonic_top_mass_mu_wrong", "leptonic top mass - mu channel (wrong selection)", 200, 0, 800);
  TH1* h_leptonicTopMassE_wrong = new TH1F("leptonic_top_mass_e_wrong", "leptonic top mass - e channel (wrong selection)", 200, 0, 800);
  TH1* h_ptSystem_wrong  = new TH1F("pt_system_wrong", "tt system pt (wrong selection)", 200, 0, 800);
  TH1* h_htFrac_wrong  = new TH1F("ht_frac_wrong", "HT frac (wrong selection)", 50, 0, 1);


  uint64_t selectedEntries = 0;
  uint64_t selectedEntriesBeforeMatching = 0;
  uint64_t selectedEntriesAfterMatching = 0;
  uint64_t selectedSemiMuEntries = 0;

  uint64_t entries = MC->GetEntries();
  //uint64_t entries = 0;

  for (uint64_t entry = 0; entry < entries; entry++) {
    MC->GetEntry(entry);
    event->GetEntry(entry);
    jets->GetEntry(entry);
    MET->GetEntry(entry);
    muons->GetEntry(entry);
    electrons->GetEntry(entry);

    if (((entry + 1) % 100000) == 0) {
      std::cout << "Processing entry " << entry + 1 << " out of " << entries << std::endl;
    }

    bool keepEvent = true;

    // Indexes
    bool isSemiMu = false;
    int leptonIndex = -1;
    int neutrinoIndex = -1;

    int leptonicBIndex = -1;
    int hadronicBIndex = -1;
    int leptonicTopIndex = -1;

    int firstJetIndex = -1;
    int secondJetIndex = -1;

    if (false) {
      std::cout << "New event" << std::endl;
      for (int i = 0; i < m_n_MCs; i++) {
        std::cout << "\tType: " << m_MC_type[i] << std::endl;
      }
    }

    for (int i = 0; i < m_n_MCs; i++) {
      // We are only interested in particles in final states:
      // - leptons / neutrinos & jets
      // All theses particles have for mother a W and for grand mother a top
      if (m_MC_imot1[i] == -1 || m_MC_imot1[m_MC_imot1[i]] == -1)
        continue;

      if (fabs(m_MC_type[m_MC_imot1[i]]) == ID_T ||
          fabs(m_MC_type[m_MC_imot1[m_MC_imot1[i]]]) == ID_T)  {


        // W? Continue
        if (fabs(m_MC_type[i]) == ID_W)
          continue;

        // Only semi-mu or semi-e events are interesting, so throw away event with a tau
        if (fabs(m_MC_type[i]) == ID_TAU) {
          keepEvent = false;
          break;
        }

        switch ((int) fabs(m_MC_type[i])) {
          case ID_E:
            if (leptonIndex != -1) {
              keepEvent = false;
              break;
            }
            isSemiMu = false;
            leptonIndex = i;
            break;

          case ID_MU:
            if (leptonIndex != -1) {
              keepEvent = false;
              break;
            }
            isSemiMu = true;
            leptonIndex = i;
            break;

          case ID_NEUTRINO_E:
          case ID_NEUTRINO_MU:
          case ID_NEUTRINO_TAU:
            if (neutrinoIndex != -1) {
              keepEvent = false;
              break;
            }
            neutrinoIndex = i;

            leptonicTopIndex = m_MC_imot1[m_MC_imot1[i]];
            break;

          case ID_B:
            if (leptonicBIndex == -1) {
              leptonicBIndex = i;
            } else {
              if (hadronicBIndex != -1) {
                keepEvent = false;
                break;
              }
              hadronicBIndex = i;
            }
            break;

          default: // Other jets
            if (firstJetIndex == -1) {
              firstJetIndex = i;
            } else {
              if (secondJetIndex != -1) {
                keepEvent = false;
                break;
              }
              secondJetIndex = i;
            }
            break;
        }

        if (! keepEvent)
          break;
      }
    }

    if (leptonIndex == -1 || neutrinoIndex == -1 || leptonicBIndex == -1 || hadronicBIndex == -1 || firstJetIndex == -1 || secondJetIndex == -1)
      keepEvent = false;

    if (! keepEvent)
      continue;

    // Reorder B jet indexes
    if (m_MC_imot1[leptonicBIndex] != leptonicTopIndex) {
      // Wrong combinaison, swap
      std::swap(leptonicBIndex, hadronicBIndex);
    }

    if (false) {
      std::cout << "Event selected" << std::endl;
      std::cout << "\tSemi-mu: " << isSemiMu << std::endl;
      std::cout << "\tLepton type: " << m_MC_type[leptonIndex] << std::endl;
      std::cout << "\tNeutrino type: " << m_MC_type[neutrinoIndex] << std::endl;
      std::cout << "\tLeptonic B type: " << m_MC_type[leptonicBIndex] << std::endl;
      std::cout << "\tHadronic B type: " << m_MC_type[hadronicBIndex] << std::endl;
      std::cout << "\tFirst jet type: " << m_MC_type[firstJetIndex] << std::endl;
      std::cout << "\tSecond jet type: " << m_MC_type[secondJetIndex] << std::endl;
    }

    if (fabs(m_MC_type[leptonicBIndex]) != ID_B ||
        fabs(m_MC_type[hadronicBIndex]) != ID_B) {
      std::cout << "ERROR: Faulty algorithm" << std::endl;
      continue;
    }

    // Fill GEN histograms
    TLorentzVector* genFirstJetP4 = static_cast<TLorentzVector*>((*m_MC_lorentzvector)[firstJetIndex]);
    TLorentzVector* genSecondJetP4 = static_cast<TLorentzVector*>((*m_MC_lorentzvector)[secondJetIndex]);
    TLorentzVector* genHadronicBJetP4 = static_cast<TLorentzVector*>((*m_MC_lorentzvector)[hadronicBIndex]);

    double genWMass = (*genFirstJetP4 + *genSecondJetP4).M();
    double genHadrTopMass = (*genFirstJetP4 + *genSecondJetP4 + *genHadronicBJetP4).M();

    h_hadronicWMassMC->Fill(genWMass);
    h_hadronicTopMassMC->Fill(genHadrTopMass);

    // MET
    TLorentzVector* neutrinoP4 = static_cast<TLorentzVector*>((*m_met_lorentzvector)[0]);
    if (neutrinoP4->Pt() < 20)
      continue;

    // Lepton
    TLorentzVector* leptonP4 = NULL;
    float ptLeptonCut = 0;
    if (isSemiMu) {
      if (n_muons != 1)
        continue;

      leptonP4 = static_cast<TLorentzVector*>((*m_muo_lorentzvector)[0]);
      ptLeptonCut = 25;
    } else {
      if (n_electrons != 1)
        continue;

      leptonP4 = static_cast<TLorentzVector*>((*m_ele_lorentzvector)[0]);
      ptLeptonCut = 30;
    }

    if (leptonP4->Pt() < ptLeptonCut || fabs(leptonP4->Eta()) > 2.1)
      continue;

    selectedEntriesBeforeMatching++;

    // Found match between gen and reco
    int recoLeptonicBIndex = -1;
    int recoHadronicBIndex = -1;
    int recoFirstJetIndex = -1;
    int recoSecondJetIndex = -1;

    for (int i = 0; i < n_jets; i++) {
      int index = m_jet_MCIndex[i];
      if (index == leptonicBIndex)
        recoLeptonicBIndex = i;
      else if (index == hadronicBIndex)
        recoHadronicBIndex = i;
      else if (index == firstJetIndex)
        recoFirstJetIndex = i;
      else if (index == secondJetIndex)
        recoSecondJetIndex = i;
    }

    if (recoFirstJetIndex == -1 || recoSecondJetIndex == -1 || recoHadronicBIndex == -1 || recoLeptonicBIndex == -1)
      continue;

    //std::cout << "Njets: " << n_jets << "; Good combinaison: " << recoLeptonicBIndex << " " << recoHadronicBIndex << " " << recoFirstJetIndex << " " << recoSecondJetIndex << std::endl;

    for (int j1 = 0; j1 < n_jets; j1++) {
      //int j1_mcIndex = m_jet_MCIndex[j1];
      TLorentzVector* firstJetP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[j1]);
      if (firstJetP4->Pt() < 30)
        continue;

      if (j1 != recoLeptonicBIndex) {
        //std::cout << "Using jet #" << j1 << " for leptonic top mass" << std::endl;
        TLorentzVector neutrino = *neutrinoP4;
        KinFit::PzNeutrino(*leptonP4, neutrino, *firstJetP4);

        double topMass = (*leptonP4 + *firstJetP4 + neutrino).M();
        if (isSemiMu)
          h_leptonicTopMassMu_wrong->Fill(topMass);
        else
          h_leptonicTopMassE_wrong->Fill(topMass);
      }

      for (int j2 = 0; j2 < n_jets; j2++) {
        if (j1 == j2)
          continue;

        //int j2_mcIndex = m_jet_MCIndex[j2];
        TLorentzVector* secondJetP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[j2]);
        if (secondJetP4->Pt() < 30)
          continue;

        if ((j1 != recoFirstJetIndex && j2 != recoSecondJetIndex) ||
            (j1 != recoSecondJetIndex && j2 != recoFirstJetIndex)) {
          //std::cout << "Using jets #" << j1 << " " << j2 << " for hadronic W mass" << std::endl;
          h_hadronicWMass_wrong->Fill((*firstJetP4 + *secondJetP4).M());
        }


        for (int j3 = 0; j3 < n_jets; j3++) {
          if (j3 == j1 || j3 == j2)
            continue;

          TLorentzVector* thirdJetP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[j3]);
          if (thirdJetP4->Pt() < 30)
            continue;

          //int j3_mcIndex = m_jet_MCIndex[j3];

          if (!(((j1 == recoFirstJetIndex) || (j1 == recoSecondJetIndex) || (j1 == recoHadronicBIndex)) &&
              ((j2 == recoFirstJetIndex) || (j2 == recoSecondJetIndex) || (j2 == recoHadronicBIndex)) &&
              ((j3 == recoFirstJetIndex) || (j3 == recoSecondJetIndex) || (j3 == recoHadronicBIndex)))) {

            //std::cout << "Using jets #" << j1 << " " << j2 << " " << j3 << " for hadronic TOP mass" << std::endl;
            h_hadronicTopMass_wrong->Fill((*firstJetP4 + *secondJetP4 + *thirdJetP4).M());
          }

          for (int j4 = 0; j4 < n_jets; j4++) {
            if (j4 == j1 || j4 == j2 || j4 == j3)
              continue;

            TLorentzVector* fourthJetP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[j4]);
            if (fourthJetP4->Pt() < 30)
              continue;

            //int j4_mcIndex = m_jet_MCIndex[j4];

            if (!(((j1 == recoFirstJetIndex) || (j1 == recoSecondJetIndex) || (j1 == recoHadronicBIndex) || (j1 == recoLeptonicBIndex)) &&
                  ((j2 == recoFirstJetIndex) || (j2 == recoSecondJetIndex) || (j2 == recoHadronicBIndex) || (j2 == recoLeptonicBIndex)) &&
                  ((j3 == recoFirstJetIndex) || (j3 == recoSecondJetIndex) || (j3 == recoHadronicBIndex) || (j3 == recoLeptonicBIndex)) &&
                  ((j4 == recoFirstJetIndex) || (j4 == recoSecondJetIndex) || (j4 == recoHadronicBIndex) || (j4 == recoLeptonicBIndex)))) {

              //std::cout << "Using jets #" << j1 << " " << j2 << " " << j3 << " " << j4 << " for pt syst mass" << std::endl;
              TLorentzVector neutrino = *neutrinoP4;
              KinFit::PzNeutrino(*leptonP4, neutrino, *firstJetP4);

              h_ptSystem_wrong->Fill((neutrino + *leptonP4 + *firstJetP4 + *secondJetP4 + *thirdJetP4 + *fourthJetP4).Pt());

              double htFrac = firstJetP4->Pt() + secondJetP4->Pt() + thirdJetP4->Pt() + fourthJetP4->Pt();
              double denom = 0;
              for (int i = 0; i < n_jets; i++) {
                TLorentzVector* p4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[i]);
                if (p4->Pt() > 30) {
                  denom += p4->Pt();
                }
              }
              htFrac /= denom;
              h_htFrac_wrong->Fill(htFrac);
            }
          }
        }
      }
    }

    selectedEntriesAfterMatching++;

    // Selection. Pt > 70 50 30 30, eta < 2.4

    TLorentzVector* firstJetP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[recoFirstJetIndex]);
    TLorentzVector* secondJetP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[recoSecondJetIndex]);
    TLorentzVector* hadronicBP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[recoHadronicBIndex]);
    TLorentzVector* leptonicBP4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[recoLeptonicBIndex]);

    std::vector<double> pts = {firstJetP4->Pt(), secondJetP4->Pt(), hadronicBP4->Pt(), leptonicBP4->Pt()};
    std::sort(pts.begin(), pts.end());

    if (pts[3] < 70 || pts[2] < 50 || pts[1] < 30 || pts[0] < 30)
      continue;

    if (fabs(firstJetP4->Eta()) > 2.4 || fabs(secondJetP4->Eta()) > 2.4 || fabs(hadronicBP4->Eta()) > 2.4 || fabs(leptonicBP4->Eta()) > 2.4)
      continue;

    h_deltaPtFirstJet->Fill((firstJetP4->Pt() - genFirstJetP4->Pt()) / firstJetP4->Pt());
    h_deltaPtSecondJet->Fill((secondJetP4->Pt() - genSecondJetP4->Pt()) / secondJetP4->Pt());

    double hadronicWMass = (*firstJetP4 + *secondJetP4).M();
    h_hadronicWMass->Fill(hadronicWMass);
    h_hadronicWMass_large->Fill(hadronicWMass);

    double hadronicTopMass = (*firstJetP4 + *secondJetP4 + *hadronicBP4).M();
    h_hadronicTopMass->Fill(hadronicTopMass);
    h_hadronicTopMass_large->Fill(hadronicTopMass);

    KinFit::PzNeutrino(*leptonP4, *neutrinoP4, *leptonicBP4);
    double leptonicTopMass = (*leptonP4 + *neutrinoP4 + *leptonicBP4).M();

    if (isSemiMu) {
      h_leptonicTopMassMu->Fill(leptonicTopMass);
      h_leptonicTopMassMu_large->Fill(leptonicTopMass);
    } else {
      h_leptonicTopMassE->Fill(leptonicTopMass);
      h_leptonicTopMassE_large->Fill(leptonicTopMass);
    }

    // Pt TTbar system
    double ttPt = (*firstJetP4 + *secondJetP4 + *hadronicBP4 + *leptonicBP4 + *leptonP4 + *neutrinoP4).Pt();
    h_ptSystem->Fill(ttPt);
    h_ptSystem_large->Fill(ttPt);

    double htFrac = firstJetP4->Pt() + secondJetP4->Pt() + hadronicBP4->Pt() + leptonicBP4->Pt();
    double denom = 0;
    for (int i = 0; i < n_jets; i++) {
      TLorentzVector* p4 = static_cast<TLorentzVector*>((*m_jet_lorentzvector)[i]);
      if (p4->Pt() > 30) {
        denom += p4->Pt();
      }
    }
    htFrac /= denom;
    h_htFrac->Fill(htFrac);

    selectedEntries++;
    if (isSemiMu)
      selectedSemiMuEntries++;
  }

  std::cout << selectedEntriesBeforeMatching << " selected entries before MC matching of out " << entries << "(" << (float) selectedEntriesBeforeMatching / entries * 100 << "%)" << std::endl;
  std::cout << selectedEntriesAfterMatching << " selected entries after MC matching of out " << entries << "(" << (float) selectedEntriesAfterMatching / entries * 100 << "%)" << std::endl;
  std::cout << selectedEntries << " selected entries of out " << entries << "(" << (float) selectedEntries / entries * 100 << "%)" << std::endl;
  std::cout << (float) selectedSemiMuEntries / selectedEntries * 100 << "% of semi-mu entries" << std::endl;

  // Fit plots
  TF1* top_gaussian = new TF1("f1", "gaus", 150, 210);
  //TF1* top_gaussian = new TF1("f1", "gaus", 0, 2000);
  //TF1* top_gaussian = new TF1("f1", bw, 0, 2000, 3);
  top_gaussian->SetParameter(0, 500);
  top_gaussian->SetParameter(1, 180);
  top_gaussian->SetParameter(2, 40);

  TF1* top_gaussian_lept_e = new TF1("f3", "gaus", 130, 200);
  //TF1* top_gaussian_lept_e = new TF1("f3", bw, 120, 240, 3);
  top_gaussian_lept_e->SetParameter(0, 500);
  top_gaussian_lept_e->SetParameter(1, 180);
  top_gaussian_lept_e->SetParameter(2, 40);

  TF1* top_gaussian_lept_mu = new TF1("f4", "gaus", 130, 200);
  top_gaussian_lept_mu->SetParameter(1, 180);

  TF1* w_gaussian = new TF1("f2", "gaus", 70, 100);
  //TF1* w_gaussian = new TF1("f2", "gaus", 0, 2000);
  //TF1* w_gaussian = new TF1("f2", bw, 70, 100, 3);
  w_gaussian->SetParameter(0, 500);
  w_gaussian->SetParameter(1, 80);
  w_gaussian->SetParameter(2, 10);

  h_hadronicWMass->Fit(w_gaussian, "QR");
  h_hadronicTopMass->Fit(top_gaussian, "QR");
  h_leptonicTopMassMu->Fit(top_gaussian_lept_mu, "RQ");
  h_leptonicTopMassE->Fit(top_gaussian_lept_e, "RQ");

  std::cout << "Report:" << std::endl;
  std::cout <<"_______________" << std::endl << std::endl;
  if (! CSV_MODE) {
    std::cout << "hadronic W mass: " << w_gaussian->GetParameter(1) << " +/- " << w_gaussian->GetParameter(2) << std::endl;
    std::cout << "hadronic top mass: " << top_gaussian->GetParameter(1) << " +/- " << top_gaussian->GetParameter(2) << std::endl;
    std::cout << "leptonic top mass (semi-mu): " << top_gaussian_lept_mu->GetParameter(1) << " +/- " << top_gaussian_lept_mu->GetParameter(2) << std::endl;
    std::cout << "leptonic top mass (semi-e): " << top_gaussian_lept_e->GetParameter(1) << " +/- " << top_gaussian_lept_e->GetParameter(2) << std::endl;
  } else {
    std::cout << w_gaussian->GetParameter(1) << "\t" << w_gaussian->GetParameter(2) << std::endl;
    std::cout << top_gaussian->GetParameter(1) << "\t" << top_gaussian->GetParameter(2) << std::endl;
    std::cout << top_gaussian_lept_mu->GetParameter(1) << "\t" << top_gaussian_lept_mu->GetParameter(2) << std::endl;
    std::cout << top_gaussian_lept_e->GetParameter(1) << "\t" << top_gaussian_lept_e->GetParameter(2) << std::endl;
  }

  TFile* output = TFile::Open(outputFile.c_str(), "recreate");
  h_hadronicWMass->Write();
  h_hadronicWMass_large->Write();
  h_hadronicTopMass->Write();
  h_hadronicTopMass_large->Write();
  h_leptonicTopMassMu->Write();
  h_leptonicTopMassMu_large->Write();
  h_leptonicTopMassE->Write();
  h_leptonicTopMassE_large->Write();

  h_hadronicWMassMC->Write();
  h_hadronicTopMassMC->Write();

  h_deltaPtFirstJet->Write();
  h_deltaPtSecondJet->Write();

  h_ptSystem->Write();
  h_ptSystem_large->Write();

  h_htFrac->Write();

  h_hadronicWMass_wrong->Write();
  h_hadronicTopMass_wrong->Write();
  h_leptonicTopMassMu_wrong->Write();
  h_leptonicTopMassE_wrong->Write();
  h_ptSystem_wrong->Write();
  h_htFrac_wrong->Write();

  output->Close();
  
  delete output;
  delete MC; delete event; delete MET; delete jets; delete muons; delete electrons;

  delete top_gaussian; delete top_gaussian_lept_e; delete top_gaussian_lept_mu; delete w_gaussian;

  delete h_hadronicWMass;
  delete h_hadronicTopMass;
  delete h_leptonicTopMassMu;
  delete h_leptonicTopMassE;
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
    TCLAP::CmdLine cmd("compute various distribution to optimize chi square", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", true, "", "string", cmd);

    /*TCLAP::ValueArg<std::string> typeArg("", "type", "current inputfile type (semie or semimu)", false, "", "string", cmd);
      TCLAP::ValueArg<std::string> pileupArg("", "pileup", "PU profile used for MC production", false, "S7", "string", cmd);*/

    TCLAP::SwitchArg csvModeArg("", "csv", "CSV mode", cmd, false);

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
