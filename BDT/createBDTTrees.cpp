#include "createBDTTrees.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <TProfile.h>
#include <fstream>
#include <TLorentzVector.h>

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

/* Compute Neutrino Pz based on W mass */
bool createBDTTrees::computeNeutrinoPz(const TLorentzVector& bJet, const TLorentzVector& lepton, TLorentzVector& neutrino) const {
  if (lepton.E() == 0)
    return false;

  const double m_w = 8.04190000E+01; // Value used in Madgraph generation
  const double m_top = 172.5; // Value used in Madgraph generation

  double x = (m_w * m_w - lepton.M() * lepton.M() + 2. * (neutrino.Px() * lepton.Px() + neutrino.Py() * lepton.Py())) / (2 * lepton.E());
  double a = 1 - (lepton.Pz() * lepton.Pz()) / (lepton.E() * lepton.E());
  double b = -2. * (lepton.Pz() / lepton.E()) * x;
  double c = neutrino.Pt() * neutrino.Pt() - x * x;

  if (!a && !b)
    return false;

  if (!a) {
    neutrino.SetPz(-1 * c / b);
    neutrino.SetE(sqrt(neutrino.Px() * neutrino.Px() + neutrino.Py() * neutrino.Py() + neutrino.Pz() * neutrino.Pz()));

    return true;
  }

  double delta = b * b - 4 * a *c;

  if (delta < 0) {   // No solution, try to correct MET
    double rat = neutrino.Py() / neutrino.Px();

    double u = 4. / (lepton.E() * lepton.E()) * ((lepton.Px() + rat * lepton.Py()) * (lepton.Px() + rat * lepton.Py()) / (1 + rat * rat)
        - (lepton.E() * lepton.E()) + (lepton.Pz() * lepton.Pz()));

    double v = 4. / (lepton.E() * lepton.E()) * (m_w * m_w - lepton.M() * lepton.M())
      * (lepton.Px() + rat * lepton.Py()) / sqrt(1 + rat * rat);

    double w = (m_w * m_w - lepton.M() * lepton.M()) * (m_w * m_w - lepton.M() * lepton.M()) / (lepton.E() * lepton.E());

    double deltan = v * v - 4 * u * w;

    if (deltan < 0)
      return false; // Hopeless, MET can't be corrected

    double pt      = 0.;
    double corfact = 0.;

    if (u == 0) {
      pt = -w / v;
      if (pt <= 0)
        return false; // There is no way out...

      corfact = pt / neutrino.Pt();
    } else { // Deltan>=0 and u!=0
      double pt2 = (v - (sqrt(deltan))) / (2 * u);
      double pt1 = (v + (sqrt(deltan))) / (2 * u);

      // Pas de correction car negative
      if (pt1 <= 0 && pt2 <= 0) return 0;

      if (pt1 > 0 && pt2 < 0) pt = pt1;
      if (pt2 > 0 && pt1 < 0) pt = pt2;
      if (pt1 > 0 && pt2 > 0) {
        (fabs(pt1 - neutrino.Pt()) <= fabs(pt2 - neutrino.Pt()))
          ? pt = pt1
          : pt = pt2;
      }

      corfact = pt / neutrino.Pt();
    }

    // Now we have the correction factor

    neutrino.SetPx(corfact * neutrino.Px());
    neutrino.SetPy(corfact * neutrino.Py());

    // Recompute the new parameters

    x = (m_w * m_w - lepton.M() * lepton.M() + 2.*(neutrino.Px() * lepton.Px() + neutrino.Py() * lepton.Py())) / (2 * lepton.E());
    a = 1 - (lepton.Pz() * lepton.Pz()) / (lepton.E() * lepton.E());
    b = -2.*(lepton.Pz() / lepton.E()) * x;
    c = neutrino.Px() * neutrino.Px() + neutrino.Py() * neutrino.Py() - x * x;

    delta = b * b - 4 * a * c;

    if (fabs(delta) < 0.000001) delta = 0.;

    if (delta != 0)
      return false; // This should not happen, but who knows...
  }

  // We can go back to the normal path:

  TLorentzVector TopCand1 = lepton + bJet;
  TLorentzVector TopCand2 = lepton + bJet;

  neutrino.SetPz((-b - (sqrt(delta))) / (2 * a));
  neutrino.SetE(sqrt(neutrino.Px()*neutrino.Px() + neutrino.Py()*neutrino.Py() + neutrino.Pz()*neutrino.Pz()));
  TopCand1 += neutrino;

  neutrino.SetPz((-b + (sqrt(delta))) / (2 * a));
  neutrino.SetE(sqrt(neutrino.Px()*neutrino.Px() + neutrino.Py()*neutrino.Py() + neutrino.Pz()*neutrino.Pz()));
  TopCand2 += neutrino;

  double mtt_1 = sqrt(std::max(0., TopCand1.M2()));
  double mtt_2 = sqrt(std::max(0., TopCand2.M2()));

  if (fabs(mtt_1 - m_top) <= fabs(mtt_2 - m_top)) { // Otherwise it's OK
    neutrino.SetPz((-b - (sqrt(delta))) / (2 * a));
    neutrino.SetE(sqrt(neutrino.Px()*neutrino.Px() + neutrino.Py()*neutrino.Py() + neutrino.Pz()*neutrino.Pz()));
  }

  return true;
}

bool createBDTTrees::jetPassSelection(uint32_t jetIndex) {
  TLorentzVector* p4 = (TLorentzVector*) (*jet_p4)[jetIndex];

  if (p4->Pt() < 30.)
    return false;

  if (fabs(p4->Eta()) > 2.4)
    return false;

  return true;
}

void createBDTTrees::Loop()
{
  // Create the output file first in order that each histogram
  // associate itself with this file
  TFile * output = TFile::Open(mOutputFile.c_str(), "recreate");
  output->cd();

  TH1::SetDefaultSumw2(true);

  TTree signal_tree("signal", "signal tree");
  signal_tree.Branch("leptonic_B_Pt", &leptonic_B_Pt, "leptonic_B_Pt/F");
  signal_tree.Branch("hadronic_B_Pt", &hadronic_B_Pt, "hadronic_B_Pt/F");
  signal_tree.Branch("lightJet1p2_Pt", &lightJet1p2_Pt, "lightJet1p2_Pt/F");

  signal_tree.Branch("leptonic_W_Pt", &leptonic_W_Pt, "leptonic_W_Pt/F");
  signal_tree.Branch("hadronic_W_Pt", &hadronic_W_Pt, "hadronic_W_Pt/F");
  signal_tree.Branch("leptonic_W_M", &leptonic_W_M, "leptonic_W_M/F");
  signal_tree.Branch("hadronic_W_M", &hadronic_W_M, "hadronic_W_M/F");

  signal_tree.Branch("leptonic_Top_Pt", &leptonic_Top_Pt, "leptonic_Top_Pt/F");
  signal_tree.Branch("hadronic_Top_Pt", &hadronic_Top_Pt, "hadronic_Top_Pt/F");
  signal_tree.Branch("leptonic_Top_M", &leptonic_Top_M, "leptonic_Top_M/F");
  signal_tree.Branch("hadronic_Top_M", &hadronic_Top_M, "hadronic_Top_M/F");

  signal_tree.Branch("delta_phi_tops", &delta_phi_tops, "delta_phi_tops/F");
  signal_tree.Branch("delta_phi_lightjets", &delta_phi_lightjets, "delta_phi_lightjets/F");
  signal_tree.Branch("delta_phi_W", &delta_phi_W, "delta_phi_W/F");

  signal_tree.Branch("delta_R_tops", &delta_R_tops, "delta_R_tops/F");
  signal_tree.Branch("delta_R_lightjets", &delta_R_lightjets, "delta_R_lightjets/F");
  signal_tree.Branch("delta_R_W", &delta_R_W, "delta_R_W/F");

  signal_tree.Branch("weight", &signal_weight, "weight/F");

  TTree* background_tree = signal_tree.CloneTree();
  background_tree->SetName("background");
  background_tree->GetBranch("weight")->SetAddress(&background_weight);

  Long64_t nentries = fMTT->GetEntries();
  //nentries = 10000;

  PUReweighter puReweighter(mIsSemiMu, mPUProfile);

  std::cout << "Processing..." << std::endl;

  signal_weight = 1;

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if (jentry % 100000 == 0)
      std::cout << "Processing entry #" << (jentry + 1) << " over " << nentries << " (" << (float) jentry / nentries * 100 << "%)" << std::endl;

    GetEntry(jentry);

    double eventWeight = 1.;
    if (mIsMC) {
      eventWeight *= puReweighter.weight(n_trueInteractions);
    }

    if (isSel != 1) {
      // Keep only event passing basic selection
      continue;
    }

    if (! eventIsAssociable) {
      // There's no RECO solution in this event, so no "signal" possible
      continue;
    }

    // Loop over all jets
    // Flag as signal the one where reco jets match gen jets,
    // and as background the rest

    TLorentzVector* leptonP4 = nullptr;
    if (mIsSemiMu) {
      leptonP4 = (TLorentzVector*) (*muon_p4)[0];
    } else {
      leptonP4 = (TLorentzVector*) (*electron_p4)[0];
    }

    // Compute background weight
    uint32_t numberOfCombinaison = 0;

    for (uint32_t i = 0; i < n_jets; i++) {
      if (! jetPassSelection(i))
        continue;

      // Hadronic B
      for (uint32_t j = 0; j < n_jets; j++) {
        if (i == j)
          continue;

        if (! jetPassSelection(j))
          continue;

        // Light jet
        for (uint32_t k = 0; k < n_jets; k++) {
          if ((j == k) || (i == k))
            continue;

          if (! jetPassSelection(k))
            continue;

          // Light jet
          for (uint32_t l = k + 1; l < n_jets; l++) {
            if ((l == i) || (l == j) || (l == k))
              continue;

            if (! jetPassSelection(l))
              continue;

            TLorentzVector* leptonicBP4 = (TLorentzVector*) (*jet_p4)[i];
            TLorentzVector neutrinoP4 = *((TLorentzVector*) (*met_p4)[0]);
            if (! computeNeutrinoPz(*leptonicBP4, *leptonP4, neutrinoP4)) {
              continue;
            }

            numberOfCombinaison++;
          }
        }
      }
    }

    background_weight = 1. / (numberOfCombinaison - 1);

    background_weight *= eventWeight;
    signal_weight = eventWeight;

    if (numberOfCombinaison <= 0)
      continue;

    // Leptonic B
    for (uint32_t i = 0; i < n_jets; i++) {
      if (! jetPassSelection(i))
        continue;

      // Hadronic B
      for (uint32_t j = 0; j < n_jets; j++) {
        if (i == j)
          continue;

        if (! jetPassSelection(j))
          continue;

        // Light jet
        for (uint32_t k = 0; k < n_jets; k++) {
          if ((j == k) || (i == k))
            continue;

          if (! jetPassSelection(k))
            continue;

          // Light jet
          for (uint32_t l = k + 1; l < n_jets; l++) {
            if ((l == i) || (l == j) || (l == k))
              continue;

            if (! jetPassSelection(l))
              continue;

            TLorentzVector* leptonicBP4 = (TLorentzVector*) (*jet_p4)[i];
            TLorentzVector* hadronicBP4 = (TLorentzVector*) (*jet_p4)[j];
            TLorentzVector* lightJet1P4 = (TLorentzVector*) (*jet_p4)[k];
            TLorentzVector* lightJet2P4 = (TLorentzVector*) (*jet_p4)[l];

            leptonic_B_Pt = leptonicBP4->Pt();
            hadronic_B_Pt = hadronicBP4->Pt();

            lightJet1p2_Pt = lightJet1P4->Pt() + lightJet2P4->Pt();

            TLorentzVector hadronicWP4 = *lightJet1P4 + *lightJet2P4;
            hadronic_W_Pt = hadronicWP4.Pt();
            hadronic_W_M = hadronicWP4.M();

            TLorentzVector hadronicTopP4 = hadronicWP4 + *hadronicBP4;
            hadronic_Top_Pt = hadronicTopP4.Pt();
            hadronic_Top_M  = hadronicTopP4.M();

            TLorentzVector neutrinoP4 = *((TLorentzVector*) (*met_p4)[0]);
            TLorentzVector correctedNeutrinoP4 = neutrinoP4;

            if (! computeNeutrinoPz(*leptonicBP4, *leptonP4, correctedNeutrinoP4)) {
              continue;
              correctedNeutrinoP4 = neutrinoP4;
            }

            TLorentzVector leptonicWP4 = *leptonP4 + correctedNeutrinoP4;
            leptonic_W_Pt = leptonicWP4.Pt();
            leptonic_W_M = leptonicWP4.M();

            TLorentzVector leptonicTopP4 = leptonicWP4 + *leptonicBP4;
            leptonic_Top_Pt = leptonicTopP4.Pt();
            leptonic_Top_M = leptonicTopP4.M();

            delta_phi_tops = fabs(leptonicTopP4.DeltaPhi(hadronicTopP4));
            delta_phi_lightjets = fabs(lightJet1P4->DeltaPhi(*lightJet2P4));
            delta_phi_W = fabs(leptonicWP4.DeltaPhi(hadronicWP4));

            delta_R_tops = leptonicTopP4.DeltaR(hadronicTopP4);
            delta_R_lightjets = lightJet1P4->DeltaR(*lightJet2P4);
            delta_R_W = leptonicWP4.DeltaR(hadronicWP4);

            if (
               jet_mc_index[i] == MC_leptonicBIndex &&
               jet_mc_index[j] == MC_hadronicBIndex &&
               (( jet_mc_index[k] == MC_hadronicFirstJetIndex && jet_mc_index[l] == MC_hadronicSecondJetIndex ) ||
                ( jet_mc_index[k] == MC_hadronicSecondJetIndex && jet_mc_index[l] == MC_hadronicFirstJetIndex))
               )
            {
              signal_tree.Fill();
            } else {
              background_tree->Fill();
            }
          }
        }
      }
    }
  }

  output->Write();
  output->Close();
  delete output;
}

void loadChain(const std::vector<std::string>& inputFiles, const std::string& treeName, TChain*& output) {

  output = new TChain(treeName.c_str());

  for (const std::string& file: inputFiles) {
    output->Add(file.c_str());
  }
}

createBDTTrees::createBDTTrees(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isSemiMu, bool isMC, int btag, PUProfile puProfile) : fMTT(0), fEvent(0)
{
  mIsSemiMu = isSemiMu;
  mIsMC = isMC;
  mOutputFile = outputFile;
  mBTag = btag;
  mPUProfile = puProfile;

  // Get trees
  loadChain(inputFiles, "Mtt", fMTT);
  loadChain(inputFiles, "event", fEvent);
  loadChain(inputFiles, "muon_PF", fMuons);
  loadChain(inputFiles, "electron_PF", fElectrons);
  loadChain(inputFiles, "jet_PF", fJet);
  loadChain(inputFiles, "MET_PF", fMET);

  Init();
}

createBDTTrees::~createBDTTrees()
{
  if (fMTT)
    delete fMTT->GetCurrentFile();
}

Int_t createBDTTrees::GetEntry(Long64_t entry)
{
  if (fMTT)
    fMTT->GetEntry(entry);

  if (fEvent)
    fEvent->GetEntry(entry);

  if (fMuons)
    fMuons->GetEntry(entry);

  if (fElectrons)
    fElectrons->GetEntry(entry);

  if (fJet)
    fJet->GetEntry(entry);

  if (fMET)
    fMET->GetEntry(entry);

  return 1;
}

bool createBDTTrees::jetComesFromTTDecay(int mcIndex) const {
  return
    mcIndex == MC_leptonicBIndex ||
    mcIndex == MC_hadronicBIndex ||
    mcIndex == MC_hadronicFirstJetIndex ||
    mcIndex == MC_hadronicSecondJetIndex;
}

void createBDTTrees::SetBranchAddress(TTree* t, const char* branchName, void* ptr, TBranch** branch) {
  t->SetBranchStatus(branchName, 1);
  t->SetBranchAddress(branchName, ptr, branch);
}

void createBDTTrees::Init()
{
  fCurrent = -1;

  fMTT->SetBranchStatus("*", 0);

  SetBranchAddress(fMTT, "MC_channel", &MC_channel, NULL);
  SetBranchAddress(fMTT, "isSel", &isSel, NULL);
  SetBranchAddress(fMTT, "eventIsAssociable", &eventIsAssociable, NULL);

  SetBranchAddress(fMTT, "MC_leptonicBIndex", &MC_leptonicBIndex, NULL);
  SetBranchAddress(fMTT, "MC_hadronicBIndex", &MC_hadronicBIndex, NULL);
  SetBranchAddress(fMTT, "MC_hadronicFirstJetIndex", &MC_hadronicFirstJetIndex, NULL);
  SetBranchAddress(fMTT, "MC_hadronicSecondJetIndex", &MC_hadronicSecondJetIndex, NULL);

  fEvent->SetMakeClass(1);
  fEvent->SetBranchAddress("nTrueInteractions", &n_trueInteractions, NULL);
  fEvent->SetBranchStatus("*", 0);
  fEvent->SetBranchStatus("nTrueInteractions", 1);

  muon_p4 = NULL;
  fMuons->SetBranchStatus("*", 0);
  fMuons->SetBranchStatus("muon_4vector", 1);
  fMuons->SetBranchStatus("n_muons", 1);
  fMuons->SetBranchAddress("muon_4vector", &muon_p4, NULL);
  fMuons->SetBranchAddress("n_muons", &n_muons, NULL);

  electron_p4 = NULL;
  fElectrons->SetBranchStatus("*", 0);
  fElectrons->SetBranchStatus("electron_4vector", 1);
  fElectrons->SetBranchStatus("n_electrons", 1);
  fElectrons->SetBranchAddress("electron_4vector", &electron_p4, NULL);
  fElectrons->SetBranchAddress("n_electrons", &n_electrons, NULL);

  jet_p4 = NULL;
  fJet->SetBranchStatus("*", 0);
  fJet->SetBranchStatus("jet_4vector", 1);

  fJet->SetBranchStatus("n_jets", 1);
  fJet->SetBranchStatus("jet_mcParticleIndex", 1);

  fJet->SetBranchAddress("jet_4vector", &jet_p4, NULL);

  fJet->SetBranchAddress("n_jets", &n_jets);
  fJet->SetBranchAddress("jet_mcParticleIndex", &jet_mc_index);

  met_p4 = NULL;
  fMET->SetBranchAddress("met_4vector", &met_p4);
}

void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

int main(int argc, char** argv) {

  try {

    TCLAP::CmdLine cmd("Convert extractor tuples to histograms", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", true, "", "string", cmd);

    TCLAP::SwitchArg dataArg("", "data", "Is this data?", false);
    TCLAP::SwitchArg mcArg("", "mc", "Is this mc?", false);

    cmd.xorAdd(dataArg, mcArg);

    TCLAP::SwitchArg semimuArg("", "semimu", "Is this semi-mu channel?", false);
    TCLAP::SwitchArg semieArg("", "semie", "Is this semi-e channel?", false);

    cmd.xorAdd(semimuArg, semieArg);

    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jet to require", true, 2, "int", cmd);

    TCLAP::ValueArg<std::string> pileupArg("", "pileup", "PU profile used for MC production", false, "S10", "string", cmd);

    TCLAP::ValueArg<std::string> pileupSystArg("", "pileup-syst", "PU profile to use for pileup reweigthing", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> pdfSystArg("", "pdf-syst", "PDF systematic to compute", false, "nominal", "string", cmd);


    TCLAP::ValueArg<int> maxEntriesArg("n", "", "Maximal number of entries to process", false, -1, "int", cmd);
    TCLAP::ValueArg<double> generatorWeightArg("", "weight", "MC generator weight", false, 1., "double", cmd);

    cmd.parse(argc, argv);

    PUProfile puProfile;
    std::string p = pileupArg.getValue();
    std::transform(p.begin(), p.end(), p.begin(), ::tolower);
    if (p == "s6")
      puProfile = PUProfile::S6;
    else if (p == "s7")
      puProfile = PUProfile::S7;
    else
      puProfile = PUProfile::S10;

    std::string puSyst = pileupSystArg.getValue();
    std::transform(puSyst.begin(), puSyst.end(), puSyst.begin(), ::tolower);
    if (puSyst != "nominal" && puSyst != "up" && puSyst != "down") {
      std::cerr << "--pilup-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    std::string pdfSyst = pdfSystArg.getValue();
    std::transform(pdfSyst.begin(), pdfSyst.end(), pdfSyst.begin(), ::tolower);
    if (pdfSyst != "nominal" && pdfSyst != "up" && pdfSyst != "down") {
      std::cerr << "--pdf-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    bool isData = dataArg.isSet();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    createBDTTrees convertor(inputFiles, outputFileArg.getValue(), semimuArg.isSet(), !isData, btagArg.getValue(), puProfile);
    convertor.Loop();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
