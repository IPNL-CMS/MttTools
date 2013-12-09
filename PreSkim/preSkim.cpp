#include "preSkim.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom2.h>
#include <vector>
#include <TProfile.h>
#include <fstream>
#include <memory>
#include <chrono>
#include <ctime>
#include <TLorentzVector.h>

#include "../PUReweighting/PUReweighter.h"
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

const int nBins = 15; 
const double bins[] = {340, 360, 380, 400, 420, 460, 500, 550, 600, 650, 750, 850, 950, 1050, 1200, 1400, 1600};

TLorentzVector* getP4(TClonesArray* array, int index) {
  return (TLorentzVector*) array->At(index);
}

void PreSkim::Loop()
{
  // Create the output file first in order that each histogram
  // associate itself with this file
  TFile * output = TFile::Open(mOutputFile.c_str(), "recreate");
  output->cd();

  TH1::SetDefaultSumw2(true);

  TH1D *hWeight = new TH1D("weight", "", 50, 0, 3);
  TH1D *hLeptonWeight = new TH1D("lepton_weight", "", 50, 0, 3);
  TH1D *hBTagWeight = new TH1D("btag_weight", "", 100, 0, 2);
  TH1D *hPUWeight = new TH1D("PU_weight", "", 50, 0, 3);
  TH1D *hGeneratorWeight = new TH1D("generator_weight", "", 100, -2, 2);

  TH1D *hIsSel = new TH1D("isSel", "", 10, 0, 10);
  TH1D *hBestSolChi2 = new TH1D("bestSolChi2", "", 400, 0, 1000);

  TH1D *hNTrueInt_nosel = new TH1D("nTrueInt_reco_nosel", "", 70, 0., 70);
  TH1D *hNTrueInt_nosel_noweight = new TH1D("nTrueInt_reco_nosel_noweight", "", 70, 0., 70);
  TH1D *hNVertex_nosel = new TH1D("nVertex_reco_nosel", "", 70, 0., 70);

  TH1D *hBoostTT_gen = new TH1D("boostTT_gen", "", 50, 0., 1.);
  TH1D *hPtTT_gen = new TH1D("ptTT_gen", "", 60, 0., 600.);
  TH1D *hEtaTT_gen = new TH1D("etaTT_gen", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hFirstJetPt_nosel = new TH1D("firstJetPt_reco_nosel", "", 100, 0., 640.);
  TH1D *hSecondJetPt_nosel = new TH1D("secondJetPt_reco_nosel", "", 100, 0., 620.);
  TH1D *hThirdJetPt_nosel = new TH1D("thirdJetPt_reco_nosel", "", 50, 0., 300.);
  TH1D *hFourthJetPt_nosel = new TH1D("fourthJetPt_reco_nosel", "", 50, 0., 300.);

  TH1D *hFirstJetEta_nosel = new TH1D("firstJetEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hSecondJetEta_nosel = new TH1D("secondJetEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hThirdJetEta_nosel = new TH1D("thirdJetEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hFourthJetEta_nosel = new TH1D("fourthJetEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hMET_nosel = new TH1D("MET_reco_nosel", "", 100, 0., 400.);

  TH1D *hNGoodJets_nosel = new TH1D("nGoodJets_reco_nosel", "", 10, 0, 10);
  TH1D *hNBtaggedJets_nosel = new TH1D("nBTaggedJets_reco_nosel", "", 6, 0, 6);

  TH1D *h_mtt_gen_nosel = new TH1D("mtt_gen_nosel", "", 300, 0., 1500.);

  TH1D *h_mtt_resolution = new TH1D("mtt_resolution", "", 100, -600., 600.);

  TH1D *hDeltaPhiTops_gen = new TH1D("deltaPhiTops_gen", "", 50, 0, M_PI);
  TH1D *hDeltaEtaTops_gen = new TH1D("deltaEtaTops_gen", "", 50, -3*M_PI, 3*M_PI);
  TH1D *hDeltaRTops_gen = new TH1D("deltaRTops_gen", "", 50, 0, 10);

  /******
   * Leptonic side
   */
  TH1* hDeltaPhiLeptonNeutrino_gen = new TH1D("deltaPhiLeptonNeutrino_gen", "", 50, 0, M_PI);
  TH1* hDeltaRLeptonNeutrino_gen = new TH1D("deltaRLeptonNeutrino_gen", "", 50, 0, 10);
  TH1* hDeltaEtaLeptonNeutrino_gen = new TH1D("deltaEtaLeptonNeutrino_gen", "", 50, -3*M_PI, 3*M_PI);

  // W 
  TH1* hLeptonicWPt_gen = new TH1D("ptWLeptonic_gen", "", 60, 0, 600.);
  TH1* hLeptonicWEta_gen = new TH1D("etaWLeptonic_gen", "", 50, -3*M_PI, 3*M_PI);

  /*****
   * Hadronic side
   */
  TH1* hDeltaPhiTwoLightJets_gen = new TH1D("deltaPhiTwoLightJets_gen", "", 50, 0, M_PI);
  TH1* hDeltaRTwoLightJets_gen = new TH1D("deltaRTwoLightJets_gen", "", 50, 0, 10);
  TH1* hDeltaEtaTwoLightJets_gen = new TH1D("deltaEtaTwoLightJets_gen", "", 50, -3*M_PI, 3.*M_PI);

  // W 
  TH1* hHadronicWPt_gen = new TH1D("ptWHadronic_gen", "", 60, 0, 600.);
  TH1* hHadronicWEta_gen = new TH1D("etaWHadronic_gen", "", 50, -3*M_PI, 3*M_PI);

  TH1* hDeltaPhiW_gen = new TH1D("deltaPhiW_gen", "", 50, 0, M_PI);
  TH1* hDeltaRW_gen = new TH1D("deltaRW_gen", "", 50, 0, 10);
  TH1* hDeltaEtaW_gen = new TH1D("deltaEtaW_gen", "", 50, -3*M_PI, 3.*M_PI);

  TH1* hHT_reco_nosel = new TH1D("HT_reco_nosel", "", 100, 0, 1000);
  TH1* hHT30_reco_nosel = new TH1D("HT30_reco_nosel", "", 100, 0, 1000);

  TH1* hHTFull_reco_nosel = new TH1D("HTFull_reco_nosel", "", 100, 0, 1000);

  hFirstJetPt_nosel->SetXTitle("1^{st} jet p_{T} [GeV/c]");

  h_mtt_gen_nosel->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  h_mtt_resolution->SetXTitle("m_{t#bar{t}}^{reco} - m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  Long64_t nentries = fMTT->GetEntries();

  //PUReweighter puReweighter(mIsSemiMu, mDataset);
  PUReweighter puReweighter(mIsSemiMu);

  std::cout << "Processing..." << std::endl;

  // Clone trees
  TTree* mtt_clone = fMTT->CloneTree(0);
  TTree* vertices_clone = fVertices->CloneTree(0);
  TTree* event_clone = fEvent->CloneTree(0);
  TTree* jet_clone = fJet->CloneTree(0);

  //TTree* electron_clone = fElectron->CloneTree(0);
  //TTree* electron_loose_clone = fElectronLoose->CloneTree(0);
  //TTree* muon_clone = fMuon->CloneTree(0);
  //TTree* muon_loose_clone = fMuonLoose->CloneTree(0);
  //TTree* mc_clone = fMC->CloneTree(0);

  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now(), end;
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if (jentry % 100000 == 0) {
      end = std::chrono::system_clock::now();
      float msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
      std::cout << "Processing entry #" << (jentry + 1) << " over " << nentries << " (" << (float) jentry / nentries * 100 << "% - " << msec / 1000. << "s)" << std::endl;
      start = std::chrono::system_clock::now();
    }

    GetEntry(jentry);

    // Compute event weight
    if (std::isnan(m_lepton_weight)) {
      std::cout << "Warning: lepton weight is NaN" << std::endl;
      m_lepton_weight = 1.;
    } else if (std::isinf(m_lepton_weight)) {
      std::cout << "Warning: lepton weight is Inf" << std::endl;
      m_lepton_weight = 1.;
    }

    if (std::isnan(m_btag_weight)) {
      std::cout << "Warning: btag weight is NaN" << std::endl;
      m_btag_weight = 1.;
    } else if (std::isinf(m_btag_weight)) {
      std::cout << "Warning: btag weight is Inf" << std::endl;
      m_btag_weight = 1.;
    }

    double eventWeight = 1.;
    if (mIsMC) {
      float puWeight = puReweighter.weight(n_trueInteractions);

      hLeptonWeight->Fill(m_lepton_weight);
      hBTagWeight->Fill(m_btag_weight);
      hGeneratorWeight->Fill(generator_weight);
      hPUWeight->Fill(puWeight);

      eventWeight *= puWeight;
      eventWeight *= generator_weight;
      eventWeight *= m_lepton_weight;
      eventWeight *= m_btag_weight;
    }

    if (std::isnan(eventWeight)) {
      std::cout << "Warning: event weight is NaN" << std::endl;
      eventWeight = 1.;
    }

    hWeight->Fill(eventWeight);

    // Fill gen value
    if (mIsMC && (MC_channel != 0)) {
      h_mtt_gen_nosel->Fill(MC_mtt, eventWeight);
      hBoostTT_gen->Fill(MC_boost_tt, eventWeight);

      hPtTT_gen->Fill(MC_pt_tt, eventWeight);
      hEtaTT_gen->Fill(MC_eta_tt, eventWeight);
      hDeltaPhiTops_gen->Fill(fabs(getP4(gen_top1_p4, 0)->DeltaPhi(*getP4(gen_top2_p4, 0))), eventWeight);
      hDeltaEtaTops_gen->Fill(getP4(gen_top1_p4, 0)->Eta() - getP4(gen_top2_p4, 0)->Eta(), eventWeight);
      hDeltaRTops_gen->Fill(getP4(gen_top1_p4, 0)->DeltaR(*getP4(gen_top2_p4, 0)), eventWeight);
    }

    if (mIsMC && (MC_channel == 1 || MC_channel == 2)) {
      TLorentzVector leptonic_W(0., 0., 0., 0.);
      if (gen_lepton_p4->GetEntriesFast() && gen_neutrino_p4->GetEntriesFast()) {
        hDeltaPhiLeptonNeutrino_gen->Fill(fabs(getP4(gen_lepton_p4, 0)->DeltaPhi(*getP4(gen_neutrino_p4, 0))), eventWeight);
        hDeltaRLeptonNeutrino_gen->Fill(getP4(gen_lepton_p4, 0)->DeltaR(*getP4(gen_neutrino_p4, 0)), eventWeight);
        hDeltaEtaLeptonNeutrino_gen->Fill(getP4(gen_lepton_p4, 0)->Eta() - getP4(gen_neutrino_p4, 0)->Eta(), eventWeight);

        leptonic_W = *getP4(gen_lepton_p4, 0) + *getP4(gen_neutrino_p4, 0);
        hLeptonicWPt_gen->Fill(leptonic_W.Pt(), eventWeight);
        hLeptonicWEta_gen->Fill(leptonic_W.Eta(), eventWeight);
      }

      TLorentzVector hadronic_W(0., 0., 0., 0.);
      if (gen_lightJet1_p4->GetEntriesFast() && gen_lightJet2_p4->GetEntriesFast()) {
        hDeltaPhiTwoLightJets_gen->Fill(fabs(getP4(gen_lightJet1_p4, 0)->DeltaPhi(*getP4(gen_lightJet2_p4, 0))), eventWeight);
        hDeltaRTwoLightJets_gen->Fill(getP4(gen_lightJet1_p4, 0)->DeltaR(*getP4(gen_lightJet2_p4, 0)), eventWeight);
        hDeltaEtaTwoLightJets_gen->Fill(getP4(gen_lightJet1_p4, 0)->Eta() - getP4(gen_lightJet2_p4, 0)->Eta(), eventWeight);

        hadronic_W = *getP4(gen_lightJet1_p4, 0) + *getP4(gen_lightJet2_p4, 0);
        hHadronicWPt_gen->Fill(leptonic_W.Pt(), eventWeight);
        hHadronicWEta_gen->Fill(leptonic_W.Eta(), eventWeight);
      }

      if (leptonic_W.Pt() != 0 && hadronic_W.Pt() != 0) {
        hDeltaPhiW_gen->Fill(fabs(leptonic_W.DeltaPhi(hadronic_W)), eventWeight);
        hDeltaRW_gen->Fill(leptonic_W.DeltaR(hadronic_W), eventWeight);
        hDeltaEtaW_gen->Fill(leptonic_W.Eta() - hadronic_W.Eta(), eventWeight);
      }
    }

    // If it's data, check that trigger's fired
    if (!mIsMC && !m_triggerPassed) {
      continue;
    }

    // Fill reco value without any selection
    hNTrueInt_nosel_noweight->Fill(n_trueInteractions);
    hNTrueInt_nosel->Fill(n_trueInteractions, eventWeight);
    hNVertex_nosel->Fill(n_vertices, eventWeight);
    hIsSel->Fill(isSel, eventWeight);
    hMET_nosel->Fill(MET, eventWeight);

    hNGoodJets_nosel->Fill(nJets, eventWeight);
    hNBtaggedJets_nosel->Fill(nBtaggedJets_CSVM, eventWeight);

    if (jet_p4->GetEntriesFast() > 0) {
      TLorentzVector* p4 = (TLorentzVector*) (*jet_p4)[0];
      hFirstJetPt_nosel->Fill(p4->Pt(), eventWeight);
      hFirstJetEta_nosel->Fill(p4->Eta(), eventWeight);
    }

    if (jet_p4->GetEntriesFast() > 1) {
      TLorentzVector* p4 = (TLorentzVector*) (*jet_p4)[1];
      hSecondJetPt_nosel->Fill(p4->Pt(), eventWeight);
      hSecondJetEta_nosel->Fill(p4->Eta(), eventWeight);
    }

    if (jet_p4->GetEntriesFast() > 2) {
      TLorentzVector* p4 = (TLorentzVector*) (*jet_p4)[2];
      hThirdJetPt_nosel->Fill(p4->Pt(), eventWeight);
      hThirdJetEta_nosel->Fill(p4->Eta(), eventWeight);
    }

    if (jet_p4->GetEntriesFast() > 3) {
      TLorentzVector* p4 = (TLorentzVector*) (*jet_p4)[3];
      hFourthJetPt_nosel->Fill(p4->Eta(), eventWeight);
      hFourthJetEta_nosel->Fill(p4->Pt(), eventWeight);
    }

    float HT = 0;
    float HT30 = 0;
    for (uint32_t i = 0; i < (uint32_t) jet_p4->GetEntriesFast(); i++) {
      float pt = ((TLorentzVector*) (*jet_p4)[i])->Pt();

      HT += pt;
      if (pt > 30)
        HT30 += pt;
    }
    hHT_reco_nosel->Fill(HT, eventWeight);
    hHT30_reco_nosel->Fill(HT30, eventWeight);

    double ptLepton = 0;
    double etaLepton = 0;
    double ptLeptonCut = 0;
    if (mIsSemiMu) {
      if (nGoodMuons > 0) {
        ptLepton = muonPt[0];
        etaLepton = muonEta[0];
        ptLeptonCut = 26.;
      }
    } else {
      if (nGoodElectrons > 0) {
        ptLepton = electronPt[0];
        etaLepton = electronEta[0];
        ptLeptonCut = 30.;
      }
    }

    float HTFull = HT + ptLepton + MET;

    hHTFull_reco_nosel->Fill(HTFull, eventWeight);

    if (isSel == 1 && numComb_chi2 > 0)
      hBestSolChi2->Fill(bestSolChi2, eventWeight);

    // We want at least on good lepton
    if (mIsSemiMu && nGoodMuons <= 0)
        continue;
    else if (!mIsSemiMu && nGoodElectrons <= 0)
        continue;

    if (!mIsSemiMu) {
      // The TOP reference selection exclude electron with
      // SuperCluster eta between 1.4442 and 1.5660
      // The TopTrigger efficiency does the same thing, but
      // using electron eta instead of SuperCluster eta.
      // Redo a cut here on electron eta
      // FIXME?
      if (fabs(etaLepton) >= 1.442 && fabs(etaLepton) < 1.5660)
        continue;
    }

    // isSel is 1 for a selected event.
    // Also check that at least one chi2 or mva combination has been done
    if (isSel != 1 || (numComb_chi2 <= 0 && numComb_MVA <= 0) || ptLepton < ptLeptonCut)
      continue;

    mtt_clone->Fill();
    event_clone->Fill();
    jet_clone->Fill();
    vertices_clone->Fill();

    //electron_loose_clone->Fill();
    //electron_clone->Fill();
    //muon_loose_clone->Fill();
    //muon_clone->Fill();
    //mc_clone->Fill();
  }

  output->Write();
  output->Close();
  delete output;
}

void PreSkim::loadChain(const std::vector<std::string>& inputFiles, const std::string& treeName, TChain*& output) {

  output = new TChain(treeName.c_str());

  for (const std::string& file: inputFiles) {
    output->Add(file.c_str());
  }

  output->SetCacheSize(30*1024*1024);
}

PreSkim::PreSkim(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isSemiMu, bool isMC) : fMTT(0), fVertices(0), fEvent(0)
{
  mIsSemiMu = isSemiMu;
  mIsMC = isMC;
  mOutputFile = outputFile;

  // Get trees
  loadChain(inputFiles, "Mtt", fMTT);
  loadChain(inputFiles, "Vertices", fVertices);
  loadChain(inputFiles, "event", fEvent);
  loadChain(inputFiles, "jet_PF", fJet);

  // Only for skimming
  //loadChain(inputFiles, "MC", fMC);
  //loadChain(inputFiles, "electron_PF", fElectron);
  //loadChain(inputFiles, "electron_loose_PF", fElectronLoose);
  //loadChain(inputFiles, "muon_PF", fMuon);
  //loadChain(inputFiles, "muon_loose_PF", fMuonLoose);

  //m_chains.insert({fElectronLoose, fElectron, fMuon, fMuonLoose, fMC});

  Init();
}

PreSkim::~PreSkim()
{
  if (fMTT)
    delete fMTT->GetCurrentFile();
}

Int_t PreSkim::GetEntry(Long64_t entry)
{

  for (auto& chain: m_chains) {
    chain->GetEntry(entry);
  }

  return 1;
}

void PreSkim::Init()
{
  fCurrent = -1;

  //fMTT->SetBranchStatus("*", 0);

  SetBranchAddress(fMTT, "MC_channel", &MC_channel);
  SetBranchAddress(fMTT, "MC_mtt", &MC_mtt);
  SetBranchAddress(fMTT, "MC_beta_tt", &MC_boost_tt);

  SetBranchAddress(fMTT, "MC_pt_tt", &MC_pt_tt);
  SetBranchAddress(fMTT, "MC_eta_tt", &MC_eta_tt);

  //SetBranchAddress(fMTT, "MC_top1Pt", &MC_top1Pt);
  //SetBranchAddress(fMTT, "MC_top2Pt", &MC_top2Pt);

  //SetBranchAddress(fMTT, "MC_top1Eta", &MC_top1Eta);
  //SetBranchAddress(fMTT, "MC_top2Eta", &MC_top2Eta);

  //SetBranchAddress(fMTT, "MC_pt_tt_com", &MC_pt_tt_com);
  //SetBranchAddress(fMTT, "MC_eta_tt_com", &MC_eta_tt_com);

  //SetBranchAddress(fMTT, "MC_top1Pt_com", &MC_top1Pt_com);
  //SetBranchAddress(fMTT, "MC_top2Pt_com", &MC_top2Pt_com);

  //SetBranchAddress(fMTT, "MC_top1Eta_com", &MC_top1Eta_com);
  //SetBranchAddress(fMTT, "MC_top2Eta_com", &MC_top2Eta_com);

  SetBranchAddress(fMTT, "nGoodMuons", &nGoodMuons);
  if (mIsSemiMu) {
    SetBranchAddress(fMTT, "muonPt", muonPt);
    SetBranchAddress(fMTT, "muonEta", muonEta);
  } else {
    SetBranchAddress(fMTT, "nGoodElectrons", &nGoodElectrons);
    SetBranchAddress(fMTT, "electronPt", &electronPt);
    SetBranchAddress(fMTT, "electronEta", &electronEta);
  }
  SetBranchAddress(fMTT, "1stjetpt", &p_1stjetpt);
  SetBranchAddress(fMTT, "2ndjetpt", &p_2ndjetpt);
  SetBranchAddress(fMTT, "3rdjetpt", &p_3rdjetpt);
  SetBranchAddress(fMTT, "4thjetpt", &p_4thjetpt);
  SetBranchAddress(fMTT, "nJets", &nJets);
  SetBranchAddress(fMTT, "jetEta", jetEta);
  //SetBranchAddress(fMTT, "jetPt", jetPt);
  SetBranchAddress(fMTT, "nBtaggedJets_CSVM", &nBtaggedJets_CSVM);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHEM", &nBtaggedJets_TCHEM);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHET", &nBtaggedJets_TCHET);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHPL", &nBtaggedJets_TCHPL);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHPM", &nBtaggedJets_TCHPM);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHPT", &nBtaggedJets_TCHPT);
  //SetBranchAddress(fMTT, "nBtaggedJets_SSVHEM", &nBtaggedJets_SSVHEM);
  //SetBranchAddress(fMTT, "nBtaggedJets_SSVHPT", &nBtaggedJets_SSVHPT);
  SetBranchAddress(fMTT, "MET", &MET);
  SetBranchAddress(fMTT, "isSel", &isSel);
  //SetBranchAddress(fMTT, "oneMatchedCombi", &oneMatchedCombi);
  SetBranchAddress(fMTT, "bestSolChi2", &bestSolChi2);
  //SetBranchAddress(fMTT, "isBestSolMatched", &isBestSolMatched);
  //SetBranchAddress(fMTT, "KFChi2", &KFChi2);
  SetBranchAddress(fMTT, "numComb_chi2", &numComb_chi2);
  SetBranchAddress(fMTT, "numComb_MVA", &numComb_MVA);
  //SetBranchAddress(fMTT, "solChi2", solChi2);

  SetBranchAddress(fMTT, "mLepTop_AfterChi2", &mLepTop_AfterChi2);
  SetBranchAddress(fMTT, "mHadTop_AfterChi2", &mHadTop_AfterChi2);

  SetBranchAddress(fMTT, "mLepTop_AfterChi2", &mLepTop_AfterChi2);
  SetBranchAddress(fMTT, "mHadTop_AfterChi2", &mHadTop_AfterChi2);
  SetBranchAddress(fMTT, "mtt_AfterChi2", &mtt_AfterChi2);
  SetBranchAddress(fMTT, "pt_tt_AfterChi2", &pt_tt_AfterChi2);
  SetBranchAddress(fMTT, "eta_tt_AfterChi2", &eta_tt_AfterChi2);
  SetBranchAddress(fMTT, "beta_tt_AfterChi2", &beta_tt_AfterChi2);
  SetBranchAddress(fMTT, "lepton_weight", &m_lepton_weight);
  SetBranchAddress(fMTT, "btag_weight", &m_btag_weight);
  SetBranchAddress(fMTT, "selectedLeptonIndex_AfterChi2", &selectedLeptonIndex_AfterChi2);

  if (fMTT->GetBranch("trigger_passed")) {
    SetBranchAddress(fMTT, "trigger_passed", &m_triggerPassed);
  } else {
    // Backward compatibilty
    m_triggerPassed = true;
  }

  gen_top1_p4 = NULL;
  gen_top2_p4 = NULL;
  SetBranchAddress(fMTT, "MC_Top1_p4", &gen_top1_p4);
  SetBranchAddress(fMTT, "MC_Top2_p4", &gen_top2_p4);

  gen_lepton_p4 = NULL;
  gen_neutrino_p4 = NULL;
  SetBranchAddress(fMTT, "MC_lepton_p4", &gen_lepton_p4);
  SetBranchAddress(fMTT, "MC_neutrino_p4", &gen_neutrino_p4);

  gen_leptonic_B_p4 = NULL;
  gen_hadronic_B_p4 = NULL;
  SetBranchAddress(fMTT, "MC_leptonic_B_p4", &gen_leptonic_B_p4);
  SetBranchAddress(fMTT, "MC_hadronic_B_p4", &gen_hadronic_B_p4);

  gen_lightJet1_p4 = NULL;
  gen_lightJet2_p4 = NULL;
  SetBranchAddress(fMTT, "MC_lightJet1_B_p4", &gen_lightJet1_p4);
  SetBranchAddress(fMTT, "MC_lightJet2_B_p4", &gen_lightJet2_p4);

  //fVertices->SetBranchStatus("*", 0);
  SetBranchAddress(fVertices, "n_vertices", &n_vertices);

  //fEvent->SetBranchStatus("*", 0);
  SetBranchAddress(fEvent, "nTrueInteractions", &n_trueInteractions);
  generator_weight = 1;
  SetBranchAddress(fEvent, "generator_weight", &generator_weight);

  jet_p4 = NULL;
  //fJet->SetBranchStatus("*", 0);
  SetBranchAddress(fJet, "jet_4vector", &jet_p4);
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

    TCLAP::CmdLine cmd("Skim extractor tuples with basic selection", ' ', "0.1");

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
    else if (p == "s10")
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
    
    PreSkim skimmer(inputFiles, outputFileArg.getValue(), semimuArg.isSet(), !isData);
    skimmer.Loop();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
