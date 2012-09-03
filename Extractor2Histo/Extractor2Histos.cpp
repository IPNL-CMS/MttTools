#include "Extractor2Histos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>

#include "../PUReweighting/PUReweighter.h"

void Extractor2Histos::Loop()
{
  TH1D *hNVtx = new TH1D("hNVtx", "", 70, 0., 70);
  TH1D *hNVtx_finalsel = new TH1D("hNVtx_finalsel", "", 70, 0., 70);

  TH1D *hNTrueInt = new TH1D("hNTrueInt", "", 70, 0., 70);

  TH1D *hLeptonPt = new TH1D("hLeptonPt", "", 50, 20., 200.);
  TH1D *hLeptonPt_finalsel = new TH1D("hLeptonPt_finalsel", "", 50, 20., 200.);

  TH1D *hMuRelIso = new TH1D("hMuRelIso", "", 50, 0., 0.15);
  TH1D *hMuRelIso_finalsel = new TH1D("hMuRelIso_finalsel", "", 50, 0., 0.15);

  TH1D *h1stjetpt = new TH1D("h1stjetpt", "", 100, 70., 640.);
  TH1D *h1stjetpt_finalsel = new TH1D("h1stjetpt_finalsel", "", 100, 70., 640.);

  TH1D *h2ndjetpt = new TH1D("h2ndjetpt", "", 100, 50., 620.);
  TH1D *h2ndjetpt_finalsel = new TH1D("h2ndjetpt_finalsel", "", 100, 50., 620.);

  TH1D *h3rdjetpt = new TH1D("h3rdjetpt", "", 50, 30., 300.);
  TH1D *h3rdjetpt_finalsel = new TH1D("h3rdjetpt_finalsel", "", 50, 30., 300.);

  TH1D *h4thjetpt = new TH1D("h4thjetpt", "", 50, 30., 300.);
  TH1D *h4thjetpt_finalsel = new TH1D("h4thjetpt_finalsel", "", 50, 30., 300.);

  TH1D *hMET = new TH1D("hMET", "", 100, 20., 400.);
  TH1D *hMET_finalsel = new TH1D("hMET_finalsel", "", 100, 20., 400.);

  TH1D *hmtlep = new TH1D("hmtlep", "", 100, 120., 240.);
  TH1D *hmthad = new TH1D("hmthad", "", 150, 120., 300.);
  TH1D *hmtt = new TH1D("hmtt", "", 100, 300., 1300.);
  TH1D *hmttSelectedall = new TH1D("hmttSelectedall", "", 100, 300., 1300.);
  TH1D *hmttSelected0b = new TH1D("hmttSelected0b", "", 100, 300., 1300.);
  TH1D *hmttSelected1b = new TH1D("hmttSelected1b", "", 100, 300., 1300.);
  TH1D *hmttSelected2b = new TH1D("hmttSelected2b", "", 100, 300., 1300.);
  TH1D *hmttSelected2bMtt500 = new TH1D("hmttSelected2bMtt500", "", 100, 300., 1300.);
  TH1D *hmtt_all = new TH1D("hmtt_all", "", 100, 300., 1300.);
  TH1D *hmtt_bestchi2 = new TH1D("hmtt_bestchi2", "", 100, 300., 1300.);
  //TH1D *hNGoodMuons = new TH1D("hNGoodMuons", "", 5, -0.5, 4.5);

  TH1D *hNGoodJets = new TH1D("hNGoodJets", "", 6, 3.5, 9.5);
  TH1D *hNGoodJets_finalsel = new TH1D("hNGoodJets_finalsel", "", 6, 3.5, 9.5);

  TH1D *hNBtaggedJets = new TH1D("hNBtaggedJets", "", 5, -0.5, 4.5);
  TH1D *hNBtaggedJets_finalsel = new TH1D("hNBtaggedJets_finalsel", "", 5, -0.5, 4.5);

  if (mIsSemiMu) {
    hLeptonPt->SetXTitle("#mu p_{T} [GeV/c]");
    hLeptonPt_finalsel->SetXTitle("#mu p_{T} [GeV/c]");
  } else {
    hLeptonPt_finalsel->SetXTitle("e p_{T} [GeV/c]");
    hLeptonPt->SetXTitle("e p_{T} [GeV/c]");
  }

  h1stjetpt->SetXTitle("1^{st} jet p_{T} [GeV/c]");
  h1stjetpt_finalsel->SetXTitle("1^{st} jet p_{T} [GeV/c]");

  h2ndjetpt->SetXTitle("2^{nd} jet p_{T} [GeV/c]");
  h2ndjetpt_finalsel->SetXTitle("2^{nd} jet p_{T} [GeV/c]");

  h3rdjetpt->SetXTitle("3^{rd} jet p_{T} [GeV/c]");
  h3rdjetpt_finalsel->SetXTitle("3^{rd} jet p_{T} [GeV/c]");

  h4thjetpt_finalsel->SetXTitle("4^{th} jet p_{T} [GeV/c]");

  hMET->SetXTitle("MET [GeV]");
  hMET_finalsel->SetXTitle("MET [GeV]");

  hmtlep->SetXTitle("leptonic m_{t} [GeV/c^{2}]");
  hmthad->SetXTitle("hadronic m_{t} [GeV/c^{2}]");
  hmtt ->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmttSelectedall->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmttSelected0b->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmttSelected1b->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmttSelected2b->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmttSelected2bMtt500->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmtt_all->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmtt_bestchi2->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  //hNGoodMuons->SetXTitle("Num good #mu");
  hNGoodJets->SetXTitle("Num good jets");
  hNGoodJets_finalsel->SetXTitle("Num good jets");

  hNBtaggedJets->SetXTitle("Num TCHEL jets");
  hNBtaggedJets_finalsel->SetXTitle("Num TCHEL jets");

  Long64_t nentries = fMTT->GetEntries();

  PUReweighter puReweighter(mIsSemiMu, mDataset);

  std::cout << "Processing " << mInputFile << " ..." << std::endl;

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if (jentry % 100000 == 0)
      std::cout << "Processing entry #" << (jentry + 1) << " over " << nentries << " (" << (float) jentry / nentries * 100 << "%)" << std::endl;

    GetEntry(jentry);

    double ptLepton = -1;
    double ptLeptonCut = -1;
    if (mIsSemiMu)
    {
      if (nGoodMuons <= 0)
        continue;

      ptLepton = muonPt[0];
      ptLeptonCut = 25.;
    }
    else
    {
      if (nGoodElectrons <= 0 || hyperTight1MC[0] != 1 || elRelIso[0] >= 0.1)
        continue;

      ptLepton = electronPt[0];
      ptLeptonCut = 30.;
    }

    //FIXME: PU reweighting
    double eventWeight = 1.;
    if (mIsMC)
    {
      eventWeight *= puReweighter.weight(n_trueInteractions);
    }

    if (ptLepton > ptLeptonCut && nJets >= 4 && nBtaggedJets_CSVM >= 2) {
      hNGoodJets_finalsel->Fill(nJets, eventWeight);
      hNBtaggedJets_finalsel->Fill(nBtaggedJets_CSVM, eventWeight);

      hLeptonPt_finalsel->Fill(ptLepton, eventWeight);
      hMuRelIso_finalsel->Fill(muRelIso[0], eventWeight);
      
      hMET_finalsel->Fill(MET, eventWeight);

      h1stjetpt_finalsel->Fill(p_1stjetpt, eventWeight);
      h2ndjetpt_finalsel->Fill(p_2ndjetpt, eventWeight);
      h3rdjetpt_finalsel->Fill(p_3rdjetpt, eventWeight);
      h4thjetpt_finalsel->Fill(p_4thjetpt, eventWeight);

      hNVtx_finalsel->Fill(n_vertices, eventWeight);
    }

    if (ptLepton > ptLeptonCut && isSel == 1)
    {
      hNGoodJets->Fill(nJets, eventWeight);
      hNBtaggedJets->Fill(nBtaggedJets_CSVM, eventWeight);
      hNVtx->Fill(n_vertices, eventWeight);
      hNTrueInt->Fill(n_trueInteractions, eventWeight);
      if (nBtaggedJets_CSVM > 1 && p_1stjetpt > 70. && p_2ndjetpt > 50)
      {
        hLeptonPt->Fill(ptLepton, eventWeight);
        h1stjetpt->Fill(p_1stjetpt, eventWeight);
        h2ndjetpt->Fill(p_2ndjetpt, eventWeight);
        h3rdjetpt->Fill(p_3rdjetpt, eventWeight);
        h4thjetpt->Fill(p_4thjetpt, eventWeight);
        hMET->Fill(MET, eventWeight);
      }

      if (mtt_AfterChi2 > 0. && mtt_AfterChi2 < 998)
      {
        hmtt_bestchi2->Fill(mtt_AfterChi2, eventWeight);
      }

      if (mtt_AfterChi2 > 0. && bestSolChi2 < 500)
      {
        hmtlep->Fill(mLepTop_AfterChi2andKF, eventWeight);
        hmthad->Fill(mHadTop_AfterChi2andKF, eventWeight);
        hmtt->Fill(mtt_AfterChi2andKF, eventWeight);
        if (p_1stjetpt > 70. && p_2ndjetpt > 50)
        {
          //if (mtt_NBtaggedJets_TCHEL>1) hmttSelected2b->Fill(mtt_AfterChi2andKF,eventWeight);
          //if (mtt_NBtaggedJets_TCHEL==1) hmttSelected1b->Fill(mtt_AfterChi2andKF,eventWeight);
          //if (mtt_NBtaggedJets_TCHEL<1) hmttSelected0b->Fill(mtt_AfterChi2andKF,eventWeight);
          //hmttSelectedall->Fill(mtt_AfterChi2andKF,eventWeight);
          if (nBtaggedJets_CSVM > 1)
          {
            hMuRelIso->Fill(muRelIso[0], eventWeight);
            hmttSelected2b->Fill(mtt_AfterChi2, eventWeight);
          }

          if (nBtaggedJets_CSVM > 1 && mtt_AfterChi2 > 500)
            hmttSelected2bMtt500->Fill(mtt_AfterChi2, eventWeight);

          if (nBtaggedJets_CSVM == 1)
            hmttSelected1b->Fill(mtt_AfterChi2, eventWeight);

          if (nBtaggedJets_CSVM < 1)
            hmttSelected0b->Fill(mtt_AfterChi2, eventWeight);

          hmttSelectedall->Fill(mtt_AfterChi2, eventWeight);
        }
      }
    }
  }

  TFile * output = TFile::Open(mOutputFile, "recreate");
  output->cd();

  hNVtx->Write();
  hNVtx_finalsel->Write();

  hNTrueInt->Write();

  hLeptonPt->Write();
  hLeptonPt_finalsel->Write();

  hMuRelIso->Write();
  hMuRelIso_finalsel->Write();

  h1stjetpt->Write();
  h1stjetpt_finalsel->Write();

  h2ndjetpt->Write();
  h2ndjetpt_finalsel->Write();

  h3rdjetpt->Write();
  h3rdjetpt_finalsel->Write();

  h4thjetpt->Write();
  h4thjetpt_finalsel->Write();

  hMET_finalsel->Write();
  hMET->Write();

  hmtlep->Write();
  hmthad->Write();
  hmtt->Write();
  hmttSelectedall->Write();
  hmttSelected0b->Write();
  hmttSelected1b->Write();
  hmttSelected2b->Write();
  hmttSelected2bMtt500->Write();
  hmtt_all->Write();
  hmtt_bestchi2->Write();
  //hNGoodMuons->Write();
  hNGoodJets->Write();
  hNGoodJets_finalsel->Write();

  hNBtaggedJets->Write();
  hNBtaggedJets_finalsel->Write();

  output->Close();
  delete output;

}

Extractor2Histos::Extractor2Histos(TString fIn, TString fOut, const std::string& dataset, bool isSemiMu, bool isMC) : fMTT(0), fVertices(0), fEvent(0)
{
  mDataset = dataset;
  mIsSemiMu = isSemiMu;
  mIsMC = isMC;
  mOutputFile = fOut;
  mInputFile = fIn;

  TFile* f = TFile::Open(fIn);
  if (! f) {
    std::cerr << "Error: can't open " << fIn << std::endl;
    return;
  }

  // Get Mtt tree
  fMTT = static_cast<TTree*>(f->Get("Mtt"));

  // Get Vertices tree
  fVertices = static_cast<TTree*>(f->Get("Vertices"));

  // Get event tree
  fEvent = static_cast<TTree*>(f->Get("event"));

  Init();
}

Extractor2Histos::~Extractor2Histos()
{
  if (fMTT)
    delete fMTT->GetCurrentFile();

  /*if (fVertices)
    delete fVertices->GetCurrentFile();

    if (fEvent)
    delete fEvent->GetCurrentFile();*/
}

Int_t Extractor2Histos::GetEntry(Long64_t entry)
{
  if (fMTT)
    fMTT->GetEntry(entry);

  if (fVertices)
    fVertices->GetEntry(entry);

  if (fEvent)
    fEvent->GetEntry(entry);

  return 1;
}

void Extractor2Histos::SetBranchAddress(TTree* t, const char* branchName, void* ptr, TBranch** branch) {
  t->SetBranchStatus(branchName, 1);
  t->SetBranchAddress(branchName, ptr, branch);
}

void Extractor2Histos::Init()
{
  fCurrent = -1;
  fMTT->SetMakeClass(1);

  fMTT->SetBranchStatus("*", 0);
  //SetBranchAddress(fMTT, "MC_channel", &MC_channel, &b_MC_channel);
  //SetBranchAddress(fMTT, "MC_mtt", &MC_mtt, &b_MC_mtt);
  //SetBranchAddress(fMTT, "MC_nPU", &MC_nPU, &b_m_nPU);
  SetBranchAddress(fMTT, "nGoodMuons", &nGoodMuons, &b_nGoodMuons);
  //SetBranchAddress(fMTT, "nLooseGoodMuons", &nLooseGoodMuons, &b_nLooseGoodMuons);
  if (mIsSemiMu) {
    SetBranchAddress(fMTT, "muonPt", muonPt, &b_muonPt);
    //SetBranchAddress(fMTT, "2DDrMin", p_2DDrMin, &b_2DDrMin);
    //SetBranchAddress(fMTT, "2DpTrel", p_2DpTrel, &b_2DpTrel);
    SetBranchAddress(fMTT, "muRelIso", muRelIso, &b_muRelIso);
  } else {
    SetBranchAddress(fMTT, "nGoodElectrons", &nGoodElectrons, &b_nGoodElectrons);
    SetBranchAddress(fMTT, "electronPt", &electronPt, &b_electronPt);
    SetBranchAddress(fMTT, "elRelIso", &elRelIso, &b_elRelIso);
    SetBranchAddress(fMTT, "hyperTight1MC", &hyperTight1MC, &b_hyperTight1MC);
  }
  SetBranchAddress(fMTT, "1stjetpt", &p_1stjetpt, &b_1stjetpt);
  SetBranchAddress(fMTT, "2ndjetpt", &p_2ndjetpt, &b_2ndjetpt);
  SetBranchAddress(fMTT, "3rdjetpt", &p_3rdjetpt, &b_3rdjetpt);
  SetBranchAddress(fMTT, "4thjetpt", &p_4thjetpt, &b_4thjetpt);
  SetBranchAddress(fMTT, "nJets", &nJets, &b_nJets);
  //SetBranchAddress(fMTT, "jetEta", jetEta, &b_jetEta);
  SetBranchAddress(fMTT, "jetPt", jetPt, &b_jetPt);
  SetBranchAddress(fMTT, "nBtaggedJets_CSVM", &nBtaggedJets_CSVM, NULL);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHEM", &nBtaggedJets_TCHEM, &b_nBtaggedJets_TCHEM);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHET", &nBtaggedJets_TCHET, &b_nBtaggedJets_TCHET);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHPL", &nBtaggedJets_TCHPL, &b_nBtaggedJets_TCHPL);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHPM", &nBtaggedJets_TCHPM, &b_nBtaggedJets_TCHPM);
  //SetBranchAddress(fMTT, "nBtaggedJets_TCHPT", &nBtaggedJets_TCHPT, &b_nBtaggedJets_TCHPT);
  //SetBranchAddress(fMTT, "nBtaggedJets_SSVHEM", &nBtaggedJets_SSVHEM, &b_nBtaggedJets_SSVHEM);
  //SetBranchAddress(fMTT, "nBtaggedJets_SSVHPT", &nBtaggedJets_SSVHPT, &b_nBtaggedJets_SSVHPT);
  SetBranchAddress(fMTT, "MET", &MET, &b_MET);
  SetBranchAddress(fMTT, "isSel", &isSel, &b_isSel);
  //SetBranchAddress(fMTT, "oneMatchedCombi", &oneMatchedCombi, &b_oneMatchedCombi);
  SetBranchAddress(fMTT, "bestSolChi2", &bestSolChi2, &b_bestSolChi2);
  //SetBranchAddress(fMTT, "isBestSolMatched", &isBestSolMatched, &b_isBestSolMatched);
  //SetBranchAddress(fMTT, "KFChi2", &KFChi2, &b_KFChi2);
  //SetBranchAddress(fMTT, "numComb", &numComb, &b_numComb);
  //SetBranchAddress(fMTT, "solChi2", solChi2, &b_solChi2);
  //SetBranchAddress(fMTT, "mLepTop_AfterChi2", &mLepTop_AfterChi2, &b_mLepTop_AfterChi2);
  //SetBranchAddress(fMTT, "mHadTop_AfterChi2", &mHadTop_AfterChi2, &b_mHadTop_AfterChi2);
  SetBranchAddress(fMTT, "mtt_AfterChi2", &mtt_AfterChi2, &b_mtt_AfterChi2);
  SetBranchAddress(fMTT, "mLepTop_AfterChi2andKF", &mLepTop_AfterChi2andKF, &b_mLepTop_AfterChi2andKF);
  SetBranchAddress(fMTT, "mHadTop_AfterChi2andKF", &mHadTop_AfterChi2andKF, &b_mHadTop_AfterChi2andKF);
  SetBranchAddress(fMTT, "mtt_AfterChi2andKF", &mtt_AfterChi2andKF, &b_mtt_AfterChi2andKF);

  fVertices->SetMakeClass(1);
  fVertices->SetBranchAddress("n_vertices", &n_vertices, NULL);
  fVertices->SetBranchStatus("*", 0);
  fVertices->SetBranchStatus("n_vertices", 1);

  fEvent->SetMakeClass(1);
  fEvent->SetBranchAddress("nTrueInteractions", &n_trueInteractions, NULL);
  fEvent->SetBranchStatus("*", 0);
  fEvent->SetBranchStatus("nTrueInteractions", 1);
}
