#include "Extractor2Histos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <TProfile.h>
#include <fstream>
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

void Extractor2Histos::Loop()
{
  // Create the output file first in order that each histogram
  // associate itself with this file
  TFile * output = TFile::Open(mOutputFile.c_str(), "recreate");
  output->cd();

  TH1::SetDefaultSumw2(true);

  TH1D *hNVtx_noweight = new TH1D("hNVtx_noweight", "", 70, 0., 70);
  TH1D *hNVtx = new TH1D("hNVtx", "", 70, 0., 70);
  TH1D *hNVtx_beforesel = new TH1D("hNVtx_beforesel", "", 70, 0., 70);

  TH1D *hIsSel = new TH1D("hIsSel", "", 10, 0, 10);

  TH1D *hNTrueInt = new TH1D("hNTrueInt", "", 70, 0., 70);
  TH1D *hNTrueInt_nosel = new TH1D("hNTrueInt_nosel", "", 70, 0., 70);

  TH1D *hLeptonPt = new TH1D("hLeptonPt", "", 50, 20., 200.);
  TH1D *hLeptonPt_beforesel = new TH1D("hLeptonPt_beforesel", "", 50, 20., 200.);
  TH1D *hLeptonPt_nosel = new TH1D("hLeptonPt_nosel", "", 100, 0., 200.);

  TH1D *hLeptonEta = new TH1D("hLeptonEta", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hLeptonEta_beforesel = new TH1D("hLeptonEta_beforesel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hLeptonEta_nosel = new TH1D("hLeptonEta_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hLeptTopPt = new TH1D("hLeptTopPt", "", 60, 20., 600.);
  TH1D *hLeptTopPt_beforesel = new TH1D("hLeptTopPt_beforesel", "", 60, 20., 600.);
  TH1D *hLeptTopPt_nosel = new TH1D("hLeptTopPt_nosel", "", 60, 20., 600.);

  //TH1D *hTopPt_gen = new TH1D("hTopPt_gen", "", 60, 20., 600.);
  //TH1D *hTopEta_gen = new TH1D("hTopEta_gen", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hBoostTT_gen = new TH1D("hBoostTT_gen", "", 50, 0., 1.);
  TH1D *hPtTT_gen = new TH1D("hPtTT_gen", "", 60, 0., 600.);
  TH1D *hEtaTT_gen = new TH1D("hEtaTT_gen", "", 50, -2*M_PI, 2*M_PI);

  //TH1D *hTopPt_com_gen = new TH1D("hTopPt_com_gen", "", 60, 20., 600.);
  //TH1D *hTopEta_com_gen = new TH1D("hTopEta_com_gen", "", 50, -2*M_PI, 2*M_PI);

  //TH1D *hPtTT_com_gen = new TH1D("hPtTT_com_gen", "", 60, 0., 600.);
  //TH1D *hEtaTT_com_gen = new TH1D("hEtaTT_com_gen", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hLeptTopEta = new TH1D("hLeptTopEta", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hLeptTopEta_beforesel = new TH1D("hLeptTopEta_beforesel", "", 50, -2*M_PI, 2*M_PI);
  TH1D *hLeptTopEta_nosel = new TH1D("hLeptTopEta_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hHadrTopPt = new TH1D("hHadrTopPt", "", 60, 20., 600.);
  TH1D *hHadrTopPt_beforesel = new TH1D("hHadrTopPt_beforesel", "", 60, 20., 600.);
  TH1D *hHadrTopPt_nosel = new TH1D("hHadrTopPt_nosel", "", 60, 20., 600.);

  TH1D *hHadrTopEta = new TH1D("hHadrTopEta", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hHadrTopEta_beforesel = new TH1D("hHadrTopEta_beforesel", "", 50, -2*M_PI, 2*M_PI);
  TH1D *hHadrTopEta_nosel = new TH1D("hHadrTopEta_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hMuRelIso = new TH1D("hMuRelIso", "", 50, 0., 0.15);
  TH1D *hMuRelIso_beforesel = new TH1D("hMuRelIso_beforesel", "", 50, 0., 0.15);
  TH1D *hMuRelIso_nosel = new TH1D("hMuRelIso_nosel", "", 70, 0., 1);

  TH1D *hElRelIso = new TH1D("hElRelIso", "", 50, 0., 0.15);
  TH1D *hElRelIso_beforesel = new TH1D("hElRelIso_beforesel", "", 50, 0., 0.15);

  TH1D *h1stjetpt = new TH1D("h1stjetpt", "", 100, 70., 640.);
  TH1D *h1stjetpt_nosel = new TH1D("h1stjetpt_nosel", "", 100, 0., 640.);
  TH1D *h1stjetpt_beforesel = new TH1D("h1stjetpt_beforesel", "", 100, 0., 640.);

  TH1D *h2ndjetpt = new TH1D("h2ndjetpt", "", 100, 50., 620.);
  TH1D *h2ndjetpt_beforesel = new TH1D("h2ndjetpt_beforesel", "", 100, 0., 620.);
  TH1D *h2ndjetpt_nosel = new TH1D("h2ndjetpt_nosel", "", 100, 0., 620.);

  TH1D *h3rdjetpt = new TH1D("h3rdjetpt", "", 50, 30., 300.);
  TH1D *h3rdjetpt_beforesel = new TH1D("h3rdjetpt_beforesel", "", 50, 0., 300.);
  TH1D *h3rdjetpt_nosel = new TH1D("h3rdjetpt_nosel", "", 50, 0., 300.);

  TH1D *h4thjetpt = new TH1D("h4thjetpt", "", 50, 30., 300.);
  TH1D *h4thjetpt_beforesel = new TH1D("h4thjetpt_beforesel", "", 50, 0., 300.);
  TH1D *h4thjetpt_nosel = new TH1D("h4thjetpt_nosel", "", 50, 0., 300.);

  TH1D *h1stjeteta = new TH1D("h1stjeteta", "", 100, -2*M_PI, 2*M_PI);
  TH1D *h1stjeteta_nosel = new TH1D("h1stjeteta_nosel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *h1stjeteta_beforesel = new TH1D("h1stjeteta_beforesel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *h2ndjeteta = new TH1D("h2ndjeteta", "", 100, -2*M_PI, 2*M_PI);
  TH1D *h2ndjeteta_beforesel = new TH1D("h2ndjeteta_beforesel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *h2ndjeteta_nosel = new TH1D("h2ndjeteta_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *h3rdjeteta = new TH1D("h3rdjeteta", "", 100, -2*M_PI, 2*M_PI);
  TH1D *h3rdjeteta_beforesel = new TH1D("h3rdjeteta_beforesel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *h3rdjeteta_nosel = new TH1D("h3rdjeteta_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *h4thjeteta = new TH1D("h4thjeteta", "", 100, -2*M_PI, 2*M_PI);
  TH1D *h4thjeteta_beforesel = new TH1D("h4thjeteta_beforesel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *h4thjeteta_nosel = new TH1D("h4thjeteta_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hMET = new TH1D("hMET", "", 100, 20., 400.);
  TH1D *hMET_nosel = new TH1D("hMET_nosel", "", 100, 0., 400.);
  TH1D *hMET_beforesel = new TH1D("hMET_beforesel", "", 100, 20., 400.);

  TH1D *hmtlep = new TH1D("hmtlep", "", 100, 120., 240.);
  TH1D *hmthad = new TH1D("hmthad", "", 150, 120., 300.);
  
  TH1D *hmttSelected_btag_sel = new TH1D("hmttSelected_btag_sel", "", 250, 0., 2500.);
  TH1D *hmttSelected_btag_sel_mass_cut = new TH1D("hmttSelected_btag_sel_mass_cut", "", 250, 0., 2500.);

  TH1D *hNGoodMuons = new TH1D("hNGoodMuons", "", 5, -0.5, 4.5);

  TH1D *hNGoodJets = new TH1D("hNGoodJets", "", 6, 3.5, 9.5);
  TH1D *hNGoodJets_beforesel = new TH1D("hNGoodJets_beforesel", "", 6, 3.5, 9.5);

  TH1D *hNBtaggedJets = new TH1D("hNBtaggedJets", "", 5, -0.5, 4.5);
  TH1D *hNBtaggedJets_beforesel = new TH1D("hNBtaggedJets_beforesel", "", 5, -0.5, 4.5);

  TH1D *h_mtt_gen_no_sel = new TH1D("h_mtt_gen_no_sel", "", 500, 0., 2500.);
  TH1D *h_mtt_gen_beforesel = new TH1D("h_mtt_gen_beforesel", "", 500, 0., 2500.);
  TH1D *h_mtt_gen = new TH1D("h_mtt_gen", "", 500, 0, 2500.);

  TH1D *h_mtt_resolution = new TH1D("h_mtt_resolution", "", 100, -600., 600.);

  TProfile *pMttResolution_btag_sel = new TProfile("pMttResolution_btag_sel", "", nBins, bins);
  pMttResolution_btag_sel->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  TH1D *hBoostTT = new TH1D("hBoostTT", "", 50, 0., 1.);
  TH1D *hPtTT = new TH1D("hPtTT", "", 60, 0., 600.);
  TH1D *hEtaTT = new TH1D("hEtaTT", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hBoostTT_beforesel = new TH1D("hBoostTT_beforesel", "", 50, 0., 1.);
  TH1D *hPtTT_beforesel = new TH1D("hPtTT_beforesel", "", 60, 0., 600.);
  TH1D *hEtaTT_beforesel = new TH1D("hEtaTT_beforesel", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hDeltaPhiTops_gen = new TH1D("hDeltaPhiTops_gen", "", 50, -4*M_PI, 4*M_PI);
  TH1D *hDeltaPhiTops_reco_fullsel = new TH1D("hDeltaPhiTops_reco_fullsel", "", 50, -4*M_PI, 4*M_PI);

  if (mIsSemiMu) {
    hLeptonPt->SetXTitle("#mu p_{T} [GeV/c]");
    hLeptonPt_beforesel->SetXTitle("#mu p_{T} [GeV/c]");
  } else {
    hLeptonPt_beforesel->SetXTitle("e p_{T} [GeV/c]");
    hLeptonPt->SetXTitle("e p_{T} [GeV/c]");
  }

  h1stjetpt->SetXTitle("1^{st} jet p_{T} [GeV/c]");
  h1stjetpt_nosel->SetXTitle("1^{st} jet p_{T} [GeV/c]");
  h1stjetpt_beforesel->SetXTitle("1^{st} jet p_{T} [GeV/c]");

  h2ndjetpt->SetXTitle("2^{nd} jet p_{T} [GeV/c]");
  h2ndjetpt_beforesel->SetXTitle("2^{nd} jet p_{T} [GeV/c]");

  h3rdjetpt->SetXTitle("3^{rd} jet p_{T} [GeV/c]");
  h3rdjetpt_beforesel->SetXTitle("3^{rd} jet p_{T} [GeV/c]");

  h4thjetpt_beforesel->SetXTitle("4^{th} jet p_{T} [GeV/c]");

  hMET->SetXTitle("MET [GeV]");
  hMET_beforesel->SetXTitle("MET [GeV]");

  hmtlep->SetXTitle("leptonic m_{t} [GeV/c^{2}]");
  hmthad->SetXTitle("hadronic m_{t} [GeV/c^{2}]");

  h_mtt_gen_no_sel->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");
  h_mtt_gen_beforesel->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");
  h_mtt_gen->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  h_mtt_resolution->SetXTitle("m_{t#bar{t}}^{reco} - m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  hmttSelected_btag_sel->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmttSelected_btag_sel_mass_cut->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hNGoodMuons->SetXTitle("Num good #mu");
  hNGoodJets->SetXTitle("Num good jets");
  hNGoodJets_beforesel->SetXTitle("Num good jets");

  hNBtaggedJets->SetXTitle("Num TCHEL jets");
  hNBtaggedJets_beforesel->SetXTitle("Num TCHEL jets");

  Long64_t nentries = fMTT->GetEntries();

  //PUReweighter puReweighter(mIsSemiMu, mDataset);
  PUReweighter puReweighter(mIsSemiMu);

  std::cout << "Processing..." << std::endl;

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if (jentry % 100000 == 0)
      std::cout << "Processing entry #" << (jentry + 1) << " over " << nentries << " (" << (float) jentry / nentries * 100 << "%)" << std::endl;

    GetEntry(jentry);

    double eventWeight = 1.;
    if (mIsMC) {
      eventWeight *= puReweighter.weight(n_trueInteractions);
    } else {
      eventWeight *= m_weight;
    }

    hNTrueInt_nosel->Fill(n_trueInteractions, eventWeight);

    if (!mIsMC && !m_triggerPassed) {
      continue;
    }

    hIsSel->Fill(isSel, eventWeight);
    h_mtt_gen_no_sel->Fill(MC_mtt, eventWeight);
    hMET_nosel->Fill(MET, eventWeight);
    if (jet_p4->GetEntriesFast() > 0) {
      TLorentzVector* p4 = (TLorentzVector*) (*jet_p4)[0];
      h1stjetpt_nosel->Fill(p4->Pt(), eventWeight);
      h1stjeteta_nosel->Fill(p4->Eta(), eventWeight);
    }

    if (jet_p4->GetEntriesFast() > 1) {
      TLorentzVector* p4 = (TLorentzVector*) (*jet_p4)[1];
      h2ndjetpt_nosel->Fill(p4->Pt(), eventWeight);
      h2ndjeteta_nosel->Fill(p4->Eta(), eventWeight);
    }

    if (jet_p4->GetEntriesFast() > 2) {
      TLorentzVector* p4 = (TLorentzVector*) (*jet_p4)[2];
      h3rdjetpt_nosel->Fill(p4->Pt(), eventWeight);
      h3rdjeteta_nosel->Fill(p4->Eta(), eventWeight);
    }

    if (jet_p4->GetEntriesFast() > 3) {
      TLorentzVector* p4 = (TLorentzVector*) (*jet_p4)[3];
      h4thjeteta_nosel->Fill(p4->Eta(), eventWeight);
      h4thjetpt_nosel->Fill(p4->Pt(), eventWeight);
    }

    //hTopPt_gen->Fill(MC_top1Pt, eventWeight);
    //hTopPt_gen->Fill(MC_top2Pt, eventWeight);
    //hTopPt_com_gen->Fill(MC_top1Pt_com, eventWeight);
    //hTopPt_com_gen->Fill(MC_top2Pt_com, eventWeight);

    //hTopEta_gen->Fill(MC_top1Eta, eventWeight);
    //hTopEta_gen->Fill(MC_top2Eta, eventWeight);
    //hTopEta_com_gen->Fill(MC_top1Eta_com, eventWeight);
    //hTopEta_com_gen->Fill(MC_top2Eta_com, eventWeight);

    hBoostTT_gen->Fill(MC_boost_tt, eventWeight);

    hPtTT_gen->Fill(MC_pt_tt, eventWeight);
    //hPtTT_com_gen->Fill(MC_pt_tt_com, eventWeight);

    hEtaTT_gen->Fill(MC_eta_tt, eventWeight);
    //hEtaTT_com_gen->Fill(MC_eta_tt_com, eventWeight);

    hDeltaPhiTops_gen->Fill(gen_top1_p4->DeltaPhi(*gen_top2_p4), eventWeight);

    if (n_muons > 0) {
      hMuRelIso_nosel->Fill(muon_relIso[0], eventWeight); 

      TLorentzVector* p4 = (TLorentzVector*) (*muon_p4)[0];
      hLeptonEta_nosel->Fill(p4->Eta(), eventWeight);
      hLeptonPt_nosel->Fill(p4->Pt(), eventWeight);
    }

    hLeptTopPt_nosel->Fill(lepTopPt_AfterChi2, eventWeight);
    hLeptTopEta_nosel->Fill(lepTopEta_AfterChi2, eventWeight);

    hHadrTopPt_nosel->Fill(hadTopPt_AfterChi2, eventWeight);
    hHadrTopEta_nosel->Fill(hadTopEta_AfterChi2, eventWeight);

    double ptLepton = -1;
    double etaLepton = -1;
    double ptLeptonCut = -1;
    if (mIsSemiMu)
    {
      if (nGoodMuons <= 0)
        continue;

      ptLepton = muonPt[0];
      etaLepton = muonEta[0];
      ptLeptonCut = 25.;
    }
    else
    {
      if (nGoodElectrons <= 0)
        continue;

      ptLepton = electronPt[0];
      etaLepton = electronEta[0];
      ptLeptonCut = 30.;
    }

    if (isSel == 1 && numComb > 0)
    {
      h_mtt_gen_beforesel->Fill(MC_mtt, eventWeight);

      hNGoodJets_beforesel->Fill(nJets, eventWeight);
      hNBtaggedJets_beforesel->Fill(nBtaggedJets_CSVM, eventWeight);

      hLeptonPt_beforesel->Fill(ptLepton, eventWeight);
      hLeptonEta_beforesel->Fill(etaLepton, eventWeight);
      if (mIsSemiMu)
        hMuRelIso_beforesel->Fill(muRelIso[0], eventWeight);
      else
        hElRelIso_beforesel->Fill(elRelIso[0], eventWeight);
      
      hMET_beforesel->Fill(MET, eventWeight);

      hLeptTopPt_beforesel->Fill(lepTopPt_AfterChi2, eventWeight);
      hLeptTopEta_beforesel->Fill(lepTopEta_AfterChi2, eventWeight);

      hHadrTopPt_beforesel->Fill(hadTopPt_AfterChi2, eventWeight);
      hHadrTopEta_beforesel->Fill(hadTopEta_AfterChi2, eventWeight);

      h1stjetpt_beforesel->Fill(p_1stjetpt, eventWeight);
      h2ndjetpt_beforesel->Fill(p_2ndjetpt, eventWeight);
      h3rdjetpt_beforesel->Fill(p_3rdjetpt, eventWeight);
      h4thjetpt_beforesel->Fill(p_4thjetpt, eventWeight);

      h1stjeteta_beforesel->Fill(jetEta[0], eventWeight);
      h2ndjeteta_beforesel->Fill(jetEta[1], eventWeight);
      h3rdjeteta_beforesel->Fill(jetEta[2], eventWeight);
      h4thjeteta_beforesel->Fill(jetEta[3], eventWeight);

      hNVtx_beforesel->Fill(n_vertices, eventWeight);

      hBoostTT_beforesel->Fill(beta_tt_AfterChi2, eventWeight);
      hPtTT_beforesel->Fill(pt_tt_AfterChi2, eventWeight);
      hEtaTT_beforesel->Fill(eta_tt_AfterChi2, eventWeight);

      bool btagSel = false;
      if (mBTag == 1)
        btagSel = nBtaggedJets_CSVM == 1;
      else if (mBTag == 2)
        btagSel = nBtaggedJets_CSVM > 1;

      if (ptLepton > ptLeptonCut && btagSel && p_1stjetpt > 70. && p_2ndjetpt > 50 && mtt_AfterChi2 > 0. && bestSolChi2 < 500)
      {

        hBoostTT->Fill(beta_tt_AfterChi2, eventWeight);
        hPtTT->Fill(pt_tt_AfterChi2, eventWeight);
        hEtaTT->Fill(eta_tt_AfterChi2, eventWeight);

        h_mtt_gen->Fill(MC_mtt, eventWeight);
        h_mtt_resolution->Fill(mtt_AfterChi2 - MC_mtt, eventWeight);

        hLeptonPt->Fill(ptLepton, eventWeight);
        hLeptonEta->Fill(etaLepton, eventWeight);
        h1stjetpt->Fill(p_1stjetpt, eventWeight);
        h2ndjetpt->Fill(p_2ndjetpt, eventWeight);
        h3rdjetpt->Fill(p_3rdjetpt, eventWeight);
        h4thjetpt->Fill(p_4thjetpt, eventWeight);
        h1stjeteta->Fill(jetEta[0], eventWeight);
        h2ndjeteta->Fill(jetEta[1], eventWeight);
        h3rdjeteta->Fill(jetEta[2], eventWeight);
        h4thjeteta->Fill(jetEta[3], eventWeight);
        hMET->Fill(MET, eventWeight);

        hLeptTopPt->Fill(lepTopPt_AfterChi2, eventWeight);
        hLeptTopEta->Fill(lepTopEta_AfterChi2, eventWeight);

        hHadrTopPt->Fill(hadTopPt_AfterChi2, eventWeight);
        hHadrTopEta->Fill(hadTopEta_AfterChi2, eventWeight);

        hNGoodJets->Fill(nJets, eventWeight);
        hNBtaggedJets->Fill(nBtaggedJets_CSVM, eventWeight);

        hNVtx_noweight->Fill(n_vertices);
        hNVtx->Fill(n_vertices, eventWeight);
        hNTrueInt->Fill(n_trueInteractions, eventWeight);

        if (mIsSemiMu)
          hMuRelIso->Fill(muRelIso[0], eventWeight);
        else
          hElRelIso->Fill(elRelIso[0], eventWeight);

        hmtlep->Fill(mLepTop_AfterChi2, eventWeight);
        hmthad->Fill(mHadTop_AfterChi2, eventWeight);

        hmttSelected_btag_sel->Fill(mtt_AfterChi2, eventWeight);
        pMttResolution_btag_sel->Fill(MC_mtt , mtt_AfterChi2, eventWeight);

        hDeltaPhiTops_reco_fullsel->Fill(lepTopP4_AfterChi2->DeltaPhi(*hadTopP4_AfterChi2), eventWeight);

        if (mtt_AfterChi2 > 500)
        {
          hmttSelected_btag_sel_mass_cut->Fill(mtt_AfterChi2, eventWeight);
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

Extractor2Histos::Extractor2Histos(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isSemiMu, bool isMC, int btag) : fMTT(0), fVertices(0), fEvent(0)
{
  mIsSemiMu = isSemiMu;
  mIsMC = isMC;
  mOutputFile = outputFile;
  mBTag = btag;

  // Get trees
  loadChain(inputFiles, "Mtt", fMTT);
  loadChain(inputFiles, "Vertices", fVertices);
  loadChain(inputFiles, "event", fEvent);
  loadChain(inputFiles, "muon_loose_PF", fLooseMuons);
  loadChain(inputFiles, "jet_PF", fJet);


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

  if (fLooseMuons)
    fLooseMuons->GetEntry(entry);

  if (fJet)
    fJet->GetEntry(entry);

  return 1;
}

void Extractor2Histos::SetBranchAddress(TTree* t, const char* branchName, void* ptr, TBranch** branch) {
  t->SetBranchStatus(branchName, 1);
  t->SetBranchAddress(branchName, ptr, branch);
}

void Extractor2Histos::Init()
{
  fCurrent = -1;

  fMTT->SetBranchStatus("*", 0);

  SetBranchAddress(fMTT, "MC_mtt", &MC_mtt, &b_MC_mtt);
  SetBranchAddress(fMTT, "MC_beta_tt", &MC_boost_tt, NULL);

  SetBranchAddress(fMTT, "MC_pt_tt", &MC_pt_tt, NULL);
  SetBranchAddress(fMTT, "MC_eta_tt", &MC_eta_tt, NULL);

  //SetBranchAddress(fMTT, "MC_top1Pt", &MC_top1Pt, NULL);
  //SetBranchAddress(fMTT, "MC_top2Pt", &MC_top2Pt, NULL);

  //SetBranchAddress(fMTT, "MC_top1Eta", &MC_top1Eta, NULL);
  //SetBranchAddress(fMTT, "MC_top2Eta", &MC_top2Eta, NULL);

  //SetBranchAddress(fMTT, "MC_pt_tt_com", &MC_pt_tt_com, NULL);
  //SetBranchAddress(fMTT, "MC_eta_tt_com", &MC_eta_tt_com, NULL);

  //SetBranchAddress(fMTT, "MC_top1Pt_com", &MC_top1Pt_com, NULL);
  //SetBranchAddress(fMTT, "MC_top2Pt_com", &MC_top2Pt_com, NULL);

  //SetBranchAddress(fMTT, "MC_top1Eta_com", &MC_top1Eta_com, NULL);
  //SetBranchAddress(fMTT, "MC_top2Eta_com", &MC_top2Eta_com, NULL);

  //SetBranchAddress(fMTT, "MC_nPU", &MC_nPU, &b_m_nPU);
  SetBranchAddress(fMTT, "nGoodMuons", &nGoodMuons, &b_nGoodMuons);
  //SetBranchAddress(fMTT, "nLooseGoodMuons", &nLooseGoodMuons, &b_nLooseGoodMuons);
  if (mIsSemiMu) {
    SetBranchAddress(fMTT, "muonPt", muonPt, &b_muonPt);
    SetBranchAddress(fMTT, "muonEta", muonEta, &b_muonPt);
    //SetBranchAddress(fMTT, "2DDrMin", p_2DDrMin, &b_2DDrMin);
    //SetBranchAddress(fMTT, "2DpTrel", p_2DpTrel, &b_2DpTrel);
    SetBranchAddress(fMTT, "muRelIso", muRelIso, &b_muRelIso);
  } else {
    SetBranchAddress(fMTT, "nGoodElectrons", &nGoodElectrons, &b_nGoodElectrons);
    SetBranchAddress(fMTT, "electronPt", &electronPt, &b_electronPt);
    SetBranchAddress(fMTT, "electronEta", &electronEta, &b_electronPt);
    SetBranchAddress(fMTT, "elRelIso", &elRelIso, &b_elRelIso);
    //SetBranchAddress(fMTT, "hyperTight1MC", &hyperTight1MC, &b_hyperTight1MC);
  }
  SetBranchAddress(fMTT, "1stjetpt", &p_1stjetpt, &b_1stjetpt);
  SetBranchAddress(fMTT, "2ndjetpt", &p_2ndjetpt, &b_2ndjetpt);
  SetBranchAddress(fMTT, "3rdjetpt", &p_3rdjetpt, &b_3rdjetpt);
  SetBranchAddress(fMTT, "4thjetpt", &p_4thjetpt, &b_4thjetpt);
  SetBranchAddress(fMTT, "nJets", &nJets, &b_nJets);
  SetBranchAddress(fMTT, "jetEta", jetEta, &b_jetEta);
  //SetBranchAddress(fMTT, "jetPt", jetPt, &b_jetPt);
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
  SetBranchAddress(fMTT, "numComb", &numComb, &b_numComb);
  //SetBranchAddress(fMTT, "solChi2", solChi2, &b_solChi2);
  SetBranchAddress(fMTT, "mLepTop_AfterChi2", &mLepTop_AfterChi2, &b_mLepTop_AfterChi2);
  SetBranchAddress(fMTT, "mHadTop_AfterChi2", &mHadTop_AfterChi2, &b_mHadTop_AfterChi2);
  SetBranchAddress(fMTT, "mtt_AfterChi2", &mtt_AfterChi2, &b_mtt_AfterChi2);
  SetBranchAddress(fMTT, "pt_tt_AfterChi2", &pt_tt_AfterChi2, NULL);
  SetBranchAddress(fMTT, "eta_tt_AfterChi2", &eta_tt_AfterChi2, NULL);
  SetBranchAddress(fMTT, "beta_tt_AfterChi2", &beta_tt_AfterChi2, NULL);

  SetBranchAddress(fMTT, "lepTopPt_AfterChi2", &lepTopPt_AfterChi2, NULL);
  SetBranchAddress(fMTT, "lepTopEta_AfterChi2", &lepTopEta_AfterChi2, NULL);
  SetBranchAddress(fMTT, "hadTopPt_AfterChi2", &hadTopPt_AfterChi2, NULL);
  SetBranchAddress(fMTT, "hadTopEta_AfterChi2", &hadTopEta_AfterChi2, NULL);
  //SetBranchAddress(fMTT, "mLepTop_AfterChi2andKF", &mLepTop_AfterChi2andKF, &b_mLepTop_AfterChi2andKF);
  //SetBranchAddress(fMTT, "mHadTop_AfterChi2andKF", &mHadTop_AfterChi2andKF, &b_mHadTop_AfterChi2andKF);
  //SetBranchAddress(fMTT, "mtt_AfterChi2andKF", &mtt_AfterChi2andKF, &b_mtt_AfterChi2andKF);
  SetBranchAddress(fMTT, "weight", &m_weight, NULL);

  if (fMTT->GetBranch("trigger_passed")) {
    SetBranchAddress(fMTT, "trigger_passed", &m_triggerPassed, NULL);
  } else {
    // Backward compatibilty
    m_triggerPassed = true;
  }

  lepTopP4_AfterChi2 = NULL;
  hadTopP4_AfterChi2 = NULL;

  fMTT->SetBranchAddress("lepTopP4_AfterChi2.", &lepTopP4_AfterChi2);
  fMTT->SetBranchAddress("hadTopP4_AfterChi2.", &hadTopP4_AfterChi2);

  fMTT->SetBranchStatus("lepTopP4_AfterChi2*", 1);
  fMTT->SetBranchStatus("hadTopP4_AfterChi2*", 1);

  gen_top1_p4 = NULL;
  gen_top2_p4 = NULL;

  fMTT->SetBranchAddress("MC_Top1_p4.", &gen_top1_p4);
  fMTT->SetBranchAddress("MC_Top2_p4.", &gen_top2_p4);

  fMTT->SetBranchStatus("MC_Top1_p4*", 1);
  fMTT->SetBranchStatus("MC_Top2_p4*", 1);

  fVertices->SetMakeClass(1);
  fVertices->SetBranchAddress("n_vertices", &n_vertices, NULL);
  fVertices->SetBranchStatus("*", 0);
  fVertices->SetBranchStatus("n_vertices", 1);

  fEvent->SetMakeClass(1);
  fEvent->SetBranchAddress("nTrueInteractions", &n_trueInteractions, NULL);
  fEvent->SetBranchStatus("*", 0);
  fEvent->SetBranchStatus("nTrueInteractions", 1);

  muon_p4 = NULL;
  fLooseMuons->SetMakeClass(1);
  fLooseMuons->SetBranchStatus("*", 0);
  fLooseMuons->SetBranchStatus("muon_4vector", 1);
  fLooseMuons->SetBranchStatus("n_muons", 1);
  fLooseMuons->SetBranchStatus("muon_deltaBetaCorrectedRelIsolation", 1);
  fLooseMuons->SetBranchAddress("muon_4vector", &muon_p4, NULL);
  fLooseMuons->SetBranchAddress("n_muons", &n_muons, NULL);
  fLooseMuons->SetBranchAddress("muon_deltaBetaCorrectedRelIsolation", &muon_relIso, NULL);

  jet_p4 = NULL;
  fJet->SetMakeClass(1);
  fJet->SetBranchStatus("*", 0);
  fJet->SetBranchStatus("jet_4vector", 1);

  fJet->SetBranchAddress("jet_4vector", &jet_p4, NULL);
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
    
    Extractor2Histos convertor(inputFiles, outputFileArg.getValue(), semimuArg.isSet(), !isData, btagArg.getValue());
    convertor.Loop();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
