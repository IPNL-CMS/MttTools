#include "Extractor2Histos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom2.h>
#include <vector>
#include <TProfile.h>
#include <fstream>
#include <memory>
#include <TLorentzVector.h>

#include "TopTriggerEfficiencyProvider.h"
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

  TH1D *hWeight = new TH1D("weight", "", 50, 0, 5);
  TH1D *hWeight_fullsel = new TH1D("weight_fullsel", "", 50, 0, 5);
  TH1D *hLeptonWeight = new TH1D("lepton_weight", "", 50, 0, 5);
  TH1D *hBTagWeight = new TH1D("btag_weight", "", 50, 0, 5);
  TH1D *hTriggerWeight = new TH1D("trigger_weight", "", 50, 0, 1);
  TH1D *hPUWeight = new TH1D("PU_weight", "", 50, 0, 5);
  TH1D *hGeneratorWeight = new TH1D("generator_weight", "", 100, -2, 2);

  TH1D *hRunPeriod = new TH1D("run_period", "", 2, 0, 2);
  hRunPeriod->GetXaxis()->SetBinLabel(1, "Run2012 A+B");
  hRunPeriod->GetXaxis()->SetBinLabel(2, "Run2012 C+D");

  TH1D *hNVtx_noweight = new TH1D("nVertex_reco_fullsel_noweight", "", 70, 0., 70);
  TH1D *hNVtx = new TH1D("nVertex_reco_fullsel", "", 70, 0., 70);
  TH1D *hNVtx_chi2sel = new TH1D("nVertex_reco_chi2sel", "", 70, 0., 70);

  TH1D *hIsSel = new TH1D("isSel", "", 10, 0, 10);

  TH1D *hNTrueInt = new TH1D("nTrueInt_reco_fullsel", "", 70, 0., 70);
  TH1D *hNTrueInt_nosel = new TH1D("nTrueInt_reco_nosel", "", 70, 0., 70);

  TH1D *hLeptonPt = new TH1D("leptonPt_reco_fullsel", "", 50, 20., 200.);
  TH1D *hLeptonPt_chi2sel = new TH1D("leptonPt_reco_chi2sel", "", 50, 20., 200.);
  TH1D *hLeptonPt_nosel = new TH1D("leptonPt_reco_nosel", "", 100, 0., 200.);

  TH1D *hLeptonEta = new TH1D("leptonEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hLeptonEta_chi2sel = new TH1D("leptonEta_reco_chi2sel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hLeptonEta_nosel = new TH1D("leptonEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hLeptTopPt = new TH1D("leptTopPt_reco_fullsel", "", 60, 20., 600.);
  TH1D *hLeptTopPt_chi2sel = new TH1D("leptTopPt_reco_chi2sel", "", 60, 20., 600.);
  TH1D *hLeptTopPt_nosel = new TH1D("leptTopPt_reco_nosel", "", 60, 20., 600.);

  //TH1D *hTopPt_gen = new TH1D("hTopPt_gen", "", 60, 20., 600.);
  //TH1D *hTopEta_gen = new TH1D("hTopEta_gen", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hBoostTT_gen = new TH1D("boostTT_gen", "", 50, 0., 1.);
  TH1D *hPtTT_gen = new TH1D("ptTT_gen", "", 60, 0., 600.);
  TH1D *hEtaTT_gen = new TH1D("etaTT_gen", "", 50, -2*M_PI, 2*M_PI);

  //TH1D *hTopPt_com_gen = new TH1D("hTopPt_com_gen", "", 60, 20., 600.);
  //TH1D *hTopEta_com_gen = new TH1D("hTopEta_com_gen", "", 50, -2*M_PI, 2*M_PI);

  //TH1D *hPtTT_com_gen = new TH1D("hPtTT_com_gen", "", 60, 0., 600.);
  //TH1D *hEtaTT_com_gen = new TH1D("hEtaTT_com_gen", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hLeptTopEta = new TH1D("leptTopEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hLeptTopEta_chi2sel = new TH1D("leptTopEta_reco_chi2sel", "", 50, -2*M_PI, 2*M_PI);
  TH1D *hLeptTopEta_nosel = new TH1D("leptTopEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hHadrTopPt = new TH1D("hadrTopPt_reco_fullsel", "", 60, 20., 600.);
  TH1D *hHadrTopPt_chi2sel = new TH1D("hadrTopPt_reco_chi2sel", "", 60, 20., 600.);
  TH1D *hHadrTopPt_nosel = new TH1D("hadrTopPt_reco_nosel", "", 60, 20., 600.);

  TH1D *hHadrTopEta = new TH1D("hadrTopEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hHadrTopEta_chi2sel = new TH1D("hadrTopEta_reco_chi2sel", "", 50, -2*M_PI, 2*M_PI);
  TH1D *hHadrTopEta_nosel = new TH1D("hadrTopEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hMuRelIso = new TH1D("muRelIso_reco_fullsel", "", 50, 0., 0.15);
  TH1D *hMuRelIso_chi2sel = new TH1D("muRelIso_reco_chi2sel", "", 50, 0., 0.15);
  TH1D *hMuRelIso_nosel = new TH1D("muRelIso_reco_nosel", "", 70, 0., 1);

  TH1D *hElRelIso = new TH1D("elRelIso_reco_fullsel", "", 50, 0., 0.15);
  TH1D *hElRelIso_chi2sel = new TH1D("elRelIso_reco_chi2sel", "", 50, 0., 0.15);

  TH1D *hFirstJetPt = new TH1D("firstJetPt_reco_fullsel", "", 100, 70., 640.);
  TH1D *hFirstJetPt_nosel = new TH1D("firstJetPt_reco_nosel", "", 100, 0., 640.);
  TH1D *hFirstJetPt_chi2sel = new TH1D("firstJetPt_reco_chi2sel", "", 100, 0., 640.);

  TH1D *hSecondJetPt = new TH1D("secondJetPt_reco_fullsel", "", 100, 50., 620.);
  TH1D *hSecondJetPt_chi2sel = new TH1D("secondJetPt_reco_chi2sel", "", 100, 0., 620.);
  TH1D *hSecondJetPt_nosel = new TH1D("secondJetPt_reco_nosel", "", 100, 0., 620.);

  TH1D *hThirdJetPt = new TH1D("thirdJetPt_reco_fullsel", "", 50, 30., 300.);
  TH1D *hThirdJetPt_chi2sel = new TH1D("thirdJetPt_reco_chi2sel", "", 50, 0., 300.);
  TH1D *hThirdJetPt_nosel = new TH1D("thirdJetPt_reco_nosel", "", 50, 0., 300.);

  TH1D *hFourthJetPt = new TH1D("fourthJetPt_reco_fullsel", "", 50, 30., 300.);
  TH1D *hFourthJetPt_chi2sel = new TH1D("fourthJetPt_reco_chi2sel", "", 50, 0., 300.);
  TH1D *hFourthJetPt_nosel = new TH1D("fourthJetPt_reco_nosel", "", 50, 0., 300.);

  TH1D *hFirstJetEta = new TH1D("firstJetEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hFirstJetEta_nosel = new TH1D("firstJetEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hFirstJetEta_chi2sel = new TH1D("firstJetEta_reco_chi2sel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hSecondJetEta = new TH1D("secondJetEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hSecondJetEta_chi2sel = new TH1D("secondJetEta_reco_chi2sel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hSecondJetEta_nosel = new TH1D("secondJetEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hThirdJetEta = new TH1D("thirdJetEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hThirdJetEta_chi2sel = new TH1D("thirdJetEta_reco_chi2sel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hThirdJetEta_nosel = new TH1D("thirdJetEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hFourthJetEta = new TH1D("fourthJetEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hFourthJetEta_chi2sel = new TH1D("fourthJetEta_reco_chi2sel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hFourthJetEta_nosel = new TH1D("fourthJetEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hMET = new TH1D("MET_reco_fullsel", "", 100, 20., 400.);
  TH1D *hMET_nosel = new TH1D("MET_reco_nosel", "", 100, 0., 400.);
  TH1D *hMET_chi2sel = new TH1D("MET_reco_chi2sel", "", 100, 20., 400.);

  TH1D *hmtlep = new TH1D("mtLep_reco_fullsel", "", 100, 120., 240.);
  TH1D *hmthad = new TH1D("mtHad_reco_fullsel", "", 150, 120., 300.);
  
  TH1D *hmttSelected_btag_sel = new TH1D("mttSelected_btag_sel_reco_fullsel", "", 250, 0., 2500.);
  TH1D *hmttSelected_btag_sel_mass_cut = new TH1D("mttSelected_btag_sel_mass_cut_reco_fullsel", "", 250, 0., 2500.);

  TH1D *hNGoodJets = new TH1D("nGoodJets_reco_fullsel", "", 6, 3.5, 9.5);
  TH1D *hNGoodJets_chi2sel = new TH1D("nGoodJets_reco_chi2sel", "", 6, 3.5, 9.5);

  TH1D *hNBtaggedJets = new TH1D("nBTaggedJets_reco_fullsel", "", 5, -0.5, 4.5);
  TH1D *hNBtaggedJets_chi2sel = new TH1D("nBTaggedJets_reco_chi2sel", "", 5, -0.5, 4.5);

  TH1D *h_mtt_gen_no_sel = new TH1D("mtt_gen_nosel", "", 500, 0., 2500.);
  TH1D *h_mtt_gen_chi2sel = new TH1D("mtt_gen_chi2sel", "", 500, 0., 2500.);
  TH1D *h_mtt_gen = new TH1D("mtt_gen_fullsel", "", 500, 0, 2500.);

  TH1D *h_mtt_resolution = new TH1D("mtt_resolution", "", 100, -600., 600.);

  TProfile *pMttResolution_btag_sel = new TProfile("pMttResolution_btag_sel", "", nBins, bins);
  pMttResolution_btag_sel->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  TH1D *hBoostTT = new TH1D("boostTT_reco_fullsel", "", 50, 0., 1.);
  TH1D *hPtTT = new TH1D("ptTT_reco_fullsel", "", 60, 0., 600.);
  TH1D *hEtaTT = new TH1D("etaTT_reco_fullsel", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hBoostTT_chi2sel = new TH1D("boostTT_reco_chi2sel", "", 50, 0., 1.);
  TH1D *hPtTT_chi2sel = new TH1D("ptTT_reco_chi2sel", "", 60, 0., 600.);
  TH1D *hEtaTT_chi2sel = new TH1D("etaTT_reco_chi2sel", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hDeltaPhiTops_gen = new TH1D("deltaPhiTops_gen", "", 50, 0, M_PI);
  TH1D *hDeltaPhiTops_reco_fullsel = new TH1D("deltaPhiTops_reco_fullsel", "", 50, 0, M_PI);

  TH1D *hDeltaEtaTops_gen = new TH1D("deltaEtaTops_gen", "", 50, -3*M_PI, 3*M_PI);
  TH1D *hDeltaEtaTops_reco_fullsel = new TH1D("deltaEtaTops_reco_fullsel", "", 50, -3*M_PI, 3*M_PI);

  TH1D *hDeltaRTops_gen = new TH1D("deltaRTops_gen", "", 50, 0, 10);
  TH1D *hDeltaRTops_reco_fullsel = new TH1D("deltaRTops_reco_fullsel", "", 50, 0, 10);

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

  if (mIsSemiMu) {
    hLeptonPt->SetXTitle("#mu p_{T} [GeV/c]");
    hLeptonPt_chi2sel->SetXTitle("#mu p_{T} [GeV/c]");
  } else {
    hLeptonPt_chi2sel->SetXTitle("e p_{T} [GeV/c]");
    hLeptonPt->SetXTitle("e p_{T} [GeV/c]");
  }

  hFirstJetPt->SetXTitle("1^{st} jet p_{T} [GeV/c]");
  hFirstJetPt_nosel->SetXTitle("1^{st} jet p_{T} [GeV/c]");
  hFirstJetPt_chi2sel->SetXTitle("1^{st} jet p_{T} [GeV/c]");

  hSecondJetPt->SetXTitle("2^{nd} jet p_{T} [GeV/c]");
  hSecondJetPt_chi2sel->SetXTitle("2^{nd} jet p_{T} [GeV/c]");

  hThirdJetPt->SetXTitle("3^{rd} jet p_{T} [GeV/c]");
  hThirdJetPt_chi2sel->SetXTitle("3^{rd} jet p_{T} [GeV/c]");

  hFourthJetPt_chi2sel->SetXTitle("4^{th} jet p_{T} [GeV/c]");

  hMET->SetXTitle("MET [GeV]");
  hMET_chi2sel->SetXTitle("MET [GeV]");

  hmtlep->SetXTitle("leptonic m_{t} [GeV/c^{2}]");
  hmthad->SetXTitle("hadronic m_{t} [GeV/c^{2}]");

  h_mtt_gen_no_sel->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");
  h_mtt_gen_chi2sel->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");
  h_mtt_gen->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  h_mtt_resolution->SetXTitle("m_{t#bar{t}}^{reco} - m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  hmttSelected_btag_sel->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmttSelected_btag_sel_mass_cut->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hNGoodJets->SetXTitle("Num good jets");
  hNGoodJets_chi2sel->SetXTitle("Num good jets");

  hNBtaggedJets->SetXTitle("Num TCHEL jets");
  hNBtaggedJets_chi2sel->SetXTitle("Num TCHEL jets");

  Long64_t nentries = fMTT->GetEntries();
  //nentries = 10000;

  //PUReweighter puReweighter(mIsSemiMu, mDataset);
  PUReweighter puReweighter(mIsSemiMu);

  std::cout << "Processing..." << std::endl;

  // 2012 luminosity
  float lumi_run2012_A = 0;
  float lumi_run2012_B = 0;
  float lumi_run2012_C = 0;
  float lumi_run2012_D = 0;

  if (mIsSemiMu) {
    lumi_run2012_A = 0.876225;
    lumi_run2012_B = 4.412;
    lumi_run2012_C = 7.044;
    lumi_run2012_D = 7.368;
  } else {
    lumi_run2012_A = 0.876225;
    lumi_run2012_B = 4.399;
    lumi_run2012_C = 7.022;
    lumi_run2012_D = 7.369;
  }

  float lumi_run2012_AB = lumi_run2012_A + lumi_run2012_B;
  float lumi_run2012_CD = lumi_run2012_C + lumi_run2012_D;
  float lumi_total = lumi_run2012_AB + lumi_run2012_CD;

  float lumi_run2012_AB_over_total = lumi_run2012_AB / lumi_total;

  std::shared_ptr<TopTriggerEfficiencyProvider> m_trigger_efficiency_provider = std::make_shared<TopTriggerEfficiencyProvider>();
  TRandom2 random_generator;

  uint64_t positive_events = 0;
  uint64_t negative_events = 0;

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if (jentry % 100000 == 0)
      std::cout << "Processing entry #" << (jentry + 1) << " over " << nentries << " (" << (float) jentry / nentries * 100 << "%)" << std::endl;

    GetEntry(jentry);

    // Choose if we are run2012 A+B, or C+D
    bool isRun2012AB = false;
    if (mIsMC) {
      double r = random_generator.Rndm();
      if (r < lumi_run2012_AB_over_total)
        isRun2012AB = true;

      if (mIsMC) {
        hRunPeriod->Fill( isRun2012AB ? 0 : 1 );
      }
    } else {
      isRun2012AB = (run <= 196531);
    }

    if (generator_weight > 0)
      positive_events++;
    else
      negative_events++;

    // First, fill gen value
    if (mIsMC && (MC_channel == 1 || MC_channel == 2)) {
      TLorentzVector leptonic_W(0., 0., 0., 0.);
      if (gen_lepton_p4->Pt() != 0 && gen_neutrino_p4->Pt() != 0) {
        hDeltaPhiLeptonNeutrino_gen->Fill(fabs(gen_lepton_p4->DeltaPhi(*gen_neutrino_p4)));
        hDeltaRLeptonNeutrino_gen->Fill(gen_lepton_p4->DeltaR(*gen_neutrino_p4));
        hDeltaEtaLeptonNeutrino_gen->Fill(gen_lepton_p4->Eta() - gen_neutrino_p4->Eta());

        leptonic_W = *gen_lepton_p4 + *gen_neutrino_p4;
        hLeptonicWPt_gen->Fill(leptonic_W.Pt());
        hLeptonicWEta_gen->Fill(leptonic_W.Eta());
      }

      TLorentzVector hadronic_W(0., 0., 0., 0.);
      if (gen_lightJet1_p4->Pt() != 0 && gen_lightJet2_p4->Pt() != 0) {
        hDeltaPhiTwoLightJets_gen->Fill(fabs(gen_lightJet1_p4->DeltaPhi(*gen_lightJet2_p4)));
        hDeltaRTwoLightJets_gen->Fill(gen_lightJet1_p4->DeltaR(*gen_lightJet2_p4));
        hDeltaEtaTwoLightJets_gen->Fill(gen_lightJet1_p4->Eta() - gen_lightJet2_p4->Eta());

        hadronic_W = *gen_lightJet1_p4 + *gen_lightJet2_p4;
        hHadronicWPt_gen->Fill(leptonic_W.Pt());
        hHadronicWEta_gen->Fill(leptonic_W.Eta());
      }

      if (leptonic_W.Pt() != 0 && hadronic_W.Pt() != 0) {
        hDeltaPhiW_gen->Fill(fabs(leptonic_W.DeltaPhi(hadronic_W)));
        hDeltaRW_gen->Fill(leptonic_W.DeltaR(hadronic_W));
        hDeltaEtaW_gen->Fill(leptonic_W.Eta() - hadronic_W.Eta());
      }
    }


    double ptLepton = 0;
    double etaLepton = 0;
    double ptLeptonCut = 0;
    double etaSCLepton = 0;
    if (mIsSemiMu) {
      if (nGoodMuons > 0) {
        ptLepton = muonPt[0];
        etaLepton = muonEta[0];
        ptLeptonCut = 25.;
        etaSCLepton = etaLepton;
      }
    } else {
      if (nGoodElectrons > 0) {
        ptLepton = electronPt[0];
        etaLepton = electronEta[0];
        ptLeptonCut = 30.;
        etaSCLepton = electron_SCEta[selectedLeptonIndex_AfterChi2];
      }
    }

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

    hNTrueInt_nosel->Fill(n_trueInteractions, eventWeight);
    hIsSel->Fill(isSel, eventWeight);
    if (mIsMC && (MC_channel != 0)) {
      h_mtt_gen_no_sel->Fill(MC_mtt, generator_weight);
    }

    hMET_nosel->Fill(MET, eventWeight);

    if (!mIsMC && !m_triggerPassed) {
      continue;
    }

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
    for (uint32_t i = 0; i < jet_p4->GetEntriesFast(); i++) {
      float pt = ((TLorentzVector*) (*jet_p4)[i])->Pt();

      HT += pt;
      if (pt > 30)
        HT30 += pt;
    }

    float HTFull = HT + ptLepton + MET;

    hHT_reco_nosel->Fill(HT, eventWeight);
    hHT30_reco_nosel->Fill(HT30, eventWeight);
    hHTFull_reco_nosel->Fill(HTFull, eventWeight);

    //hTopPt_gen->Fill(MC_top1Pt, eventWeight);
    //hTopPt_gen->Fill(MC_top2Pt, eventWeight);
    //hTopPt_com_gen->Fill(MC_top1Pt_com, eventWeight);
    //hTopPt_com_gen->Fill(MC_top2Pt_com, eventWeight);

    //hTopEta_gen->Fill(MC_top1Eta, eventWeight);
    //hTopEta_gen->Fill(MC_top2Eta, eventWeight);
    //hTopEta_com_gen->Fill(MC_top1Eta_com, eventWeight);
    //hTopEta_com_gen->Fill(MC_top2Eta_com, eventWeight);

    if (mIsMC && MC_channel != 0)
      hBoostTT_gen->Fill(MC_boost_tt);

    if (mIsMC && MC_channel != 0) {
      hPtTT_gen->Fill(MC_pt_tt);
      //hPtTT_com_gen->Fill(MC_pt_tt_com, eventWeight);

      hEtaTT_gen->Fill(MC_eta_tt);
      //hEtaTT_com_gen->Fill(MC_eta_tt_com, eventWeight);

      hDeltaPhiTops_gen->Fill(fabs(gen_top1_p4->DeltaPhi(*gen_top2_p4)));
      hDeltaEtaTops_gen->Fill(gen_top1_p4->Eta() - gen_top2_p4->Eta());
      hDeltaRTops_gen->Fill(gen_top1_p4->DeltaR(*gen_top2_p4));
    }

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

    if (mIsSemiMu)
    {
      if (nGoodMuons <= 0)
        continue;
    }
    else
    {
      if (nGoodElectrons <= 0)
        continue;
    }

    if (isSel == 1 && numComb > 0)
    {
      if (mIsMC && MC_channel != 0)
        h_mtt_gen_chi2sel->Fill(MC_mtt, generator_weight);

      hNGoodJets_chi2sel->Fill(nJets, eventWeight);
      hNBtaggedJets_chi2sel->Fill(nBtaggedJets_CSVM, eventWeight);

      hLeptonPt_chi2sel->Fill(ptLepton, eventWeight);
      hLeptonEta_chi2sel->Fill(etaLepton, eventWeight);
      if (mIsSemiMu)
        hMuRelIso_chi2sel->Fill(muRelIso[0], eventWeight);
      else
        hElRelIso_chi2sel->Fill(elRelIso[0], eventWeight);
      
      hMET_chi2sel->Fill(MET, eventWeight);

      hLeptTopPt_chi2sel->Fill(lepTopPt_AfterChi2, eventWeight);
      hLeptTopEta_chi2sel->Fill(lepTopEta_AfterChi2, eventWeight);

      hHadrTopPt_chi2sel->Fill(hadTopPt_AfterChi2, eventWeight);
      hHadrTopEta_chi2sel->Fill(hadTopEta_AfterChi2, eventWeight);

      hFirstJetPt_chi2sel->Fill(p_1stjetpt, eventWeight);
      hSecondJetPt_chi2sel->Fill(p_2ndjetpt, eventWeight);
      hThirdJetPt_chi2sel->Fill(p_3rdjetpt, eventWeight);
      hFourthJetPt_chi2sel->Fill(p_4thjetpt, eventWeight);

      hFirstJetEta_chi2sel->Fill(jetEta[0], eventWeight);
      hSecondJetEta_chi2sel->Fill(jetEta[1], eventWeight);
      hThirdJetEta_chi2sel->Fill(jetEta[2], eventWeight);
      hFourthJetEta_chi2sel->Fill(jetEta[3], eventWeight);

      hNVtx_chi2sel->Fill(n_vertices, eventWeight);

      hBoostTT_chi2sel->Fill(beta_tt_AfterChi2, eventWeight);
      hPtTT_chi2sel->Fill(pt_tt_AfterChi2, eventWeight);
      hEtaTT_chi2sel->Fill(eta_tt_AfterChi2, eventWeight);

      bool btagSel = false;
      if (mBTag == 1)
        btagSel = nBtaggedJets_CSVM == 1;
      else if (mBTag == 2)
        btagSel = nBtaggedJets_CSVM > 1;

      float firstJetCut = 0, secondJetCut = 0, thirdJetCut = 0;
      if (isRun2012AB) {
        firstJetCut = 45;
        secondJetCut = 45;
        thirdJetCut = 45;
      } else {
        firstJetCut = 55;
        secondJetCut = 45;
        thirdJetCut = 35;
      }

      if (ptLepton > ptLeptonCut && btagSel && p_1stjetpt > firstJetCut && p_2ndjetpt > secondJetCut && p_3rdjetpt > thirdJetCut && mtt_AfterChi2 > 0. && (isinf(bestSolChi2) || bestSolChi2 < 500))
      {

        if (mIsMC) {
          // Compute trigger weight
          double triggerWeight = m_trigger_efficiency_provider->get_weight(ptLepton, etaSCLepton, p_4thjetpt, jetEta[3], n_vertices, nJets, mIsSemiMu, TopTriggerEfficiencyProvider::NOMINAL)[0];
          hTriggerWeight->Fill(triggerWeight);

          eventWeight *= triggerWeight;

          hWeight_fullsel->Fill(eventWeight);
        }

        hBoostTT->Fill(beta_tt_AfterChi2, eventWeight);
        hPtTT->Fill(pt_tt_AfterChi2, eventWeight);
        hEtaTT->Fill(eta_tt_AfterChi2, eventWeight);

        if (mIsMC && MC_channel != 0) {
          h_mtt_gen->Fill(MC_mtt, generator_weight);
          h_mtt_resolution->Fill(mtt_AfterChi2 - MC_mtt, eventWeight);
        }

        hLeptonPt->Fill(ptLepton, eventWeight);
        hLeptonEta->Fill(etaLepton, eventWeight);
        hFirstJetPt->Fill(p_1stjetpt, eventWeight);
        hSecondJetPt->Fill(p_2ndjetpt, eventWeight);
        hThirdJetPt->Fill(p_3rdjetpt, eventWeight);
        hFourthJetPt->Fill(p_4thjetpt, eventWeight);
        hFirstJetEta->Fill(jetEta[0], eventWeight);
        hSecondJetEta->Fill(jetEta[1], eventWeight);
        hThirdJetEta->Fill(jetEta[2], eventWeight);
        hFourthJetEta->Fill(jetEta[3], eventWeight);
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

        hDeltaPhiTops_reco_fullsel->Fill(fabs(lepTopP4_AfterChi2->DeltaPhi(*hadTopP4_AfterChi2)), eventWeight);
        hDeltaEtaTops_reco_fullsel->Fill(lepTopP4_AfterChi2->Eta() - hadTopP4_AfterChi2->Eta(), eventWeight);
        hDeltaRTops_reco_fullsel->Fill(lepTopP4_AfterChi2->DeltaR(*hadTopP4_AfterChi2), eventWeight);

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
  loadChain(inputFiles, "electron_PF", fElectrons);


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

  if (fElectrons)
    fElectrons->GetEntry(entry);

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

  SetBranchAddress(fMTT, "MC_channel", &MC_channel, NULL);
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
  SetBranchAddress(fMTT, "numComb_chi2", &numComb, &b_numComb);
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
  SetBranchAddress(fMTT, "lepton_weight", &m_lepton_weight, NULL);
  SetBranchAddress(fMTT, "btag_weight", &m_btag_weight, NULL);
  SetBranchAddress(fMTT, "selectedLeptonIndex_AfterChi2", &selectedLeptonIndex_AfterChi2, NULL);

  if (fMTT->GetBranch("trigger_passed")) {
    SetBranchAddress(fMTT, "trigger_passed", &m_triggerPassed, NULL);
  } else {
    // Backward compatibilty
    m_triggerPassed = true;
  }

  lepTopP4_AfterChi2 = NULL;
  hadTopP4_AfterChi2 = NULL;

  fMTT->SetBranchStatus("lepTopP4_AfterChi2*", 1);
  fMTT->SetBranchStatus("hadTopP4_AfterChi2*", 1);
  fMTT->SetBranchAddress("lepTopP4_AfterChi2.", &lepTopP4_AfterChi2);
  fMTT->SetBranchAddress("hadTopP4_AfterChi2.", &hadTopP4_AfterChi2);

  gen_top1_p4 = NULL;
  gen_top2_p4 = NULL;

  fMTT->SetBranchStatus("MC_Top1_p4*", 1);
  fMTT->SetBranchStatus("MC_Top2_p4*", 1);
  fMTT->SetBranchAddress("MC_Top1_p4.", &gen_top1_p4);
  fMTT->SetBranchAddress("MC_Top2_p4.", &gen_top2_p4);

  gen_lepton_p4 = NULL;
  gen_neutrino_p4 = NULL;

  fMTT->SetBranchStatus("MC_lepton_p4*", 1);
  fMTT->SetBranchStatus("MC_neutrino_p4*", 1);
  fMTT->SetBranchAddress("MC_lepton_p4.", &gen_lepton_p4);
  fMTT->SetBranchAddress("MC_neutrino_p4.", &gen_neutrino_p4);

  gen_leptonic_B_p4 = NULL;
  gen_hadronic_B_p4 = NULL;

  fMTT->SetBranchStatus("MC_leptonic_B_p4*", 1);
  fMTT->SetBranchStatus("MC_hadronic_B_p4*", 1);
  fMTT->SetBranchAddress("MC_leptonic_B_p4.", &gen_leptonic_B_p4);
  fMTT->SetBranchAddress("MC_hadronic_B_p4.", &gen_hadronic_B_p4);

  gen_lightJet1_p4 = NULL;
  gen_lightJet2_p4 = NULL;

  fMTT->SetBranchStatus("MC_lightJet1_B_p4*", 1);
  fMTT->SetBranchStatus("MC_lightJet2_B_p4*", 1);
  fMTT->SetBranchAddress("MC_lightJet1_B_p4.", &gen_lightJet1_p4);
  fMTT->SetBranchAddress("MC_lightJet2_B_p4.", &gen_lightJet2_p4);

  fVertices->SetMakeClass(1);
  fVertices->SetBranchStatus("*", 0);
  fVertices->SetBranchStatus("n_vertices", 1);
  fVertices->SetBranchAddress("n_vertices", &n_vertices, NULL);

  fEvent->SetMakeClass(1);
  fEvent->SetBranchStatus("*", 0);
  fEvent->SetBranchStatus("nTrueInteractions", 1);
  fEvent->SetBranchAddress("nTrueInteractions", &n_trueInteractions, NULL);

  if (fEvent->GetBranch("generator_weight")) {
    fEvent->SetBranchStatus("generator_weight", 1);
    fEvent->SetBranchAddress("generator_weight", &generator_weight, NULL);
  } else {
    generator_weight = 1;
  }

  run = 0;
  fEvent->SetBranchStatus("run", 1);
  fEvent->SetBranchAddress("run", &run, NULL);

  muon_p4 = NULL;
  fLooseMuons->SetMakeClass(1);
  fLooseMuons->SetBranchStatus("*", 0);
  fLooseMuons->SetBranchStatus("muon_4vector", 1);
  fLooseMuons->SetBranchStatus("n_muons", 1);
  fLooseMuons->SetBranchStatus("muon_deltaBetaCorrectedRelIsolation", 1);
  fLooseMuons->SetBranchAddress("muon_4vector", &muon_p4, NULL);
  fLooseMuons->SetBranchAddress("n_muons", &n_muons, NULL);
  fLooseMuons->SetBranchAddress("muon_deltaBetaCorrectedRelIsolation", &muon_relIso, NULL);

  fElectrons->SetBranchStatus("*", 0);
  fElectrons->SetBranchStatus("electron_SCEta", 1);
  fElectrons->SetBranchAddress("electron_SCEta", &electron_SCEta);

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
