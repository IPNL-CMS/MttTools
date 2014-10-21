#include "Extractor2Histos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <THStack.h>
#include <TRandom2.h>
#include <vector>
#include <TProfile.h>
#include <TMath.h>
#include <fstream>
#include <memory>
#include <chrono>
#include <ctime>
#include <TLorentzVector.h>
#include <Math/GenVector/VectorUtil.h>

#include <boost/filesystem.hpp>

#include "GaussianProfile.h"
#include "TopTriggerEfficiencyProvider.h"

#include <PUReweighting/PUReweighter.h>
#include <BkgVsTTBDTReader.h>

#include "tclap/CmdLine.h"

#include "TMVA/Reader.h"

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

const int nBins = 12;
const double bins[] = {340, 400, 450, 500, 550, 600, 650, 750, 850, 950, 1050, 1200, 1400};

//const int nBins = 10;
//const double bins[] = {340, 400, 450, 500, 600, 700, 800, 900, 1000, 1200, 1400};

TLorentzVector* getP4(TClonesArray* array, int index) {
  return (TLorentzVector*) array->At(index);
}

void Extractor2Histos::Loop()
{
  // Create the output file first in order that each histogram
  // associate itself with this file
  TFile * output = TFile::Open(mOutputFile.c_str(), "recreate");
  output->cd();

  TH1::SetDefaultSumw2(true);

  TH1D *hWeight = new TH1D("weight", "", 100, 0, 5);
  TH1D *hWeight_fullsel = new TH1D("weight_fullsel", "", 100, 0, 2);
  TH1D *hLeptonWeight = new TH1D("lepton_weight", "", 100, 0, 2);
  TH1D *hLeptonWeight_fullsel = new TH1D("lepton_weight_fullsel", "", 100, 0.5, 2);
  TH1D *hBTagWeight = new TH1D("btag_weight", "", 100, 0, 2);
  TH1D *hBTagWeight_fullsel = new TH1D("btag_weight_fullsel", "", 100, 0.5, 2);
  TH1D *hTriggerWeight_fullsel = new TH1D("trigger_weight_fullsel", "", 100, 0, 1);
  TH1D *hPUWeight = new TH1D("PU_weight", "", 100, 0, 2);
  TH1D *hPUWeight_fullsel = new TH1D("PU_weight_fullsel", "", 100, 0, 2);
  TH1D *hTopPtWeight = new TH1D("top_pt_weight", "", 100, 0, 2);
  TH1D *hTopPtWeight_fullsel = new TH1D("top_pt_weight_fullsel", "", 100, 0, 2);
  TH1D *hGeneratorWeight = new TH1D("generator_weight", "", 100, -2, 2);
  TH1D *hGeneratorWeight_fullsel = new TH1D("generator_weight_fullsel", "", 100, -2, 2);

  TH1D *hRunPeriod = new TH1D("run_period", "", 2, 0, 2);
  hRunPeriod->GetXaxis()->SetBinLabel(1, "Run2012 A+B");
  hRunPeriod->GetXaxis()->SetBinLabel(2, "Run2012 C+D");

  TH1D *hNVtx_noweight = new TH1D("nVertex_reco_fullsel_noweight", "", 70, 0., 70);
  TH1D *hNVtx = new TH1D("nVertex_reco_fullsel", "", 70, 0., 70);
  TH1D *hNVtx_chi2sel = new TH1D("nVertex_reco_chi2sel", "", 70, 0., 70);
  TH1D *hNVtx_nosel = new TH1D("nVertex_reco_nosel", "", 70, 0., 70);

  TH1D *hIsSel = new TH1D("isSel", "", 10, 0, 10);
  TH1D *hBestSolChi2 = new TH1D("bestSolChi2", "", 400, 0, 1000);
  TH1D *hBestSolChi2_fullsel = new TH1D("bestSolChi2_fullsel", "", 400, 0, 500);
  TH1D *hBestSolChi2Exp_fullsel = new TH1D("bestSolChi2_exp_fullsel", "", 100, 0., 1.);
  hBestSolChi2Exp_fullsel->SetXTitle("exp(-#chi^{2})");
  TH1D *hBestSolChi2Proba_fullsel = new TH1D("bestSolChi2Proba_fullsel", "", 100, 0., 1.);
  hBestSolChi2Proba_fullsel->SetXTitle("#chi^{2} prob");

  TH1D *hNTrueInt = new TH1D("nTrueInt_reco_fullsel", "", 70, 0., 70);
  TH1D *hNTrueInt_nosel = new TH1D("nTrueInt_reco_nosel", "", 70, 0., 70);

  TH1D *hLeptonPt = new TH1D("leptonPt_reco_fullsel", "", 50, 20., 200.);
  TH1D *hLeptonPt_chi2sel = new TH1D("leptonPt_reco_chi2sel", "", 50, 20., 200.);
  TH1D *hLeptonPt_nosel = new TH1D("leptonPt_reco_nosel", "", 100, 0., 200.);

  TH1D *hLeptonEta = new TH1D("leptonEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hLeptonEta_chi2sel = new TH1D("leptonEta_reco_chi2sel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hLeptonEta_nosel = new TH1D("leptonEta_reco_nosel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hLeptTopPt = new TH1D("leptTopPt_reco_fullsel", "", 60, 20., 600.);
  TH1D *hLeptTopMt = new TH1D("leptTopMt_reco_fullsel", "", 60, 20., 600.);
  TH1D *hLeptTopPt_chi2sel = new TH1D("leptTopPt_reco_chi2sel", "", 60, 20., 600.);
  TH1D *hLeptTopPz = new TH1D("leptTopPz_reco_fullsel", "", 60, 20., 600.);
  TH1D *hLeptTopE = new TH1D("leptTopE_reco_fullsel", "", 60, 20., 600.);
  hLeptTopE->SetXTitle("leptonic top E [GeV/c]");
  hLeptTopMt->SetXTitle("leptonic top M_{t} [GeV/c]");
  TH1D *hLeptTopPhi = new TH1D("leptTopPhi_reco_fullsel", "", 200, -4., 4.);
  hLeptTopPhi->SetXTitle("leptonic top #phi");

  TH1D *hBoostTT_gen = new TH1D("boostTT_gen", "", 50, 0., 1.);
  TH1D *hPtTT_gen = new TH1D("ptTT_gen", "", 60, 0., 600.);
  TH1D *hEtaTT_gen = new TH1D("etaTT_gen", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hLeptTopEta = new TH1D("leptTopEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hLeptTopEta_chi2sel = new TH1D("leptTopEta_reco_chi2sel", "", 50, -2*M_PI, 2*M_PI);

  TH1D *hHadrTopPt = new TH1D("hadrTopPt_reco_fullsel", "", 60, 20., 600.);
  TH1D *hHadrTopMt = new TH1D("hadrTopMt_reco_fullsel", "", 60, 20., 600.);
  TH1D *hHadrTopPt_chi2sel = new TH1D("hadrTopPt_reco_chi2sel", "", 60, 20., 600.);
  TH1D *hHadrTopPz = new TH1D("hadrTopPz_reco_fullsel", "", 60, 20., 600.);
  TH1D *hHadrTopE = new TH1D("hadrTopE_reco_fullsel", "", 60, 20., 600.);
  hHadrTopE->SetXTitle("hadronic top E [GeV/C]");
  hHadrTopMt->SetXTitle("hadronic top M_{t} [GeV/C]");
  TH1D *hHadrTopPhi = new TH1D("hadrTopPhi_reco_fullsel", "", 200, -4., 4.);
  hHadrTopPhi->SetXTitle("hadronic top #phi");

  TH1D *hTopsDeltaY = new TH1D("TopsDeltaPhi_reco_fullsel", "", 200, -400., 400.);
  hTopsDeltaY->SetXTitle("Top #Delta y");

  TH1D *hHadrTopEta = new TH1D("hadrTopEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hHadrTopEta_chi2sel = new TH1D("hadrTopEta_reco_chi2sel", "", 50, -2*M_PI, 2*M_PI);

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
  TH1D *hMETx = new TH1D("METx_reco_fullsel", "", 100, -200., 200.);
  TH1D *hMETy = new TH1D("METy_reco_fullsel", "", 100, -200., 200.);
  TH1D *hMETPhi = new TH1D("METPhi_reco_fullsel", "", 800, -4., 4.);
  TH1D *hMET_nosel = new TH1D("MET_reco_nosel", "", 100, 0., 400.);
  TH1D *hMET_chi2sel = new TH1D("MET_reco_chi2sel", "", 100, 20., 400.);

  TH1D *hmtlep = new TH1D("mtLep_reco_fullsel", "", 100, 120., 240.);
  TH1D *hmthad = new TH1D("mtHad_reco_fullsel", "", 150, 120., 300.);

  TH1D *hmttSelected_btag_sel_positive = new TH1D("mttSelected_btag_sel_reco_fullsel_positive", "", 150, 0., 2000.);
  TH1D *hmttSelected_btag_sel_negative = new TH1D("mttSelected_btag_sel_reco_fullsel_negative", "", 150, 0., 2000.);

  TH1D *hmtt_AfterChi2 = new TH1D("mtt_AfterChi2", "", 150, 0., 2000.);
  TH1D *hmtt_AfterKF = new TH1D("mtt_AfterKF", "", 150, 0., 2000.);
  TH1D *hmttSelected_btag_sel = new TH1D("mttSelected_btag_sel_reco_fullsel", "", 150, 0., 2000.);
  TH1D *hmttSelected_btag_sel_no_gen_weight = new TH1D("mttSelected_btag_sel_reco_fullsel_no_gen_weight", "", 150, 0., 2000.);
  TH1D *hmttSelected_btag_sel_mass_cut = new TH1D("mttSelected_btag_sel_mass_cut_reco_fullsel", "", 150, 0., 2000.);

  TH1D *hSelectedFirstJetPt = new TH1D("selectedFirstJetPt_reco_fullsel", "", 100, 70., 640.);
  TH1D *hSelectedSecondJetPt = new TH1D("selectedSecondJetPt_reco_fullsel", "", 50, 30., 300.);
  TH1D *hSelectedHadronicBPt = new TH1D("selectedHadronicBPt_reco_fullsel", "", 100, 70., 640.);
  TH1D *hSelectedLeptonicBPt = new TH1D("selectedLeptonicBPt_reco_fullsel", "", 100, 70., 640.);
  TH1D *hSelectedLeptonPt = new TH1D("selectedLeptonPt_reco_fullsel", "", 50, 20., 200.);
  TH1D *hSelectedNeutrinoPt = new TH1D("selectedNeutrinoPt_reco_fullsel", "", 50, 20., 200.);
  TH1D *hSelectedNeutrinoPz = new TH1D("selectedNeutrinoPz_reco_fullsel", "", 50, 20., 200.);

  TH1D *hSelectedFirstJetEta = new TH1D("selectedFirstJetEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hSelectedSecondJetEta = new TH1D("selectedSecondJetEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hSelectedHadronicBEta = new TH1D("selectedHadronicBEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hSelectedLeptonicBEta = new TH1D("selectedLeptonicBEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hSelectedLeptonEta = new TH1D("selectedLeptonEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);
  TH1D *hSelectedNeutrinoEta = new TH1D("selectedNeutrinoEta_reco_fullsel", "", 100, -2*M_PI, 2*M_PI);

  TH1D *hSelectedFirstJetPhi = new TH1D("selectedFirstJetPhi_reco_fullsel", "", 200, -4., 4.);
  TH1D *hSelectedSecondJetPhi = new TH1D("selectedSecondJetPhi_reco_fullsel", "", 200, -4., 4.);
  TH1D *hSelectedHadronicBPhi = new TH1D("selectedHadronicBPhi_reco_fullsel", "", 200, -4., 4.);
  TH1D *hSelectedLeptonicBPhi = new TH1D("selectedLeptonicBPhi_reco_fullsel", "", 200, -4., 4.);
  TH1D *hSelectedLeptonPhi = new TH1D("selectedLeptonPhi_reco_fullsel", "", 200, -4., 4.);
  TH1D *hSelectedNeutrinoPhi = new TH1D("selectedNeutrinoPhi_reco_fullsel", "", 200, -4., 4.);

  TH1D *hSelectedFirstJetE = new TH1D("selectedFirstJetE_reco_fullsel", "", 100, 70., 640.);
  TH1D *hSelectedSecondJetE = new TH1D("selectedSecondJetE_reco_fullsel", "", 50, 30., 300.);
  TH1D *hSelectedHadronicBE = new TH1D("selectedHadronicBE_reco_fullsel", "", 100, 70., 640.);
  TH1D *hSelectedLeptonicBE = new TH1D("selectedLeptonicBE_reco_fullsel", "", 100, 70., 640.);
  TH1D *hSelectedLeptonE = new TH1D("selectedLeptonE_reco_fullsel", "", 50, 20., 200.);
  TH1D *hSelectedNeutrinoE = new TH1D("selectedNeutrinoE_reco_fullsel", "", 50, 20., 200.);

  TH1* hBDTDiscriminant = new TH1D("bdt_discriminant", "", 50, -1, 1);

  TH1D *hkf_chisquare = new TH1D("kf_chisquare_fullsel", "", 1800, 0., 450.);
  TH1D *hkf_proba = new TH1D("kf_proba_fullsel", "", 200, 0., 1.);
  TH1D *hkf_proba_zoom = new TH1D("kf_proba_fullsel_zoom", "", 100, 0., 0.1);

  // Shape variables
  TH1D *hAplanarity = new TH1D("aplanarity_reco_fullsel", "", 50, 0., 0.5);
  TH1D *hCircularity = new TH1D("circularity_reco_fullsel", "", 50, 0., 1.);
  TH1D *hSphericity = new TH1D("sphericity_reco_fullsel", "", 50, 0., 1.);
  TH1D *hIsotropy = new TH1D("isotropy_reco_fullsel", "", 50, 0., 1.);
  TH1D *hD = new TH1D("D_reco_fullsel", "", 50, 0., 1.);
  TH1D *hC = new TH1D("C_reco_fullsel", "", 50, 0., 1.);

  hAplanarity->SetXTitle("Aplanarity");
  hCircularity->SetXTitle("Circularity");
  hSphericity->SetXTitle("Spericity");
  hIsotropy->SetXTitle("Isotropy");
  hD->SetXTitle("D");
  hC->SetXTitle("C");

  hkf_chisquare->SetXTitle("#Chi^{2}_{KF}");
  hkf_proba->SetXTitle("Proba_{KF}");
  hkf_proba_zoom->SetXTitle("Proba_{KF}");

  if (mIsSemiMu) {
    hSelectedLeptonPt->SetXTitle("selected #mu p_{T} [GeV/c]");
    hSelectedLeptonEta->SetXTitle("selected #mu #eta");
    hSelectedLeptonPhi->SetXTitle("selected #mu #phi");
    hSelectedLeptonE->SetXTitle("selected #mu E [GeV/c]");
  } else {
    hSelectedLeptonPt->SetXTitle("selected e p_{T} [GeV/c]");
    hSelectedLeptonEta->SetXTitle("selected e #eta");
    hSelectedLeptonPhi->SetXTitle("selected e #phi");
    hSelectedLeptonE->SetXTitle("selected e E [GeV/c]");
  }
  hSelectedFirstJetPt->SetXTitle("selected 1^{st} jet p_{T} [GeV/c]");
  hSelectedSecondJetPt->SetXTitle("selected 2^{nd} jet p_{T} [GeV/c]");
  hSelectedHadronicBPt->SetXTitle("selected had B p_{T} [GeV/c]");
  hSelectedLeptonicBPt->SetXTitle("selected lep B p_{T} [GeV/c]");
  hSelectedNeutrinoPt->SetXTitle("selected #nu p_{T} [GeV/c]");
  hSelectedNeutrinoPz->SetXTitle("selected #nu p_{Z} [GeV/c]");

  hSelectedFirstJetEta->SetXTitle("selected 1^{st} jet #eta");
  hSelectedSecondJetEta->SetXTitle("selected 2^{nd} jet #eta");
  hSelectedHadronicBEta->SetXTitle("selected had B #eta");
  hSelectedLeptonicBEta->SetXTitle("selected lep B #eta");
  hSelectedNeutrinoEta->SetXTitle("selected #nu #eta");

  hSelectedFirstJetPhi->SetXTitle("selected 1^{st} jet #phi");
  hSelectedSecondJetPhi->SetXTitle("selected 2^{nd} jet #phi");
  hSelectedHadronicBPhi->SetXTitle("selected had B #phi");
  hSelectedLeptonicBPhi->SetXTitle("selected lep B #phi");
  hSelectedNeutrinoPhi->SetXTitle("selected #nu #phi");

  hSelectedFirstJetE->SetXTitle("selected 1^{st} jet E [GeV/c]");
  hSelectedSecondJetE->SetXTitle("selected 2^{nd} jet E [GeV/c]");
  hSelectedHadronicBE->SetXTitle("selected had B E [GeV/c]");
  hSelectedLeptonicBE->SetXTitle("selected lep B E [GeV/c]");
  hSelectedNeutrinoE->SetXTitle("selected #nu E [GeV/c]");

  hmtt_AfterChi2->SetXTitle("m_{t#bar{t}}^{Chi2} [GeV/c^{2}]");
  hmtt_AfterKF->SetXTitle("m_{t#bar{t}}^{KF} [GeV/c^{2}]");

  TH1D *hNGoodJets = new TH1D("nGoodJets_reco_fullsel", "", 6, 3.5, 9.5);
  TH1D *hNGoodJets_chi2sel = new TH1D("nGoodJets_reco_chi2sel", "", 6, 3.5, 9.5);

  TH1D *hNBtaggedJets = new TH1D("nBTaggedJets_reco_fullsel", "", 5, -0.5, 4.5);
  TH1D *hNBtaggedJets_chi2sel = new TH1D("nBTaggedJets_reco_chi2sel", "", 5, -0.5, 4.5);

  TH1D *h_mtt_gen_nosel = new TH1D("mtt_gen_nosel", "", 300, 0., 1500.);

  TH1D *h_mtt_resolution_AfterChi2 = new TH1D("mtt_resolution_AfterChi2", "", 100, -600., 600.);
  TH1D *h_mtt_resolution_AfterKF = new TH1D("mtt_resolution_AfterKF", "", 100, -600., 600.);
  TH1D *h_mtt_resolution = new TH1D("mtt_resolution", "", 100, -600., 600.);
  TH1D *h_mtt_resolution_four_jets = new TH1D("mtt_resolution_four_jets", "", 100, -600., 600.);

  GaussianProfile mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel("mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel", nBins, bins);
  GaussianProfile mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel("mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel", nBins, bins, 100, -1.2, 1.2);
  GaussianProfile mtt_gen_vs_mtt_reco_resolution_rmsVSmttreco_sortingAlgoSel("mtt_gen_vs_mtt_reco_resolution_rmsVSmttreco_sortingAlgoSel", nBins, bins, 100, -1.2, 1.2);

  GaussianProfile mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel_lowKfProb("mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel_lowKfProb", nBins, bins);
  GaussianProfile mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel_lowKfProb("mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel_lowKfProb", nBins, bins, 100, -1.2, 1.2);

  GaussianProfile mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel_highKfProb("mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel_highKfProb", nBins, bins);
  GaussianProfile mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel_highKfProb("mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel_highKfProb", nBins, bins, 100, -1.2, 1.2);

  GaussianProfile mtt_gen_vs_mtt_reco_linearity_fullSel("mtt_gen_vs_mtt_reco_linearity_fullSel", nBins, bins);
  GaussianProfile mtt_gen_vs_mtt_reco_resolution_fullSel("mtt_gen_vs_mtt_reco_resolution_fullSel", nBins, bins, 100, -1.2, 1.2);

  //GaussianProfile mtt_gen_vs_mtt_reco_chi2sel("mtt_gen_vs_mtt_reco_chi2sel", nBins, bins);
  //GaussianProfile mtt_gen_vs_mtt_reso_chi2sel("mtt_gen_vs_mtt_reso_chi2sel", nBins, bins, 100, -1.2, 1.2);

  GaussianProfile mtt_gen_vs_mtt_reco_linearity_fourJets_sortingAlgoSel("mtt_gen_vs_mtt_reco_linearity_fourJets_sortingAlgoSel", nBins, bins);
  GaussianProfile mtt_gen_vs_mtt_reco_resolution_fourJets_sortingAlgoSel("mtt_gen_vs_mtt_reco_resolution_fourJets_sortingAlgoSel", nBins, bins, 100, -1.2, 1.2);

  GaussianProfile mtt_gen_vs_mtt_reco_linearity_fourJets_fullSel("mtt_gen_vs_mtt_reco_linearity_fourJets_fullSel", nBins, bins);
  GaussianProfile mtt_gen_vs_mtt_reco_resolution_fourJets_fullSel("mtt_gen_vs_mtt_reco_resolution_fourJets_fullSel", nBins, bins, 100, -1.2, 1.2);

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

  TH1D *hDeltaThetaTops_reco_fullsel = new TH1D("deltaThetaTops_reco_fullsel", "", 50, 0, M_PI);

  /******
   * Leptonic side
   */
  TH1* hDeltaPhiLeptonNeutrino_gen = new TH1D("deltaPhiLeptonNeutrino_gen", "", 50, 0, M_PI);

  TH1* hDeltaRLeptonNeutrino_gen = new TH1D("deltaRLeptonNeutrino_gen", "", 50, 0, 10);

  TH1* hDeltaEtaLeptonNeutrino_gen = new TH1D("deltaEtaLeptonNeutrino_gen", "", 50, -3*M_PI, 3*M_PI);

  // W
  TH1* hLeptonicWPt_gen = new TH1D("ptWLeptonic_gen", "", 60, 0, 600.);
  TH1* hLeptonicWMt = new TH1D("mtWLeptonic_reco_fullsel", "", 60, 0, 600.);
  hLeptonicWMt->SetXTitle("Leptonic W M_{t} [GeV/c]");
  TH1* hLeptonicWEta_gen = new TH1D("etaWLeptonic_gen", "", 50, -3*M_PI, 3*M_PI);

  /*****
   * Hadronic side
   */
  TH1* hDeltaPhiTwoLightJets_gen = new TH1D("deltaPhiTwoLightJets_gen", "", 50, 0, M_PI);
  TH1* hDeltaRTwoLightJets_gen = new TH1D("deltaRTwoLightJets_gen", "", 50, 0, 10);
  TH1* hDeltaEtaTwoLightJets_gen = new TH1D("deltaEtaTwoLightJets_gen", "", 50, -3*M_PI, 3.*M_PI);

  // W
  TH1* hHadronicWPt_gen = new TH1D("ptWHadronic_gen", "", 60, 0, 600.);
  TH1* hHadronicWMt = new TH1D("mtWHadronic_reco_fullsel", "", 60, 0, 600.);
  hHadronicWMt->SetXTitle("Hadronic W M_{t} [GeV/c]");
  TH1* hHadronicWEta_gen = new TH1D("etaWHadronic_gen", "", 50, -3*M_PI, 3*M_PI);


  TH1* hDeltaPhiW_gen = new TH1D("deltaPhiW_gen", "", 50, 0, M_PI);
  TH1* hDeltaRW_gen = new TH1D("deltaRW_gen", "", 50, 0, 10);
  TH1* hDeltaEtaW_gen = new TH1D("deltaEtaW_gen", "", 50, -3*M_PI, 3.*M_PI);


  TH1* hHT_reco_nosel = new TH1D("HT_reco_nosel", "", 100, 0, 1000);
  TH1* hHT30_reco_nosel = new TH1D("HT30_reco_nosel", "", 100, 0, 1000);

  TH1* hHTFull_reco_nosel = new TH1D("HTFull_reco_nosel", "", 100, 0, 1000);

  TH1* hHTFrac = new TH1D("HTFrac", "", 100, 0, 1);

  if (mIsSemiMu) {
    hLeptonPt->SetXTitle("#mu p_{T} [GeV/c]");
    hLeptonPt_chi2sel->SetXTitle("#mu p_{T} [GeV/c]");
  } else {
    hLeptonPt_chi2sel->SetXTitle("e p_{T} [GeV/c]");
    hLeptonPt->SetXTitle("e p_{T} [GeV/c]");
  }

  hLeptTopPz->SetXTitle("leptonic top p_{Z} [GeV/c]");
  hHadrTopPz->SetXTitle("hadronic top p_{Z} [GeV/c]");

  hFirstJetPt->SetXTitle("1^{st} jet p_{T} [GeV/c]");
  hFirstJetPt_nosel->SetXTitle("1^{st} jet p_{T} [GeV/c]");
  hFirstJetPt_chi2sel->SetXTitle("1^{st} jet p_{T} [GeV/c]");

  hSecondJetPt->SetXTitle("2^{nd} jet p_{T} [GeV/c]");
  hSecondJetPt_chi2sel->SetXTitle("2^{nd} jet p_{T} [GeV/c]");

  hThirdJetPt->SetXTitle("3^{rd} jet p_{T} [GeV/c]");
  hThirdJetPt_chi2sel->SetXTitle("3^{rd} jet p_{T} [GeV/c]");

  hFourthJetPt_chi2sel->SetXTitle("4^{th} jet p_{T} [GeV/c]");

  hMET->SetXTitle("MET [GeV]");
  hMETx->SetXTitle("MET_{X} [GeV]");
  hMETy->SetXTitle("MET_{Y} [GeV]");
  hMETPhi->SetXTitle("MET #phi [GeV]");
  hMET_chi2sel->SetXTitle("MET [GeV]");

  hmtlep->SetXTitle("leptonic m_{t} [GeV/c^{2}]");
  hmthad->SetXTitle("hadronic m_{t} [GeV/c^{2}]");

  h_mtt_gen_nosel->SetXTitle("m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  h_mtt_resolution->SetXTitle("m_{t#bar{t}}^{reco} - m_{t#bar{t}}^{gen} [GeV/c^{2}]");

  hmttSelected_btag_sel->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hmttSelected_btag_sel_mass_cut->SetXTitle("m_{t#bar{t}} [GeV/c^{2}]");
  hNGoodJets->SetXTitle("Num good jets");
  hNGoodJets_chi2sel->SetXTitle("Num good jets");

  hNBtaggedJets->SetXTitle("Num CSVM jets");
  hNBtaggedJets_chi2sel->SetXTitle("Num CSVM jets");

  Long64_t nentries = fMTT->GetEntries();
  //nentries = 10000;

  //PUReweighter puReweighter(mIsSemiMu, mDataset);

  Systematic puSyst = Systematic::NOMINAL;
  if (mPUSyst == "up")
    puSyst = Systematic::UP;
  else if (mPUSyst == "down")
    puSyst = Systematic::DOWN;

  PUReweighter puReweighter(mIsSemiMu, PUProfile::S10, puSyst);

  TopTriggerEfficiencyProvider::JES triggerJESSyst = TopTriggerEfficiencyProvider::NOMINAL;
  if (mJECSyst == "up")
    triggerJESSyst = TopTriggerEfficiencyProvider::UP;
  else if (mJECSyst == "down")
    triggerJESSyst = TopTriggerEfficiencyProvider::DOWN;

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

  BkgVsTTBDTReader bkgVsTTBDTReader(m_inputFiles);
  bkgVsTTBDTReader.initMVA(m_bdtWeights);

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

    if (generator_weight > 0) {
      positive_events++;
      generator_weight = 1;
    } else {
      negative_events++;
      generator_weight = -1;
    }

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
    double eventWeightNoGenWeight = 1.;
    float puWeight = 1.;
    float topPtWeight = 1.;
    if (mIsMC) {
      puWeight = puReweighter.weight(n_trueInteractions);

      if (mLeptonSyst == "up")
        m_lepton_weight += m_lepton_weight_error;
      else if (mLeptonSyst == "down")
        m_lepton_weight -= m_lepton_weight_error;

      if (mBTagSyst == "up")
        m_btag_weight += m_btag_weight_error;
      else if (mBTagSyst == "down")
        m_btag_weight -= m_btag_weight_error;

      hLeptonWeight->Fill(m_lepton_weight);
      hBTagWeight->Fill(m_btag_weight);
      hGeneratorWeight->Fill(generator_weight);
      hPUWeight->Fill(puWeight);

      eventWeight *= puWeight;
      eventWeight *= m_lepton_weight;
      eventWeight *= m_btag_weight;

      // Top pt reweighting: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
      if (MC_channel != 0) {
        auto SF = [](TLorentzVector* top) {
          if (top->Pt() > 400)
            return 1.;

          return std::exp(0.159 - 0.00141 * top->Pt());
        };

        topPtWeight = std::sqrt(SF(getP4(gen_top1_p4, 0)) * SF(getP4(gen_top2_p4, 0)));
        //hTopPtWeight->Fill(topPtWeight);

        //eventWeight *= topPtWeight;
      }

      eventWeightNoGenWeight = eventWeight;
      eventWeight *= generator_weight;
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

    if (!mIsMC && !m_triggerPassed) {
      continue;
    }

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

    float HT = 0;
    float HT30 = 0;
    for (uint32_t i = 0; i < (uint32_t) jet_p4->GetEntriesFast(); i++) {
      float pt = ((TLorentzVector*) (*jet_p4)[i])->Pt();

      HT += pt;
      if (pt > 30)
        HT30 += pt;
    }

    float HTFull = HT + ptLepton + MET;

    float mtt_four_leading_jets = 0.;
    if (! mSkim) {
      hNTrueInt_nosel->Fill(n_trueInteractions, eventWeight);
      hNVtx_nosel->Fill(n_vertices, eventWeight);
      hIsSel->Fill(isSel, eventWeight);

      hMET_nosel->Fill(MET, eventWeight);

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

      hHT_reco_nosel->Fill(HT, eventWeight);
      hHT30_reco_nosel->Fill(HT30, eventWeight);
      hHTFull_reco_nosel->Fill(HTFull, eventWeight);

      if (n_muons > 0) {
        hMuRelIso_nosel->Fill(muon_relIso[0], eventWeight);

        TLorentzVector* p4 = (TLorentzVector*) (*muon_p4)[0];
        hLeptonEta_nosel->Fill(p4->Eta(), eventWeight);
        hLeptonPt_nosel->Fill(p4->Pt(), eventWeight);
      }

      if (mIsSemiMu)
      {
        if (nGoodMuons <= 0)
          continue;
      }
      else
      {
        if (nGoodElectrons <= 0)
          continue;

        // The TOP reference selection exclude electron with
        // SuperCluster eta between 1.4442 and 1.5660
        // The TopTrigger efficiency does the same thing, but
        // using electron eta instead of SuperCluster eta.
        // Redo a cut here on electron eta
        // FIXME?
        if (fabs(etaLepton) >= 1.442 && fabs(etaLepton) < 1.5660)
          continue;
      }
    }

    if (mUseHybrid) {
      mUseChi2 = false;
      mUseKF = false;    
    }

    if (mUseMVA) {
      if (!mSkim && (isSel != 1 || numComb <= 0))
        continue;
      mtt_AfterReco = mtt_AfterMVA;
    } else if(mUseKF) {
      if (! kf_converged)
        continue;
      Init(true);
    } else if (mUseChi2) {
      if (numComb_chi2 == 0)
        continue;
      Init(false);
    } else if (mUseHybrid) {
      if (numComb_chi2 == 0) {
        if (kf_converged) {
          mUseKF = true;
          mUseChi2 = false;
        } else
          continue;
      } else {
        if (! kf_converged) {
          mUseChi2 = true;
          mUseKF = false;
        } else {
          if (kf_proba > 0.9) {
            mUseKF = true;
            mUseChi2 = false;          
          } else {
            mUseChi2 = true;
            mUseKF = false;          
          }
        }
      } 

      if (! mUseChi2 && ! mUseKF) continue;
      if (!mSkim && isSel != 1)
        continue;
      Init(mUseKF);
    }

    hBestSolChi2->Fill(bestSolChi2, eventWeight);
    hNGoodJets_chi2sel->Fill(nJets, eventWeight);
    hNBtaggedJets_chi2sel->Fill(nBtaggedJets_CSVM, eventWeight);

    hLeptonPt_chi2sel->Fill(ptLepton, eventWeight);
    hLeptonEta_chi2sel->Fill(etaLepton, eventWeight);
    if (mIsSemiMu)
      hMuRelIso_chi2sel->Fill(muRelIso[0], eventWeight);
    else
      hElRelIso_chi2sel->Fill(elRelIso[0], eventWeight);

    hMET_chi2sel->Fill(MET, eventWeight);

    hLeptTopPt_chi2sel->Fill(lepTopPt_AfterReco, eventWeight);
    hLeptTopEta_chi2sel->Fill(lepTopEta_AfterReco, eventWeight);

    hHadrTopPt_chi2sel->Fill(hadTopPt_AfterReco, eventWeight);
    hHadrTopEta_chi2sel->Fill(hadTopEta_AfterReco, eventWeight);

    hFirstJetPt_chi2sel->Fill(p_1stjetpt, eventWeight);
    hSecondJetPt_chi2sel->Fill(p_2ndjetpt, eventWeight);
    hThirdJetPt_chi2sel->Fill(p_3rdjetpt, eventWeight);
    hFourthJetPt_chi2sel->Fill(p_4thjetpt, eventWeight);

    hFirstJetEta_chi2sel->Fill(jetEta[0], eventWeight);
    hSecondJetEta_chi2sel->Fill(jetEta[1], eventWeight);
    hThirdJetEta_chi2sel->Fill(jetEta[2], eventWeight);
    hFourthJetEta_chi2sel->Fill(jetEta[3], eventWeight);

    hNVtx_chi2sel->Fill(n_vertices, eventWeight);

    hBoostTT_chi2sel->Fill(beta_tt_AfterReco, eventWeight);
    hPtTT_chi2sel->Fill(pt_tt_AfterReco, eventWeight);
    hEtaTT_chi2sel->Fill(eta_tt_AfterReco, eventWeight);

    LorentzVector selectedLeptonP4_LV_AfterReco;
    if (mUseKF) {
      selectedLeptonP4_LV_AfterReco = *selectedLeptonP4_AfterKF;
    } else {
      selectedLeptonP4_LV_AfterReco = LorentzVector(getP4(selectedLeptonP4_AfterReco, 0)->Pt(), getP4(selectedLeptonP4_AfterReco, 0)->Eta(), getP4(selectedLeptonP4_AfterReco, 0)->Phi(), getP4(selectedLeptonP4_AfterReco, 0)->E());
    }

    if (mIsMC && MC_channel != 0) {
      mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel.fill(MC_mtt, mtt_AfterReco, eventWeight);
      mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel.fill(MC_mtt, (mtt_AfterReco - MC_mtt) / MC_mtt, eventWeight);
      mtt_gen_vs_mtt_reco_resolution_rmsVSmttreco_sortingAlgoSel.fill(mtt_AfterReco, (mtt_AfterReco - MC_mtt) / MC_mtt, eventWeight);

      if (mUseKF) {
        if (kf_proba < 0.2) {
          mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel_lowKfProb.fill(MC_mtt, mtt_AfterReco, eventWeight);
          mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel_lowKfProb.fill(MC_mtt, (mtt_AfterReco - MC_mtt) / MC_mtt, eventWeight);

        } else {
          mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel_highKfProb.fill(MC_mtt, mtt_AfterReco, eventWeight);
          mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel_highKfProb.fill(MC_mtt, (mtt_AfterReco - MC_mtt) / MC_mtt, eventWeight);
        }
      }

      if (true) {

        // Reconstruct mtt using the four leading jet

        LorentzVector tt_system = selectedLeptonP4_LV_AfterReco + *selectedNeutrinoP4_AfterReco;
        int njets = 0;
        for (int i = 0; i < jet_p4->GetEntriesFast(); i++) {

          if (njets > 3)
            break;

          TLorentzVector* p = (TLorentzVector*) (*jet_p4)[i];

          if (p->Pt() < 30 || fabs(p->Eta()) > 2.4)
            continue;

          LorentzVector p_LV(p->Pt(), p->Eta(), p->Phi(), p->E());

          tt_system += p_LV;
          njets++;
        }

        mtt_four_leading_jets = tt_system.M();

        mtt_gen_vs_mtt_reco_linearity_fourJets_sortingAlgoSel.fill(MC_mtt, mtt_four_leading_jets, eventWeight);
        mtt_gen_vs_mtt_reco_resolution_fourJets_sortingAlgoSel.fill(MC_mtt, (mtt_four_leading_jets - MC_mtt) / MC_mtt, eventWeight);

        h_mtt_resolution_four_jets->Fill(mtt_four_leading_jets - MC_mtt, eventWeight);
      }
    }

    bool btagSel = false;
    if (mBTag < 0)
      btagSel = true;
    else if (mBTag == 0)
      btagSel = nBtaggedJets_CSVM == 0;
    else if (mBTag == 1)
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

    if (!mSkim && ptLepton <= ptLeptonCut)
      continue;

    if (btagSel && p_1stjetpt > firstJetCut && p_2ndjetpt > secondJetCut && p_3rdjetpt > thirdJetCut) {
      if (mIsMC) {
        if (isRun2012AB) {
          m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunA, lumi_run2012_A);
          m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunB, lumi_run2012_B);
          m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunC, 0);
          m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunD, 0);
        } else {
          m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunA, 0);
          m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunB, 0);
          m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunC, lumi_run2012_C);
          m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunD, lumi_run2012_D);
        }

        // Compute trigger weight
        std::vector<double> triggerWeights = m_trigger_efficiency_provider->get_weight(ptLepton, etaLepton, p_4thjetpt, jetEta[3], n_vertices, nJets, mIsSemiMu, triggerJESSyst);
        double triggerWeight = triggerWeights[0];
        if (mTriggerSyst == "up")
          triggerWeight = triggerWeight + triggerWeights[1];
        else if (mTriggerSyst == "down")
          triggerWeight = triggerWeight - triggerWeights[1];

        //eventWeight *= triggerWeight;

        hWeight_fullsel->Fill(eventWeight);
        hLeptonWeight_fullsel->Fill(m_lepton_weight);
        hBTagWeight_fullsel->Fill(m_btag_weight);
        hGeneratorWeight_fullsel->Fill(generator_weight);
        hPUWeight_fullsel->Fill(puWeight);
        //hTopPtWeight_fullsel->Fill(topPtWeight);
        hTriggerWeight_fullsel->Fill(triggerWeight);

        mtt_gen_vs_mtt_reco_linearity_fullSel.fill(MC_mtt, mtt_AfterReco, eventWeight);
        mtt_gen_vs_mtt_reco_resolution_fullSel.fill(MC_mtt, (mtt_AfterReco - MC_mtt) / MC_mtt, eventWeight);
        mtt_gen_vs_mtt_reco_linearity_fourJets_fullSel.fill(MC_mtt, mtt_four_leading_jets, eventWeight);
        mtt_gen_vs_mtt_reco_resolution_fourJets_fullSel.fill(MC_mtt, (mtt_four_leading_jets - MC_mtt) / MC_mtt, eventWeight);
      }

      hBestSolChi2_fullsel->Fill(bestSolChi2, eventWeight);
      hBestSolChi2Exp_fullsel->Fill(exp(-bestSolChi2), eventWeight);
      hBestSolChi2Proba_fullsel->Fill(TMath::Prob(bestSolChi2, 3), eventWeight); // number of dof = n-1 with n = 4 (terms in chi2)
      hBoostTT->Fill(beta_tt_AfterReco, eventWeight);
      hPtTT->Fill(pt_tt_AfterReco, eventWeight);
      hEtaTT->Fill(eta_tt_AfterReco, eventWeight);

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

      if (met_p4 && met_p4->GetEntriesFast() > 0) {
        TLorentzVector* p4 = (TLorentzVector*) (*met_p4)[0];
        hMETPhi->Fill(p4->Phi(), eventWeight);
        hMETx->Fill(p4->Px(), eventWeight);
        hMETy->Fill(p4->Py(), eventWeight);
      }

      hLeptTopPt->Fill(lepTopPt_AfterReco, eventWeight);
      hLeptTopEta->Fill(lepTopEta_AfterReco, eventWeight);

      hHadrTopPt->Fill(hadTopPt_AfterReco, eventWeight);
      hHadrTopEta->Fill(hadTopEta_AfterReco, eventWeight);

      hNGoodJets->Fill(nJets, eventWeight);
      hNBtaggedJets->Fill(nBtaggedJets_CSVM, eventWeight);

      hNVtx_noweight->Fill(n_vertices);
      hNVtx->Fill(n_vertices, eventWeight);
      hNTrueInt->Fill(n_trueInteractions, eventWeight);

      if (mIsSemiMu)
        hMuRelIso->Fill(muRelIso[0], eventWeight);
      else
        hElRelIso->Fill(elRelIso[0], eventWeight);

      hmtlep->Fill(mLepTop_AfterReco, eventWeight);
      hmthad->Fill(mHadTop_AfterReco, eventWeight);

      if (mUseHybrid) {
        if (mUseChi2) {
          hmtt_AfterChi2->Fill(mtt_AfterReco, eventWeight);
          if (mIsMC && MC_channel != 0) {
            h_mtt_resolution_AfterChi2->Fill(mtt_AfterReco - MC_mtt, eventWeight);
          }
        }

        else if (mUseKF) {
          hmtt_AfterKF->Fill(mtt_AfterReco, eventWeight);
          if (mIsMC && MC_channel != 0) {
             h_mtt_resolution_AfterKF->Fill(mtt_AfterReco - MC_mtt, eventWeight);
          }
        }
      }

      hmttSelected_btag_sel->Fill(mtt_AfterReco, eventWeight);
      hmttSelected_btag_sel_no_gen_weight->Fill(mtt_AfterReco, eventWeightNoGenWeight);

      if (eventWeight > 0)
        hmttSelected_btag_sel_positive->Fill(mtt_AfterReco, eventWeight);
      else
        hmttSelected_btag_sel_negative->Fill(mtt_AfterReco, eventWeight);

      pMttResolution_btag_sel->Fill(MC_mtt , mtt_AfterReco, eventWeight);

      hDeltaPhiTops_reco_fullsel->Fill(fabs(ROOT::Math::VectorUtil::DeltaPhi(*lepTopP4_AfterReco, *hadTopP4_AfterReco)), eventWeight);
      hDeltaEtaTops_reco_fullsel->Fill(lepTopP4_AfterReco->Eta() - hadTopP4_AfterReco->Eta(), eventWeight);
      hDeltaRTops_reco_fullsel->Fill(ROOT::Math::VectorUtil::DeltaR(*lepTopP4_AfterReco, *hadTopP4_AfterReco), eventWeight);
      hDeltaThetaTops_reco_fullsel->Fill(fabs(lepTopP4_AfterReco->Theta() - hadTopP4_AfterReco->Theta()), eventWeight);

      LorentzVector selectedHadronicWP4 = *selectedFirstJetP4_AfterReco + *selectedSecondJetP4_AfterReco + *selectedHadronicBP4_AfterReco;
      LorentzVector selectedLeptonicWP4 = *selectedLeptonicBP4_AfterReco + *selectedNeutrinoP4_AfterReco + selectedLeptonP4_LV_AfterReco;

      hLeptonicWMt->Fill(selectedLeptonicWP4.Mt(), eventWeight);
      hHadronicWMt->Fill(selectedHadronicWP4.Mt(), eventWeight);

      hLeptTopPz->Fill(lepTopP4_AfterReco->Pz(), eventWeight);
      hHadrTopPz->Fill(hadTopP4_AfterReco->Pz(), eventWeight);

      hLeptTopPhi->Fill(lepTopP4_AfterReco->Phi(), eventWeight);
      hHadrTopPhi->Fill(hadTopP4_AfterReco->Phi(), eventWeight);

      hTopsDeltaY->Fill(hadTopP4_AfterReco->Rapidity() - lepTopP4_AfterReco->Rapidity(), eventWeight);

      hLeptTopE->Fill(lepTopP4_AfterReco->E(), eventWeight);
      hLeptTopMt->Fill(lepTopP4_AfterReco->Mt(), eventWeight);
      hHadrTopE->Fill(hadTopP4_AfterReco->E(), eventWeight);
      hHadrTopMt->Fill(hadTopP4_AfterReco->Mt(), eventWeight);

      selectedFirstJetPt_AfterReco = selectedFirstJetP4_AfterReco->Pt();
      selectedSecondJetPt_AfterReco = selectedSecondJetP4_AfterReco->Pt();
      selectedHadronicBPt_AfterReco = selectedHadronicBP4_AfterReco->Pt();
      selectedLeptonicBPt_AfterReco = selectedLeptonicBP4_AfterReco->Pt();

      hSelectedFirstJetPt->Fill(selectedFirstJetPt_AfterReco, eventWeight);
      hSelectedSecondJetPt->Fill(selectedSecondJetPt_AfterReco, eventWeight);
      hSelectedHadronicBPt->Fill(selectedHadronicBPt_AfterReco, eventWeight);
      hSelectedLeptonicBPt->Fill(selectedLeptonicBPt_AfterReco, eventWeight);
      hSelectedNeutrinoPt->Fill(selectedNeutrinoP4_AfterReco->Pt(), eventWeight);
      hSelectedNeutrinoPz->Fill(selectedNeutrinoP4_AfterReco->Pz(), eventWeight);

      hSelectedFirstJetEta->Fill(selectedFirstJetP4_AfterReco->Eta(), eventWeight);
      hSelectedSecondJetEta->Fill(selectedSecondJetP4_AfterReco->Eta(), eventWeight);
      hSelectedHadronicBEta->Fill(selectedHadronicBP4_AfterReco->Eta(), eventWeight);
      hSelectedLeptonicBEta->Fill(selectedLeptonicBP4_AfterReco->Eta(), eventWeight);
      hSelectedNeutrinoEta->Fill(selectedNeutrinoP4_AfterReco->Eta(), eventWeight);

      hSelectedFirstJetPhi->Fill(selectedFirstJetP4_AfterReco->Phi(), eventWeight);
      hSelectedSecondJetPhi->Fill(selectedSecondJetP4_AfterReco->Phi(), eventWeight);
      hSelectedHadronicBPhi->Fill(selectedHadronicBP4_AfterReco->Phi(), eventWeight);
      hSelectedLeptonicBPhi->Fill(selectedLeptonicBP4_AfterReco->Phi(), eventWeight);
      hSelectedNeutrinoPhi->Fill(selectedNeutrinoP4_AfterReco->Phi(), eventWeight);

      hSelectedFirstJetE->Fill(selectedFirstJetP4_AfterReco->E(), eventWeight);
      hSelectedSecondJetE->Fill(selectedSecondJetP4_AfterReco->E(), eventWeight);
      hSelectedHadronicBE->Fill(selectedHadronicBP4_AfterReco->E(), eventWeight);
      hSelectedLeptonicBE->Fill(selectedLeptonicBP4_AfterReco->E(), eventWeight);
      hSelectedNeutrinoE->Fill(selectedNeutrinoP4_AfterReco->E(), eventWeight);
      if (mUseKF) {
        hSelectedLeptonPt->Fill(selectedLeptonP4_AfterKF->Pt(), eventWeight);
        hSelectedLeptonEta->Fill(selectedLeptonP4_AfterKF->Eta(), eventWeight);
        hSelectedLeptonPhi->Fill(selectedLeptonP4_AfterKF->Phi(), eventWeight);
        hSelectedLeptonE->Fill(selectedLeptonP4_AfterKF->E(), eventWeight);
      } else {
        hSelectedLeptonPt->Fill(getP4(selectedLeptonP4_AfterReco, 0)->Pt(), eventWeight);
        hSelectedLeptonEta->Fill(getP4(selectedLeptonP4_AfterReco, 0)->Eta(), eventWeight);
        hSelectedLeptonPhi->Fill(getP4(selectedLeptonP4_AfterReco, 0)->Phi(), eventWeight);
        hSelectedLeptonE->Fill(getP4(selectedLeptonP4_AfterReco, 0)->E(), eventWeight);
      }

      sumPt4JetsSel = selectedFirstJetPt_AfterReco + selectedSecondJetPt_AfterReco + selectedHadronicBPt_AfterReco + selectedLeptonicBPt_AfterReco;
      sumPtJetsInEvent = 0.;
      for (int i = 0; i < nJets; i++) {
        sumPtJetsInEvent += jetPt[i];
      }
      HTFracValue = sumPt4JetsSel/sumPtJetsInEvent;
      hHTFrac->Fill(HTFracValue, eventWeight);

      if (mtt_AfterReco > 500)
        hmttSelected_btag_sel_mass_cut->Fill(mtt_AfterReco, eventWeight);

      h_mtt_resolution->Fill(mtt_AfterReco - MC_mtt, eventWeight);
      hkf_chisquare->Fill(kf_chisquare, eventWeight);
      hkf_proba->Fill(kf_proba, eventWeight);
      hkf_proba_zoom->Fill(kf_proba, eventWeight);
      hAplanarity->Fill(p_aplanarity, eventWeight);
      hCircularity->Fill(p_circularity, eventWeight);
      hSphericity->Fill(p_sphericity, eventWeight);
      hIsotropy->Fill(p_isotropy, eventWeight);
      hD->Fill(p_D, eventWeight);
      hC->Fill(p_C, eventWeight);

      float discriminant = bkgVsTTBDTReader.evaluate(jentry);
      hBDTDiscriminant->Fill(discriminant, eventWeightNoGenWeight);
    }
  }

  output->cd();

  THStack* hsmtt_chi2_kf_fraction = new THStack("mtt_chi2_kf_fraction", "mtt_chi2_kf_fraction");
  hmtt_AfterChi2->SetLineColor(kRed-7);
  hmtt_AfterChi2->SetFillColor(kRed-7);
  hmtt_AfterKF->SetLineColor(kGreen-8);
  hmtt_AfterKF->SetFillColor(kGreen-8);
  hsmtt_chi2_kf_fraction->Add(hmtt_AfterChi2);
  hsmtt_chi2_kf_fraction->Add(hmtt_AfterKF);
  hsmtt_chi2_kf_fraction->Write();

  THStack* hs_mtt_resolution_chi2_kf_fraction = new THStack("mtt_resolution_chi2_kf_fraction", "mtt_resolution_chi2_kf_fraction");
  h_mtt_resolution_AfterChi2->SetLineColor(kRed-7);
  h_mtt_resolution_AfterChi2->SetFillColor(kRed-7);
  h_mtt_resolution_AfterKF->SetLineColor(kGreen-8);
  h_mtt_resolution_AfterKF->SetFillColor(kGreen-8);
  hs_mtt_resolution_chi2_kf_fraction->Add(h_mtt_resolution_AfterChi2);
  hs_mtt_resolution_chi2_kf_fraction->Add(h_mtt_resolution_AfterKF);
  hs_mtt_resolution_chi2_kf_fraction->Write();

  if (mIsMC && hTopPtWeight_fullsel->GetEntries() > 0) {
    std::cout << "Top pt reweighting mean weight: " << hTopPtWeight_fullsel->GetMean() << std::endl;
    std::cout << "Use this value to compute the new number of generated events for this sample" << std::endl;

    // Save some useful informations
    boost::filesystem::path p(mOutputFile);
    p.replace_extension("info");

    std::ofstream f(p.string());
    //f << hTopPtWeight_fullsel->GetMean() << std::endl;
    f << 1 << std::endl;
    f.close();
  }

  output->Write();

  if (mIsMC) {
    mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel.write(output);
    mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel.write(output);
    mtt_gen_vs_mtt_reco_resolution_rmsVSmttreco_sortingAlgoSel.write(output);
    mtt_gen_vs_mtt_reco_linearity_fourJets_sortingAlgoSel.write(output);
    mtt_gen_vs_mtt_reco_resolution_fourJets_sortingAlgoSel.write(output);

    mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel_lowKfProb.write(output);
    mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel_lowKfProb.write(output);
    mtt_gen_vs_mtt_reco_resolution_sortingAlgoSel_highKfProb.write(output);
    mtt_gen_vs_mtt_reco_linearity_sortingAlgoSel_highKfProb.write(output);

    mtt_gen_vs_mtt_reco_linearity_fullSel.write(output);
    mtt_gen_vs_mtt_reco_resolution_fullSel.write(output);
    mtt_gen_vs_mtt_reco_linearity_fourJets_fullSel.write(output);
    mtt_gen_vs_mtt_reco_resolution_fourJets_fullSel.write(output);
  }

  output->Close();
  delete output;
}

void loadChain(const std::vector<std::string>& inputFiles, const std::string& treeName, TChain*& output) {

  output = new TChain(treeName.c_str());

  for (const std::string& file: inputFiles) {
    output->Add(file.c_str());
  }

  output->SetCacheSize(30*1024*1024);
}

Extractor2Histos::Extractor2Histos(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isSemiMu, bool isMC, int btag, bool skim, bool mva, bool chi2, bool kf, bool hybrid, const std::string& triggerSyst, const std::string& jecSyst, const std::string& puSyst, const std::string& pdfSyst, const std::string& leptonSyst, const std::string& btagSyst, const std::string& bdtWeights) : fMTT(0), fVertices(0), fEvent(0)
{
  mIsSemiMu = isSemiMu;
  mIsMC = isMC;
  mOutputFile = outputFile;
  mBTag = btag;
  mSkim = skim;
  mUseMVA = mva;
  mUseChi2 = chi2;
  mUseKF = kf;
  mUseHybrid = hybrid;
  mTriggerSyst = triggerSyst;
  mJECSyst = jecSyst;
  mPUSyst = puSyst;
  mPDFSyst = pdfSyst;
  mLeptonSyst = leptonSyst;
  mBTagSyst = btagSyst;
  m_bdtWeights = bdtWeights;
  m_inputFiles = inputFiles;

  fLooseMuons = nullptr;
  fJet = nullptr;
  fMET = nullptr;

  // Get trees
  loadChain(inputFiles, "Mtt", fMTT);
  loadChain(inputFiles, "Vertices", fVertices);
  loadChain(inputFiles, "event", fEvent);

  if (! mSkim) {
    loadChain(inputFiles, "muon_loose_PF", fLooseMuons);
  }

  loadChain(inputFiles, "jet_PF", fJet);
  //loadChain(inputFiles, "MET_PF", fMET);

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

  if (fMET)
    fMET->GetEntry(entry);

  //if (fElectrons)
  //fElectrons->GetEntry(entry);

  return 1;
}

void ActivateBranch(TBranch* branch) {
  branch->SetStatus(1);
  TObjArray* objArray = branch->GetListOfBranches();
  for (int i = 0; i < objArray->GetEntries(); i++) {
    TBranch* b = static_cast<TBranch*> (objArray->At(i));
    if (b) {
      ActivateBranch(b);
    }
  }
}

void Extractor2Histos::SetBranchAddress(TTree* t, const char* branchName, void* ptr) {
  TBranch* branch = t->FindBranch(branchName);
  if (!branch)
    return;

  ActivateBranch(branch);
  t->SetBranchAddress(branchName, ptr, NULL);
}

void Extractor2Histos::Init(bool useKF) 
{
  if (useKF) {
    mtt_AfterReco = mtt_AfterKF;
    mLepTop_AfterReco = mLepTop_AfterKF;
    mHadTop_AfterReco = mHadTop_AfterKF;
    pt_tt_AfterReco = pt_tt_AfterKF;
    eta_tt_AfterReco = eta_tt_AfterKF;
    beta_tt_AfterReco = beta_tt_AfterKF;

    lepTopPt_AfterReco = lepTopPt_AfterKF;
    lepTopEta_AfterReco = lepTopEta_AfterKF;
    hadTopPt_AfterReco = hadTopPt_AfterKF;
    hadTopEta_AfterReco = hadTopEta_AfterKF;

    lepTopP4_AfterReco = lepTopP4_AfterKF;
    hadTopP4_AfterReco = hadTopP4_AfterKF;
    selectedFirstJetP4_AfterReco = selectedFirstJetP4_AfterKF;
    selectedSecondJetP4_AfterReco = selectedSecondJetP4_AfterKF;
    selectedNeutrinoP4_AfterReco = selectedNeutrinoP4_AfterKF;
    selectedHadronicBP4_AfterReco = selectedHadronicBP4_AfterKF;
    selectedLeptonicBP4_AfterReco = selectedLeptonicBP4_AfterKF;
  } else {
    mtt_AfterReco = mtt_AfterChi2;
    mLepTop_AfterReco = mLepTop_AfterChi2;
    mHadTop_AfterReco = mHadTop_AfterChi2;
    pt_tt_AfterReco = pt_tt_AfterChi2;
    eta_tt_AfterReco = eta_tt_AfterChi2;
    beta_tt_AfterReco = beta_tt_AfterChi2;

    lepTopPt_AfterReco = lepTopPt_AfterChi2;
    lepTopEta_AfterReco = lepTopEta_AfterChi2;
    hadTopPt_AfterReco = hadTopPt_AfterChi2;
    hadTopEta_AfterReco = hadTopEta_AfterChi2;

    lepTopP4_AfterReco = lepTopP4_AfterChi2;
    hadTopP4_AfterReco = hadTopP4_AfterChi2;
    selectedFirstJetP4_AfterReco = selectedFirstJetP4_AfterChi2;
    selectedSecondJetP4_AfterReco = selectedSecondJetP4_AfterChi2;
    selectedNeutrinoP4_AfterReco = selectedNeutrinoP4_AfterChi2;
    selectedHadronicBP4_AfterReco = selectedHadronicBP4_AfterChi2;
    selectedLeptonicBP4_AfterReco = selectedLeptonicBP4_AfterChi2;
  }
}

void Extractor2Histos::Init()
{
  fCurrent = -1;

  fMTT->SetBranchStatus("*", 0);

  lepTopP4_AfterReco = NULL;
  hadTopP4_AfterReco = NULL;
  selectedFirstJetP4_AfterReco = NULL;
  selectedSecondJetP4_AfterReco = NULL;
  selectedHadronicBP4_AfterReco = NULL;
  selectedLeptonicBP4_AfterReco = NULL;
  selectedNeutrinoP4_AfterReco = NULL;
  selectedLeptonP4_AfterReco = NULL;
  selectedLeptonP4_AfterKF = NULL;
  lepTopP4_AfterChi2 = NULL;
  hadTopP4_AfterChi2 = NULL;
  selectedFirstJetP4_AfterChi2 = NULL;
  selectedSecondJetP4_AfterChi2 = NULL;
  selectedHadronicBP4_AfterChi2 = NULL;
  selectedLeptonicBP4_AfterChi2 = NULL;
  selectedNeutrinoP4_AfterChi2 = NULL;
  lepTopP4_AfterKF = NULL;
  hadTopP4_AfterKF = NULL;
  selectedFirstJetP4_AfterKF = NULL;
  selectedSecondJetP4_AfterKF = NULL;
  selectedHadronicBP4_AfterKF = NULL;
  selectedLeptonicBP4_AfterKF = NULL;
  selectedNeutrinoP4_AfterKF = NULL;
  selectedLeptonP4_AfterKF = NULL;

  if (! mUseMVA) { // if we don't use MVA algo, variables to choose between Chi2 and KF
    SetBranchAddress(fMTT, "kf_converged", &kf_converged);
    SetBranchAddress(fMTT, "mtt_AfterKF", &mtt_AfterKF);
    SetBranchAddress(fMTT, "kf_proba", &kf_proba);
    SetBranchAddress(fMTT, "kf_chisquare", &kf_chisquare);
    SetBranchAddress(fMTT, "mtt_AfterChi2", &mtt_AfterChi2);
    SetBranchAddress(fMTT, "numComb_chi2", &numComb_chi2);

    SetBranchAddress(fMTT, "mLepTop_AfterKF", &mLepTop_AfterKF);
    SetBranchAddress(fMTT, "mHadTop_AfterKF", &mHadTop_AfterKF);
    SetBranchAddress(fMTT, "pt_tt_AfterKF", &pt_tt_AfterKF);
    SetBranchAddress(fMTT, "eta_tt_AfterKF", &eta_tt_AfterKF);
    SetBranchAddress(fMTT, "beta_tt_AfterKF", &beta_tt_AfterKF);

    SetBranchAddress(fMTT, "lepTopPt_AfterKF", &lepTopPt_AfterKF);
    SetBranchAddress(fMTT, "lepTopEta_AfterKF", &lepTopEta_AfterKF);
    SetBranchAddress(fMTT, "hadTopPt_AfterKF", &hadTopPt_AfterKF);
    SetBranchAddress(fMTT, "hadTopEta_AfterKF", &hadTopEta_AfterKF);

    SetBranchAddress(fMTT, "lepTopP4_AfterKF", &lepTopP4_AfterKF);
    SetBranchAddress(fMTT, "hadTopP4_AfterKF", &hadTopP4_AfterKF);
    SetBranchAddress(fMTT, "selectedFirstJetP4_AfterKF", &selectedFirstJetP4_AfterKF);
    SetBranchAddress(fMTT, "selectedSecondJetP4_AfterKF", &selectedSecondJetP4_AfterKF);
    SetBranchAddress(fMTT, "selectedNeutrinoP4_AfterKF", &selectedNeutrinoP4_AfterKF);
    SetBranchAddress(fMTT, "selectedHadronicBP4_AfterKF", &selectedHadronicBP4_AfterKF);
    SetBranchAddress(fMTT, "selectedLeptonicBP4_AfterKF", &selectedLeptonicBP4_AfterKF);

    SetBranchAddress(fMTT, "selectedLeptonP4_AfterKF", &selectedLeptonP4_AfterKF);
    SetBranchAddress(fMTT, "selectedLeptonP4", &selectedLeptonP4_AfterReco);

    SetBranchAddress(fMTT, "mLepTop_AfterChi2", &mLepTop_AfterChi2);
    SetBranchAddress(fMTT, "mHadTop_AfterChi2", &mHadTop_AfterChi2);
    SetBranchAddress(fMTT, "mtt_AfterChi2", &mtt_AfterChi2);
    SetBranchAddress(fMTT, "pt_tt_AfterChi2", &pt_tt_AfterChi2);
    SetBranchAddress(fMTT, "eta_tt_AfterChi2", &eta_tt_AfterChi2);
    SetBranchAddress(fMTT, "beta_tt_AfterChi2", &beta_tt_AfterChi2);

    SetBranchAddress(fMTT, "lepTopPt_AfterChi2", &lepTopPt_AfterChi2);
    SetBranchAddress(fMTT, "lepTopEta_AfterChi2", &lepTopEta_AfterChi2);
    SetBranchAddress(fMTT, "hadTopPt_AfterChi2", &hadTopPt_AfterChi2);
    SetBranchAddress(fMTT, "hadTopEta_AfterChi2", &hadTopEta_AfterChi2);

    SetBranchAddress(fMTT, "lepTopP4_AfterChi2", &lepTopP4_AfterChi2);
    SetBranchAddress(fMTT, "hadTopP4_AfterChi2", &hadTopP4_AfterChi2);
    SetBranchAddress(fMTT, "selectedFirstJetP4_AfterChi2", &selectedFirstJetP4_AfterChi2);
    SetBranchAddress(fMTT, "selectedSecondJetP4_AfterChi2", &selectedSecondJetP4_AfterChi2);
    SetBranchAddress(fMTT, "selectedNeutrinoP4_AfterChi2", &selectedNeutrinoP4_AfterChi2);
    SetBranchAddress(fMTT, "selectedHadronicBP4_AfterChi2", &selectedHadronicBP4_AfterChi2);
    SetBranchAddress(fMTT, "selectedLeptonicBP4_AfterChi2", &selectedLeptonicBP4_AfterChi2);
  }

  SetBranchAddress(fMTT, "MC_channel", &MC_channel);
  SetBranchAddress(fMTT, "MC_mtt", &MC_mtt);

  SetBranchAddress(fMTT, "MC_beta_tt", &MC_boost_tt);

  SetBranchAddress(fMTT, "MC_pt_tt", &MC_pt_tt);
  SetBranchAddress(fMTT, "MC_eta_tt", &MC_eta_tt);

  SetBranchAddress(fMTT, "nGoodMuons", &nGoodMuons);
  //SetBranchAddress(fMTT, "nLooseGoodMuons", &nLooseGoodMuons);
  if (mIsSemiMu) {
    SetBranchAddress(fMTT, "muonPt", muonPt);
    SetBranchAddress(fMTT, "muonEta", muonEta);
    SetBranchAddress(fMTT, "muRelIso", muRelIso);
  } else {
    SetBranchAddress(fMTT, "nGoodElectrons", &nGoodElectrons);
    SetBranchAddress(fMTT, "electronPt", &electronPt);
    SetBranchAddress(fMTT, "electronEta", &electronEta);
    SetBranchAddress(fMTT, "elRelIso", &elRelIso);
  }

  SetBranchAddress(fMTT, "1stjetpt", &p_1stjetpt);
  SetBranchAddress(fMTT, "2ndjetpt", &p_2ndjetpt);
  SetBranchAddress(fMTT, "3rdjetpt", &p_3rdjetpt);
  SetBranchAddress(fMTT, "4thjetpt", &p_4thjetpt);
  SetBranchAddress(fMTT, "nJets", &nJets);
  SetBranchAddress(fMTT, "jetEta", jetEta);
  SetBranchAddress(fMTT, "jetPt", jetPt);
  SetBranchAddress(fMTT, "nBtaggedJets_CSVM", &nBtaggedJets_CSVM);
  SetBranchAddress(fMTT, "MET", &MET);
  SetBranchAddress(fMTT, "isSel", &isSel);
  SetBranchAddress(fMTT, "bestSolChi2", &bestSolChi2);
  SetBranchAddress(fMTT, "aplanarity", &p_aplanarity);
  SetBranchAddress(fMTT, "circularity", &p_circularity);
  SetBranchAddress(fMTT, "sphericity", &p_sphericity);
  SetBranchAddress(fMTT, "isotropy", &p_isotropy);
  SetBranchAddress(fMTT, "D", &p_D);
  SetBranchAddress(fMTT, "C", &p_C);

  if (mUseMVA) {
    SetBranchAddress(fMTT, "numComb_MVA", &numComb);
    SetBranchAddress(fMTT, "mLepTop_AfterMVA", &mLepTop_AfterReco);
    SetBranchAddress(fMTT, "mHadTop_AfterMVA", &mHadTop_AfterReco);
    SetBranchAddress(fMTT, "mtt_AfterMVA", &mtt_AfterMVA);
    SetBranchAddress(fMTT, "pt_tt_AfterMVA", &pt_tt_AfterReco);
    SetBranchAddress(fMTT, "eta_tt_AfterMVA", &eta_tt_AfterReco);
    SetBranchAddress(fMTT, "beta_tt_AfterMVA", &beta_tt_AfterReco);

    SetBranchAddress(fMTT, "lepTopPt_AfterMVA", &lepTopPt_AfterReco);
    SetBranchAddress(fMTT, "lepTopEta_AfterMVA", &lepTopEta_AfterReco);
    SetBranchAddress(fMTT, "hadTopPt_AfterMVA", &hadTopPt_AfterReco);
    SetBranchAddress(fMTT, "hadTopEta_AfterMVA", &hadTopEta_AfterReco);
  } 
  SetBranchAddress(fMTT, "lepton_weight", &m_lepton_weight);
  SetBranchAddress(fMTT, "btag_weight", &m_btag_weight);

  if (mLeptonSyst == "up")
    SetBranchAddress(fMTT, "lepton_weight_error_high", &m_lepton_weight_error);

  if (mLeptonSyst == "down")
    SetBranchAddress(fMTT, "lepton_weight_error_low", &m_lepton_weight_error);

  if (mBTagSyst == "up")
    SetBranchAddress(fMTT, "btag_weight_error_high", &m_btag_weight_error);

  if (mBTagSyst == "down")
    SetBranchAddress(fMTT, "btag_weight_error_low", &m_btag_weight_error);

  if (fMTT->GetBranch("trigger_passed")) {
    SetBranchAddress(fMTT, "trigger_passed", &m_triggerPassed);
  } else {
    // Backward compatibilty
    m_triggerPassed = true;
  }

  if (mUseMVA) {
    SetBranchAddress(fMTT, "lepTopP4_AfterMVA", &lepTopP4_AfterReco);
    SetBranchAddress(fMTT, "hadTopP4_AfterMVA", &hadTopP4_AfterReco);
    SetBranchAddress(fMTT, "selectedFirstJetP4_AfterMVA", &selectedFirstJetP4_AfterReco);
    SetBranchAddress(fMTT, "selectedSecondJetP4_AfterMVA", &selectedSecondJetP4_AfterReco);
    SetBranchAddress(fMTT, "selectedLeptonP4", &selectedLeptonP4_AfterReco);
    SetBranchAddress(fMTT, "selectedNeutrinoP4_AfterMVA", &selectedNeutrinoP4_AfterReco);
    SetBranchAddress(fMTT, "selectedHadronicBP4_AfterMVA", &selectedHadronicBP4_AfterReco);
    SetBranchAddress(fMTT, "selectedLeptonicBP4_AfterMVA", &selectedLeptonicBP4_AfterReco);
  } 
  gen_top1_p4 = NULL;
  gen_top2_p4 = NULL;

  if (mIsMC) {
    SetBranchAddress(fMTT, "MC_Top1_p4", &gen_top1_p4);
    SetBranchAddress(fMTT, "MC_Top2_p4", &gen_top2_p4);
  }

  gen_lepton_p4 = NULL;
  gen_neutrino_p4 = NULL;

  if (mIsMC) {
    SetBranchAddress(fMTT, "MC_lepton_p4", &gen_lepton_p4);
    SetBranchAddress(fMTT, "MC_neutrino_p4", &gen_neutrino_p4);
  }

  gen_leptonic_B_p4 = NULL;
  gen_hadronic_B_p4 = NULL;

  if (mIsMC) {
    SetBranchAddress(fMTT, "MC_leptonic_B_p4", &gen_leptonic_B_p4);
    SetBranchAddress(fMTT, "MC_hadronic_B_p4", &gen_hadronic_B_p4);
  }

  gen_lightJet1_p4 = NULL;
  gen_lightJet2_p4 = NULL;

  if (mIsMC) {
    SetBranchAddress(fMTT, "MC_lightJet1_B_p4", &gen_lightJet1_p4);
    SetBranchAddress(fMTT, "MC_lightJet2_B_p4", &gen_lightJet2_p4);
  }

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

  if (! mSkim) {
    muon_p4 = NULL;
    fLooseMuons->SetMakeClass(1);
    fLooseMuons->SetBranchStatus("*", 0);
    fLooseMuons->SetBranchStatus("muon_4vector", 1);
    fLooseMuons->SetBranchStatus("n_muons", 1);
    fLooseMuons->SetBranchStatus("muon_deltaBetaCorrectedRelIsolation", 1);
    fLooseMuons->SetBranchAddress("muon_4vector", &muon_p4, NULL);
    fLooseMuons->SetBranchAddress("n_muons", &n_muons, NULL);
    fLooseMuons->SetBranchAddress("muon_deltaBetaCorrectedRelIsolation", &muon_relIso, NULL);
  }

  jet_p4 = NULL;
  fJet->SetMakeClass(1);
  fJet->SetBranchStatus("*", 0);
  fJet->SetBranchStatus("jet_4vector", 1);

  fJet->SetBranchAddress("jet_4vector", &jet_p4, NULL);

  met_p4 = NULL;
  if (fMET) {
    fMET->SetMakeClass(1);
    fMET->SetBranchStatus("*", 0);
    fMET->SetBranchStatus("met_4vector", 1);

    fMET->SetBranchAddress("met_4vector", &met_p4, NULL);
  }
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

    TCLAP::SwitchArg skimArg("", "skim", "Run over a skimmed file", cmd, false);

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

    TCLAP::ValueArg<std::string> bdtWeightsArg("", "bdt-weights", "Full path to the BDT weights (the XML file)", true, "", "string", cmd);

    TCLAP::ValueArg<std::string> pileupArg("", "pileup", "PU profile used for MC production", false, "S10", "string", cmd);

    TCLAP::ValueArg<std::string> jecSystArg("", "jec-syst", "Computing trigger weight for this JEC up / down", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> triggerSystArg("", "trigger-syst", "Computing trigger weight systematic", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> pileupSystArg("", "pileup-syst", "PU profile to use for pileup reweigthing", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> pdfSystArg("", "pdf-syst", "PDF systematic to compute", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> btagSystArg("", "btag-syst", "Compute btag weight systematic", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> leptonSystArg("", "lepton-syst", "Compute lepton weight systematic", false, "nominal", "string", cmd);


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

    std::string triggerSyst = triggerSystArg.getValue();
    std::transform(triggerSyst.begin(), triggerSyst.end(), triggerSyst.begin(), ::tolower);
    if (triggerSyst != "nominal" && triggerSyst != "up" && triggerSyst != "down") {
      std::cerr << "--trigger-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    std::string jecSyst = jecSystArg.getValue();
    std::transform(jecSyst.begin(), jecSyst.end(), jecSyst.begin(), ::tolower);
    if (jecSyst != "nominal" && jecSyst != "up" && jecSyst != "down") {
      std::cerr << "--jec-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

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

    std::string leptonSyst = leptonSystArg.getValue();
    std::transform(leptonSyst.begin(), leptonSyst.end(), leptonSyst.begin(), ::tolower);
    if (leptonSyst != "nominal" && leptonSyst != "up" && leptonSyst != "down") {
      std::cerr << "--lepton-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    std::string btagSyst = btagSystArg.getValue();
    std::transform(btagSyst.begin(), btagSyst.end(), btagSyst.begin(), ::tolower);
    if (btagSyst != "nominal" && btagSyst != "up" && btagSyst != "down") {
      std::cerr << "--btag-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    bool isData = dataArg.isSet();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    Extractor2Histos convertor(inputFiles, outputFileArg.getValue(), semimuArg.isSet(), !isData, btagArg.getValue(), skimArg.getValue(), mvaArg.getValue(), chi2Arg.getValue(), kfArg.getValue(), hybridArg.getValue(), triggerSyst, jecSyst, puSyst, pdfSyst, leptonSyst, btagSyst, bdtWeightsArg.getValue());
    convertor.Loop();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
