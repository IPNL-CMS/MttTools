#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <sys/wait.h>
#include <cmath>
#include <regex>
#include <cstdio>

#include <TString.h>
#include <TRandom3.h>
#include <TDatime.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLatex.h>
#include <TChain.h>
#include <TVectorD.h>

#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooExponential.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooSimultaneous.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooMinuit.h>
#include <RooWorkspace.h>
#include <RooSuperCategory.h>
#include <Roo1DTable.h>
#include <RooDataHist.h>
#include <RooChi2Var.h>
#include <RooExtendPdf.h>
#include <RooIntegralMorph.h>
#include <RooMomentMorph.h>
#include <RooHistPdf.h>

using namespace RooFit;

#include <list>
#include "tdrstyle.C"

#include "Utils.h"
#include "Functions.h"

#include <tclap/CmdLine.h>
#include <json/json.h>

#include <TString.h>

void drawPdfs(RooRealVar& observable, RooAbsPdf& interpolated_pdf, RooAbsPdf& lowPdf, RooAbsPdf& highPdf, const TString& filename, int btag, int mass) {

  const int padWidth = 900;
  const int padHeight = 900;

  const int canvasWidth = padWidth * 1;
  const int canvasHeight = padHeight * 1;

  TCanvas *canvas = new TCanvas("canvas", "mTT fit", canvasWidth, canvasHeight);

  RooPlot* plot = observable.frame();
  plot->SetTitle(""); //FIXME

  lowPdf.plotOn(plot, LineStyle(kDashed));
  highPdf.plotOn(plot, LineStyle(kDashed));
  interpolated_pdf.plotOn(plot, LineColor(kRed));

  std::string leptonName = filename.Contains("muon", TString::kIgnoreCase) ? "muons" : "electrons";
  float binningSize = (observable.getBinning().highBound() - observable.getBinning().lowBound()) / (float)observable.getBins();

  plot->SetXTitle(TString::Format("#font[132]{#font[12]{m_{TT}} (GeV/#font[12]{c}^{2}), %s}", leptonName.c_str()));
  plot->SetYTitle(TString::Format("#font[132]{Events/(%0.2f GeV/#font[12]{c}^{2})}", binningSize));
  plot->SetTitleOffset(1.42, "Y");
  plot->SetNdivisions(505);

  TLatex t;
  t.SetNDC();
  t.SetTextSize(0.04);

  gPad->SetTopMargin(0.05); 
  gPad->SetBottomMargin(0.12); 
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.050);

  plot->Draw();

  TString btagLabel = "";
  if (btag == 2)
    btagLabel = "#geq 2 b-tags";
  else
    btagLabel = TString::Format("%d b-tag", btag);

  std::string leptonShortcutName = filename.Contains("muon", TString::kIgnoreCase) ? "#mu" : "e";

  TString legendLabel = TString::Format("#font[42]{%s, #geq 4 jets, %s}", leptonShortcutName.c_str(), btagLabel.Data());

  t.DrawLatex(0.53, 0.88, "#font[42]{CMS preliminary}");
  t.DrawLatex(0.53, 0.84, "#font[42]{Z' interpolation}");
  t.DrawLatex(0.53, 0.80, TString::Format("#font[42]{m = %d GeV}", mass));
  t.DrawLatex(0.53, 0.76, legendLabel);


  canvas->Update();
  canvas->SaveAs(filename);

  std::cout << "Plot saved as '" << filename << "'" << std::endl;

  delete canvas;
  delete plot;
}

void doInterpolation(int mass, const std::string& jec, int btag) {

  std::string analysisName = getAnalysisName();
  std::string base_path = "analysis/" + getAnalysisUUID();

  TString prefix = TString::Format("%s-%s%d_%s_%d_btag", jec.c_str(), getAnalysisPrefix(), mass, analysisName.c_str(), btag);

  int lowMass = 0;
  int highMass = 0;

  if (mass > 750 && mass < 1000) {
    lowMass = 750;
    highMass = 1000;
  } else if (mass > 1000 && mass < 1250) {
    lowMass = 1000;
    highMass = 1250;
  } else if (mass > 1250 && mass < 1500) {
    lowMass = 1250;
    highMass = 1500;
  }/* else if (mass > 1500 && mass < 2000) {
    lowMass = 1500;
    highMass = 2000;
  }*/

  if (lowMass == 0 || highMass == 0) {
    std::cout << "Error: please use a mass between 750 and 1500 GeV." << std::endl;
    return;
  }

  // Open workspaces and retrieve PDFs
  TString lowMass_prefix = TString::Format("%s-%s%d_%s_%d_btag", jec.c_str(), getAnalysisPrefix(), lowMass, analysisName.c_str(), btag);
  TString lowMass_workspaceFile = base_path + "/frit/" + lowMass_prefix + "_workspace.root";

  TString highMass_prefix = TString::Format("%s-%s%d_%s_%d_btag", jec.c_str(), getAnalysisPrefix(), highMass, analysisName.c_str(), btag);
  TString highMass_workspaceFile = base_path + "/frit/" + highMass_prefix + "_workspace.root";

  TFile lowMass_file(lowMass_workspaceFile);
  RooAbsPdf* lowMass_pdf_e = static_cast<RooWorkspace*>(lowMass_file.Get("w"))->pdf("signal_electron");
  lowMass_pdf_e->SetName("lowMass_signal_electron");
  RooAbsPdf* lowMass_pdf_mu = static_cast<RooWorkspace*>(lowMass_file.Get("w"))->pdf("signal_muon");
  lowMass_pdf_mu->SetName("lowMass_signal_muon");
  lowMass_file.Close();

  TFile highMass_file(highMass_workspaceFile);
  RooAbsPdf* highMass_pdf_e = static_cast<RooWorkspace*>(highMass_file.Get("w"))->pdf("signal_electron");
  highMass_pdf_e->SetName("highMass_signal_electron");
  RooAbsPdf* highMass_pdf_mu = static_cast<RooWorkspace*>(highMass_file.Get("w"))->pdf("signal_muon");
  highMass_pdf_mu->SetName("highMass_signal_muon");
  highMass_file.Close();

  //double alpha = 1. - (double) (mass - lowMass) / (double) (highMass - lowMass);
  double alpha = (double) (mass - lowMass) / (double) (highMass - lowMass);

  RooRealVar mtt("mtt", "mtt", 500, 2000, "GeV/c^2");
  RooRealVar rAlpha("alpha", "alpha", alpha, 0, 1);

  mtt.setBins(1000, "cache");
  rAlpha.setBins(1000, "cache");

  // Interpolate
  std::cout << "Interpolate ..." << std::endl;
  //RooIntegralMorph interpolation_muon("signal_muon", "signal_muon", *lowMass_pdf_mu, *highMass_pdf_mu, mtt, rAlpha, true);
  //RooIntegralMorph interpolation_e("signal_electron", "signal_electron", *lowMass_pdf_e, *highMass_pdf_e, mtt, rAlpha, true);
  
  TVectorD hypoMass(2);
  hypoMass(0) = 0; 
  hypoMass(1) = 1;

  RooMomentMorph interpolation_muon("interpolated_signal_muon", "signal_muon", rAlpha, RooArgList(mtt), RooArgList(*lowMass_pdf_mu, *highMass_pdf_mu), hypoMass, RooMomentMorph::Linear);
  RooMomentMorph interpolation_e("interpolated_signal_electron", "signal_electron", rAlpha, RooArgList(mtt), RooArgList(*lowMass_pdf_e, *highMass_pdf_e), hypoMass, RooMomentMorph::Linear);
  
  std::cout << "Done." << std::endl;

  drawPdfs(mtt, interpolation_muon, *lowMass_pdf_mu, *highMass_pdf_mu, base_path + "/frit/" + prefix + "_interpolation_muon.pdf", btag, mass);
  drawPdfs(mtt, interpolation_e, *lowMass_pdf_e, *highMass_pdf_e, base_path + "/frit/" + prefix + "_interpolation_electron.pdf", btag, mass);

  // Transform PDF to RooHistPdf
  // Binning: 1 GeV/bin
  mtt.setBins((mtt.getMax() - mtt.getMin()) / 1.);
  RooDataHist binned_dataset_muon("binned_dataset_muon", "binned_dataset_muon", RooArgSet(mtt));
  interpolation_muon.fillDataHist(&binned_dataset_muon, NULL, 1.);

  RooDataHist binned_dataset_electron("binned_dataset_electron", "binned_dataset_electron", RooArgSet(mtt));
  interpolation_e.fillDataHist(&binned_dataset_electron, NULL, 1.);

  RooHistPdf hist_pdf_muon("signal_muon", "hist_pdf_muon", RooArgSet(mtt), binned_dataset_muon);
  RooHistPdf hist_pdf_electron("signal_electron", "hist_pdf_electron", RooArgSet(mtt), binned_dataset_electron);

  RooWorkspace workspace("w", "Interpolation signal workspace");
//  workspace.import(interpolation_muon);
//  workspace.import(interpolation_e);
  workspace.import(hist_pdf_muon);
  workspace.import(hist_pdf_electron);

  TString workspaceFile = base_path + "/frit/" + prefix + "_workspace.root";
  workspace.writeToFile(workspaceFile, true);
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Create interpolation between two pdf", ' ', "0.1");

    TCLAP::ValueArg<std::string> jecArg("", "jec", "Interpolate using this JEC.", false, "nominal", "string");
    TCLAP::ValueArg<int> massArg("m", "mass", "Zprime mass", true, 750, "integer");
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", true, 2, "int");

    cmd.add(jecArg);
    cmd.add(massArg);
    cmd.add(btagArg);

    cmd.parse(argc, argv);

    std::string jec = jecArg.getValue();

    // configure root
    gROOT->Clear();
    gROOT->SetStyle("Plain");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    setTDRStyle();

    // Set RooFit verbosity
    RooMsgService::instance().setStreamStatus(0,false);
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().addStream(RooFit::ERROR); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
    RooMsgService::instance().setSilentMode(true);

    // Disable stderr
    freopen("/dev/null", "w", stderr);
    doInterpolation(massArg.getValue(), jec, btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
  }
}
