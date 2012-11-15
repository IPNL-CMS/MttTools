#include <math.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <sstream>

#include <TROOT.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>

#include "tdrstyle.C"

#include "Utils.h"

#include <tclap/CmdLine.h>
#include <json/json.h>

void loadObservedLimits(int btag, const std::vector<double>& masses, std::vector<double>& obs) {
  Json::Reader reader;
  Json::Value root;
  std::ifstream file("likelihood_scan.json");
  bool success = reader.parse(file, root);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "likelihood_scan.json" << "'. Exiting." << std::endl;
    std::cerr << "Please run the likelihood scan using 'runLikelihoodScan'" << std::endl;
    exit(1);
  }

  std::stringstream ss2;
  ss2 << btag;
  std::string btagStr = ss2.str();

  for (std::vector<double>::const_iterator it = masses.begin(); it != masses.end(); ++it) {
    std::stringstream ss;
    ss << (int) *it;
    std::string strMass = ss.str();

    if (! root.isMember(strMass)) {
      std::cerr << "ERROR: mass '" << *it << "' not found in JSON file. Exiting." << std::endl;
      exit(1);
    }

    double observedLimit = root[strMass][btagStr]["scan_wsyst_cut_limit"].asDouble();
    obs.push_back(observedLimit);
  }
}

void loadExpectedLimits(int btag, const std::vector<double>& masses, std::vector<double>& exp, std::vector<double>& error_h_95, std::vector<double>& error_l_95, std::vector<double>& error_h_68, std::vector<double>& error_l_68) {
  Json::Reader reader;
  Json::Value root;
  std::ifstream file("expected_limits.json");
  bool success = reader.parse(file, root);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "expected_limits.json" << "'. Exiting." << std::endl;
    std::cerr << "Please run 'treatToyStuff'" << std::endl;
    exit(1);
  }

  std::stringstream ss2;
  ss2 << btag;
  std::string btagStr = ss2.str();

  for (std::vector<double>::const_iterator it = masses.begin(); it != masses.end(); ++it) {
    std::stringstream ss;
    ss << (int) *it;
    std::string strMass = ss.str();

    if (! root.isMember(strMass)) {
      std::cerr << "ERROR: mass '" << *it << "' not found in JSON file. Exiting." << std::endl;
      exit(1);
    }

    double expectedLimit = root[strMass][btagStr]["median"].asDouble();
    double eyl68 = root[strMass][btagStr]["widthM68"].asDouble();
    double eyh68 = root[strMass][btagStr]["widthP68"].asDouble();
    double eyl95 = root[strMass][btagStr]["widthM95"].asDouble();
    double eyh95 = root[strMass][btagStr]["widthP95"].asDouble();

    exp.push_back(expectedLimit);
    error_h_95.push_back(eyh95);
    error_l_95.push_back(eyl95);
    error_h_68.push_back(eyh68);
    error_l_68.push_back(eyl68);
  }
}

void drawLimitCurve(int btag) {

  gROOT->Clear();
  gStyle->SetOptStat(0);
  setTDRStyle();

  std::vector<double> masses;
  masses.push_back(750);
  masses.push_back(1000);
  masses.push_back(1250);
  masses.push_back(1500);

  std::vector<double> observed; // Observed limits
  std::vector<double> expected; // Expected limits

  std::vector<double> error_high_95; // 95%
  std::vector<double> error_low_95; // 95%

  std::vector<double> error_high_68;
  std::vector<double> error_low_68;

  double no_error[4] = {0};

  double xsections[4] = {
    2.74 * 1.3,
    0.753 * 1.3,
    0.236 * 1.3,
    0.082 * 1.3};

  loadObservedLimits(btag, masses, observed);
  loadExpectedLimits(btag, masses, expected, error_high_95, error_low_95, error_high_68, error_low_68);

  TCanvas* c1 = new TCanvas("canvas", "canvas", 600, 600);

  TGraphAsymmErrors* gra95 = new TGraphAsymmErrors(4, &masses[0], &expected[0], no_error, no_error, &error_low_95[0], &error_high_95[0]);//the 95% band of the 95% expected limits
  gra95->SetMinimum(0);
  gra95->SetMaximum(20);
  gra95->SetTitle("95% C.L. upper limit on the cross section");
  gra95->GetXaxis()->SetTitle("m(t#bar{t}) (GeV/c^{2})")	  ;
  gra95->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z') #times BR(Z' #rightarrow t#bar{t}) (pb)");
  gra95->SetFillColor(15);
  gra95->Draw("AE3");

  TGraphAsymmErrors* gra = new TGraphAsymmErrors(4, &masses[0], &expected[0], no_error, no_error, &error_low_68[0], &error_high_68[0]);//the 68% band of the 95% expected limits
  gra->SetMinimum(0);
  gra->SetMaximum(20);
  gra->SetTitle("95% C.L. upper limit on the cross section");
  gra->GetXaxis()->SetTitle("m(t#bar{t}) (GeV/c^{2})")	  ;
  gra->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z') #times BR(Z' #rightarrow t#bar{t}) (pb)");
  gra->SetFillColor(18);
  gra->Draw("E3,same");

  TGraph* gra2 = new TGraph(4, &masses[0], &expected[0]); //the expected limits
  gra2->SetMarkerStyle(20);
  gra2->SetMarkerSize(1.2);
  gra2->Draw("P,same");

  /*TGraph* gra3 = new TGraph(4, &masses[0], xsections); //
  gra3->SetMarkerStyle(20);
  gra3->SetMarkerSize(1.2);
  gra3->SetMarkerColor(6);
  gra3->SetLineWidth(3);
  gra3->SetLineColor(6);
  gra3->Draw("C,same");*/

  TGraphErrors* graobs = new TGraphErrors(4, &masses[0], &observed[0], no_error, no_error);//the observed limits
  graobs->SetMarkerStyle(20);
  graobs->SetMarkerColor(2);
  graobs->SetMarkerSize(1.2);
  graobs->Draw("P,same");

  TLegend* legendLimit = new TLegend(0.30,0.72,0.89,0.89);
  legendLimit->AddEntry(gra2,"Expected Upper Limit (95% CL)","p");
  legendLimit->AddEntry(graobs,"Observed Upper  Limit (95% CL)","p");
  legendLimit->AddEntry(gra,"95% CL exclusion: 68% band","f");
  legendLimit->AddEntry(gra95,"95% CL exclusion: 95% band","f");
  //legendLimit->AddEntry(gra3, "Z' leptophobic", "l");
  legendLimit->SetFillStyle(0);
  legendLimit->Draw("SAME");

  TLatex t;
  t.SetTextSize(0.04);
  t.SetNDC();
  t.DrawLatex(0.63, 0.67, "#font[42]{CMS preliminary}");
  t.DrawLatex(0.63, 0.62, "#font[42]{5 fb^{-1} at #sqrt{s}=8 TeV}");

  TString prefix = TString::Format("limitCurve_2012_%s_%s_%dbtag", getSignalPdfName().c_str(), getFitBackgroundPdfName().c_str(), btag);

  c1->Print(prefix + ".png");
  c1->SaveAs(prefix + ".pdf");

  delete c1;
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Draw limit curve", ' ', "0.1");

    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", true, 2, "int", cmd);

    cmd.parse(argc, argv);

    drawLimitCurve(btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
