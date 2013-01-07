#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include "tdrstyle.C"

#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>

#include <TROOT.h>
#include <TH1.h>

#include "Utils.h"
#include <tclap/CmdLine.h>
#include <json/json.h>

#include <sys/types.h>
#include <sys/stat.h>

std::string INPUT_PATH;

void saveToyLimits(const int mass, const int btag, const double mean, const double rms, const double median, const double m68, const double p68,
    const double m95, const double p95) {

  const double widthM68 = fabs(m68 - median);
  const double widthP68 = fabs(p68 - median);
  const double widthM95 = fabs(m95 - median);
  const double widthP95 = fabs(p95 - median);

  //const double widthM68 = m68 - median;
  //const double widthP68 = p68 - median;
  //const double widthM95 = m95 - median;
  //const double widthP95 = p95 - median;

  Json::Reader reader;
  Json::Value root;
  std::ifstream file("expected_limits.json");
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str("");
  ss << btag;
  std::string btagStr = ss.str();

  std::string uuid = getAnalysisUUID();
  
  root[uuid][strMass][btagStr]["mean"] = mean;
  root[uuid][strMass][btagStr]["rms"] = rms;
  root[uuid][strMass][btagStr]["median"] = median;
  root[uuid][strMass][btagStr]["m68"] = m68;
  root[uuid][strMass][btagStr]["p68"] = p68;
  root[uuid][strMass][btagStr]["m95"] = m95;
  root[uuid][strMass][btagStr]["p95"] = p95;
  root[uuid][strMass][btagStr]["widthM68"] = widthM68;
  root[uuid][strMass][btagStr]["widthP68"] = widthP68;
  root[uuid][strMass][btagStr]["widthM95"] = widthM95;
  root[uuid][strMass][btagStr]["widthP95"] = widthP95;

  std::ofstream ofile;
  ofile.open("expected_limits.json", std::ios::out | std::ios::trunc);
  Json::StyledWriter writer;
  ofile << writer.write(root);
  ofile.close();
}


void treatToyStuff(bool writeTxtFile, bool savePsGifFiles, int massZprime, int btag) {


  gROOT->Clear();  
  setTDRStyle();

  std::string analysisName = getAnalysisName();

  TString str = TString::Format("%s/%d-btag/data_2012_nominal_%d_toylimit_%s.root", INPUT_PATH.c_str(), btag, massZprime, analysisName.c_str());
  TFile* f0 = TFile::Open(str);

  TCanvas* cToy = new TCanvas("cToy", "Toy MC pulls", 800, 800);
  cToy->Divide(1, 2);
  cToy->cd(1);
  TH1* hSigFrac_mu = (TH1*) f0->Get("hSigFrac_mu");
  hSigFrac_mu->Draw();
  hSigFrac_mu->Fit("gaus");
  cToy->cd(2);
  //hLimit_mu->Draw();
  TH1* hResidual_mu = (TH1*) f0->Get("hResidual_mu");
  hResidual_mu->Draw();
  hResidual_mu->Fit("gaus");
  cToy->Draw();

  TCanvas* cToy2 = new TCanvas("cToy2", "Toy MC pulls", 400, 400);
  //    cToy2->cd();   
  TH1* hSigma = (TH1*) f0->Get("hSigma");
  hSigma->Draw();
  hSigma->Fit("gaus");

  TCanvas* cerr = new TCanvas("cerr", "Toy MC errors", 400, 400);
  cerr->Divide(1,2);
  cerr->cd(1);
  TH1* hLoErr_mu = (TH1*) f0->Get("hLoErr_mu");
  hLoErr_mu->Draw();
  cerr->cd(2);
  TH1* hHiErr_mu = (TH1*) f0->Get("hHiErr_mu");
  hHiErr_mu->Draw();

  TCanvas* cLik = new TCanvas("cLik","Likelihoods", 400,400); 
  //   cLik->cd();
  TH1* hMinNll = (TH1*) f0->Get("hMinNll");
  hMinNll->Draw();
  /*TArrow* arrow = new TArrow(minNllFit, nToyExp/50., minNllFit, 0., 0.03, ">");
    arrow->SetLineWidth(2.);
    arrow->SetLineColor(2);
    arrow->Draw();*/
  TCanvas* cZ = new TCanvas("cZ","Z cross section", 400,400); 
  //    cZ->cd();
  TH1* hLimit_Z = (TH1*) f0->Get("hLimit_Z");
  //hLimit_Z->GetXaxis()->SetTitle("95% C.L. upper limit on #sigma(pp #rightarrow Z') #times BR(Z' #rightarrow t#bar{t}) (pb)");
  hLimit_Z->GetXaxis()->SetRangeUser(0.,10.);
  hLimit_Z->Draw();
  std::cout << "Mean of the Z limit distribution " << hLimit_Z->GetMean() << std::endl;
  std::cout << "RMS of the Z limit distribution " << hLimit_Z->GetRMS() << std::endl;

  double areaq = 0.50;
  double medianq = 0.;
  hLimit_Z->GetQuantiles(1, &medianq, &areaq);
  std::cout << "Median of the Z limit distribution " << medianq << std::endl;

  double areaqP68 = 0.84;
  double P68band = 0.;
  hLimit_Z->GetQuantiles(1, &P68band, &areaqP68);

  double areaqM68 = 0.16;
  double M68band = 0.;
  hLimit_Z->GetQuantiles(1, &M68band, &areaqM68);

  double areaqP95 = 0.975;
  double P95band = 0.;
  hLimit_Z->GetQuantiles(1, &P95band, &areaqP95);    

  double areaqM95 = 0.025;
  double M95band = 0.;
  hLimit_Z->GetQuantiles(1, &M95band, &areaqM95);

  std::cout << "68% band around the median " << M68band << " " << P68band << std::endl;
  std::cout << "95% band around the median " << M95band << " " << P95band << std::endl;
  std::cout << "width of the 68% band around the median " << M68band - medianq << " " << P68band - medianq << std::endl;
  std::cout << "width of the 95% band around the median " << M95band - medianq << " " << P95band - medianq << std::endl;


  if (savePsGifFiles) {
    mkdir(TString::Format("toys/plots/%d-btag", btag).Data(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    TString prefix = TString::Format("toys/plots/%d-btag/data_2012_Zprime%d_%s", btag, massZprime, analysisName.c_str());

    cToy->Print(TString::Format("%s_LimitPlotsEMu.gif", prefix.Data()));
    cToy2->Print(TString::Format("%s_LimitPlotsSigma.gif", prefix.Data()));
    cZ->Print(TString::Format("%s_LimitPlotZ.gif", prefix.Data()));
    cerr->Print(TString::Format("%s_LimitErrors.gif", prefix.Data()));
    cLik->Print(TString::Format("%s_LimitNLLToyExp.gif", prefix.Data()));
    cToy->Print(TString::Format("%s_LimitPlotsEMu.pdf", prefix.Data()));
    cToy2->Print(TString::Format("%s_pull.pdf", prefix.Data()));
    cZ->Print(TString::Format("%s_LimitPlotZ.pdf", prefix.Data()));
    cerr->Print(TString::Format("%s_LimitErrors.pdf", prefix.Data()));
    cLik->Print(TString::Format("%s_LimitNLLToyExp.pdf", prefix.Data()));
  }

  /*
  if (writeTxtFile) {
    ofstream outParFileLimit(TString::Format("data_2012_Zprime%d_fitRes_ToyLimit_median_%s.txt", massZprime, pdfSignalName.c_str()));
    outParFileLimit << "Mean of the Z limit distribution " << hLimit_Z->GetMean() << std::endl;
    outParFileLimit << "RMS of the Z limit distribution " << hLimit_Z->GetRMS() << std::endl;
    outParFileLimit << "Median of the Z limit distribution " << medianq << std::endl;
    outParFileLimit << "68% band around the median " << M68band << " " << P68band << std::endl;
    outParFileLimit << "95% band around the median " << M95band << " " << P95band << std::endl;      
    outParFileLimit << "width of the 68% band around the median " << M68band - medianq << " " << P68band - medianq << std::endl;
    outParFileLimit << "width of the 95% band around the median " << M95band - medianq << " " << P95band - medianq << std::endl;

    outParFileLimit.close();
  }
  */

  saveToyLimits(massZprime, btag, hLimit_Z->GetMean(), hLimit_Z->GetRMS(), medianq, M68band, P68band, M95band, P95band);

  f0->Close();
  delete f0;
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Treat toy infos", ' ', "0.1");

    std::string defaultInputPath = "analysis/" + getAnalysisUUID() + "/toys/results/";

    TCLAP::SwitchArg writeArg("", "dont-write-root", "Don't save root files", cmd, true);
    TCLAP::SwitchArg saveArg("", "dont-save", "Don't save images", cmd, true);
    TCLAP::ValueArg<std::string> inputPathArg("", "input-path", "Where loading files", false, defaultInputPath, "string", cmd);
    TCLAP::MultiArg<int> massArg("m", "mass", "Zprime mass", false, "integer", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", true, 2, "int", cmd);

    cmd.parse(argc, argv);

    std::vector<int> masses = massArg.getValue();
    if (masses.size() == 0) {
      masses.push_back(750);
      masses.push_back(1000);
      masses.push_back(1250);
      masses.push_back(1500);
    }

    INPUT_PATH = inputPathArg.getValue();

    //FIXME: Parallelize
    for (std::vector<int>::iterator mass = masses.begin(); mass != masses.end(); ++mass) {
      treatToyStuff(writeArg.getValue(), saveArg.getValue(), *mass, btagArg.getValue());
    }

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
