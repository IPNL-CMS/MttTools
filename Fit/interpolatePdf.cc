#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <sys/wait.h>
#include <cmath>
#include <regex>

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

using namespace RooFit;

#include <list>
#include "tdrstyle.C"

#include "Utils.h"
#include "Functions.h"

#include <tclap/CmdLine.h>
#include <json/json.h>

#include <TString.h>

void doInterpolation(int mass, const std::string& jec, int btag) {

  std::string analysisName = getAnalysisName();
  TString prefix = TString::Format("%s-Zprime%d_%s_%d_btag", jec.c_str(), mass, analysisName.c_str(), btag);

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
    std::cerr << "Error: please use a mass between 750 and 1500 GeV." << std::endl;
    return;
  }

  // Open workspaces and retrieve PDFs
  TString lowMass_prefix = TString::Format("%s-Zprime%d_%s_%d_btag", jec.c_str(), lowMass, analysisName.c_str(), btag);
  TString lowMass_workspaceFile = "frit/" + lowMass_prefix + "_workspace.root";

  TString highMass_prefix = TString::Format("%s-Zprime%d_%s_%d_btag", jec.c_str(), highMass, analysisName.c_str(), btag);
  TString highMass_workspaceFile = "frit/" + highMass_prefix + "_workspace.root";

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

  double alpha = 1. - (double) (mass - lowMass) / (double) (highMass - lowMass);

  RooRealVar mtt("mtt", "mtt", 750, 2000, "GeV/c^2");
  RooRealVar rAlpha("alpha", "alpha", alpha, 0, 1);

  mtt.setBins(10, "cache");
  rAlpha.setBins(10, "cache");

  // Interpolate
  std::cout << "Interpolate ..." << std::endl;
  RooIntegralMorph interpolation_muon("signal_muon", "signal_muon", *lowMass_pdf_mu, *highMass_pdf_mu, mtt, rAlpha, true);
  RooIntegralMorph interpolation_e("signal_electron", "signal_electron", *lowMass_pdf_e, *highMass_pdf_e, mtt, rAlpha, true);
  std::cout << "Done." << std::endl;

  RooWorkspace workspace("w", "Interpolation signal workspace");
  workspace.import(interpolation_muon);
  workspace.import(interpolation_e);

  TString workspaceFile = "frit/" + prefix + "_workspace.root";
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
    doInterpolation(massArg.getValue(), jec, btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
  }
}

void saveParameter(int mass, const std::string& jecType, int btag, const std::string& name, double nSig, double nSig_err, double chi2) {

  FILE* lock = fopen("frit_efficiencies.lock", "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  if (fileExists("frit_efficiencies.json")) {
    std::ifstream file("frit_efficiencies.json");
    reader.parse(file, root);
    file.close();
  }

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  Json::Value& node = root[getAnalysisUUID()][strMass][btagStr][jecType][name];
  node["events"] = nSig;
  node["error"] = nSig_err;
  node["chi2"] = chi2;

  FILE* fd = fopen("frit_efficiencies.json", "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

void saveCombinedChiSquare(int mass, const std::string& jecType, int btag, double chi2) {

  FILE* lock = fopen("frit_efficiencies.lock", "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  if (fileExists("frit_efficiencies.json")) {
    std::ifstream file("frit_efficiencies.json");
    reader.parse(file, root);
    file.close();
  }

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  Json::Value& node = root[getAnalysisUUID()][strMass][btagStr][jecType];
  node["chi2"] = chi2;


  FILE* fd = fopen("frit_efficiencies.json", "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

void drawHistograms(RooAbsCategoryLValue& categories, RooRealVar& observable, int nBins, RooDataSet& dataset, RooSimultaneous& simPdfs, std::map<std::string, std::shared_ptr<BaseFunction>>& backgroundPdfs, int btag, const std::string& prefix, bool log) {

  std::cout << std::endl;
  std::cout << "Drawing histograms..." << std::endl;


  std::vector<std::shared_ptr<RooPlot>> plots;

  int n = categories.numTypes();
  int x = 2, y = (n + 0.5) / x;
  std::cout << "\t" << n << " categories. Dividing into " << x << "; " << y << std::endl;

  const int padWidth = 900;
  const int padHeight = 900;
  const float LUMI = 5.1;

  const int canvasWidth = padWidth * x;
  const int canvasHeight = padHeight * y;

  TCanvas *canvas = new TCanvas("canvas", "mTT fit", canvasWidth, canvasHeight);
  canvas->Divide(x, y);
  int currentPad = 1;

  TIterator* it = categories.typeIterator();
  RooCatType* type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    canvas->cd(currentPad++);

    RooPlot* plot = observable.frame(nBins);
    plot->SetTitle(""); //FIXME

    std::string category = type->GetName();
    categories = category.c_str();

    std::stringstream ss;

    RooSuperCategory* foo = dynamic_cast<RooSuperCategory*>(&categories);
    if (foo) {
      const RooArgSet& parentCategories = foo->inputCatList();
      TIterator *it2 = parentCategories.createIterator();
      RooAbsCategory* parentCategory = nullptr;
      bool first = true;
      while ((parentCategory = static_cast<RooAbsCategory*>(it2->Next()))) {
        ss << (first ? "" : " && ") << parentCategory->GetName() << " == " << parentCategory->GetName() << "::" << parentCategory->getLabel();
        first = false;
      }
    } else {
      ss << categories.GetName() << " == " << categories.GetName() << "::" << category;
    }

    std::string cut = ss.str();
    //std::cout << "Cut: " << cut << std::endl;

    dataset.plotOn(plot, Cut(cut.c_str()));

    if (backgroundPdfs[category].get())
      simPdfs.plotOn(plot, Slice(categories), ProjWData(categories, dataset), Components(backgroundPdfs[category]->getPdf()), LineStyle(kDashed), LineColor(kRed), LineWidth(2));
    simPdfs.plotOn(plot, Slice(categories), ProjWData(categories, dataset), LineColor(kBlue), LineWidth(2));

    //simPdfs.plotOn(plot, Slice(categories, category.c_str()), ProjWData(categories, dataset), Components(backgroundPdfs[category]->getPdf()), LineStyle(kDashed), LineColor(kRed), LineWidth(2));
    //simPdfs.plotOn(plot, Slice(categories, category.c_str()), ProjWData(categories, dataset), LineColor(kBlue), LineWidth(2));

    std::string leptonName = TString(category).Contains("muon", TString::kIgnoreCase) ? "muons" : "electrons";
    float binningSize = (observable.getBinning().highBound() - observable.getBinning().lowBound()) / (float) nBins;

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

    gPad->SetLogy(log);
    canvas->SetLogy(log);

    plot->Draw();

    // Find number of b-tag
    std::regex regex("([0-9]+)-btag", std::regex_constants::basic);
    std::smatch regexResults;

    if (std::regex_search(category, regexResults, regex)) {
      std::cout << "Found " << regexResults[1] << " b-tag" << std::endl;
      std::string result(regexResults[1].first, regexResults[2].second);
      btag = atoi(result.c_str());
    } else {
      std::cout << "No b-tag found." << std::endl;
    }

    TString btagLabel = "";
    if (btag == 2)
      btagLabel = "#geq 2 b-tags";
    else
      btagLabel = TString::Format("%d b-tag", btag);

    std::string leptonShortcutName = TString(category).Contains("muon", TString::kIgnoreCase) ? "#mu" : "e";

    TString legendLabel = TString::Format("#font[42]{%s, #geq 4 jets, %s}", leptonShortcutName.c_str(), btagLabel.Data());

    t.DrawLatex(0.53, 0.88, "#font[42]{CMS preliminary}");
    t.DrawLatex(0.53, 0.84, TString::Format("#font[42]{%0.2f fb^{-1} at #sqrt{s}=8 TeV}", LUMI));
    t.DrawLatex(0.53, 0.80, legendLabel);

    plots.push_back(std::shared_ptr<RooPlot>(plot));

    canvas->Update();
  }

  if (! log) {
    canvas->Print((prefix + "_fitCB.pdf").c_str());
    canvas->Print((prefix + "_fitCB.png").c_str());
  } else {
    canvas->Print((prefix + "_fitCB_log.pdf").c_str());
    canvas->Print((prefix + "_fitCB_log.png").c_str());
  }
  
  delete canvas;

  std::cout << std::endl;
}

void fritSignal(TChain* chain, const std::string& jecType, const std::string& configFile, int massZprime, int btag) {

  std::cout << "[" << getpid() << "] Processing for " << jecType << std::endl;

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

  //fit region
  Float_t minmTT = 500;
  Float_t maxmTT = 2000;
  Int_t nBins = 50;

  RooRealVar mtt("mtt", "mtt", minmTT, maxmTT, "GeV/c^2");
  RooRealVar weight("weight", "weight", 0, 10000);

  RooCategory whichLepton("lepton_type", "lepton_type");
  whichLepton.defineType("electron", 11);
  whichLepton.defineType("muon", 13);

  // Get dataset from external file
  //RooDataSet *dataset = RooDataSet::read(file.c_str(), RooArgList(Mtt_KF_reco, whichLepton));
  RooDataSet *dataset = new RooDataSet("dataset", "dataset", RooArgSet(mtt, whichLepton, weight), Import(*chain)/*, WeightVar(weight)*/);

  // Reduce data set
  //RooDataSet* RedData = static_cast<RooDataSet*>(dataOrig->reduce("Mtt_KF_reco > 500 && Mtt_KF_reco < 2045"));

  if (dataset->numEntries() == 0) {
    std::cerr << "[" << getpid() << "] ERROR: No entry in dataset" << std::endl;
    return;
  } else {
    std::cout << "[" << getpid() << "] " << dataset->numEntries() << " entries loaded" << std::endl;
  }

  std::string analysisName = getAnalysisName();

  std::map<std::string, std::shared_ptr<BaseFunction>> backgroundPdfs = getCategoriesPdf("./fit_configuration", configFile, mtt, NULL, massZprime, "background", whichLepton, nullptr);

  // It takes some time, inform the user of what's going on
  std::cout << "Loading signal pdfs... (may takes some time)" << std::endl;
  std::map<std::string, std::shared_ptr<BaseFunction>> signalPdfs = getCategoriesPdf("./fit_configuration", configFile, mtt, dataset, massZprime, "signal", whichLepton, NULL);
  std::cout << "Done." << std::endl;

  for (auto& pdf: backgroundPdfs) {
    std::cout << "Background Pdf: " << pdf.first << " -> ";
    pdf.second->getPdf().Print();
  }

  for (auto& pdf: signalPdfs) {
    std::cout << "Signal Pdf: " << pdf.first << " -> ";
    pdf.second->getPdf().Print();
  }

  TString prefix = TString::Format("%s-Zprime%d_%s_%d_btag", jecType.c_str(), massZprime, analysisName.c_str(), btag);

  std::map<std::string, std::shared_ptr<RooAbsPdf>> globalPdfs;
  std::map<std::string, std::shared_ptr<RooRealVar>> globalPdfsEvents;
  std::map<std::string, std::vector<std::shared_ptr<RooRealVar>>> globalPdfsParameters;

  TIterator* it = whichLepton.typeIterator();
  RooCatType * type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {

    // Create parameters
    std::string name = type->GetName();

    std::shared_ptr<RooRealVar> numberOfEvents { std::shared_ptr<RooRealVar>(new RooRealVar((name + "_nSig").c_str(), "number of sig events", 5000., 0, 200000.)) };

    globalPdfsEvents[name] = numberOfEvents;

    if (signalPdfs[name]->isExtended()) {
      // Create Extended pdf for fit
      globalPdfs[name] = std::shared_ptr<RooAbsPdf>(new RooExtendPdf((name + "_global_pdf").c_str(), "extended signal pdf", signalPdfs[name]->getPdf(), *numberOfEvents));
      backgroundPdfs[name].reset();
    } else {
      std::vector<std::shared_ptr<RooRealVar>> parameters {
        std::shared_ptr<RooRealVar>(new RooRealVar((name + "_nSig").c_str(), "number of sig events", 5000., -100., 20000.)),
        std::shared_ptr<RooRealVar>(new RooRealVar((name + "_nBkg").c_str(), "number of bkg events", 300., 0., 50000))
      };

      globalPdfsParameters[name] = parameters;
      globalPdfs[name] = std::shared_ptr<RooAbsPdf>(new RooAddPdf((name + "_global_pdf").c_str(), "signal + background", RooArgList(signalPdfs[name]->getPdf(), backgroundPdfs[name]->getPdf()), RooArgList(*parameters[0], *parameters[1])));
    }
  }

  RooSimultaneous simPdf("simPdf","simultaneous pdf",whichLepton) ;

  it = whichLepton.typeIterator();
  type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    std::string name = type->GetName();

    const RooAbsPdf& pdf = *globalPdfs[name];
    //const RooAbsPdf& pdf = signalPdfs[name]->getPdf();
    std::cout << "Adding pdf ";
    pdf.Print();
    std::cout << " for category " << name << std::endl;
      
    simPdf.addPdf(pdf, name.c_str());
  }

  Double_t minmTTFit = minmTT + 0.0;
  Double_t maxmTTFit = maxmTT - 0.0;
  mtt.setRange(minmTTFit, maxmTTFit);

  RooFitResult* fitResult = NULL;
  bool fitIsGood = false;
  int fitIterations = 5;
  do {
    delete fitResult;

    fitResult = simPdf.fitTo(*dataset, Save(), Optimize(1), RooFit::Strategy(1));

    fitIsGood = !std::isnan(fitResult->edm());
    fitIterations--;

    if (! fitIsGood) {
      std::cout << "[" << getpid() << "] Warning: Fit return a 'NaN' estimated distance to minimum. Re-running fit." << std::endl;
    }

  } while (! fitIsGood && fitIterations > 0);

  drawHistograms(whichLepton, mtt, nBins, *dataset, simPdf, backgroundPdfs, btag, std::string("frit/" + prefix), false);
  drawHistograms(whichLepton, mtt, nBins, *dataset, simPdf, backgroundPdfs, btag, std::string("frit/" + prefix), true);

  fitResult->Print("v");

  // Compute Chi2
  // We need to transform our unbinned dataset to a binned one
  // Set a resolution of 10 GeV, which is much greater than our detector resolution
  const float resolution = 5.;
  const int nBinsForChi2 = (mtt.getMax() - mtt.getMin() + 0.5) / resolution;
  mtt.setBins(nBinsForChi2);
  std::cout << "Binning dataset with " << nBinsForChi2 << " bins for chi2 computation (" << mtt.getMin() << " -> " << mtt.getMax() << " ; resolution: " << resolution << " GeV)" << std::endl;

  RooArgSet* floatingParameters = static_cast<RooArgSet*>(simPdf.getParameters(RooArgSet(mtt))->selectByAttrib("Constant", false));
  int numberOfFloatingParams = floatingParameters->getSize() - 1;
  delete floatingParameters;

  Roo1DTable * table = dataset->table(whichLepton);

  // Chi2 computation inspired from http://cmssdt.cern.ch/SDT/doxygen/CMSSW_5_1_1/doc/html/d4/d3a/GenericTnPFitter_8h_source.html, line 291

  std::map<std::string, float> chiSquares;
  float combinedChi2 = 0.;
  it = whichLepton.typeIterator();
  type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    std::cout << "For category '" << type->GetName() << "':" << std::endl;
    if (signalPdfs[type->GetName()]->isExtended()) {
      std::cout << "\tNumber of signal events: " << globalPdfsEvents[type->GetName()]->getVal() << " out of " << table->get(type->GetName()) << std::endl;
    } else {
      std::cout << "\tNumber of signal events: " << globalPdfsParameters[type->GetName()][0]->getVal() << " out of " << table->get(type->GetName()) << std::endl;
      std::cout << "\tNumber of background events: " << globalPdfsParameters[type->GetName()][1]->getVal() << " out of " << table->get(type->GetName()) << std::endl;

    }

    std::stringstream cut;
    cut << whichLepton.GetName() << " == " << type->getVal();

    RooDataHist* binnedDataset = new RooDataHist("binnedDataset", "binned dataset", RooArgSet(mtt), *dataset->reduce(Cut(cut.str().c_str())));

    float chi2 = RooChi2Var("chi2", "chi2", *globalPdfs[type->GetName()], *binnedDataset, DataError(RooAbsData::Poisson)).getVal();
    std::cout << "\tChi2: " << (chi2 / (mtt.getBins() - numberOfFloatingParams * .5)) << std::endl;
    chiSquares[type->GetName()] = (chi2 / (mtt.getBins() - numberOfFloatingParams * .5));

    combinedChi2 += chi2;

    delete binnedDataset;
  }

  std::cout << std::endl;

  delete table;

  combinedChi2 /= (2 * mtt.getBins() - numberOfFloatingParams);
  std::cout << "Combined Chi2: " << combinedChi2 << std::endl;

  //
  // pour savoir les noms des objets dans un RooPlot: MttPlot_mu->Print();
  // Attention la définition ci dessous depend de l'ordre dans lequel on trace dans le rooplot  MttPlot_mu
  // De plus le nombre de degre de liberté n'est pas facilement calculable car c'est une fit simultané
  // on met nfloatparam/2 , on s'en contente pour l'instant.
  /*
  Double_t chi2mu = MttPlot_mu->chiSquare(MttPlot_mu->nameOf(1),MttPlot_mu->nameOf(0),nfloatparam/2);
  Double_t chi2e = MttPlot_e->chiSquare(MttPlot_e->nameOf(1),MttPlot_e->nameOf(0),nfloatparam/2); // 3 for falt and faltB
  cout << "[" << getpid() << "] Chi2/Ndf for muons = " << chi2mu << "  |  Chi2/Ndf for electrons = " << chi2e << endl;
  Double_t chi2tot=nBins/((float)2*nBins-nfloatparam)*(MttPlot_mu->chiSquare(MttPlot_mu->nameOf(1),MttPlot_mu->nameOf(0))+MttPlot_e->chiSquare(MttPlot_e->nameOf(1),MttPlot_e->nameOf(0)));
  cout << "[" << getpid() << "] Chi2/Ndf total = " << chi2tot << endl;
  */

  // write efficiencies to file
  /*
  ofstream outEffFile(prefix+"_efficiencies.txt");
  outEffFile << "#Found nZprime, semi-mu " << nSig_mu.getVal() 
    << " +- " << nSig_mu.getError() << " out of n.mu = " << nEventsRed_mu.getVal() << std::endl;
  outEffFile << "#Found nZprime, semi-e  " << nSig_e.getVal() 
    << " +- " << nSig_e.getError() << " out of n.e = " << nEventsRed_e.getVal() << std::endl;
  outEffFile << "#Chi2/Ndf for muons = " << chi2mu << std::endl;
  outEffFile << "#Chi2/Ndf for electrons = " << chi2e << std::endl;
  outEffFile << "#Chi2/Ndf total = " << chi2tot << std::endl;
  outEffFile << "# The following values are for computeEff" << std::endl;
  outEffFile << nSig_mu.getVal() << " " << nSig_mu.getError() << " " << nSig_e.getVal() << " " << nSig_e.getError() << std::endl;
  outEffFile.close();
  */

  //FIXME
  it = whichLepton.typeIterator();
  type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    std::string name = type->GetName();
    
    if (signalPdfs[name]->isExtended()) {
      saveParameter(massZprime, jecType, btag, name, globalPdfsEvents[name]->getVal(), globalPdfsEvents[name]->getError(), chiSquares[name]);
    } else {
      saveParameter(massZprime, jecType, btag, name, globalPdfsParameters[name][0]->getVal(), globalPdfsParameters[name][0]->getError(), chiSquares[name]);
    }
  }
  saveCombinedChiSquare(massZprime, jecType, btag, combinedChi2);

  // Save our signal functions into the workspace
  // We will need it for the fit on data
  TString workspaceFile = "frit/" + prefix + "_workspace.root";
  RooWorkspace workspace("w", "Frit signal workspace");

  for (auto& pdf: signalPdfs) {
    workspace.import(pdf.second->getPdf());
  }

  workspace.import(*fitResult, "fitResult");
  workspace.writeToFile(workspaceFile, true);

  delete fitResult;
}
