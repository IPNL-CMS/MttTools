#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <boost/regex.hpp>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include "tdrstyle.C"

#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TDatime.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TApplication.h>
#include <TLatex.h>
#include <TChain.h>
#include <TGaxis.h>
#include <TVectorD.h>

#include <RooGlobalFunc.h>
#include <RooMsgService.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooFitResult.h>
#include <TMatrixDSym.h>
#include <RooCBShape.h>
#include <RooGenericPdf.h>
#include <RooAddPdf.h>
#include <RooSimultaneous.h>
#include <RooFormulaVar.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooRandom.h>
#include <RooExtendPdf.h>
#include <RooExponential.h>
#include <RooHist.h>
#include <RooPlotable.h>
#include <RooFitResult.h>
#include <RooMinimizer.h>
#include <RooMinuit.h>
#include <RooWorkspace.h>
#include <RooSuperCategory.h>
#include <Roo1DTable.h>
#include <RooDataHist.h>
#include <RooChi2Var.h>
#include <RooUniform.h>
#include <RooIntegralMorph.h>
#include <RooMomentMorph.h>
#include <RooHistPdf.h>

#include "Utils.h"
#include <tclap/CmdLine.h>
#include <json/json.h>

#include "SignalFunctions.h"
#include "Functions.h"

#define SAFE_DELETE(p) { delete p; p = NULL; }

using namespace RooFit;

// Nice output on Bash
namespace Bash {
  enum Color {
    NONE = 0, BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE
  };

  std::string set_color(Color foreground = NONE, Color background = NONE) {
    std::stringstream num_s;
    std::string s = "\033[";

    if (!foreground && ! background) s += "0"; // reset colors if no params

    if (foreground) {
      num_s << (29 + foreground);
      s += num_s.str();

      if (background) s += ";";
    }

    if (background) {
      num_s << (39 + background);
      s += num_s.str();
    }

    return s + "m";
  }

}

void fitMtt(std::map<int, TChain*> eventChain, int massZprime, string fitConfigurationFile, int btag, const std::string& customWorkspaceFile, const std::string& outputFile);

std::string BASE_PATH;
std::string EFF_FILE;
bool VERBOSE = false;

void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

TChain* loadChain(const std::vector<std::string>& files, const std::string& name) {
  TChain* c = new TChain(name.c_str());

  for (const std::string& file: files) {
    c->Add(file.c_str());
  }

  return c;
}

int main(int argc, char** argv)
{
  try
  {
    TCLAP::CmdLine cmd("Fit tt mass spectrum", ' ', "0.1");

    TCLAP::ValueArg<int> massArg("m", "mass", "Zprime mass", true, 750, "integer");
    cmd.add(massArg);

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    std::string default_base_path = "./analysis/" + getAnalysisUUID() + "/";
    TCLAP::ValueArg<std::string> pathArg("", "path", "Folder where to load files", false, default_base_path, "string", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);
    TCLAP::SwitchArg verboseArg("v", "verbose", "Verbose mode", cmd);
    TCLAP::ValueArg<std::string> fitConfigFileArg("", "config-file", "Configuration file name containing fit parameters", false, "fit_pdf_faltb.json", "string", cmd);
    TCLAP::ValueArg<std::string> effArg("", "eff-file", "File where efficiences are stored (JSON format)", false, "efficiencies.json", "string", cmd);

    // Workspace
    TCLAP::ValueArg<std::string> workspaceArg("", "workspace", "Use a custom workspace containing the signal pdf to use", false, "", "string", cmd);

    TCLAP::ValueArg<std::string> outputArg("o", "output-file", "Use a custom workspace containing the signal pdf to use", false, "", "string", cmd);

    cmd.parse(argc, argv);

    BASE_PATH = pathArg.getValue();
    BASE_PATH.erase(std::remove(BASE_PATH.begin(), BASE_PATH.end(), '\"'), BASE_PATH.end());

    EFF_FILE = effArg.getValue();

    VERBOSE = verboseArg.getValue();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    std::map<int, TChain*> chains;
    if (btagArg.getValue() == 3) {
      for (int i = 0; i < 3; i++) {
        // Format tree name
        TString treeName = TString::Format("dataset_%dbtag", i);
        chains[i] = loadChain(inputFiles, treeName.Data());
      }

    } else {
      // Format tree name
      TString treeName = TString::Format("dataset_%dbtag", btagArg.getValue());
      chains[btagArg.getValue()] = loadChain(inputFiles, treeName.Data());
    }

    fitMtt(chains, massArg.getValue(), fitConfigFileArg.getValue(), btagArg.getValue(), workspaceArg.getValue(), outputArg.getValue());

    for (auto& chain: chains)
      delete chain.second;

  }
  catch (TCLAP::ArgException& e)
  {
    std::cerr << e.error() << std::endl;
  }
}

void loadEfficiencies(int mass, const std::string& jecType, int btag, double& a, double& b, double& c, double& d, double& e,
    double& f, double& g, double& h)
{

  std::string path = (BASE_PATH + EFF_FILE);
  Json::Reader reader;
  Json::Value root;
  std::ifstream file(path.c_str());
  bool success = reader.parse(file, root);
  file.close();
  if (! success)
  {
    std::cerr << "ERROR: Failed to parse '" << path << "'. Exiting." << std::endl;
    exit(1);
  }

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  root = root[getAnalysisUUID(BASE_PATH)];

  if (! root.isMember(strMass))
  {
    std::cerr << "ERROR: mass '" << mass << "' not found in JSON file. Exiting." << std::endl;
    exit(1);
  }

  Json::Value massNode = root[strMass][btagStr];

  if (! massNode.isMember(jecType))
  {
    std::cerr << "ERROR: '" << jecType << "' not found for m=" << mass << " in efficiencies JSON file. Exiting." << std::endl;
    exit(1);
  }

  Json::Value effNode = massNode[jecType];

  if (! effNode.isArray() || effNode.size() != 8)
  {
    std::cerr << "ERROR: malformated JSON file. Exiting." << std::endl;
    exit(1);
  }

  a = effNode[0].asDouble();
  b = effNode[1].asDouble();
  c = effNode[2].asDouble();
  d = effNode[3].asDouble();
  e = effNode[4].asDouble();
  f = effNode[5].asDouble();
  g = effNode[6].asDouble();
  h = effNode[7].asDouble();

  std::cout << "Efficiencies loaded successfully for " << mass << ":" << jecType << std::endl;
}

void loadSigmaRef(int mass, int btag, double& sigma)
{

  std::string path = (BASE_PATH + "sigma_reference.json");
  Json::Reader reader;
  Json::Value root;
  std::ifstream file(path.c_str());
  bool success = reader.parse(file, root);
  file.close();
  if (! success)
  {
    //std::cerr << "ERROR: Failed to parse '" << path << "'. Exiting." << std::endl;
    //std::cerr << "Are you sure you have run me without --scan first?" << std::endl;
    //exit(1);
    std::cout << "WARNING: No sigma_ref found. Setting to 0" << std::endl;
    return;
  }

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  root = root[getAnalysisUUID(BASE_PATH)];

  if (! root.isMember(strMass))
  {
    //std::cerr << "ERROR: mass '" << mass << "' not found in JSON file. Exiting." << std::endl;
    //std::cerr << "Are you sure you have run me without --scan first?" << std::endl;
    //exit(1);

    std::cout << "WARNING: No sigma_ref found. Setting to 0" << std::endl;
    return;
  }

  Json::Value massNode = root[strMass][btagStr];
  if (! massNode.isMember("sigma"))
  {
    //std::cerr << "ERROR: malformated JSON file. Exiting." << std::endl;
    //std::cerr << "Are you sure you have run me without --scan first?" << std::endl;
    //exit(1);

    std::cout << "WARNING: No sigma_ref found. Setting to 0" << std::endl;
    return;
  }

  sigma = fabs(massNode["sigma"].asDouble());
}

void loadSystematics(int mass, int btag, double& jec, double& pdf, double& pdf_cb)
{

  std::string path = (BASE_PATH + "systematics.json");
  Json::Reader reader;
  Json::Value root;
  std::ifstream file(path.c_str());
  bool success = reader.parse(file, root);
  file.close();
  if (! success)
  {
    //std::cerr << "ERROR: Failed to parse '" << path << "'. Exiting." << std::endl;
    //exit(1);
    std::cout << "WARNING: Systematics error not found. Setting to 0" << std::endl;
    return;
  }

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  root = root[getAnalysisUUID(BASE_PATH)];

  if (! root.isMember(strMass))
  {
    //std::cerr << "ERROR: mass '" << mass << "' not found in JSON file. Exiting." << std::endl;
    //exit(1);
    std::cout << "WARNING: Systematics error not found. Setting to 0" << std::endl;
    return;
  }

  Json::Value massNode = root[strMass][btagStr];
  if (! analysisUseJECSyst(BASE_PATH)) {
    jec = 0;
  } else {
    if (! massNode.isMember("jec"))
    {
      //std::cerr << "ERROR: malformated JSON file, JEC systematic not found. Exiting." << std::endl;
      //exit(1);
      std::cout << "WARNING: JEC systematics error not found. Setting to 0" << std::endl;
      jec = 0.;
    }
    else
    {
      jec = fabs(massNode["jec"].asDouble());
    }
  }

  if (! analysisUseBkgSyst(BASE_PATH)) {
    pdf = 0;
  } else {
    if (! massNode.isMember("background_pdf"))
    {
      //std::cerr << "ERROR: malformated JSON file, PDF background systematic not found. Exiting." << std::endl;
      //exit(1);
      std::cout << "WARNING: Background PDF systematic error not found. Setting to 0" << std::endl;
      pdf = 0;
    }
    else
    {
      pdf = fabs(massNode["background_pdf"].asDouble());
    }
  }

  if (! analysisUseSignalSyst(BASE_PATH)) {
    pdf_cb = 0;
  } else {
    if (! massNode.isMember("signal_pdf"))
    {
      //std::cerr << "ERROR: malformated JSON file, PDF CB systematic not found. Exiting." << std::endl;
      //std::cerr << "Maybe you just need to run <nom_utilitaire> first?" << std::endl;
      //exit(1);
      std::cout << "WARNING: Signal PDF systematic error not found. Setting to 0" << std::endl;
      pdf_cb = 0.;
    }
    else
    {
      pdf_cb = fabs(massNode["signal_pdf"].asDouble());
    }
  }
}

void renameAndSetPdfParametersConst(const RooArgSet& observables, const RooAbsPdf& pdf, const std::string& prefix)
{
  RooArgSet* params = pdf.getParameters(observables);
  TIterator* iter = params->createIterator();
  RooRealVar* var = NULL;
  while ((var = static_cast<RooRealVar*>(iter->Next())))
  {
    TString cleanedName = TString(var->GetName()).ReplaceAll("electron", "").ReplaceAll("muon", "");
    TString newName = prefix + cleanedName;

    var->setConstant(true);
    var->SetNameTitle(newName, newName);
  }

  delete iter;
  delete params;

  params = pdf.getComponents();
  iter = params->createIterator();
  var = NULL;
  while ((var = static_cast<RooRealVar*>(iter->Next())))
  {
    if ((void *) var == (const void*) &pdf)
      continue;

    TString cleanedName = TString(var->GetName()).ReplaceAll("electron", "").ReplaceAll("muon", "").ReplaceAll("signal_", "");
    TString newName = cleanedName + prefix;

    var->SetNameTitle(newName, newName);
  }

  delete iter;
  delete params;
}

void setPdfParametersConst(const RooArgSet& observables, const RooAbsPdf& pdf, bool constant)
{
  RooArgSet* params = pdf.getParameters(observables);
  TIterator* iter = params->createIterator();
  RooRealVar* var = NULL;
  while ((var = static_cast<RooRealVar*>(iter->Next())))
  {
    var->setConstant(constant);
  }

  delete iter;
  delete params;
}

void setPdfParametersRange(const RooArgSet& observables, const RooAbsPdf& pdf, double sigma)
{
  RooArgSet* params = pdf.getParameters(observables);
  TIterator* iter = params->createIterator();
  RooRealVar* var = NULL;
  while ((var = static_cast<RooRealVar*>(iter->Next())))
  {
    var->setRange(var->getVal() - sigma * var->getError(), var->getVal() + sigma * var->getError());
  }

  delete iter;
  delete params;
}

void mixEfficiencies(const std::map<int, double>& effs1, const std::map<int, double>& effs2, const std::string& prefix1, const std::string& prefix2, RooWorkspace& workspace) {
  for (const auto& eff1: effs1) {
    for (const auto& eff2: effs2) {
      if (prefix1 == prefix2 && eff1.first == eff2.first)
        continue;

      std::stringstream ss;
      ss << "eff_ratio_" << prefix1 << "_" << eff1.first << "b_" << prefix2 << "_" << eff2.first << "b[" << eff1.second / eff2.second << "]";

      std::cout << ss.str() << std::endl;

      workspace.factory(ss.str().c_str());
    }
  }
}

void parseConfigFile(const std::string& filename, /*RooAbsCategoryLValue& categories,*/RooWorkspace& workspace) {

  std::fstream configFile((BASE_PATH + "/" + filename).c_str(), std::ios::in);
  std::string line;
  while (std::getline(configFile, line)) {
    if (line[0] == '#' || line.length() == 0)
      continue;

    workspace.factory(line.c_str());
  }

  workspace.Print("v");
}

RooDataSet* binnedToUnbinnedDataset(RooDataHist& in, RooRealVar& observable, RooRealVar& weight, const std::string& name) {

  RooDataSet* ds = new RooDataSet(name.c_str(), name.c_str(), RooArgSet(observable, weight), RooFit::WeightVar(weight));
  for (int i = 0; i < in.numEntries(); i++) {
    const RooArgSet* set = in.get(i);
    ds->add(*set, in.weight(*set));
  }

  return ds;
}

RooAbsPdf* getInterpolatedPdf(RooRealVar& observable, double massZprime, const std::string& jec, int btag, const std::string& categoryName, const std::string& suffix = "") {

  // Interpolation
  int lowMass = 0;
  int highMass = 0;

  const int ALGO_MOMENT_MORPH = 1;
  const int ALGO_INTEGRAL_MORPH = 2;

  int algo;

  if (massZprime > 500 && massZprime < 750) {
    lowMass = 500;
    highMass = 750;
    algo = ALGO_INTEGRAL_MORPH;
  } else if (massZprime > 750 && massZprime < 1000) {
    lowMass = 750;
    highMass = 1000;
    algo = ALGO_MOMENT_MORPH;
  } else if (massZprime > 1000 && massZprime < 1250) {
    lowMass = 1000;
    highMass = 1250;
    algo = ALGO_MOMENT_MORPH;
  } else if (massZprime > 1250 && massZprime < 1500) {
    lowMass = 1250;
    highMass = 1500;
    algo = ALGO_MOMENT_MORPH;
  }/* else if (massZprime > 1500 && massZprime < 2000) {
      lowMass = 1500;
      highMass = 2000;
      }*/

  std::string cleanedCategory = TString(categoryName.c_str()).ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", "").Data(); // For workspace names

  algo = ALGO_INTEGRAL_MORPH;

  if (lowMass == 0 || highMass == 0) {
    std::cout << "Error: please use a mass between 500 and 1500 GeV." << std::endl;
    return NULL;
  }

  std::string analysisName = getAnalysisName(BASE_PATH);
  std::string pdfName = "signal_" + std::string((TString(categoryName.c_str()).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron"));
  std::string goodSuffix = (suffix.length() == 0) ? cleanedCategory : suffix;

  // Open workspaces and retrieve PDFs
  TString lowMass_prefix = TString::Format("%s-Zprime%d_%s_%d_btag", jec.c_str(), lowMass, analysisName.c_str(), btag);
  TString lowMass_workspaceFile = BASE_PATH + "/frit/" + lowMass_prefix + "_workspace.root";

  TString highMass_prefix = TString::Format("%s-Zprime%d_%s_%d_btag", jec.c_str(), highMass, analysisName.c_str(), btag);
  TString highMass_workspaceFile = BASE_PATH + "/frit/" + highMass_prefix + "_workspace.root";

  TFile lowMass_file(lowMass_workspaceFile);
  RooAbsPdf* lowMass_pdf = static_cast<RooWorkspace*>(lowMass_file.Get("w"))->pdf(pdfName.c_str());
  lowMass_pdf->SetName(TString::Format("lowMass_%s", pdfName.c_str()));
  lowMass_file.Close();

  TFile highMass_file(highMass_workspaceFile);
  RooAbsPdf* highMass_pdf = static_cast<RooWorkspace*>(highMass_file.Get("w"))->pdf(pdfName.c_str());
  highMass_pdf->SetName(TString::Format("highMass_%s", pdfName.c_str()));
  highMass_file.Close();

  double alpha = 0.;
  if (algo == ALGO_MOMENT_MORPH) {
    alpha = (double) (massZprime - lowMass) / (double) (highMass - lowMass);
  } else {
    alpha = 1. - (double) (massZprime - lowMass) / (double) (highMass - lowMass);
  }

  RooRealVar* rAlpha = new RooRealVar("alpha", "alpha", alpha, 0, 1);

  //mtt.setBins(5000, "cache");
  //rAlpha->setBins(1000, "cache");
  observable.setBins(5000, "cache");

  // Interpolate
  std::cout << "Interpolate ..." << std::endl;

  RooAbsPdf* interpolation;

  TString temporaryName = TString::Format("interpolation_%s", goodSuffix.c_str());

  if (algo == ALGO_INTEGRAL_MORPH) {

    std::cout << "Using RooIntergralMorph for interpolation" << std::endl;
    interpolation = new RooIntegralMorph(temporaryName, temporaryName, *lowMass_pdf, *highMass_pdf, observable, *rAlpha, false);

  } else {

    std::cout << "Using RooMomentMorph for interpolation" << std::endl;
    TVectorD hypoMass(2);
    hypoMass(0) = 0; 
    hypoMass(1) = 1;

    interpolation = new RooMomentMorph(temporaryName, temporaryName, *rAlpha, RooArgList(observable), RooArgList(*lowMass_pdf, *highMass_pdf), hypoMass, RooMomentMorph::Linear);
  }

  std::cout << "Done." << std::endl;

  int oldBinning = observable.getBins();
  observable.setBins(3000);

  RooDataHist * binnedInterpolatedDataset = new RooDataHist(std::string("binned_dataset_signal_" + goodSuffix).c_str(), "", RooArgSet(observable));
  interpolation->fillDataHist(binnedInterpolatedDataset, NULL, 1.);

  RooRealVar weight("__weight__", "__weight__", 0, 0, 100000);
  RooDataSet* unbinned_dataset = binnedToUnbinnedDataset(*binnedInterpolatedDataset, observable, weight, "test");

  RooKeysPdf* keys_pdf = new RooKeysPdf(std::string("signal_" + categoryName).c_str(), "Keys pdf for signal", observable, *unbinned_dataset, RooKeysPdf::MirrorBoth, 1);

  delete unbinned_dataset;

  RooPlot* p = observable.frame();
  interpolation->plotOn(p);
  keys_pdf->plotOn(p, RooFit::LineColor(kRed));

  TCanvas c("c");
  p->Draw();
  c.Print("test.root");

  //RooHistPdf* hist_pdf = new RooHistPdf(std::string("signal_" + categoryName).c_str(), std::string("signal_" + categoryName).c_str(), RooArgSet(observable), *binnedInterpolatedDataset);

  //interpolation->SetName(std::string("signal_" + categoryName).c_str());
  //renameAndSetPdfParametersConst(RooArgSet(observable), *interpolation, (suffix.length() == 0) ? cleanedCategory : suffix);
  renameAndSetPdfParametersConst(RooArgSet(observable), *keys_pdf, goodSuffix);

  observable.setBins(oldBinning);

  return keys_pdf;
}

void fitMtt(std::map<int, TChain*> eventChain, int massZprime, string fitConfigurationFile, int btag, const std::string& customWorkspaceFile, const std::string& outputFile)
{

  if (fitConfigurationFile == "auto" || fitConfigurationFile.empty())
    fitConfigurationFile = "fit_pdf_faltb.json";

  std::cout << "Loading fit configuration from '" << fitConfigurationFile << "'" << std::endl;

  const bool combine = (btag > 2);

  int minBTag = 0;
  int maxBTag = 0;
  std::string scriptNamePrefix;
  switch (btag) {
    case 3:
      minBTag = 1;
      maxBTag = 2;
      scriptNamePrefix = "1+2btag";
      break;

    case 4:
      minBTag = 0;
      maxBTag = 2;
      scriptNamePrefix = "0+1+2btag";
      break;

    case 5:
      minBTag = 0;
      maxBTag = 1;
      scriptNamePrefix = "0+1btag";
      break;

    default:
      minBTag = maxBTag = btag;
  }

  int nCombinedBTag = maxBTag - minBTag + 1;

  // configure root
  gROOT->Clear();
  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  setTDRStyle();

  // Set RooFit verbosity
  RooMsgService::instance().setStreamStatus(0, false);
  RooMsgService::instance().setStreamStatus(1, false);
  if (! VERBOSE) {
    RooMsgService::instance().addStream(RooFit::FATAL); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
    RooMsgService::instance().setSilentMode(true);
  } else {
    RooMsgService::instance().addStream(RooFit::INFO); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
  }

  //fit region
  Float_t minmTT = 550;
  Float_t maxmTT = 2000;
  
  // Create main workspace for global pdf
  RooWorkspace mainWorkspace("mainWorkspace", "main workspace");

  RooRealVar mtt("mtt", "m_{t#bar{t}}", minmTT, maxmTT, "GeV/c^2");
  Double_t minmTTFit = minmTT + 0.0;
  Double_t maxmTTFit = maxmTT - 0.0;
  mtt.setRange(minmTTFit, maxmTTFit);

  // Set binning to 1 GeV
  mtt.setBins((maxmTTFit - minmTTFit) / 4.);

  RooRealVar weight("weight", "weight", 0, 100000);
  mainWorkspace.import(mtt);
  mainWorkspace.import(weight);

  RooCategory lepton_type("lepton_type", "lepton_type");
  lepton_type.defineType("muon", 13); 
  lepton_type.defineType("electron", 11);
  mainWorkspace.import(lepton_type);

  RooCategory btagCategory("btag", "btag");

  if (minBTag == 0)
    btagCategory.defineType("0-btag");

  if (minBTag <= 1 && maxBTag >= 1)
    btagCategory.defineType("1-btag");

  if (minBTag <= 2 && maxBTag >= 2)
    btagCategory.defineType("2-btag");

  RooSuperCategory superCategory("superCat", "superCat", RooArgList(lepton_type, btagCategory));

  RooAbsCategoryLValue& mainCategory = combine ? static_cast<RooAbsCategoryLValue&>(superCategory) : static_cast<RooAbsCategoryLValue&>(lepton_type);

  std::string analysisName = getAnalysisName(BASE_PATH);

  std::map<std::string, std::shared_ptr<BaseFunction>> backgroundPdfs = getCategoriesPdf(BASE_PATH + "/fit_configuration", fitConfigurationFile, mtt, NULL, massZprime, "background", mainCategory, NULL);

  std::map<std::string, RooAbsPdf*> backgroundPdfsFromWorkspace;

  for (auto& pdf: backgroundPdfs) {
    std::cout << "Background pdf: " << pdf.first << " ";
    pdf.second->getPdf().Print();
    mainWorkspace.import(pdf.second->getPdf());

    backgroundPdfsFromWorkspace[pdf.first] = mainWorkspace.pdf(pdf.second->getPdf().GetName());
  }

  double lumi_mu = 19576.993;
  double lumi_e  = 19578.948;

  double s_lumi_mu_percent  = 4.4 / 100.; // 2.2%

  Double_t br_semil = 1.0; //0.14815; now included in the efficiency


  /*
   * Load efficiencies.
   * If we are working in combined mode, load for each b-tag categories.
   * Otherwise, load only the requested b-tag category
   */
  std::map<int, double> sel_eff_mu; // Selection efficiency in semi-mu channel for each b-tag
  std::map<int, double> sel_eff_e; // Selection efficiency in semi-e channel for each b-tag

  std::map<int, double> hlt_eff_mu; // HLT efficiency in semi-mu channel for each b-tag
  std::map<int, double> hlt_eff_e; // HLT efficiency in semi-e channel for each b-tag

  std::map<int, double> s_sel_eff_mu; // = Δ(nSig_mu) / nSig_mu ; relative error
  std::map<int, double> s_sel_eff_e;

  std::map<int, double> s_hlt_eff_mu; // ?
  std::map<int, double> s_hlt_eff_e; // ?

  double s_sys_JEC; // JEC systematic error
  double s_sys_PDF; // Signal PDF systematic error
  double s_sys_PDF_CB; // Background PDF systematic

  double sigma_ref = 0;; // Reference cross-section, obtained when fitting once
  loadSigmaRef(massZprime, btag, sigma_ref);

  loadSystematics(massZprime, btag, s_sys_JEC, s_sys_PDF, s_sys_PDF_CB);

  if (! combine) {
    loadEfficiencies(massZprime, "nominal", btag, sel_eff_mu[btag], sel_eff_e[btag], hlt_eff_mu[btag], hlt_eff_e[btag], s_sel_eff_mu[btag], s_sel_eff_e[btag], s_hlt_eff_mu[btag], s_hlt_eff_e[btag]);
  } else {
    // Load for 0, 1 and 2 btag
    for (int i = minBTag; i <= maxBTag; i++) {
      loadEfficiencies(massZprime, "nominal", i, sel_eff_mu[i], sel_eff_e[i], hlt_eff_mu[i], hlt_eff_e[i], s_sel_eff_mu[i], s_sel_eff_e[i], s_hlt_eff_mu[i], s_hlt_eff_e[i]);
    }
  }

  /*
     std::cout << "sel_eff_mu = " << sel_eff_mu << std::endl;
     std::cout << "sel_eff_e = " << sel_eff_e << std::endl;
     std::cout << "hlt_eff_mu = " << hlt_eff_mu << std::endl;
     std::cout << "hlt_eff_e = " << hlt_eff_e << std::endl;
     std::cout << "s_sel_eff_mu = " << s_sel_eff_mu << std::endl;
     std::cout << "s_sel_eff_e = " << s_sel_eff_e << std::endl;
     std::cout << "s_hlt_eff_mu = " << s_hlt_eff_mu << std::endl;
     std::cout << "s_hlt_eff_e = " << s_hlt_eff_e << std::endl;
     */

  std::cout << "Loading signal pdf..." << std::endl;

  std::map<std::string, RooAbsPdf*> signalPdfsFromWorkspace;

  // Read signal PDF from the workspace
  TIterator* it = mainCategory.typeIterator();
  RooCatType * type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    // Format correct filename
    
    int categoryBTag = btag;
    std::string category = type->GetName();

    boost::regex regex("([0-9]+)-btag");
    boost::smatch regexResults;

    if (boost::regex_search(category, regexResults, regex)) {
      std::string result(regexResults[1].first, regexResults[2].second);
      categoryBTag = atoi(result.c_str());
    }

    std::string pdfName = "signal_" + std::string((TString(type->GetName()).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron"));

    if (massZprime == 500 || massZprime == 750 || massZprime == 1000 || massZprime == 1250 || massZprime == 1500) {

      TString workspaceFile = TString::Format("%s/frit/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), "nominal", massZprime, analysisName.c_str(), categoryBTag);
      if (customWorkspaceFile.length() > 0) {
        workspaceFile = TString::Format(customWorkspaceFile.c_str(), categoryBTag);
        std::cout << Bash::set_color(Bash::Color::BLUE) << "Using custom workspace '" << workspaceFile << "'" << Bash::set_color() << std::endl;
      }

      std::shared_ptr<TFile> file(TFile::Open(workspaceFile));
      if (! file.get()) {
        std::cerr << "ERROR: MC signal parameters file not found (Was trying to open '" << workspaceFile << "')" << std::endl;
        exit(1);
      }

      RooWorkspace* workspace = static_cast<RooWorkspace*>(file->Get("w"));
      if (! workspace) {
        std::cerr << "ERROR: Workspace not found!" << std::endl;
        exit(1);
      }

      RooAbsPdf* pdf = workspace->pdf(pdfName.c_str());
      if (! pdf) {
        std::cerr << "ERROR: Signal pdf " << pdfName << " not found inside workspace." << std::endl;
        exit(1);
      }    

      pdf->SetName(std::string("signal_" + std::string(type->GetName())).c_str());
      renameAndSetPdfParametersConst(RooArgSet(mtt), *pdf, type->GetName());
      mainWorkspace.import(*pdf);

      signalPdfsFromWorkspace[category] = mainWorkspace.pdf(pdf->GetName());

    } else {

      RooAbsPdf* interpolation = getInterpolatedPdf(mtt, massZprime, "nominal", categoryBTag, type->GetName());
      mainWorkspace.import(*interpolation);

      signalPdfsFromWorkspace[category] = mainWorkspace.pdf(interpolation->GetName());
    }
  }

  std::cout << "Done." << std::endl;

  // Compute total systematic error
  std::map<int, double> eff_mu;
  std::map<int, double> eff_e;

  double total_efficiency = 0;
  
  if (! combine) {

    if (btag == 2) {
      eff_mu[2] = computeEfficiencyMuons_2btag(sel_eff_mu[2], hlt_eff_mu[2]);
      eff_e[2]  = computeEfficiencyElectrons_2btag(sel_eff_e[2], hlt_eff_e[2]);
      total_efficiency = eff_mu[2];
    } else if (btag == 1) {
      eff_mu[1] = computeEfficiencyMuons_1btag(sel_eff_mu[1], hlt_eff_mu[1]);
      eff_e[1]  = computeEfficiencyElectrons_1btag(sel_eff_e[1], hlt_eff_e[1]);
      total_efficiency = eff_mu[1];
    }

  } else {

    eff_mu[1] = computeEfficiencyMuons_1btag(sel_eff_mu[1], hlt_eff_mu[1]);
    eff_e[1]  = computeEfficiencyElectrons_1btag(sel_eff_e[1], hlt_eff_e[1]);

    eff_mu[2] = computeEfficiencyMuons_2btag(sel_eff_mu[2], hlt_eff_mu[2]);
    eff_e[2]  = computeEfficiencyElectrons_2btag(sel_eff_e[2], hlt_eff_e[2]);

    total_efficiency = eff_mu[2]; // Our parameters is nSig_mu for 2 btag. Use its efficiency for sigma computation
  }

  std::cout << std::endl;
  std::cout << "Analysis systematics and efficiencies" << std::endl;
  std::cout << "-------------------------------------" << std::endl << std::endl;

  std::map<int, double> s_eff_mu_percent;
  std::map<int, double> s_eff_mu_pb;

  std::map<int, double> s_eff_e_percent;
  std::map<int, double> s_eff_e_pb;

  std::map<int, double> s_yield_percent;
  std::map<int, double> s_yield_pb;

  double systematics_error_pb           = 0.;
  double systematics_error_events       = 0.;
  double systematics_error_percent      = 0.;

  // Compute systematics errors
  // Three part, lumi, b-tagging, efficiencies
  // See doc/systematics_errors.pdf for more details
  
  // Lumi
  const double lumi_systematic_relative = s_lumi_mu_percent;

  // B-tagging
  double b_tagging_systematic_relative = b_tagging_scale_factor_error / b_tagging_scale_factor;

  // Efficiency
  const double trigger_scale_factor_muons_relative = trigger_scale_factor_muons_error / trigger_scale_factor_muons;
  const double trigger_scale_factor_electrons_relative = trigger_scale_factor_electrons_error / trigger_scale_factor_electrons;

  const double muonID_scale_factor_relative = muonID_scale_factor_error / muonID_scale_factor;
  const double muonIso_scale_factor_relative = muonIso_scale_factor_error / muonIso_scale_factor;

  const double electron_scale_factor_relative = electron_scale_factor_error / electron_scale_factor;

  std::map<int, double> selection_systematic_relative_square_muons;
  std::map<int, double> selection_systematic_relative_square_electrons;

  std::map<int, double> selection_systematic_relative_muons;
  std::map<int, double> selection_systematic_relative_electrons;

  for (int i = minBTag; i <= maxBTag; i++) {

    selection_systematic_relative_square_muons[i] = trigger_scale_factor_muons_relative * trigger_scale_factor_muons_relative +
      muonID_scale_factor_relative * muonID_scale_factor_relative +
      muonIso_scale_factor_relative * muonIso_scale_factor_relative + 
      s_hlt_eff_mu[i] * s_hlt_eff_mu[i] +
      s_sel_eff_mu[i] * s_sel_eff_mu[i];

    selection_systematic_relative_square_electrons[i] = trigger_scale_factor_electrons_relative * trigger_scale_factor_electrons_relative +
      electron_scale_factor_relative * electron_scale_factor_relative +
      s_hlt_eff_e[i] * s_hlt_eff_e[i] +
      s_sel_eff_e[i] * s_sel_eff_e[i];

    if (i == 1) {

      selection_systematic_relative_square_muons[i] += pow((1 - 2 * b_tagging_efficiency * b_tagging_scale_factor) / (1 - b_tagging_efficiency * b_tagging_scale_factor), 2) * b_tagging_systematic_relative * b_tagging_systematic_relative;
      selection_systematic_relative_square_electrons[i] += pow((1 - 2 * b_tagging_efficiency * b_tagging_scale_factor) / (1 - b_tagging_efficiency * b_tagging_scale_factor), 2) * b_tagging_systematic_relative * b_tagging_systematic_relative;

    } else if (i == 2) {

      selection_systematic_relative_square_muons[i] += 2 * b_tagging_systematic_relative * b_tagging_systematic_relative;
      selection_systematic_relative_square_electrons[i] += 2 * b_tagging_systematic_relative * b_tagging_systematic_relative;

    }

    selection_systematic_relative_muons[i] = sqrt(selection_systematic_relative_square_muons[i]);
    selection_systematic_relative_electrons[i] = sqrt(selection_systematic_relative_square_electrons[i]);
  }

  std::cout << "Luminosity: " << Bash::set_color(Bash::Color::MAGENTA) << lumi_mu << " /pb" << Bash::set_color() << std::endl;
  std::cout << "Reference cross-section: " << Bash::set_color(Bash::Color::MAGENTA) << sigma_ref << " pb" << Bash::set_color() << std::endl;
  std::cout << "Total efficiency (selection * trigger * SFs): " << std::endl;

  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << "\tMuonic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << eff_mu[i] * 100 << " %" << Bash::set_color() << std::endl;
    std::cout << "\tElectronic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << eff_e[i] * 100 << " %" << Bash::set_color() << std::endl;
  }
  std::cout << "Expected number of signal events: " << std::endl;

  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << "\tMuonic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << lumi_mu * eff_mu[i] << Bash::set_color() << std::endl;
    std::cout << "\tElectronic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << lumi_e * eff_e[i] << Bash::set_color() << std::endl;
  }

  std::cout << "Efficiency systematic: " << std::endl;

  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << "\tMuonic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << selection_systematic_relative_muons[i] * 100 << " % (" << 1 + selection_systematic_relative_muons[i] << ")" << Bash::set_color() << std::endl;
    std::cout << "\tElectronic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << selection_systematic_relative_electrons[i] * 100 << " % (" << 1 + selection_systematic_relative_electrons[i] << ")" << Bash::set_color() << std::endl;
  }

  std::cout << "Efficiency systematic details: " << std::endl;

  int oldPre = std::cout.precision();
  std::cout.precision(3);
  std::cout << "Efficiency (\\%) & log normal";
  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << " & " << selection_systematic_relative_muons[i] * 100 << " & " << selection_systematic_relative_electrons[i] * 100;
  }
  std::cout << "\\\\" << std::endl << "\\hline" << std::endl;

  std::cout << "~~\\textbullet~Muon ID + iso &";
  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << " & " << sqrt(muonID_scale_factor_relative * muonID_scale_factor_relative + muonIso_scale_factor_relative * muonIso_scale_factor_relative) * 100  << " & -";
  }
  std::cout << "\\\\" << std::endl << "\\hline" << std::endl;
  std::cout << "~~\\textbullet~Electron ID + iso &";
  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << " & - & " << sqrt(electron_scale_factor_relative * electron_scale_factor_relative) * 100;
  }
  std::cout << "\\\\" << std::endl << "\\hline" << std::endl;
  std::cout << "~~\\textbullet~HLT scale factor &";
  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << " & " << trigger_scale_factor_muons_relative * 100 << " & " << trigger_scale_factor_electrons_relative * 100;
  }
  std::cout << "\\\\" << std::endl << "\\hline" << std::endl;
  std::cout << "~~\\textbullet~HLT efficiency &";
  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << " & " << s_hlt_eff_mu[i] * 100 << " & " << s_hlt_eff_e[i] * 100;
  }
  std::cout << "\\\\" << std::endl << "\\hline" << std::endl;
  std::cout << "~~\\textbullet~sel. efficiency &";
  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << " & " << s_sel_eff_mu[i] * 100<< " & " << s_sel_eff_e[i] * 100;
  }
  std::cout << "\\\\" << std::endl << "\\hline" << std::endl;
  std::cout << "~~\\textbullet~b-tagging scale factor &";
  for (int i = minBTag; i <= maxBTag; i++) {
    double a, b;
    if (i == 1) {

      a = pow((1 - 2 * b_tagging_efficiency * b_tagging_scale_factor) / (1 - b_tagging_efficiency * b_tagging_scale_factor), 2) * b_tagging_systematic_relative * b_tagging_systematic_relative;
      b = pow((1 - 2 * b_tagging_efficiency * b_tagging_scale_factor) / (1 - b_tagging_efficiency * b_tagging_scale_factor), 2) * b_tagging_systematic_relative * b_tagging_systematic_relative;

    } else if (i == 2) {

      a = 2 * b_tagging_systematic_relative * b_tagging_systematic_relative;
      b = 2 * b_tagging_systematic_relative * b_tagging_systematic_relative;

    }

    std::cout << " & " << sqrt(a) * 100<< " & " << sqrt(b) * 100;
  }
  std::cout << "\\\\" << std::endl << "\\hline" << std::endl;

  std::cout.precision(oldPre);

  std::cout << "Total systematic errors: ";

  std::cout << Bash::set_color(Bash::Color::MAGENTA);

  std::cout << systematics_error_events << " events ; " << systematics_error_pb << " pb ; " << systematics_error_percent * 100 << " %";

  std::cout << Bash::set_color();

  std::cout << std::endl << std::endl;

  /*
  // Output systematics error
  if (! combine) {
    std::cout << std::endl;

    double err_base = eff_mu[btag] * lumi_mu * br_semil;
    std::cout << "Efficiencies: " << std::endl;
    std::cout << " - Eff: " << total_efficiency * 100 << " %" << std::endl;
    std::cout << " - Lumi mu: " << lumi_mu << " /pb" << std::endl;
    std::cout << std::endl;
    std::cout << "Total syst errors: ";

    if (! useSystematics)
      std::cout << Bash::set_color(Bash::Color::RED);
      
    std::cout << err_sys_events << " events ; " << err_sys_pb << " pb ; " << err_sys_percent << " %";

    if (! useSystematics)
      std::cout << Bash::set_color();

    std::cout << std::endl;


    std::cout << " - Eff mu: " << err_base * s_eff_mu_pb[btag] << " events ; " << s_eff_mu_pb[btag] << " pb ; " << s_eff_mu_percent[btag] * 100 << " %" << std::endl;
    std::cout << " - Eff e: " << err_base * s_eff_e_pb[btag] << " events ; " << s_eff_e_pb[btag] << " pb ; " << s_eff_e_percent[btag] * 100 << " %" << std::endl;
    std::cout << " - Lumi mu: " << err_base * s_lumi_mu_pb << " events ; " << s_lumi_mu_pb << " pb ; " << s_lumi_mu_percent * 100 << " %" << std::endl;
    std::cout << " - B-tagging: " << err_base * b_tagging_corr_error_pb << " events ; " << b_tagging_corr_error_pb << " pb ; " << b_tagging_scale_factor_error * 100 << " %" << std::endl;
    std::cout << " - Yield: " << std::endl;
    std::cout << "    - JEC: " << err_base * sigma_ref * s_sys_JEC[btag] << " events ; " << sigma_ref * s_sys_JEC[btag] << " pb ; " << s_sys_JEC[btag]  * 100 << " %" << std::endl;
    std::cout << "    - Bkg PDF: " << err_base * sigma_ref * s_sys_PDF[btag] << " events ; " << sigma_ref * s_sys_PDF[btag] << " pb ; " << s_sys_PDF[btag] * 100 << " %" << std::endl;
    std::cout << "    - Signal PDF: " << err_base * sigma_ref * s_sys_PDF_CB[btag] << " events ; " << sigma_ref * s_sys_PDF_CB[btag] << " pb ; " << s_sys_PDF_CB[btag] * 100 << " %" << std::endl;
    std::cout << std::endl;
  } else {
    std::cout << std::endl << Bash::set_color(Bash::Color::RED) << "Systematics are not yet implemented for combined analysis !!" << Bash::set_color() << std::endl << std::endl;
    err_sys_events = err_sys_pb = 0;
  }
  */

  //exit(0);

  /*
  if (! muonsOnly && doLimitCurve) {
    // Try to correct the correct muonsOnly flag.
    TString file = TString::Format("%s/%s_fitRes_%s.root", BASE_PATH.c_str(), prefix.Data(), suffix.Data());
    TFile *fitFile = TFile::Open(file);
    if (! fitFile)
      fileFile = TFile::Open(TString::Format("%s/%s_fitRes_%s.root", OUTPUT_PATH.c_str(), prefix.Data(), suffix.Data()));

    RooRealVar *nEventsRed_e = static_cast<RooRealVar*>(fitFile->Get("nEventsRed_e"));

    if (! nEventsRed_e) {
      std::cout << Bash::set_color(Bash::GREEN) << "You ask me to do some toys for semi-mu and semi-e. However, \"nEventsRed_e\" is missing from \"" << file.Data() << "\"." << std::endl;
      std::cout << "Running toys only for semi-mu data." << Bash::set_color() << std::endl;
      muonsOnly = true;
    }
    fitFile->Close();
    delete fitFile;
  }

  if (muonsOnly) {
    std::cout << Bash::set_color(Bash::BLUE) << "Running using only semi-mu data." << Bash::set_color() << std::endl;
  }
  */

  //Background::Function* bkgFct = Background::createFunction(bkgfit_str, mtt);

  if (combine) {
    mixEfficiencies(eff_e, eff_mu, "e", "mu", mainWorkspace);
    mixEfficiencies(eff_e, eff_e, "e", "e", mainWorkspace);
    mixEfficiencies(eff_mu, eff_mu, "mu", "mu", mainWorkspace);
  } else {
    std::stringstream ss;
    ss << "eff_ratio_e_mu[" << eff_e[btag] / eff_mu[btag]<< "]";
    mainWorkspace.factory(ss.str().c_str());
  }

  RooRealVar lumiRatio("lumiRatio", "luminosity ratio", lumi_e / lumi_mu, "");
  mainWorkspace.import(lumiRatio);

  // Read config file for global pdf
  std::string configFile = (combine)
    ? TString::Format("combined_%s_global_pdf.script", scriptNamePrefix.c_str()).Data()
    : "individual_btag_global_pdf.script";

  parseConfigFile(configFile, /*mainCategory,*/ mainWorkspace);

  RooSimultaneous simPdfBackgroundOnly("simPdfBackgroundOnly", "simultaneous pdf with background only pdfs", mainCategory);

  it = mainCategory.typeIterator();
  type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    std::string name = type->GetName();

    std::string cleanedCategory = TString(type->GetName()).ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", "").Data(); // For workspace names

    std::string workspaceName = "global_pdf_" + cleanedCategory;
    simPdfBackgroundOnly.addPdf(*backgroundPdfsFromWorkspace[name], name.c_str());
  }

  RooDataSet *dataOrig = NULL;
  if (combine) {
    // Combine b-tags
    std::cout << "Combining b-tag dataset" << std::endl;

    std::map<int, RooDataSet*> datasets;

    if (minBTag == 0)
      datasets[0] = new RooDataSet("dataset", "dataset", RooArgSet(mtt, lepton_type, weight), Import(*(eventChain[0])), WeightVar(weight));

    if (minBTag <= 1 && maxBTag >= 1)
      datasets[1] = new RooDataSet("dataset", "dataset", RooArgSet(mtt, lepton_type, weight), Import(*(eventChain[1])), WeightVar(weight));

    if (minBTag <= 2 && maxBTag >= 2)
      datasets[2] = new RooDataSet("dataset", "dataset", RooArgSet(mtt, lepton_type, weight), Import(*(eventChain[2])), WeightVar(weight));

    switch (btag) {
      case 3:
          dataOrig = new RooDataSet("combData", "combined data", RooArgSet(mtt, lepton_type, weight), Index(btagCategory), Import("1-btag", *datasets[1]), Import("2-btag", *datasets[2]), WeightVar(weight));
        break;

      case 4:
          dataOrig = new RooDataSet("combData", "combined data", RooArgSet(mtt, lepton_type/*, weight*/), Index(btagCategory), Import("0-btag", *datasets[0]), Import("1-btag", *datasets[1]), Import("2-btag", *datasets[2])/*, WeightVar(weight)*/);
        break;

      case 5:
          dataOrig = new RooDataSet("combData", "combined data", RooArgSet(mtt, lepton_type/*, weight*/), Index(btagCategory), Import("0-btag", *datasets[0]), Import("1-btag", *datasets[1])/*, WeightVar(weight)*/);
        break;
    }

    dataOrig->Print("");
    dataOrig->table(superCategory)->Print("v");

  } else {
    // Dataset is inside eventChain. Load RooDataSet from tree
    dataOrig = new RooDataSet("dataset", "dataset", RooArgSet(mtt, lepton_type, weight), Import(*(eventChain[btag])), WeightVar(weight));
    dataOrig->table(lepton_type)->Print("v");
  }

  // Create a binned dataset
  std::shared_ptr<RooDataHist> binnedDataset = std::shared_ptr<RooDataHist>(dataOrig->binnedClone());

  std::shared_ptr<RooAbsData> datasetToFit = binnedDataset; // Binned likelihood
  //std::shared_ptr<RooAbsData> datasetToFit(dataOrig); // Unbinned likelihood

  std::cout << "Dataset entries: " << dataOrig->numEntries() << std::endl;

  // First, fit with background only pdfs
  RooFitResult* fitResult = simPdfBackgroundOnly.fitTo(*datasetToFit, Save(), Optimize(0));
  fitResult->Print("v");
  delete fitResult;

  RooWorkspace higgsWorkspace("w");

  it = mainCategory.typeIterator();
  type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {

    std::string category = type->GetName();
    std::string cleanedCategory = TString(type->GetName()).ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", "").Data(); // For workspace names

    mainCategory = category.c_str();
    std::string cut = buildCutFormula(mainCategory);
    std::cout << "Cut: " << cut << std::endl;

    std::string workspaceName = "global_pdf_" + cleanedCategory;

    std::string leptonName = TString(category).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron";
    std::string leptonNameShort = TString(category).Contains("muon", TString::kIgnoreCase) ? "mu" : "e";
    static boost::regex regex("([0-9]+)-btag");
    boost::smatch regexResults;

    int extractedBTag = -1;
    if (boost::regex_search(category, regexResults, regex)) {
      std::string result(regexResults[1].first, regexResults[2].second);
      extractedBTag = atoi(result.c_str());
    }
    if (extractedBTag < 0)
      extractedBTag = btag;

    TString workspace_suffix = TString::Format("%s_%db", leptonNameShort.c_str(), extractedBTag);
    
    RooAbsData* dataset_reduced = datasetToFit->reduce(RooArgSet(mtt), cut.c_str());
    higgsWorkspace.import(*dataset_reduced, RooFit::Rename(TString::Format("data_obs_%s", workspace_suffix.Data())));

    // Import background
    TString name = TString::Format("background_%s", category.c_str());
    RooAbsPdf* pdf = mainWorkspace.pdf(name);
    setPdfParametersRange(RooArgSet(mtt), *pdf, 10);
    higgsWorkspace.import(
        *pdf,
        RooFit::RenameVariable(name, TString::Format("background_%s", workspace_suffix.Data()))
        );

    // Import signal
    name = TString::Format("signal_%s", category.c_str());
    pdf = mainWorkspace.pdf(name);
    setPdfParametersRange(RooArgSet(mtt), *pdf, 10);
    higgsWorkspace.import(
        *pdf,
        RooFit::RenameVariable(name, TString::Format("signal_%s", workspace_suffix.Data()))
        );

    if (massZprime == 500 || massZprime == 750 || massZprime == 1000 || massZprime == 1250 || massZprime == 1500) {

      name = TString::Format("signal_%s", leptonName.c_str());

      // Import JEC up
      TString workspaceFile = TString::Format("%s/frit/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), "JECup", massZprime, analysisName.c_str(), extractedBTag);
      std::shared_ptr<TFile> f(TFile::Open(workspaceFile.Data()));
      pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(name);
      pdf->SetName(TString::Format("signal_%s_jecUp", workspace_suffix.Data()));
      higgsWorkspace.import(*pdf);
      
      // Import JEC down
      workspaceFile = TString::Format("%s/frit/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), "JECdown", massZprime, analysisName.c_str(), extractedBTag);
      f.reset(TFile::Open(workspaceFile.Data()));
      pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(name);
      pdf->SetName(TString::Format("signal_%s_jecDown", workspace_suffix.Data()));
      higgsWorkspace.import(*pdf);

      // Import PU up
      workspaceFile = TString::Format("%s/frit/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), "puUp", massZprime, analysisName.c_str(), extractedBTag);
      f.reset(TFile::Open(workspaceFile.Data()));
      pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(name);
      pdf->SetName(TString::Format("signal_%s_puUp", workspace_suffix.Data()));
      higgsWorkspace.import(*pdf);

      // Import PU down
      workspaceFile = TString::Format("%s/frit/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), "puDown", massZprime, analysisName.c_str(), extractedBTag);
      f.reset(TFile::Open(workspaceFile.Data()));
      pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(name);
      pdf->SetName(TString::Format("signal_%s_puDown", workspace_suffix.Data()));
      higgsWorkspace.import(*pdf);

      // Import JER up
      workspaceFile = TString::Format("%s/frit/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), "JERup", massZprime, analysisName.c_str(), extractedBTag);
      f.reset(TFile::Open(workspaceFile.Data()));
      pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(name);
      pdf->SetName(TString::Format("signal_%s_jerUp", workspace_suffix.Data()));
      higgsWorkspace.import(*pdf);

      // Import JER down
      workspaceFile = TString::Format("%s/frit/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), "JERdown", massZprime, analysisName.c_str(), extractedBTag);
      f.reset(TFile::Open(workspaceFile.Data()));
      pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(name);
      pdf->SetName(TString::Format("signal_%s_jerDown", workspace_suffix.Data()));
      higgsWorkspace.import(*pdf);
      

    } else {

      // Interpolation
      name = TString::Format("signal_%s_jecUp", workspace_suffix.Data());
      pdf = getInterpolatedPdf(mtt, massZprime, "JECup", extractedBTag, category, name.Data());
      pdf->SetName(name);
      higgsWorkspace.import(*pdf);

      name = TString::Format("signal_%s_jecDown", workspace_suffix.Data());
      pdf = getInterpolatedPdf(mtt, massZprime, "JECdown", extractedBTag, category, name.Data());
      pdf->SetName(name);
      higgsWorkspace.import(*pdf);
    }
  }

  TString outputFileName = TString::Format("zprime_%d_workspace.root", massZprime);
  if (outputFile.length() > 0)
    outputFileName = outputFile;

  higgsWorkspace.writeToFile(outputFileName);
}

