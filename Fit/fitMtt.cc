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
#include <RooIntegralMorph.h>
#include <RooMomentMorph.h>
#include <RooHistPdf.h>
#include <RooConstVar.h>

#include "Utils.h"
#include <tclap/CmdLine.h>
#include <json/json.h>

#include "SignalFunctions.h"
#include "Functions.h"

#define SAFE_DELETE(p) { delete p; p = NULL; }

using namespace RooFit;

/**
 * Profiling method
 */
double get_time_usec()
{
  static struct timeval _t;
  static struct timezone tz;
  gettimeofday(&_t, &tz);
  return (double) _t.tv_sec + (double)_t.tv_usec / (1000 * 1000);
}

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

void fitMtt(std::map<int, TChain*> eventChain, int massParticle, bool fit, string fitConfigurationFile, bool writeRootFile, bool saveFigures, bool doBiasTest, int nToyExp, int index, string syst_str, bool bkgOnly, bool muonsOnly, int btag, const std::string& customWorkspaceFile, bool fixBackground, bool saveWorkspace, bool createCombineWorkspace, const std::string& combineWorkspaceFilename);

std::string BASE_PATH;
std::string OUTPUT_PATH;
std::string EFF_FILE;
bool VERBOSE = false;
bool BATCH_MODE = false;

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

    TCLAP::ValueArg<int> massArg("m", "mass", "Particle mass", true, 750, "integer");
    TCLAP::SwitchArg fitArg("", "no-fit", "Don't do the reference fit", true);
    TCLAP::SwitchArg writeRootArg("", "no-root-files", "Don't write root files", true);
    TCLAP::SwitchArg saveFiguresArg("", "no-figs", "Don't save figures", true);
    TCLAP::SwitchArg doBiasTestArg("", "bias-test", "Do the bias test", false);
    TCLAP::ValueArg<int> nToyArg("", "toys", "Number of toys exp.", false, 1, "integer");
    TCLAP::ValueArg<int> indexArg("", "index", "Index", false, 1, "integer");

    std::vector<std::string> jec;
    jec.push_back("nominal");
    jec.push_back("JECup");
    jec.push_back("JECdown");
    TCLAP::ValuesConstraint<std::string> allowedVals(jec);
    TCLAP::ValueArg<std::string> systArg("", "syst", "Systematic", false, "nominal", &allowedVals);

    cmd.add(massArg);
    cmd.add(fitArg);
    cmd.add(writeRootArg);
    cmd.add(saveFiguresArg);
    cmd.add(doBiasTestArg);
    cmd.add(nToyArg);
    cmd.add(indexArg);
    cmd.add(systArg);

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    std::string default_base_path = "./analysis/" + getAnalysisUUID() + "/";
    TCLAP::ValueArg<std::string> pathArg("", "path", "Folder where to load files", false, default_base_path, "string", cmd);
    TCLAP::ValueArg<std::string> outputPathArg("", "output-path", "Folder where output files are stored", false, default_base_path, "string", cmd);
    TCLAP::ValueArg<std::string> effArg("", "eff-file", "File where efficiences are stored (JSON format)", false, "efficiencies.json", "string", cmd);
    TCLAP::SwitchArg bkgOnlyArg("", "bkg-only", "Fit using background model only", cmd);
    TCLAP::SwitchArg onlyMuonArg("", "muons-only", "Compute limits using only semi-mu data", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);
    TCLAP::SwitchArg verboseArg("v", "verbose", "Verbose mode", cmd);
    TCLAP::SwitchArg batchArg("", "batch", "Run in batch mode", cmd);
    TCLAP::ValueArg<std::string> fitConfigFileArg("", "config-file", "Configuration file name containing fit parameters", false, "fit_pdf_falt.json", "string", cmd);

    // Workspace
    TCLAP::ValueArg<std::string> workspaceArg("", "workspace", "Use a custom workspace containing the signal pdf to use", false, "", "string", cmd);

    // Fix background
    TCLAP::SwitchArg fixBackgroundArg("", "fix-background", "Fit once with background only, then fix the background and refit background + signal", cmd);

    TCLAP::SwitchArg saveWorkspaceArg("", "save-workspace", "Save the workspace for redoing plot after", cmd);

    TCLAP::SwitchArg createCombineWorkspaceArg("", "combine", "Create a workspace compatible with combine for limits computation", cmd);
    TCLAP::ValueArg<std::string> combineWorkspaceFilenameArg("o", "output-file", "The output filename of the workspace to produce", false, "", "string", cmd);

    cmd.parse(argc, argv);

    BASE_PATH = pathArg.getValue();
    OUTPUT_PATH = outputPathArg.getValue();

    BASE_PATH.erase(std::remove(BASE_PATH.begin(), BASE_PATH.end(), '\"'), BASE_PATH.end());
    OUTPUT_PATH.erase(std::remove(OUTPUT_PATH.begin(), OUTPUT_PATH.end(), '\"'), OUTPUT_PATH.end());

    EFF_FILE = effArg.getValue();

    VERBOSE = verboseArg.getValue();
    BATCH_MODE = batchArg.getValue();

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

    fitMtt(chains, massArg.getValue(), fitArg.getValue(), fitConfigFileArg.getValue(), writeRootArg.getValue(), saveFiguresArg.getValue(), doBiasTestArg.getValue(), nToyArg.getValue(), indexArg.getValue(), systArg.getValue(), bkgOnlyArg.getValue(), onlyMuonArg.getValue(), btagArg.getValue(), workspaceArg.getValue(), fixBackgroundArg.getValue(), saveWorkspaceArg.getValue(), createCombineWorkspaceArg.getValue(), combineWorkspaceFilenameArg.getValue());

    for (auto& chain: chains)
      delete chain.second;

  }
  catch (TCLAP::ArgException& e)
  {
    std::cerr << e.error() << std::endl;
  }
}

Json::Value getEfficiencyJSONNode(int mass, const std::string& syst, int btag) {

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

  if (! massNode.isMember(syst))
  {
    std::cerr << "ERROR: '" << syst << "' not found for m=" << mass << " in efficiencies JSON file. Exiting." << std::endl;
    exit(1);
  }

  Json::Value effNode = massNode[syst];

  if (! effNode.isObject())
  {
    std::cerr << "ERROR: malformated JSON file. Exiting." << std::endl;
    exit(1);
  }

  return effNode;
}

void loadSelectionEfficiencies(int mass, const std::string& syst, int btag, double& eff_mu, double& eff_e, double& error_eff_mu, double& error_eff_e, bool silent = false)
{

  Json::Value effNode = getEfficiencyJSONNode(mass, syst, btag);

  eff_mu = effNode["eff_mu"].asDouble();
  eff_e = effNode["eff_e"].asDouble();
  error_eff_mu = effNode["error_eff_mu"].asDouble();
  error_eff_e = effNode["error_eff_e"].asDouble();

  if (! silent)
    std::cout << "Selection efficiencies loaded successfully for " << mass << " and systematic " << syst << std::endl;
}

void loadTriggerEfficiencies(int mass, const std::string& syst, int btag, double& eff_mu, double& eff_e, double& error_eff_mu, double& error_eff_e, bool silent = false)
{

  Json::Value effNode = getEfficiencyJSONNode(mass, syst, btag);

  eff_mu = effNode["trigger_eff_mu"].asDouble();
  eff_e = effNode["trigger_eff_e"].asDouble();
  error_eff_mu = effNode["trigger_error_eff_mu"].asDouble();
  error_eff_e = effNode["trigger_error_eff_e"].asDouble();

  if (! silent)
    std::cout << "Trigger efficiencies loaded successfully for " << mass << " and systematic " << syst << std::endl;
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

void drawHistograms(RooAbsCategoryLValue& categories, RooRealVar& observable, RooAbsData& dataset, RooSimultaneous& simPdfs, std::map<std::string, RooAbsPdf*>& backgroundPdfs, int btag, bool savePlots, const std::string& prefix, const std::string& suffix, bool drawSignal, bool logToo, TFile* outputFile, bool drawOnlyData = false) {

  if (outputFile == nullptr && ! savePlots)
    return;

  std::vector<std::shared_ptr<RooPlot>> plots;
  std::vector<std::shared_ptr<TPad>> pads;

  int n = categories.numTypes();
  int x = std::min(n, 2), y = (int) ceil((float) n / (float) x);
  std::cout << n << " categories. Dividing into " << x << "; " << y << std::endl;

  //x = 4; y = 1;

  const float resolution = 25.0;
  const int nBinsForHisto = (observable.getMax() - observable.getMin() + 0.5) / resolution;

  const int padWidth = 900;
  const int padHeight = 900;
  const float LUMI = 19.7;

  const int canvasWidth = padWidth * x;
  const int canvasHeight = padHeight * y;

  TCanvas *canvas = new TCanvas("canvas", "mTT fit", canvasWidth, canvasHeight);
  canvas->Divide(x, y);
  int currentPad = 1;

  Roo1DTable* table = nullptr;
  if (outputFile) {
    table = dataset.table(categories);
  }

  TIterator* it = categories.typeIterator();
  RooCatType* type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    canvas->cd(currentPad);

    RooPlot* plot = observable.frame(nBinsForHisto);
    RooPlot* zoom_plot = observable.frame(nBinsForHisto);

    plot->SetTitle(""); //FIXME
    zoom_plot->SetTitle("");

    std::string category = type->GetName();
    std::string cleanedCategory = TString(category).ReplaceAll(";", "_").Data(); // Root does not like ';' in names
    TString signalCategory = TString(type->GetName()).ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", ""); // For workspace names

    categories = category.c_str();

    std::string cut = buildCutFormula(categories);

    RooDataSet* subData = static_cast<RooDataSet*>(dataset.reduce(cut.c_str()));
    subData->plotOn(plot);
    subData->plotOn(zoom_plot);

    if (! drawOnlyData) {

      if (drawSignal) {
        simPdfs.plotOn(plot, Slice(categories), ProjWData(*subData), LineColor(kBlue), LineWidth(1), Range("FULL"));
        simPdfs.plotOn(zoom_plot, Slice(categories), ProjWData(*subData), LineColor(kBlue), LineWidth(1), Range("FULL"));
      }

      simPdfs.plotOn(plot, Slice(categories), ProjWData(*subData), Components(*backgroundPdfs[category]), LineStyle(kDashed), LineColor(kRed), LineWidth(1), Range("FULL"));
      simPdfs.plotOn(plot, Slice(categories), ProjWData(*subData), Components(TString::Format("positive_signal_%s", signalCategory.Data())), LineStyle(kDashed), LineColor(kRed), LineWidth(1));
      simPdfs.plotOn(zoom_plot, Slice(categories), ProjWData(*subData), Components(*backgroundPdfs[category]), LineStyle(kDashed), LineColor(kRed), LineWidth(2), Range("FULL"));
    }

    delete subData;

    std::string leptonName = TString(category).Contains("muon", TString::kIgnoreCase) ? "muons" : "electrons";
    float binningSize = (observable.getBinning().highBound() - observable.getBinning().lowBound()) / (float) nBinsForHisto;

    plot->SetXTitle(TString::Format("#font[132]{#font[12]{M_{t#bar{t}}} (GeV), %s}", leptonName.c_str()));
    plot->SetYTitle(TString::Format("#font[132]{Events / (%0.0f GeV)}", binningSize));
    plot->SetTitleOffset(1.42, "Y");

    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.04);

    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.050);

    plot->Draw();

    // Find number of b-tag
    boost::regex regex("([0-9]+)-btag");
    boost::smatch regexResults;

    if (boost::regex_search(category, regexResults, regex)) {
      std::string result(regexResults[1].first, regexResults[2].second);
      btag = atoi(result.c_str());
    }

    TString padName = TString::Format("pad_zoom_%d", currentPad);

    double x1 = 0.18, y1 = 0.14; 

    /*
    if (btag == 2) {
      x1 += 0.337;
      x2 += 0.337;
      y1 += 0.42;
      y2 += 0.42;
    }
    */

    /*
    TPad* pad = new TPad(padName, "zoom", x1, y1, x2, y2);
    pad->Draw();
    pad->cd();
    pad->SetLogy(true);
    pad->SetRightMargin(0);

    zoom_plot->SetAxisRange(550, 900, "X");
    zoom_plot->SetMinimum(1000);
    zoom_plot->SetTitleSize(0, "X");
    zoom_plot->SetTitleSize(0, "Y");
    zoom_plot->Draw();
    */

    canvas->cd(currentPad);


    TString btagLabel = "";
    if (btag == 2)
      btagLabel = "#geq 2 b-tags";
    else
      btagLabel = TString::Format("%d b-tag", btag);

    std::string leptonShortcutName = TString(category).Contains("muon", TString::kIgnoreCase) ? "#mu" : "e";

    TString legendLabel = TString::Format("#font[42]{%s, #geq 4 jets, %s}", leptonShortcutName.c_str(), btagLabel.Data());

    x1 = 0.53;
    y1 = 0.88;

    /*
    if (btag == 2) {
      x1 -= 0.31;
      y1 -= 0.60;
    }
    */

    t.DrawLatex(x1, y1, "#font[42]{CMS preliminary}");
    t.DrawLatex(x1, y1 - 0.04, TString::Format("#font[42]{%0.2f fb^{-1} at #sqrt{s}=8 TeV}", LUMI));
    t.DrawLatex(x1, y1 - 2 * 0.04, legendLabel);

    if (outputFile) {
      if (drawSignal && !drawOnlyData) {
        RooHist* residual = plot->residHist(plot->nameOf(0), plot->nameOf(2));
        RooHist* pull     = plot->pullHist(plot->nameOf(0), plot->nameOf(2));

        RooPlot* residual_frame = observable.frame();
        residual_frame->SetXTitle(TString::Format("#font[132]{#font[12]{m_{TT}} (GeV/#font[12]{c}^{2}), %s}", leptonName.c_str()));
        residual_frame->SetYTitle("#font[132]{Data-Fit}");
        residual_frame->addPlotable(dynamic_cast<RooPlotable*>(residual), "P");

        RooPlot* pull_frame = observable.frame();
        pull_frame->SetXTitle(TString::Format("#font[132]{#font[12]{m_{TT}} (GeV/#font[12]{c}^{2}), %s}", leptonName.c_str()));
        pull_frame->SetYTitle("#font[132]{Pull}");
        pull_frame->addPlotable(dynamic_cast<RooPlotable*>(pull), "P");

        outputFile->cd();
        residual_frame->Write(TString::Format("%s_residual", cleanedCategory.c_str()));
        pull_frame->Write(TString::Format("%s_pull", cleanedCategory.c_str()));

        delete residual_frame;
        delete pull_frame;
      }

      outputFile->cd();
      plot->Write(TString::Format("%s_plot", cleanedCategory.c_str()));

      RooRealVar nEvents("nEvents", "number of events", table->get(category.c_str()));
      nEvents.Write(TString::Format("nEvents_%s", cleanedCategory.c_str()));
    }

    plots.push_back(std::shared_ptr<RooPlot>(plot));

    currentPad++;
  }

  delete table;

  if (savePlots) {
    canvas->Print((OUTPUT_PATH + prefix + "_fitRes_" + suffix + ".pdf").c_str());
    canvas->Print((OUTPUT_PATH + prefix + "_fitRes_" + suffix + ".png").c_str());
  }

  if (!logToo) {
    delete canvas;
    return;
  }

  currentPad = 1;
  it = categories.typeIterator();
  type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    canvas->cd(currentPad++);

    gPad->SetLogy(true);
    canvas->SetLogy(true);
    gPad->Modified();
    gPad->Update();
  }

  if (savePlots) {
    canvas->Print((OUTPUT_PATH + prefix + "_fitRes_" + suffix + "_log.pdf").c_str());
    canvas->Print((OUTPUT_PATH + prefix + "_fitRes_" + suffix + "_log.png").c_str());
  }

  delete canvas;
}

std::map<std::string, float> computeChi2(RooRealVar& observable, const RooSimultaneous& simPdf, RooAbsCategoryLValue& categories, RooAbsData& dataset, RooWorkspace& workspace, bool bkgOnly = false) {

  std::map<std::string, float> results;

  const float resolution = 4.;
  const int nBinsForChi2 = (observable.getMax() - observable.getMin() + 0.5) / resolution;
  const int oldBinning = observable.getBins();
  observable.setBins(nBinsForChi2);
  std::cout << "Binning dataset with " << nBinsForChi2 << " bins for chi2 computation (" << observable.getMin() << " -> " << observable.getMax() << " ; resolution: " << resolution << " GeV)" << std::endl;

  RooArgSet* floatingParameters = static_cast<RooArgSet*>(simPdf.getParameters(RooArgSet(observable))->selectByAttrib("Constant", false));
  int numberOfFloatingParams = floatingParameters->getSize() - 2;
  delete floatingParameters;

  Roo1DTable* table = dataset.table(categories);

  // Chi2 computation inspired from http://cmssdt.cern.ch/SDT/doxygen/CMSSW_5_1_1/doc/html/d4/d3a/GenericTnPFitter_8h_source.html, line 291

  const int numberOfCategories = categories.numTypes();
  float combinedChi2 = 0.;
  TIterator* it = categories.typeIterator();
  RooCatType* type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    std::cout << "For category '" << type->GetName() << "':" << std::endl;

    std::string cleanedCategory = TString(type->GetName())
      .ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", "").Data(); // For workspace names

    std::string workspaceName = "global_pdf_" + cleanedCategory;
    RooAddPdf* pdf = static_cast<RooAddPdf*>(workspace.pdf(workspaceName.c_str()));

    if (! bkgOnly) {
      const RooArgList& coefs = pdf->coefList();

      std::cout << "\tNumber of signal events: " << ((RooAbsReal&) coefs[0]).getVal() << " out of " << table->get(type->GetName()) << std::endl;
      std::cout << "\tNumber of background events: " << ((RooAbsReal&) coefs[1]).getVal() << " out of " << table->get(type->GetName()) << std::endl;
    }

    categories = type->GetName();
    std::string cut = buildCutFormula(categories);

    std::shared_ptr<RooAbsData> reducedDataset(dataset.reduce(Cut(cut.c_str())));
    std::shared_ptr<RooDataHist> binnedDataset(new RooDataHist("binnedDataset", "binned dataset", RooArgSet(observable), *reducedDataset));

    auto errorType = RooAbsData::Poisson;
    if (getAnalysisType() == HIGGS)
      errorType = RooAbsData::SumW2;

    float chi2 = RooChi2Var("chi2", "chi2", *pdf, *binnedDataset, DataError(errorType)).getVal();
    int NDF = (observable.getBins() - numberOfFloatingParams / numberOfCategories);
    float chi2NDF = (chi2 / NDF);
    std::cout << "\tChi2: " << chi2 << " / NDF: " << NDF << std::endl;
    std::cout << "\tChi2/NDF: " << chi2NDF << std::endl;

    results[type->GetName()] = chi2NDF;

    combinedChi2 += chi2;
  }

  std::cout << std::endl;
  delete table;

  combinedChi2 /= (numberOfCategories * observable.getBins() - numberOfFloatingParams);
  std::cout << "Combined Chi2: " << combinedChi2 << std::endl;
  results["combined"] = combinedChi2;

  observable.setBins(oldBinning);

  return results;
}

void addEfficienciesToWorkspace(
    const std::map<int, double>& effs1,
    const std::map<int, double>& effs2,
    const std::map<int, double>& err_effs1,
    const std::map<int, double>& err_effs2,
    const std::string& prefix1,
    const std::string& prefix2,
    RooWorkspace& workspace) {

  for (const auto& err_eff1: err_effs1) {

    std::stringstream ss;

    ss.str(std::string());
    ss << "error_eff_" << prefix1 << "_" << err_eff1.first << "b[" << err_eff1.second << "]";
    workspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "error_eff_" << prefix2 << "_" << err_eff1.first << "b[" << err_effs2.at(err_eff1.first) << "]";
    workspace.factory(ss.str().c_str());
  }

  for (const auto& eff1: effs1) {

    std::stringstream ss;

    ss.str(std::string());
    ss << "eff_" << prefix1 << "_" << eff1.first << "b[" << eff1.second << "]";
    workspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "eff_" << prefix2 << "_" << eff1.first << "b[" << effs2.at(eff1.first) << "]";
    workspace.factory(ss.str().c_str());
  }
}

void parseConfigFile(const std::string& filename, /*RooAbsCategoryLValue& categories,*/RooWorkspace& workspace) {

  std::cout << "Loading fit configuration file from " << filename << std::endl;
  std::fstream configFile((BASE_PATH + "/" + filename).c_str(), std::ios::in);
  std::string line;
  while (std::getline(configFile, line)) {
    if (line[0] == '#' || line.length() == 0)
      continue;

    workspace.factory(line.c_str());
  }

  std::cout << "Content of workspace" << std::endl;
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

RooAbsPdf* getInterpolatedPdf(RooRealVar& observable, double massParticle, const std::string& jec, int btag, const std::string& categoryName, const std::string& suffix = "") {

  // Interpolation
  int lowMass = 0;
  int highMass = 0;

  const int ALGO_MOMENT_MORPH = 1;
  const int ALGO_INTEGRAL_MORPH = 2;

  int algo;

  if (massParticle > 500 && massParticle < 750) {
    lowMass = 500;
    highMass = 750;
    algo = ALGO_INTEGRAL_MORPH;
  } else if (massParticle > 750 && massParticle < 1000) {
    lowMass = 750;
    highMass = 1000;
    algo = ALGO_MOMENT_MORPH;
  } else if (massParticle > 1000 && massParticle < 1250) {
    lowMass = 1000;
    highMass = 1250;
    algo = ALGO_MOMENT_MORPH;
  } else if (massParticle > 1250 && massParticle < 1500) {
    lowMass = 1250;
    highMass = 1500;
    algo = ALGO_MOMENT_MORPH;
  }/* else if (massParticle > 1500 && massParticle < 2000) {
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
  TString lowMass_prefix = TString::Format("%s-%s%d_%s_%d_btag", jec.c_str(), getAnalysisPrefix(), lowMass, analysisName.c_str(), btag);
  TString lowMass_workspaceFile = BASE_PATH + "/frit/" + lowMass_prefix + "_workspace.root";

  TString highMass_prefix = TString::Format("%s-%s%d_%s_%d_btag", jec.c_str(), getAnalysisPrefix(), highMass, analysisName.c_str(), btag);
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
    alpha = (double) (massParticle - lowMass) / (double) (highMass - lowMass);
  } else {
    alpha = 1. - (double) (massParticle - lowMass) / (double) (highMass - lowMass);
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

  //RooHistPdf* hist_pdf = new RooHistPdf(std::string("signal_" + categoryName).c_str(), std::string("signal_" + categoryName).c_str(), RooArgSet(observable), *binnedInterpolatedDataset);

  //interpolation->SetName(std::string("signal_" + categoryName).c_str());
  //renameAndSetPdfParametersConst(RooArgSet(observable), *interpolation, (suffix.length() == 0) ? cleanedCategory : suffix);
  renameAndSetPdfParametersConst(RooArgSet(observable), *keys_pdf, goodSuffix);

  observable.setBins(oldBinning);

  return keys_pdf;
}

bool isInterpolated(AnalysisType type, int mass) {
  if (type == HIGGS) {
    return !(mass == 400 || mass == 500 || mass == 600 || mass == 700 || mass == 800);
  } else {
    return !(mass == 500 || mass == 750 || mass == 1000 || mass == 1250 || mass == 1500 || mass == 2000);
  }
}

void fitMtt(std::map<int, TChain*> eventChain, int massParticle, bool fit, string fitConfigurationFile, bool writeRootFile, bool saveFigures, bool doBiasTest, int nToyExp, int index, string syst_str, bool bkgOnly, bool muonsOnly, int btag, const std::string& customWorkspaceFile, bool fixBackground, bool saveWorkspace, bool createCombineWorkspace, const std::string& combineWorkspaceFilename)
{

  if ((syst_str != "nominal") && (syst_str != "JECup") && (syst_str != "JECdown"))
  {
    cout << "ERROR : Badly defined JEC configuration" << endl;
    return;
  }

  if (fitConfigurationFile == "auto" || fitConfigurationFile.empty())
    fitConfigurationFile = "fit_pdf_falt.json";

  std::cout << "Loading fit configuration from '" << fitConfigurationFile << "'" << std::endl;

  const bool combine = (btag > 2);

  // Hard code some analysis
  // THIS HAS TO BE FIXME

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

  AnalysisType analysisType = getAnalysisType();

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
    RooMsgService::instance().addStream(RooFit::ERROR); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
    RooMsgService::instance().setSilentMode(true);
  } else {
    RooMsgService::instance().addStream(RooFit::INFO); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
  }

  //fit region
  Float_t minmTT = 325;
  Float_t maxmTT = 1000;

  if (analysisType == ZPRIME) {
    minmTT = 550;
    maxmTT = 2000;
  }

  RooRealVar mtt("mtt", "M_{t#bar{t}}", minmTT, maxmTT, "GeV");
  RooRealVar weight("weight", "weight", 0, 100000);

  Double_t minmTTFit = minmTT + 0.0;
  Double_t maxmTTFit = maxmTT - 0.0;

  mtt.setRange(minmTTFit, maxmTTFit);

  mtt.setRange("R1", minmTTFit, massParticle - 50);
  mtt.setRange("R2", massParticle + 50, maxmTTFit);

  // Set binning to 1 GeV
  std::cout << mtt.getBins() << std::endl;
  mtt.setBins((maxmTTFit - minmTTFit) / 4.);
  mtt.setBins(5000, "cache");

  RooCategory lepton_type("lepton_type", "lepton_type");
  lepton_type.defineType("muon", 13);
  if (! muonsOnly) {
    lepton_type.defineType("electron", 11);
  }

  RooCategory btagCategory("btag", "btag");

  if (minBTag == 0)
    btagCategory.defineType("0-btag");

  if (minBTag <= 1 && maxBTag >= 1)
    btagCategory.defineType("1-btag");

  if (minBTag <= 2 && maxBTag >= 2)
    btagCategory.defineType("2-btag");

  RooSuperCategory superCategory("superCat", "superCat", RooArgList(lepton_type, btagCategory));

  RooAbsCategoryLValue& mainCategory = combine ? static_cast<RooAbsCategoryLValue&>(superCategory) : static_cast<RooAbsCategoryLValue&>(lepton_type);

  // Create main workspace for global pdf
  RooWorkspace mainWorkspace("mainWorkspace", "main workspace");

  std::string analysisName = getAnalysisName(BASE_PATH);

  std::map<std::string, std::shared_ptr<BaseFunction>> backgroundPdfs = getCategoriesPdf(BASE_PATH + "/fit_configuration", fitConfigurationFile, mtt, NULL, massParticle, "background", mainCategory, NULL);

  std::map<std::string, RooAbsPdf*> backgroundPdfsFromWorkspace;

  for (auto& pdf: backgroundPdfs) {
    std::cout << "Background pdf: " << pdf.first << " ";
    pdf.second->getPdf().Print();
    mainWorkspace.import(pdf.second->getPdf());

    backgroundPdfsFromWorkspace[pdf.first] = mainWorkspace.pdf(pdf.second->getPdf().GetName());
  }

  TString prefix = TString::Format("data_2012_%s_%d", syst_str.c_str(), massParticle);
  //TString suffix = TString::Format("%s_%s", pdfSignalName.c_str(), bkgfit_str.c_str());
  TString suffix = TString::Format("%s", analysisName.c_str());
  TString indexJob = TString::Format("job%d", index);


  // Luminosity updated for full 2012 run, with PIXEL calculation, and Winter13 rereco
  double lumi_mu = 19700;
  double lumi_e  = 19667;

  double s_lumi_percent  = 2.6 / 100.; // 2.6%

  // Define theoretical cross section
  RooConstVar sigma("sigma", "theoretical cross-section", 1);

  /*
   * Load efficiencies.
   * If we are working in combined mode, load for each b-tag categories.
   * Otherwise, load only the requested b-tag category
   */
  std::map<int, double> sel_eff_mu; // Selection efficiency in semi-mu channel for each b-tag
  std::map<int, double> sel_eff_e; // Selection efficiency in semi-e channel for each b-tag

  std::map<int, double> hlt_eff_mu; // HLT efficiency in semi-mu channel for each b-tag
  std::map<int, double> hlt_eff_e; // HLT efficiency in semi-e channel for each b-tag

  std::map<int, double> s_sel_eff_mu; // = Δ(eff_mu)
  std::map<int, double> s_sel_eff_e; // = Δ(eff_e)

  std::map<int, double> s_hlt_eff_mu; // Error on HLT efficiency
  std::map<int, double> s_hlt_eff_e; // Error on HLT efficiency

  if (! combine) {
    loadSelectionEfficiencies(massParticle, syst_str, btag, sel_eff_mu[btag], sel_eff_e[btag], s_sel_eff_mu[btag], s_sel_eff_e[btag]);
    if (analysisType == ZPRIME)
      loadTriggerEfficiencies(massParticle, syst_str, btag, hlt_eff_mu[btag], hlt_eff_e[btag], s_hlt_eff_mu[btag], s_hlt_eff_e[btag]);
  } else {
    // Load for 0, 1 and 2 btag
    for (int i = minBTag; i <= maxBTag; i++) {
      loadSelectionEfficiencies(massParticle, syst_str, i, sel_eff_mu[i], sel_eff_e[i], s_sel_eff_mu[i], s_sel_eff_e[i]);
      if (analysisType == ZPRIME)
        loadTriggerEfficiencies(massParticle, syst_str, i, hlt_eff_mu[i], hlt_eff_e[i], s_hlt_eff_mu[i], s_hlt_eff_e[i]);
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

  struct SignalPdf {
    RooAbsPdf* positive_pdf;
    RooRealVar* positive_integral;

    RooAbsPdf* negative_pdf;
    RooRealVar* negative_integral;
  };

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

    std::string pdfName = std::string((TString(type->GetName()).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron")) + "_%s_signal_pdf";
    std::string integralName = std::string((TString(type->GetName()).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron")) + "_%s_integral";

    if (! isInterpolated(analysisType, massParticle)) {

      TString workspaceFile = TString::Format("%s/frit/%s-%s%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), syst_str.c_str(), getAnalysisPrefix(), massParticle, analysisName.c_str(), categoryBTag);
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

      TString cleanedCategory = TString(type->GetName()).ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", ""); // For workspace names
      if (analysisType == HIGGS) {

        RooAbsPdf* positive_pdf = workspace->pdf(TString::Format(pdfName.c_str(), "positive"));
        positive_pdf->SetName(std::string("positive_signal_" + std::string(cleanedCategory.Data())).c_str());
        renameAndSetPdfParametersConst(RooArgSet(mtt), *positive_pdf, type->GetName());
        mainWorkspace.import(*positive_pdf);

        RooAbsPdf* negative_pdf = workspace->pdf(TString::Format(pdfName.c_str(), "negative"));
        if (! positive_pdf) {
          std::cerr << "ERROR: Signal pdf " << pdfName << " not found inside workspace." << std::endl;
          exit(1);
        }

        negative_pdf->SetName(std::string("negative_signal_" + std::string(cleanedCategory.Data())).c_str());
        renameAndSetPdfParametersConst(RooArgSet(mtt), *negative_pdf, type->GetName());
        mainWorkspace.import(*negative_pdf);

        RooRealVar* positive_integral = workspace->var(TString::Format(integralName.c_str(), "positive"));
        positive_integral->SetName(TString::Format("%s_positive_integral", cleanedCategory.Data()));
        mainWorkspace.import(*positive_integral);

        RooRealVar* negative_integral = workspace->var(TString::Format(integralName.c_str(), "negative"));
        negative_integral->SetName(TString::Format("%s_negative_integral", cleanedCategory.Data()));
        mainWorkspace.import(*negative_integral);

      } else {
        std::string zprimePdfName = "signal_" + std::string((TString(type->GetName()).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron"));

        RooAbsPdf* pdf = workspace->pdf(zprimePdfName.c_str());
        pdf->SetName(std::string("signal_" + std::string(cleanedCategory)).c_str());
        renameAndSetPdfParametersConst(RooArgSet(mtt), *pdf, cleanedCategory.Data());
        mainWorkspace.import(*pdf);
      }
    } else {
      RooAbsPdf* interpolation = getInterpolatedPdf(mtt, massParticle, "nominal", categoryBTag, type->GetName());
      mainWorkspace.import(*interpolation);
    }
  }

  std::cout << "Done." << std::endl;

  //if (systCB != "none") // Signal systematics mode
  //{
    //RooAbsArg* arg = mainWorkspace.arg(systCB.c_str());
    //if (! arg)
    //{
      //std::cerr << "ERROR: parameter " << systCB << " not found inside workspace!" << std::endl;
      //exit(1);
    //}

    //RooRealVar* var = static_cast<RooRealVar*>(arg);
    //double oldValue = var->getVal();
    //if (systCBsign == "up")
    //{
      //var->setVal(var->getVal() + var->getErrorHi());
    //}
    //else
    //{
      //// Note: getErrorLo() returns something negative.
      //var->setVal(var->getVal() + var->getErrorLo());
    //}

    //std::cout << systCB << " value changed from " << oldValue << " to " << var->getVal() << std::endl;
  //}

  // In Z' analysis, total_eff is the product of selection efficiency and trigger efficiency
  // For the Higgs analysis, trigger efficiency is already included in selection efficiency
  std::map<int, double> total_eff_mu;
  std::map<int, double> total_eff_e;

  std::map<int, double> total_err_eff_mu;
  std::map<int, double> total_err_eff_e;

  auto combineEfficiencies = [](float sel, float hlt) {
    return sel * hlt;
  };

  auto combineEfficiencyErrors = [combineEfficiencies](float sel, float error_sel, float hlt, float error_hlt) {
    float sel_ratio = (error_sel * error_sel) / (sel * sel);
    float hlt_ratio = (error_hlt * error_hlt) / (hlt * hlt);
    float eff = combineEfficiencies(sel, hlt);
    return sqrt(eff * eff * (sel_ratio + hlt_ratio));
  };

  if (analysisType == HIGGS) {
    total_eff_mu = sel_eff_mu;
    total_eff_e = sel_eff_e;

    total_err_eff_mu = s_sel_eff_mu;
    total_err_eff_e = s_sel_eff_e;
  } else {

    if (! combine) {
      total_eff_mu[btag] = combineEfficiencies(sel_eff_mu[btag], hlt_eff_mu[btag]);
      total_eff_e[btag] = combineEfficiencies(sel_eff_e[btag], hlt_eff_e[btag]);

      total_err_eff_mu[btag] = combineEfficiencyErrors(sel_eff_mu[btag], s_sel_eff_mu[btag], hlt_eff_mu[btag], s_hlt_eff_mu[btag]);
      total_err_eff_e[btag] = combineEfficiencyErrors(sel_eff_e[btag], s_sel_eff_e[btag], hlt_eff_e[btag], s_hlt_eff_e[btag]);
    } else {
      // Load for 0, 1 and 2 btag
      for (int i = minBTag; i <= maxBTag; i++) {
        total_eff_mu[i] = combineEfficiencies(sel_eff_mu[i], hlt_eff_mu[i]);
        total_eff_e[i] = combineEfficiencies(sel_eff_e[i], hlt_eff_e[i]);

        total_err_eff_mu[i] = combineEfficiencyErrors(sel_eff_mu[i], s_sel_eff_mu[i], hlt_eff_mu[i], s_hlt_eff_mu[i]);
        total_err_eff_e[i] = combineEfficiencyErrors(sel_eff_e[i], s_sel_eff_e[i], hlt_eff_e[i], s_hlt_eff_e[i]);
      }
    }
  }

  std::cout << std::endl;
  std::cout << "Analysis efficiencies" << std::endl;
  std::cout << "-------------------------------------" << std::endl << std::endl;

  std::cout << "Luminosity: " << Bash::set_color(Bash::Color::MAGENTA) << lumi_mu << " /pb" << Bash::set_color() << std::endl;
  std::cout << "Total efficiency (selection * SFs): " << std::endl;

  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << "\tMuonic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << total_eff_mu[i] * 100 << " % +/- " << total_err_eff_mu[i] * 100 << " %" << Bash::set_color() << std::endl;
    std::cout << "\tElectronic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << total_eff_e[i] * 100 << " % +/- " << total_err_eff_e[i] * 100 << " %" << Bash::set_color() << std::endl;
  }
  std::cout << "Expected number of signal events: " << std::endl;

  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << "\tMuonic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << sigma.getVal() * lumi_mu * total_eff_mu[i] << Bash::set_color() << std::endl;
    std::cout << "\tElectronic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << sigma.getVal() * lumi_e * total_eff_e[i] << Bash::set_color() << std::endl;
  }

  if (combine) {
    addEfficienciesToWorkspace(total_eff_e, total_eff_mu, total_err_eff_e, total_err_eff_mu, "e", "mu", mainWorkspace);
  } else {
    std::stringstream ss;
    ss << "eff_ratio_e_mu[" << total_eff_e[btag] / total_eff_mu[btag]<< "]";
    mainWorkspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "eff_e[" << total_eff_e[btag] << "]";
    mainWorkspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "eff_mu[" << total_eff_mu[btag] << "]";
    mainWorkspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "error_eff_e[" << total_err_eff_e[btag] << "]";
    mainWorkspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "error_eff_mu[" << total_err_eff_mu[btag] << "]";
    mainWorkspace.factory(ss.str().c_str());
  }

  RooRealVar lumiRatio("lumiRatio", "luminosity ratio", lumi_e / lumi_mu, "");
  mainWorkspace.import(lumiRatio);

  // Import luminosity inside workspace
  mainWorkspace.factory(TString::Format("lumi_e[%.5f]", lumi_e));
  mainWorkspace.factory(TString::Format("lumi_mu[%.5f]", lumi_mu));
  mainWorkspace.factory(TString::Format("error_lumi_e[%.5f]", lumi_e * s_lumi_percent));
  mainWorkspace.factory(TString::Format("error_lumi_mu[%.5f]", lumi_mu * s_lumi_percent));

  // Create constraints
  // First, efficiency
  if (combine) {
    for (auto& eff_e: total_eff_e) {
      TString d = TString::Format("eff_e_%db_constrained[%.8f,0,1]", eff_e.first, eff_e.second);
      mainWorkspace.factory(d);
    }

    for (auto& eff_mu: total_eff_mu) {
      TString d = TString::Format("eff_mu_%db_constrained[%.8f,0,1]", eff_mu.first, eff_mu.second);
      mainWorkspace.factory(d);
    }
  } else {
      TString d = TString::Format("eff_e_constrained[%.8f,0,1]", total_eff_e[btag]);
      mainWorkspace.factory(d);

      d = TString::Format("eff_mu_constrained[%.8f,0,1]", total_eff_mu[btag]);
      mainWorkspace.factory(d);
  }

  // Second, luminosity
  {
    TString d = TString::Format("lumi_e_constrained[%.5f, %.5f, %.5f]", lumi_e, lumi_e - 5 * (lumi_e * s_lumi_percent), lumi_e + 5 * (lumi_e * s_lumi_percent));
    mainWorkspace.factory(d);

    d = TString::Format("lumi_mu_constrained[%.5f, %.5f, %.5f]", lumi_mu, lumi_mu - 5 * (lumi_mu * s_lumi_percent), lumi_mu + 5 * (lumi_mu * s_lumi_percent));
    mainWorkspace.factory(d);
  }

  if (! combine) {
    std::cout << "Efficiencies ratio x lumis ratio: " << total_eff_e[btag] / total_eff_mu[btag] * lumi_e / lumi_mu << std::endl;
  }

  // Read config file for global pdf
  std::string configFile = (combine)
    ? TString::Format("combined_%s_global_pdf_%s_analysis.script", scriptNamePrefix.c_str(), getAnalysisPrefix()).Data()
    : TString::Format("individual_btag_global_pdf_%s_analysis.script", getAnalysisPrefix()).Data();

  mainWorkspace.import(sigma);

  std::cout << std::endl;
  parseConfigFile(configFile, /*mainCategory,*/ mainWorkspace);

  // mu is the signal strength, sigma_signal / sigma_theoretical
  RooRealVar& mu = *mainWorkspace.var("mu");

  // Constraints
  std::vector<RooRealVar*> constrained_variables;

  RooArgSet variables = mainWorkspace.allVars();
  it = variables.createIterator();
  for (RooRealVar* var = (RooRealVar*) it->Next(); var != nullptr; var = (RooRealVar*) it->Next()) {
    if (TString(var->GetName()).Contains("constrained")) {
      constrained_variables.push_back(var);
    }
  }

  RooArgSet constraint_pdfs;

  RooArgSet pdfs = mainWorkspace.allPdfs();
  it = pdfs.createIterator();
  for (RooAbsPdf* pdf = (RooAbsPdf*) it->Next(); pdf != nullptr; pdf = (RooAbsPdf*) it->Next()) {
    if (TString(pdf->GetName()).Contains("constraint")) {
      constraint_pdfs.add(*pdf);
    }
  }

  //FIXME
  bool useConstraints = false;
  if (! useConstraints) {
    for (RooRealVar* var: constrained_variables) {
      var->setConstant(true);
    }
  }

  RooSimultaneous simPdf("simPdf", "simultaneous pdf", mainCategory);
  RooSimultaneous simPdfBackgroundOnly("simPdfBackgroundOnly", "simultaneous pdf with background only pdfs", mainCategory);

  if (bkgOnly)
    std::cout << "Warning: fitting with background model only!" << std::endl;

  // NO shared_ptr<> here.
  std::map<std::string, RooAbsPdf*> globalPdfs;

  it = mainCategory.typeIterator();
  type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    std::string name = type->GetName();

    std::string cleanedCategory = TString(type->GetName()).ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", "").Data(); // For workspace names

    std::string workspaceName = "global_pdf_" + cleanedCategory;
    std::cout << "Looking for " << workspaceName << " inside workspace" << std::endl;
    globalPdfs[name] = mainWorkspace.pdf(workspaceName.c_str());

    const RooAbsPdf& pdf = *globalPdfs[name];

    std::cout << "Adding pdf ";
    pdf.Print();
    std::cout << " for category " << name << std::endl;

    simPdf.addPdf(pdf, name.c_str());
    simPdfBackgroundOnly.addPdf(*backgroundPdfsFromWorkspace[name], name.c_str());
  }

  if (bkgOnly) {
    mu.setVal(0);
    mu.setConstant(true);
  }

  // Create output folder
  if (! BATCH_MODE) {
    TString folderName = TString::Format("%s/%s_%s_%d_btag/", OUTPUT_PATH.c_str(), prefix.Data(), suffix.Data(), btag);
    mkdir(folderName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    OUTPUT_PATH = folderName;
  }

  if (fit)
  {
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

    if (createCombineWorkspace) {

      // Bash code for datacards
      Roo1DTable* table = dataOrig->table(superCategory);
      for (int i = minBTag; i <= maxBTag; i++) {

        int index = (i == 1) ? 0 : 2;
        TString filename = TString::Format("datacard_%d_%dbtag.txt", massParticle, i);

        superCategory.setIndex(index);

        std::cout << "sed -i \"6s/.*/observation ";
        int nMu = table->get(superCategory.getLabel());
        std::cout << nMu;

        std::cout << "     ";

        superCategory.setIndex(index + 1);
        int nE = table->get(superCategory.getLabel());
        std::cout << nE;

        std::cout << "/\" " << filename.Data() << std::endl;

        std::cout << "sed -i \"13s/.*/rate        ";
        std::cout << sigma.getVal() * lumi_mu * total_eff_mu[i] << "   " << nMu << "         " << sigma.getVal() * lumi_e * total_eff_e[i] << "   " << nE << "/\" " << filename.Data() << std::endl;

        std::cout << "sed -i \"s/eff.*/eff  lnN    ";
        int oldPre = std::cout.precision();
        std::cout.precision(8);
        std::cout << std::fixed << 1 + total_err_eff_mu[i] << "   -             " << std::fixed << 1 + total_err_eff_e[i] << "   -/\" " << filename.Data() << std::endl;
        std::cout.precision(oldPre);
      }

      std::cout << std::endl << "--- Bash code for systematics datacards ---" << std::endl;

      std::vector<std::string> systematics = {"jecUp", "jecDown", "jerUp", "jerDown", "puUp", "puDown"};
      std::vector<std::string> systematicsName = {"JECup", "JECdown", "JERup", "JERdown", "puUp", "puDown"};
      //std::vector<std::string> systematics = {"jecUp", "jecDown", "jerUp", "jerDown", "puUp", "puDown", "pdfUp", "pdfDown"};
      //std::vector<std::string> systematicsName = {"JECup", "JECdown", "JERup", "JERdown", "puUp", "puDown", "pdfUp", "pdfDown"};

      int index_syst = -1;
      for (const std::string& syst: systematics) {
        index_syst++;

        std::map<int, double> sel_eff_mu_syst; // Selection efficiency in semi-mu channel for each b-tag
        std::map<int, double> sel_eff_e_syst; // Selection efficiency in semi-e channel for each b-tag

        std::map<int, double> hlt_eff_mu_syst; // HLT efficiency in semi-mu channel for each b-tag
        std::map<int, double> hlt_eff_e_syst; // HLT efficiency in semi-e channel for each b-tag

        std::map<int, double> s_sel_eff_mu_syst; // = Δ(nSig_mu) / nSig_mu ; relative error
        std::map<int, double> s_sel_eff_e_syst;

        std::map<int, double> s_hlt_eff_mu_syst; // ?
        std::map<int, double> s_hlt_eff_e_syst; // ?

        // Get efficiency for systematic
        std::map<int, double> total_eff_mu_syst;
        std::map<int, double> total_eff_e_syst;

        for (int i = minBTag; i <= maxBTag; i++) {
          loadSelectionEfficiencies(massParticle, systematicsName[index_syst], i, sel_eff_mu_syst[i], sel_eff_e_syst[i], s_sel_eff_mu_syst[i], s_sel_eff_e_syst[i], true);
          if (analysisType == ZPRIME)
            loadTriggerEfficiencies(massParticle, systematicsName[index_syst], i, hlt_eff_mu_syst[i], hlt_eff_e_syst[i], s_hlt_eff_mu_syst[i], s_hlt_eff_e_syst[i], true);
        }

        if (analysisType == HIGGS) {
          total_eff_mu = sel_eff_mu_syst;
          total_eff_e = sel_eff_e_syst;

        } else {

          if (! combine) {
            total_eff_mu_syst[btag] = combineEfficiencies(sel_eff_mu_syst[btag], hlt_eff_mu_syst[btag]);
            total_eff_e_syst[btag] = combineEfficiencies(sel_eff_e[btag], hlt_eff_e[btag]);
          } else {
            // Load for 0, 1 and 2 btag
            for (int i = minBTag; i <= maxBTag; i++) {
              total_eff_mu_syst[i] = combineEfficiencies(sel_eff_mu_syst[i], hlt_eff_mu_syst[i]);
              total_eff_e_syst[i] = combineEfficiencies(sel_eff_e_syst[i], hlt_eff_e_syst[i]);
            }
          }
        }

        std::vector<double> rates_signal;
        std::vector<double> rates_bkg;

        for (int i = minBTag; i <= maxBTag; i++) {

          int index = (i == 1) ? 0 : 2;

          superCategory.setIndex(index);
          int nMu = table->get(superCategory.getLabel());
          rates_bkg.push_back(nMu);

          superCategory.setIndex(index + 1);
          int nE = table->get(superCategory.getLabel());
          rates_bkg.push_back(nE);

          rates_signal.push_back(sigma.getVal() * lumi_mu * total_eff_mu_syst[i]);
          rates_signal.push_back(sigma.getVal() * lumi_e * total_eff_e_syst[i]);
        }

        TString filename = TString::Format("datacard_%d_1+2btag_%s.txt", massParticle, syst.c_str());
        std::cout << "sed -i \"25s/.*/rate        " << rates_signal[0] << "    " << rates_bkg[0]
          << "  " << rates_signal[1] << "    " << rates_bkg[1]
          << "  " << rates_signal[2] << "    " << rates_bkg[2]
          << "  " << rates_signal[3] << "    " << rates_bkg[3] << "/\" " << filename.Data() << std::endl;

      }

      std::cout << std::endl;
      delete table;

      std::cout << "Efficiency systematic details: " << std::endl;

      /*
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
      */

      std::cout << std::endl << std::endl;
    }

    std::cout << "Fitting..." << std::endl;

    TFile* outputFile = nullptr;
    if (writeRootFile) {
      outputFile = TFile::Open(OUTPUT_PATH + prefix + "_fitRes_" + suffix + ".root", "RECREATE");
    }


    RooFitResult *fitResult = nullptr;

    // First, background only fit
    std::cout << "Background only ..." << std::endl;

    fitResult = simPdfBackgroundOnly.fitTo(*datasetToFit, Save(),/*, Optimize(0),*/ Strategy(1), Minimizer(
          "Minuit", "migrad")
        );

    fitResult->Print("v");
    SAFE_DELETE(fitResult);

    if (! createCombineWorkspace)
      drawHistograms(mainCategory, mtt, *dataOrig, simPdfBackgroundOnly, backgroundPdfsFromWorkspace, btag, saveFigures, std::string(prefix), std::string(suffix) + "_bkg_only", false, true, nullptr);

    auto createWorkspace = [&](const std::string& filename, bool afterFit = false) {
      // This function create a workspace named 'filename' containing all the pdfs

      RooWorkspace wspace("w");

      if (afterFit) {
        fitResult->SetName("fit_results");
        wspace.import(*fitResult);
        wspace.import(*dataOrig);
      }

      // Iterates over all category of the analysis
      it = mainCategory.typeIterator();
      type = nullptr;
      while ((type = static_cast<RooCatType*>(it->Next()))) {

        std::string category = type->GetName();
        std::string cleanedCategory = TString(type->GetName()).ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", "").Data(); // For workspace names

        mainCategory = category.c_str();
        std::string cut = buildCutFormula(mainCategory);

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

        // Import the data, reduced based on the current categogry
        wspace.import(*dataset_reduced, RooFit::Rename(TString::Format("data_obs_%s", workspace_suffix.Data())));

        auto importPdf = [&](const char* name, const char* newName) {
          RooAbsPdf* pdf = mainWorkspace.pdf(name);
          //setPdfParametersRange(RooArgSet(mtt), *pdf, 10);
          wspace.import(*pdf, RooFit::RenameVariable(name, newName));
        };

        // Import background
        TString name = TString::Format("background_%s", cleanedCategory.c_str());
        importPdf(name, TString::Format("background_%s", workspace_suffix.Data()));

        // Import signal
        if (analysisType == HIGGS) {
          name = TString::Format("positive_signal_%s", cleanedCategory.c_str());
          importPdf(name, TString::Format("positive_signal_%s", workspace_suffix.Data()));

          name = TString::Format("negative_signal_%s", cleanedCategory.c_str());
          importPdf(name, TString::Format("positive_signal_%s", workspace_suffix.Data()));
        } else {
          name = TString::Format("signal_%s", cleanedCategory.c_str());
          importPdf(name, TString::Format("signal_%s", workspace_suffix.Data()));
        }

        if (afterFit) {
          // Import global pdf
          RooAddPdf* globalPdf = dynamic_cast<RooAddPdf*>(mainWorkspace.pdf(workspaceName.c_str()));
          wspace.import(*globalPdf, RooFit::RecycleConflictNodes(), RooFit::RenameVariable(workspaceName.c_str(), TString::Format("global_pdf_%s", workspace_suffix.Data())));
        }

        // Import signal systematics
        std::function<void (const char*, const char*)> importSystPdf;

        if (! isInterpolated(analysisType, massParticle)) {

          importSystPdf = [&](const char* oldSyst, const char* newSyst) {

            TString workspaceFile = TString::Format("%s/frit/%s-%s%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), oldSyst, getAnalysisPrefix(), massParticle, analysisName.c_str(), extractedBTag);
            std::shared_ptr<TFile> f(TFile::Open(workspaceFile.Data()));

            if (analysisType == HIGGS) {
              TString n = TString::Format("%s_positive_signal_pdf", leptonName.c_str());
              RooAbsPdf* pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(n);
              pdf->SetName(TString::Format("positive_signal_%s_%s", workspace_suffix.Data(), newSyst));
              wspace.import(*pdf);

              n = TString::Format("%s_negative_signal_pdf", leptonName.c_str());
              pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(n);
              pdf->SetName(TString::Format("negative_signal_%s_%s", workspace_suffix.Data(), newSyst));
              wspace.import(*pdf);
            } else {
              TString n = TString::Format("signal_%s", leptonName.c_str());
              RooAbsPdf* pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(n);
              pdf->SetName(TString::Format("signal_%s_%s", workspace_suffix.Data(), newSyst));
              wspace.import(*pdf);
            }
          };

        } else {

          importSystPdf = [&](const char* oldSyst, const char* newSyst) {
            TString n = TString::Format("signal_%s_%s", workspace_suffix.Data(), newSyst);
            RooAbsPdf* pdf = getInterpolatedPdf(mtt, massParticle, oldSyst, extractedBTag, category, n.Data());
            pdf->SetName(n);
            wspace.import(*pdf);
          };

        }

        importSystPdf("JECup", "jecUp");
        importSystPdf("JECdown", "jecDown");

        importSystPdf("JERup", "jerUp");
        importSystPdf("JERdown", "jerDown");

        importSystPdf("puUp", "puUp");
        importSystPdf("puUp", "puDown");

        //importSystPdf("pdfUp", "pdfUp");
        //importSystPdf("pdfUp", "pdfDown");
      }

      wspace.writeToFile(filename.c_str());
    };

    // Set parameter range to +/- 10 sigmas for background + signal fit
    it = mainCategory.typeIterator();
    type = nullptr;
    while ((type = static_cast<RooCatType*>(it->Next()))) {
      setPdfParametersRange(RooArgSet(mtt), *simPdfBackgroundOnly.getPdf(type->GetName()), 10);
    }

    if (createCombineWorkspace) {
      TString outputFileName = TString::Format("%s_%d_workspace.root", getAnalysisPrefix(), massParticle);
      if (combineWorkspaceFilename.length() > 0)
        outputFileName = combineWorkspaceFilename;

      createWorkspace(outputFileName.Data(), false);

      // If create combine workspace, exit as soon as it's done
      return;
    }

    if (bkgOnly) {
      computeChi2(mtt, simPdfBackgroundOnly, mainCategory, *dataOrig, mainWorkspace, true);
    } else {
      std::cout << "Background (floating) + signal ..." << std::endl;

      if (useConstraints) {
        fitResult = simPdf.fitTo(*datasetToFit, Save(), Strategy(1), Minimizer("Minuit", "migrad"), ExternalConstraints(constraint_pdfs));
      } else {
        fitResult = simPdf.fitTo(*datasetToFit, Save(), Strategy(1), Minimizer("Minuit", "migrad"));
      }

      fitResult->Print("v");

      std::cout << "Done." << std::endl;

      if (saveWorkspace) {
        // Save fitted pdf and datasets in order to redo some plots
        TString outputFileName = TString::Format("%s_%d_workspace_after_fit.root", getAnalysisPrefix(), massParticle);
        createWorkspace(outputFileName.Data(), true);
      }

      drawHistograms(mainCategory, mtt, *dataOrig, simPdf, backgroundPdfsFromWorkspace, btag, saveFigures, std::string(prefix), std::string(suffix), !bkgOnly, true, outputFile);

      std::map<std::string, float> chi2 = computeChi2(mtt, simPdf, mainCategory, *dataOrig, mainWorkspace);

      if (outputFile) {
        outputFile->cd();
        fitResult->Write("fitResult");
      }

      if (outputFile) {
        outputFile->Close();
        delete outputFile;
      }

      if (fixBackground) {
        // Set background pdf parameters as constant for likelihood scan
        it = mainCategory.typeIterator();
        type = nullptr;
        while ((type = static_cast<RooCatType*>(it->Next()))) {
          setPdfParametersConst(mtt, *backgroundPdfsFromWorkspace[type->GetName()], true);
        }
      }

      double sigma_signal = mu.getVal() * sigma.getVal();
      double sigma_error = mu.getError() / mu.getVal() * sigma_signal;

      double Limit_Z_obs_pb = sigma_signal + 2. * sigma_error;

      cout << "The signal cross section is " << sigma_signal << " +/- " << sigma_error << " pb" << endl;
      cout << "The estimated 95% C.L. upper limit on the signal cross section is " << Limit_Z_obs_pb << " pb" << endl;

      SAFE_DELETE(fitResult);
    }

    SAFE_DELETE(fitResult);

    // WORKHERE
    if (doBiasTest)
    {

      std::cout << "Starting bias analysis..." << std::endl;

      // 0 means a random seed
      RooRandom::randomGenerator()->SetSeed(0);

      // Create PDF for toys generation
      std::map<std::string, std::shared_ptr<BaseFunction>> backgroundPdfsForToys = getCategoriesPdf(BASE_PATH + "/fit_configuration", fitConfigurationFile, mtt, NULL, massParticle, "background", mainCategory, nullptr, "toy");

      std::map<std::string, std::shared_ptr<RooAbsPdf>> globalPdfsForToys;

      // Create PDF for fitting
      RooSimultaneous simPdfToyFit("simPdfToyFit", "simultaneous pdf for toys fitting", mainCategory);
      RooSimultaneous simPdfToyFitBackgroundOnly("simPdfToyFitBackgroundOnly", "simultaneous pdf for toys fitting (background only)", mainCategory);

      TFile* toyResFile = NULL;
      //if (writeRootFile) {
      toyResFile = new TFile(OUTPUT_PATH + prefix + "_biastest_" + suffix + "_" + indexJob + ".root", "RECREATE");
      toyResFile->cd();
      //}

      // Book historams for toy MC results
      TH1F* hNToyEvents = new TH1F("hNToyEvents", "Number of events of toy experiment", 100, 200000, 300000);
      TH1F* hToyChi2 = new TH1F("hToyChi2", "chi2 of toys experiments", 50, 0, 5);
      TH1F* hToyMu = new TH1F("hToyMu", "mu of toys experiments", 100, -1, 1);
      TH1F* hToyMuErr = new TH1F("hToyMuErr", "mu error of toys experiments", 100, 0, 1);
      TH1F* hToyMuPull = new TH1F("hToyPull", "pull of toys experiments", 100, -3, 3);

      std::cout << "Preparing to draw from data" << std::endl;
      std::map<int, std::shared_ptr<RooDataHist>> muon_datahist;
      std::map<int, std::shared_ptr<RooDataHist>> electron_datahist;

      std::map<int, std::shared_ptr<RooAbsPdf>> muon_pdf_from_data;
      std::map<int, std::shared_ptr<RooAbsPdf>> electron_pdf_from_data;

      std::map<int, int> nMuons;
      std::map<int, int> nElectrons;

      it = mainCategory.typeIterator();
      type = nullptr;
      while ((type = static_cast<RooCatType*>(it->Next()))) {
        std::string category = type->GetName();

        mainCategory = category.c_str();
        std::string cut = buildCutFormula(mainCategory);
        //std::cout << "Cut: " << cut << std::endl;

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


        if (leptonNameShort=="e") {
          electron_datahist[extractedBTag].reset(static_cast<RooDataHist*>(binnedDataset->reduce(RooFit::SelectVars(RooArgSet(mtt, weight)), RooFit::Cut(cut.c_str()), RooFit::Name(TString::Format("electron_binned_dataset_%d", extractedBTag)))));
          electron_pdf_from_data[extractedBTag].reset(new RooHistPdf(TString::Format("electron_hist_pdf_%d", extractedBTag), "electron_hist_pdf", RooArgSet(mtt), *electron_datahist[extractedBTag]));
          nElectrons[extractedBTag]=electron_datahist[extractedBTag]->sumEntries();
        }
        if (leptonNameShort=="mu") {
          muon_datahist[extractedBTag].reset(static_cast<RooDataHist*>(binnedDataset->reduce(RooFit::SelectVars(RooArgSet(mtt, weight)), RooFit::Cut(cut.c_str()), RooFit::Name(TString::Format("muon_binned_dataset_%d", extractedBTag)))));
          muon_pdf_from_data[extractedBTag].reset(new RooHistPdf(TString::Format("muon_hist_pdf_%d", extractedBTag), "muon_hist_pdf", RooArgSet(mtt), *muon_datahist[extractedBTag]));
          nMuons[extractedBTag]=muon_datahist[extractedBTag]->sumEntries();
        }
      }


      std::map<int, std::shared_ptr<RooAbsPdf::GenSpec>> genSpecsMuon;
      std::map<int, std::shared_ptr<RooAbsPdf::GenSpec>> genSpecsElectron;

      if (! combine) {
        genSpecsMuon[btag].reset(muon_pdf_from_data[btag]->prepareMultiGen(RooArgSet(mtt), NumEvents(nMuons[btag]), Extended(true),RooFit::AutoBinned(false)/*, Verbose(true)*/));
        genSpecsElectron[btag].reset(electron_pdf_from_data[btag]->prepareMultiGen(RooArgSet(mtt), NumEvents(nElectrons[btag]), Extended(true),RooFit::AutoBinned(false)/*, Verbose(true)*/));
      } else {
        for (int i = minBTag; i <= maxBTag; i++) {
          genSpecsMuon[i].reset(muon_pdf_from_data[i]->prepareMultiGen(RooArgSet(mtt), NumEvents(nMuons[i]), Extended(true),RooFit::AutoBinned(false)/*, Verbose(true)*/));
          genSpecsElectron[i].reset(electron_pdf_from_data[i]->prepareMultiGen(RooArgSet(mtt), NumEvents(nElectrons[i]), Extended(true),RooFit::AutoBinned(false)/*, Verbose(true)*/));
        }
      }


      for (int itoys=0 ; itoys<nToyExp ; itoys++) {
        std::cout << "---> Toy #" << itoys << std::endl;
        std::map<int, std::shared_ptr<RooDataSet>> generated_dataset_muon;
        std::map<int, std::shared_ptr<RooDataSet>> generated_dataset_electron;

        if (!combine) {
          std::cout << "Generating #" << nMuons[btag] << " muons and #" << nElectrons[btag] << " electrons..." << std::endl;
          generated_dataset_muon[btag].reset(muon_pdf_from_data[btag]->generate(*genSpecsMuon[btag]));
          generated_dataset_electron[btag].reset(electron_pdf_from_data[btag]->generate(*genSpecsMuon[btag]));
          std::cout << "done." << std::endl;
        } else {
          for (int i = minBTag; i <= maxBTag; i++) {
            std::cout << "Generating #" << nMuons[i] << " muons and #" << nElectrons[i] << " electrons..." << std::endl;
            generated_dataset_muon[i].reset(muon_pdf_from_data[i]->generate(*genSpecsMuon[i]));
            generated_dataset_electron[i].reset(electron_pdf_from_data[i]->generate(*genSpecsMuon[i]));
            std::cout << "done." << std::endl;
          }
        }

        std::map<std::string, RooDataSet*> all_generated_datasets_electron;
        std::map<std::string, RooDataSet*> all_generated_datasets_muon;
        it = mainCategory.typeIterator();
        type = nullptr;
        while ((type = static_cast<RooCatType*>(it->Next()))) {
          std::string category_str = type->GetName();
          std::cout <<  category_str << std::endl;
          std::string leptonName = TString(category_str).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron";
          static boost::regex regex("([0-9]+)-btag");
          boost::smatch regexResults;
          int extractedBTag = -1;
          if (boost::regex_search(category_str, regexResults, regex)) {
            std::string result(regexResults[1].first, regexResults[2].second);
            extractedBTag = atoi(result.c_str());
            if (leptonName=="muon")
              all_generated_datasets_muon[result]=generated_dataset_muon[extractedBTag].get();
            else if (leptonName=="electron")
              all_generated_datasets_electron[result]=generated_dataset_electron[extractedBTag].get();
          }
          if (extractedBTag < 0) {
            stringstream ss;
            ss << btag;
            std::string result = ss.str()+"-btag";
            if (leptonName=="muon")
              all_generated_datasets_muon[result]=generated_dataset_muon[btag].get();
            else if (leptonName=="electron")
              all_generated_datasets_electron[result]=generated_dataset_electron[btag].get();
          }
        }

        std::shared_ptr<RooDataSet> toyDataset_muon(new RooDataSet("combData_muon", "combined data muon", RooArgSet(mtt), Index(btagCategory), Import(all_generated_datasets_muon)));
        std::shared_ptr<RooDataSet> toyDataset_electron(new RooDataSet("combData_electron", "combined data electron", RooArgSet(mtt), Index(btagCategory), Import(all_generated_datasets_electron)));

        std::shared_ptr<RooDataSet> toyDataset(new RooDataSet("combData", "combined data", RooArgSet(mtt, btagCategory), Index(lepton_type), Import("muon",*toyDataset_muon), Import("electron",*toyDataset_electron)));

        //toyDataset->Print("v");
        //toyDataset->table(superCategory)->Print("v");

        std::shared_ptr<RooDataHist> toyBinnedDataset = std::shared_ptr<RooDataHist>(toyDataset->binnedClone());

        std::shared_ptr<RooAbsData> toyDatasetToFit = toyBinnedDataset; // Binned likelihood

        std::cout << " fit background + signal ..." << std::endl;
        fitResult = simPdf.fitTo(*toyDatasetToFit, Save(),/*, Optimize(0),*/ Strategy(1), Minimizer("Minuit2", "Migrad"));
        //fitResult->Print("v");

        std::map<std::string, float> chi2_toy = computeChi2(mtt, simPdf, mainCategory, *toyDataset, mainWorkspace);
        hNToyEvents->Fill(toyDataset->numEntries());
        std::cout << toyDataset->numEntries() << std::endl;
        hToyChi2->Fill(chi2_toy["combined"]);
        hToyMu->Fill(mu.getVal());;
        hToyMuErr->Fill(mu.getError());;
        hToyMuPull->Fill(mu.getVal()/mu.getError());;

        delete fitResult;
      } // ntoys
      toyResFile->cd();
      hNToyEvents->Write();
      hToyChi2->Write();
      hToyMu->Write();
      hToyMuErr->Write();
      hToyMuPull->Write();
      toyResFile->Close();
    } // doBiasTest
  }
}

