#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <sys/wait.h>
#include <cmath>
#include <regex>
#include <unistd.h>

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
#include <TH1.h>

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

using namespace RooFit;

#include <list>
#include "tdrstyle.C"

#include "Utils.h"
#include "Functions.h"

#include <tclap/CmdLine.h>
#include <json/json.h>

#include <TString.h>

std::vector<std::string> JEC;
bool FIT_WARNINGS = false;
int FIT_ERROR_LEVEL = -1;

std::string base_path = "";

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void fritSignal(TChain* chain, const std::string& jecType, const std::string& jer, const std::string& pu, const std::string& pdf, const std::string& configFile, int massZprime, int btag, bool saveWorkspace);

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

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Frit Zprime Signal", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> jecArg("", "jec", "Run the frit for this specific jec.", false, "nominal", "string");
    TCLAP::ValueArg<std::string> jerArg("", "jer", "Run the frit for this specific jer.", false, "nominal", "string");
    TCLAP::ValueArg<std::string> puArg("", "pileup", "Run the frit for this specific pileup syst.", false, "nominal", "string");
    TCLAP::ValueArg<std::string> pdfArg("", "pdf", "Run the frit for this specific pdf syst.", false, "nominal", "string");
    TCLAP::ValueArg<int> massArg("m", "mass", "Zprime mass", true, 750, "integer");
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", true, 2, "int");

    cmd.add(jecArg);
    cmd.add(jerArg);
    cmd.add(puArg);
    cmd.add(pdfArg);
    cmd.add(massArg);
    cmd.add(btagArg);

    TCLAP::ValueArg<std::string> configFileArg("", "config-file", "Fit configuration file", false, "frit_pdf.json", "string", cmd);

    TCLAP::SwitchArg saveWorkspaceArg("", "save-workspace", "Save the workspace for redoing plot after", cmd);

    cmd.parse(argc, argv);

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    // Format tree name
    TString treeName = TString::Format("dataset_%dbtag", btagArg.getValue());

    TChain * chain = loadChain(inputFiles, treeName.Data());

    std::string jec = jecArg.getValue();
    if (jec != "nominal" && jec != "JECup" && jec != "JECdown") {
      std::cerr << "--jec can only be 'nominal', 'JECup' or 'JECdown'" << std::endl;
      exit(1);
    }

    std::string pu = puArg.getValue();
    if (pu != "nominal" && pu != "up" && pu != "down") {
      std::cerr << "--pileup can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    std::string pdf = pdfArg.getValue();
    if (pdf != "nominal" && pdf != "up" && pdf != "down") {
      std::cerr << "--pdf can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    std::string jer = jerArg.getValue();
    if (jer != "nominal" && jer != "up" && jer != "down") {
      std::cerr << "--jer can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    if (pu != "nominal" && jec != "nominal" && pdf != "nominal") {
      std::cerr << "Please set --pileup, --jec and --pdf separately" << std::endl;
    }

    fritSignal(chain, jec, jer, pu, pdf, configFileArg.getValue(), massArg.getValue(), btagArg.getValue(), saveWorkspaceArg.getValue());

    delete chain;

  } catch (TCLAP::ArgException& e) {
  }
}

void saveParameter(int mass, const std::string& systType, int btag, const std::string& name, double integral, double integral_error) {

  FILE* lock = fopen((base_path + "/frit_efficiencies.lock").c_str(), "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  if (fileExists(base_path + "/frit_efficiencies.json")) {
    std::ifstream file(base_path + "/frit_efficiencies.json");
    reader.parse(file, root);
    file.close();
  }

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  Json::Value& node = root[getAnalysisUUID()][strMass][btagStr][systType][name];
  node["integral"] = integral;
  node["error"] = integral_error;

  FILE* fd = fopen((base_path + "/frit_efficiencies.json").c_str(), "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

void saveCombinedChiSquare(int mass, const std::string& jecType, int btag, double chi2) {

  FILE* lock = fopen((base_path + "/frit_efficiencies.lock").c_str(), "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  if (fileExists(base_path + "/frit_efficiencies.json")) {
    std::ifstream file(base_path + "/frit_efficiencies.json");
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


  FILE* fd = fopen((base_path + "/frit_efficiencies.json").c_str(), "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

void drawHistograms(RooAbsCategoryLValue& categories, RooRealVar& observable, int mass, int nBins, RooWorkspace& workspace, int btag, const std::string& prefix, bool log) {

  std::cout << std::endl;
  std::cout << "Drawing histograms..." << std::endl;

  std::vector<std::shared_ptr<RooPlot>> plots;

  int n = categories.numTypes();
  int x = 2, y = (n + 0.5) / x;
  std::cout << "\t" << n << " categories. Dividing into " << x << "; " << y << std::endl;

  const int padWidth = 900;
  const int padHeight = 900;

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

    RooAbsData* dataset_positive = workspace.data(TString::Format("%s_positive_signal_data", category.c_str()));
    RooAbsData* dataset_negative = workspace.data(TString::Format("%s_negative_signal_data", category.c_str()));

    RooAbsPdf* pdf_positive = workspace.pdf(TString::Format("%s_positive_signal_pdf", category.c_str()));
    RooAbsPdf* pdf_negative = workspace.pdf(TString::Format("%s_negative_signal_pdf", category.c_str()));

    if (!pdf_positive || !pdf_negative)
      continue;

    dataset_positive->plotOn(plot, RooFit::Invisible());
    pdf_positive->plotOn(plot, RooFit::LineWidth(1));
    dataset_positive->plotOn(plot, RooFit::IgnoreEmptyBins());

    //dataset_negative->plotOn(plot, RooFit::Invisible());
    //pdf_negative->plotOn(plot, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), RooFit::LineWidth(1));

    dataset_negative->plotOn(plot, RooFit::Rescale(-1), RooFit::RefreshNorm(), RooFit::Invisible());
    pdf_negative->plotOn(plot, RooFit::Normalization(-1., RooAbsPdf::ScaleType::Relative), RooFit::LineColor(kRed), RooFit::LineWidth(1));
    dataset_negative->plotOn(plot, RooFit::Rescale(-1), RooFit::IgnoreEmptyBins());

    std::string leptonName = TString(category).Contains("muon", TString::kIgnoreCase) ? "muons" : "electrons";
    float binningSize = (observable.getBinning().highBound() - observable.getBinning().lowBound()) / (float) nBins;

    plot->SetXTitle(TString::Format("#font[132]{#font[12]{M_{t#bar{t}}} (GeV), %s}", leptonName.c_str()));
    plot->SetYTitle(TString::Format("#font[132]{Events/(%0.0f GeV)}", binningSize));
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

    double x_ = 0.53;
    double y_ = 0.88;

    if (mass > 1200) {
      x_ = 0.22;
      y_ = 0.25;
    }

    t.DrawLatex(x_, y_, "#font[42]{CMS preliminary}");
    t.DrawLatex(x_, y_ - 0.04, "#font[42]{simulation at #sqrt{s}=8 TeV}");
    t.DrawLatex(x_, y_ - 0.08, legendLabel);

    plots.push_back(std::shared_ptr<RooPlot>(plot));

    canvas->Update();
  }

  if (! log) {
    canvas->Print((prefix + "_fit.pdf").c_str());
    canvas->Print((prefix + "_fit.png").c_str());
  } else {
    canvas->Print((prefix + "_fit_log.pdf").c_str());
    canvas->Print((prefix + "_fit_log.png").c_str());
  }

  delete canvas;

  std::cout << std::endl;
}

std::vector<std::shared_ptr<TH1>> splitHistogram(TH1* h) {

  std::shared_ptr<TH1> positive_h(static_cast<TH1*>(h->Clone()));
  positive_h->SetName("positive_histo");
  std::shared_ptr<TH1> negative_h(static_cast<TH1*>(h->Clone()));
  negative_h->SetName("negative_histo");
  for (uint16_t i = 1; i <= (uint16_t) positive_h->GetNbinsX(); i++) {
    if (positive_h->GetBinContent(i) < 0) {
      positive_h->SetBinContent(i, 0);
      positive_h->SetBinError(i, 0);
    }

    if (negative_h->GetBinContent(i) < 0) {
      negative_h->SetBinContent(i, -1. * negative_h->GetBinContent(i));
    } else {
      negative_h->SetBinContent(i, 0);
      negative_h->SetBinError(i, 0);
    }
  }

  return {positive_h, negative_h};
}

RooDataSet* binnedToUnbinnedDataset(RooDataHist& in, RooRealVar& observable, RooRealVar& weight, const std::string& name) {

  RooDataSet* ds = new RooDataSet(name.c_str(), name.c_str(), RooArgSet(observable, weight), RooFit::WeightVar(weight));
  for (int i = 0; i < in.numEntries(); i++) {
    const RooArgSet* set = in.get(i);
    ds->add(*set, in.weight(*set));
  }

  return ds;
}

float shiftUp(TH1* h) {

  // Shift the histogram bins up to make sure every bins are positive

  int minBin = h->GetMinimumBin();

  float minError = h->GetBinError(minBin);
  float minValue = h->GetBinContent(minBin);

  if (minValue > 0)
    return 0;

  // Substract twice the error to leave room for PDF conversion
  float min = minValue - 2 * minError;

  min *= -1;

  for (uint16_t i = 1; i <= (uint16_t) h->GetNbinsX(); i++) {
    float error = h->GetBinError(i);
    h->SetBinContent(i, h->GetBinContent(i) + min);

    // Don't change error on each bins
    h->SetBinError(i, error);
  }

  return min;
}

void fritSignal(TChain* chain, const std::string& jecType, const std::string& jer, const std::string& pu, const std::string& pdfSyst, __attribute__((unused)) const std::string& configFile, int massZprime, int btag, __attribute__((unused)) bool saveWorkspace) {

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
  Float_t minmTT = 325;
  Float_t maxmTT = 1000;
  Int_t nBins = (maxmTT - minmTT) / (25.);

  RooRealVar mtt("mtt", "mtt", minmTT, maxmTT, "GeV/c^2");
  mtt.setBins(nBins);

  RooRealVar weight("weight", "weight", -10, 10);

  RooCategory whichLepton("lepton_type", "lepton_type");
  whichLepton.defineType("electron", 11);
  whichLepton.defineType("muon", 13);

  // Get dataset from external file
  //RooDataSet *dataset = RooDataSet::read(file.c_str(), RooArgList(Mtt_KF_reco, whichLepton));
  RooDataSet *dataset = new RooDataSet("dataset", "dataset", RooArgSet(mtt, whichLepton, weight), Import(*chain), WeightVar(weight));

  if (dataset->numEntries() == 0) {
    std::cerr << "[" << getpid() << "] ERROR: No entry in dataset" << std::endl;
    return;
  } else {
    std::cout << "[" << getpid() << "] " << dataset->numEntries() << " entries loaded" << std::endl;
  }

  std::string analysisName = getAnalysisName();
  std::string analysisUUID = getAnalysisUUID();

  base_path = "./analysis/" + analysisUUID;

  TString systPrefix;
  if (jecType == "nominal" && pu == "nominal" && jer == "nominal" && pdfSyst == "nominal")
    systPrefix = "nominal";
  else if (jecType != "nominal")
    systPrefix = jecType.c_str();
  else if (jer != "nominal") {
    systPrefix = TString::Format("JER%s", jer.c_str());
  } else if (pu != "nominal") {
    if (pu == "up")
      systPrefix = "puUp";
    else
      systPrefix = "puDown";
  } else if (pdfSyst != "nominal") {
    if (pdfSyst == "up")
      systPrefix = "pdfUp";
    else
      systPrefix = "pdfDown";
  }

  TString prefix = TString::Format("%s-%s%d_%s_%d_btag", systPrefix.Data(), getAnalysisPrefix(), massZprime, analysisName.c_str(), btag);

  // Save our signal functions into the workspace
  // We will need it for the fit on data
  TString workspaceFile = base_path + "/frit/" + prefix + "_workspace.root";
  RooWorkspace workspace("w", "Frit signal workspace");

  TIterator* it = whichLepton.typeIterator();
  const RooCatType* catType = nullptr;

  // Store integral for each category
  std::map<
    std::string, // Name
    std::pair<float, float> // Integral + error
  > integrals;

  while ((catType = static_cast<RooCatType*>(it->Next()))) {
    if (! catType)
      continue;

    std::string name = catType->GetName();

    RooDataSet* reducedDataset = nullptr;
    // Reduce dataset to keep only data for this category
    whichLepton.setIndex(catType->getVal());
    std::string cutFormula = buildCutFormula(whichLepton);
    reducedDataset = static_cast<RooDataSet*>(dataset->reduce(cutFormula.c_str()));

    // Bin dataset to create histograms
    std::shared_ptr<RooDataHist> binned_dataset(reducedDataset->binnedClone());

    // Create histogram from binned dataset
    TH1* h = binned_dataset->createHistogram(name.c_str(), mtt);

    workspace.factory(TString::Format("%s_integral[%.5f]", name.c_str(), h->Integral()));

    // Split histogram in two part: one positive, and one negative
    auto splitted_histograms = splitHistogram(h);

    workspace.factory(TString::Format("%s_positive_integral[%.5f]", name.c_str(), splitted_histograms[0]->Integral()));
    workspace.factory(TString::Format("%s_negative_integral[%.5f]", name.c_str(), splitted_histograms[1]->Integral()));

    double positive_integral_error = 0;
    double positive_integral = splitted_histograms[0]->IntegralAndError(splitted_histograms[0]->GetXaxis()->GetFirst(), splitted_histograms[0]->GetXaxis()->GetLast(), positive_integral_error);

    double negative_integral_error = 0;
    double negative_integral = splitted_histograms[1]->IntegralAndError(splitted_histograms[1]->GetXaxis()->GetFirst(), splitted_histograms[1]->GetXaxis()->GetLast(), negative_integral_error);

    double integral = positive_integral + negative_integral;
    double integral_error = std::sqrt(positive_integral_error * positive_integral_error + negative_integral_error * negative_integral_error);

    integrals[name] = std::make_pair(integral, integral_error);

    if (false) {
      // This histogram constains negative events.
      // Add a given number to each bin to make sure
      // all bins content are positive
      float factor = shiftUp(h);
      workspace.factory(TString::Format("%s_shift_factor[%.5f]", name.c_str(), factor));

      std::shared_ptr<RooDataHist> data = std::make_shared<RooDataHist>(TString::Format("%s_signal_data", name.c_str()), "", mtt, RooFit::Import(*h));

      RooRealVar weight_("__weight__", "__weight__", 0, 0, 100000);
      std::shared_ptr<RooDataSet> unbinnedData(binnedToUnbinnedDataset(*data.get(), mtt, weight_, TString::Format("%s_signal_unbinned_data", name.c_str()).Data()));

      std::shared_ptr<RooKeysPdf> keys_pdf = std::make_shared<RooKeysPdf>(std::string("signal_" + name).c_str(), "Keys pdf for signal", mtt, *unbinnedData, RooKeysPdf::MirrorRight, 0.4);

      workspace.import(*data);
      workspace.import(*unbinnedData);
      workspace.import(*keys_pdf);

    } else {

      auto TH1ToPdf = [&](TH1* histogram, const std::string& p) -> std::pair<std::shared_ptr<RooHistPdf>, std::shared_ptr<RooDataHist>> {
        // First, create a RooDataHist from the TH1
        std::shared_ptr<RooDataHist> data = std::make_shared<RooDataHist>(TString::Format("%s_signal_data", p.c_str()), "", mtt, RooFit::Import(*histogram));
        std::shared_ptr<RooHistPdf> pdf = std::make_shared<RooHistPdf>(TString::Format("%s_signal_pdf", p.c_str()), "", RooArgSet(mtt), *data, 2);

        return std::make_pair(pdf, data);
      };

      TString p = TString::Format("%s_positive", name.c_str());
      auto pdf = TH1ToPdf(splitted_histograms[0].get(), p.Data());

      workspace.import(*(pdf.first));
      workspace.import(*(pdf.second));

      p = TString::Format("%s_negative", name.c_str());
      pdf = TH1ToPdf(splitted_histograms[1].get(), p.Data());

      workspace.import(*(pdf.first));
      workspace.import(*(pdf.second));
    }

    delete h;
  }

  workspace.writeToFile(workspaceFile);

  drawHistograms(whichLepton, mtt, massZprime, nBins, workspace, btag, std::string(base_path + "/frit/" + prefix), false);
  //drawHistograms(whichLepton, mtt, massZprime, true, nBins, workspace, btag, std::string(base_path + "/frit/" + prefix), true);

  std::cout << "Workspace saved in " << workspaceFile << std::endl;

  // Save number of events for each category for efficiency computations
  it = whichLepton.typeIterator();
  catType = nullptr;
  while ((catType = static_cast<RooCatType*>(it->Next()))) {
    std::string name = catType->GetName();
    saveParameter(massZprime, systPrefix.Data(), btag, name, integrals[name].first, integrals[name].second);
  }

  //if (saveWorkspace) {
    //// Save fitted pdf and datasets in order to redo some plots
    //RooWorkspace workspace("w");

    //fitResult->SetName("fit_results");
    //workspace.import(*fitResult);

    //it = whichLepton.typeIterator();
    //type = nullptr;
    //while ((type = static_cast<RooCatType*>(it->Next()))) {

      //std::string category = type->GetName();
      //std::string cleanedCategory = TString(type->GetName()).ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", "").Data(); // For workspace names

      //whichLepton.setLabel(category.c_str());
      //std::string cut = buildCutFormula(whichLepton);

      //std::string workspaceName = "global_pdf_" + cleanedCategory;

      //std::string leptonName = TString(category).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron";
      //std::string leptonNameShort = TString(category).Contains("muon", TString::kIgnoreCase) ? "mu" : "e";

      //TString workspace_suffix = TString::Format("%s", leptonNameShort.c_str());

      //RooAbsData* dataset_reduced = dataset->reduce(RooArgSet(mtt), cut.c_str());
      //workspace.import(*dataset_reduced, RooFit::Rename(TString::Format("data_obs_%s", workspace_suffix.Data())));

      //TString name = TString::Format("n_signal_%s", category.c_str());

      //// Get number of fitted signal events
      //RooAbsPdf* globalPdf = globalPdfs[leptonName].get();
      //workspace.import(
          //*globalPdf,
          //RooFit::RecycleConflictNodes(),
          //RooFit::RenameVariable(leptonName.c_str(), TString::Format("global_pdf_%s", workspace_suffix.Data()))
          //);

      //TString outputFileName = TString::Format("zprime_signal_%d_workspace_after_frit_%d_btag.root", massZprime, btag);
      //workspace.writeToFile(outputFileName);
    //}
  //}
}
