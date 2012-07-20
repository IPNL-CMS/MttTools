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

#include "Utils.h"
#include <tclap/CmdLine.h>
#include <json/json.h>

#include "SignalFunctions.h"
#include "Functions.h"

#define SAFE_DELETE(p) { delete p; p = NULL; }

using namespace RooFit;

#define NO_SYST // No systematics

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

  template<typename T>
void deleteVector(std::vector<T*>& vector)
{
  for (typename std::vector<T*>::iterator it = vector.begin(); it != vector.end(); ++it)
  {
    delete *it;
  }
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

void fitMtt(int massZprime, bool fit, string bkgfit_str, bool doLikScan, bool writeRootFile, bool writeTxtFile, bool saveFigures, bool doLimitCurve, int nToyExp, bool doLikScanInToys, int index, bool doDiscCurve, string syst_str, string systCBsign, string systCB, bool bkgOnly, bool muonsOnly, int btag);

std::string BASE_PATH;
std::string OUTPUT_PATH;
std::string EFF_FILE;
bool DO_SYST_COMPUTATION = false;
bool SAVE_SIGMA = false;
bool ONLY_LUMI_SYST = false;
bool VERBOSE = false;
bool BATCH_MODE = false;

int main(int argc, char** argv)
{
  try
  {
    TCLAP::CmdLine cmd("Fit tt mass spectrum", ' ', "0.1");

    TCLAP::ValueArg<int> massArg("m", "mass", "Zprime mass", true, 750, "integer");
    TCLAP::SwitchArg fitArg("", "no-fit", "Don't do the reference fit", true);
    TCLAP::SwitchArg doLikScanArg("", "scan", "Do the likelihood scan", false);
    TCLAP::SwitchArg writeRootArg("", "no-root-files", "Don't write root files", true);
    TCLAP::SwitchArg writeTxtArg("", "no-text-files", "Don't write text files", true);
    TCLAP::SwitchArg saveFiguresArg("", "no-figs", "Don't save figures", true);
    TCLAP::SwitchArg doLimitCurveArg("", "limit-curve", "Do the limit curve", false);
    TCLAP::ValueArg<int> nToyArg("", "toys", "Number of toys exp.", false, 1, "integer");
    TCLAP::SwitchArg doLikInToyArg("", "no-scan-in-toys", "Don't do the likelihood scan in toys", true);
    TCLAP::ValueArg<int> indexArg("", "index", "Index", false, 1, "integer");
    TCLAP::SwitchArg doDiscCurveArg("", "disc-curve", "Do the disc curve", false);

    std::vector<std::string> jec;
    jec.push_back("nominal");
    jec.push_back("JECup");
    jec.push_back("JECdown");
    TCLAP::ValuesConstraint<std::string> allowedVals(jec);
    TCLAP::ValueArg<std::string> systArg("", "syst", "Systematic", false, "nominal", &allowedVals);

    std::vector<std::string> sign;
    sign.push_back("up");
    sign.push_back("down");
    TCLAP::ValuesConstraint<std::string> allowedSign(sign);

    TCLAP::ValueArg<std::string> systSignArg("", "syst-sign", "?", false, "up", &allowedSign);
    TCLAP::ValueArg<std::string> systCBArg("", "systCB", "?", false, "none", "string");

    cmd.add(massArg);
    cmd.add(fitArg);
    cmd.add(doLikScanArg);
    cmd.add(writeRootArg);
    cmd.add(writeTxtArg);
    cmd.add(saveFiguresArg);
    cmd.add(doLimitCurveArg);
    cmd.add(nToyArg);
    cmd.add(doLikInToyArg);
    cmd.add(indexArg);
    cmd.add(doDiscCurveArg);
    cmd.add(systArg);
    cmd.add(systSignArg);
    cmd.add(systCBArg);

    TCLAP::ValueArg<std::string> pathArg("", "path", "Folder where to load files", false, "./", "string", cmd);
    TCLAP::ValueArg<std::string> outputPathArg("", "output-path", "Folder where output files are stored", false, "./", "string", cmd);
    TCLAP::ValueArg<std::string> effArg("", "eff-file", "File where efficiences are stored (JSON format)", false, "efficiencies.json", "string", cmd);
    TCLAP::SwitchArg doSystComputationArg("", "syst-computation", "Set this flag if you are running me to compute systematic", cmd);
    TCLAP::SwitchArg saveSigmaArg("", "save-sigma-ref", "Set this flag if you want to save the computed sigma as the ref one", cmd);
    TCLAP::SwitchArg bkgOnlyArg("", "bkg-only", "Fit using background model only", cmd);
    TCLAP::SwitchArg onlyMuonArg("", "muons-only", "Compute limits using only semi-mu data", cmd);
    TCLAP::SwitchArg onlyLumiSystArg("", "only-lumi-syst", "Only use luminosity error for systematics", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);
    TCLAP::SwitchArg verboseArg("v", "verbose", "Verbose mode", cmd);
    TCLAP::SwitchArg batchArg("", "batch", "Run in batch mode", cmd);
    TCLAP::ValueArg<std::string> fitConfigFileArg("", "config-file", "Configuration file name containing fit parameters", false, "fit_pdf_faltb.json", "string", cmd);

    cmd.parse(argc, argv);

    BASE_PATH = pathArg.getValue();
    OUTPUT_PATH = outputPathArg.getValue();

    BASE_PATH.erase(std::remove(BASE_PATH.begin(), BASE_PATH.end(), '\"'), BASE_PATH.end());
    OUTPUT_PATH.erase(std::remove(OUTPUT_PATH.begin(), OUTPUT_PATH.end(), '\"'), OUTPUT_PATH.end());

    EFF_FILE = effArg.getValue();
    DO_SYST_COMPUTATION = doSystComputationArg.getValue();
    SAVE_SIGMA = saveSigmaArg.getValue();

    ONLY_LUMI_SYST = onlyLumiSystArg.getValue();

    VERBOSE = verboseArg.getValue();
    BATCH_MODE = batchArg.getValue();


    fitMtt(massArg.getValue(), fitArg.getValue(), fitConfigFileArg.getValue(), doLikScanArg.getValue(), writeRootArg.getValue(), writeTxtArg.getValue(),
        saveFiguresArg.getValue(), doLimitCurveArg.getValue(), nToyArg.getValue(), doLikInToyArg.getValue(), indexArg.getValue(),
        doDiscCurveArg.getValue(), systArg.getValue(), systSignArg.getValue(), systCBArg.getValue(), bkgOnlyArg.getValue(), onlyMuonArg.getValue(), btagArg.getValue());

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

  if (! root.isMember(strMass))
  {
    //std::cerr << "ERROR: mass '" << mass << "' not found in JSON file. Exiting." << std::endl;
    //exit(1);
    std::cout << "WARNING: Systematics error not found. Setting to 0" << std::endl;
    return;
  }

  Json::Value massNode = root[strMass][btagStr];
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

void saveSigma(int mass, int btag, double sigma, double chi, double events, RooFitResult* fitRes)
{

  FILE* lock = fopen("sigma_reference.lock", "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  std::ifstream file("sigma_reference.json");
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  root[strMass][btagStr]["sigma"] = sigma;
  root[strMass][btagStr]["chi2"] = chi;
  root[strMass][btagStr]["events"] = events;
  root[strMass][btagStr]["fit_covQual"] = fitRes->covQual();
  root[strMass][btagStr]["fit_status"] = fitRes->status();


  std::ofstream ofile;
  ofile.open("sigma_reference.json", std::ios::out | std::ios::trunc);
  Json::StyledWriter writer;
  ofile << writer.write(root);
  ofile.close();

  fclose(lock);
}

void saveSystematicParameter(int mass, int btag, const std::string& type, const std::string& param, const std::string& subparam, double nEvents, double sigma, double chi2, RooFitResult* fitRes)
{

  FILE* lock = fopen("systematics_parameters.lock", "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  std::ifstream file("systematics_parameters.json");
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  if (subparam.length() > 0)
  {
    root[strMass][btagStr][type][param]["sigma"][subparam] = sigma;
    root[strMass][btagStr][type][param]["events"][subparam] = nEvents;
    root[strMass][btagStr][type][param]["chi2"][subparam] = chi2;
    root[strMass][btagStr][type][param]["fit_covQual"][subparam] = fitRes->covQual();
    root[strMass][btagStr][type][param]["fit_status"][subparam] = fitRes->status();
  }
  else
  {
    root[strMass][btagStr][type][param]["sigma"] = sigma;
    root[strMass][btagStr][type][param]["events"] = nEvents;
    root[strMass][btagStr][type][param]["chi2"] = chi2;
    root[strMass][btagStr][type][param]["fit_covQual"] = fitRes->covQual();
    root[strMass][btagStr][type][param]["fit_status"] = fitRes->status();
  }

  //TODO: reswitch to ofstream
  FILE* fd = fopen("systematics_parameters.json", "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

#define MIN_LIKELIHOOD -200.
#define MAX_LIKELIHOOD 300.

/**
 * Compute the likelihood range
 */
void computeLikeloodRange(const RooRealVar& v, double& xLow, double& xHigh)
{

  double meanValue = v.getVal();
  xLow = meanValue - 300;
  xHigh = meanValue + 600;
  // Ensure xLow is -300 at least
  xLow = std::min(MIN_LIKELIHOOD, xLow);
  // Ensure xHigh is 300 at least
  xHigh = std::max(MAX_LIKELIHOOD, xHigh);
}

double getXForBin(const int bin, const double center, const int steps)
{
  return center + bin * steps;
}

struct LikelihoodResults
{
  TH1* likscan;
  TH1* pdfscan;
  TH1* pdfscan_wsyst;
  TH1* pdfscan_wsyst_cut;

  double scan_limit;
  double scan_wsyst_limit;
  double scan_wsyst_cut_limit;

  LikelihoodResults():
    likscan(NULL), pdfscan(NULL), pdfscan_wsyst(NULL), pdfscan_wsyst_cut(NULL),
    scan_limit(0.), scan_wsyst_limit(0.), scan_wsyst_cut_limit(0.)
  {}

  void release()
  {
    SAFE_DELETE(likscan);
    SAFE_DELETE(pdfscan);
    SAFE_DELETE(pdfscan_wsyst);
    SAFE_DELETE(pdfscan_wsyst_cut);
  }

  ~LikelihoodResults()
  {
    release();
  }
};

void saveLikelihoodResults(const int mass, int btag, const LikelihoodResults& results, const double denominator)
{
  FILE* lock = fopen("likelihood_scan.lock", "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  std::ifstream file("likelihood_scan.json");
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  root[strMass][btagStr]["scan_limit"] = results.scan_limit / denominator;
  root[strMass][btagStr]["scan_wsyst_limit"] = results.scan_wsyst_limit / denominator;
  root[strMass][btagStr]["scan_wsyst_cut_limit"] = results.scan_wsyst_cut_limit / denominator;

  std::ofstream ofile;
  ofile.open("likelihood_scan.json", std::ios::out | std::ios::trunc);
  Json::StyledWriter writer;
  ofile << writer.write(root);
  ofile.close();

  fclose(lock);
}

void doLikelihoodScan(RooDataSet& dataset, RooAbsPdf& pdf, RooRealVar& observable, const double mass, const double minNLL, const int steps, const double systError, LikelihoodResults& results)
{

  RooFitResult* fitResults = NULL;

  double xLow = 0, xHigh = 0;
  computeLikeloodRange(observable, xLow, xHigh);
  int nLHBins = (int)(xHigh - xLow) / steps;

  std::cout << "Likehood from " << xLow << " to " << xHigh << "; bins: " << nLHBins << std::endl;

  results.likscan = new TH1D("likscan", "likscan", nLHBins, xLow, xHigh);
  results.pdfscan = new TH1D("pdfscan", "pdfscan", nLHBins, xLow, xHigh);

  // Take into account gaussian smearing
  double xSmearedHigh = xHigh + 3 * systError;
  int nLHSmearedBins = (int)(xSmearedHigh - xLow) / steps;
  int nLHCuttedBins = (int) xSmearedHigh / steps;
  std::cout << "Likehood smearing from " << xLow << " to " << xSmearedHigh << "; bins: " << nLHSmearedBins << std::endl;
  results.pdfscan_wsyst = new TH1D("pdfscan_wsyst", "pdfscan_wsyst", nLHSmearedBins, xLow, xSmearedHigh);
  results.pdfscan_wsyst_cut = new TH1D("pdfscan_wsyst_cut", "pdfscan_wsyst_cut", nLHCuttedBins, 0, xSmearedHigh);

  double center = observable.getVal();

  int nLHLowerBins = (int)(center - xLow) / steps;
  int nLHUpperBins = (int)(xHigh - center) / steps;

  // std::cout << "Lower: " << nLHLowerBins << std::endl;
  // std::cout << "Upper: " << nLHUpperBins << std::endl;

  //const double maxNLL = -1 * log(1e-4);
  double DNll = 0.;
  bool fill = false;
  bool done = false;
  const double minExpNllValue = 3.72665e-06;

  // First, create the nLL. We only need to create it once, so do it right now
  RooAbsReal* nll = pdf.createNLL(dataset, RooFit::Optimize(1)/*, RooFit::NumCPU(2)*/);

  std::cout << "[M=" << mass << "] Scanning negative values from " << center << " to " << getXForBin(-1 * nLHLowerBins, center, steps) << std::endl;

  for (int s = -1; s >= (-1 * nLHLowerBins); s--)
  {
    if (abs(s) % 10 == 0)
    {
      cout << "[M=" << mass << "] Likelihood scan step " << (-1 * s) << " over " << nLHBins << endl;
    }
    double x = getXForBin(s, center, steps);

    // std::cout << x << std::endl;

    observable.setVal(x);
    observable.setConstant(true);

    RooMinuit* minimizer = new RooMinuit(*nll);
    minimizer->setEvalErrorWall(1);
    minimizer->optimizeConst(1);
    minimizer->setStrategy(2);
    minimizer->migrad();
    fitResults = minimizer->save();

    if ((fitResults->minNll() - minNLL) < -20. || (fitResults->minNll() - minNLL) > 10000.)
    {
      cout << "[M=" << mass << "] WARNING : Delta(NLL) = " << fitResults->minNll() - minNLL << endl;
      DNll = 0.;
      fill = false;
    }
    else
    {
      DNll = fitResults->minNll() - minNLL;
      // std::cout << "nll: "  << DNll << std::endl;
      fill = true;
    }

    bool converged = fitResults->status() == 0 && fitResults->covQual() == 3;
    if (! converged)
    {
      std::cout << "[M=" << mass << "] WARNING: Fit has not converged!" << std::endl;
    }

    if (fill)
    {
      const double e = exp(-1. * DNll);
      results.likscan->SetBinContent(results.likscan->FindBin(x), DNll);
      results.pdfscan->SetBinContent(results.pdfscan->FindBin(x), e);
      // Negative part: exp shoud be lesser than 5e-6 AND x lesser than MIN_LIKELIHOOD
      done = (e <= minExpNllValue && x <= MIN_LIKELIHOOD);
    }

    SAFE_DELETE(fitResults);
    SAFE_DELETE(minimizer);

    if (done)
    {
      std::cout << "Iterations stopped because likelihood is big enough" << std::endl;
      break;
    }
  }

  std::cout << "[M=" << mass << "] Scanning positive values from " << center << " to " << getXForBin(nLHUpperBins, center, steps) << std::endl;

  done = false;
  for (int s = 0 ; s < nLHUpperBins; s++)
  {
    if (s % 10 == 0)
    {
      std::cout << "[M=" << mass << "] Likelihood scan step " << s + nLHLowerBins << " over " << nLHBins << std::endl;
    }

    double x = getXForBin(s, center, steps);

    observable.setVal(x);
    observable.setConstant(true);

    RooMinuit* minimizer = new RooMinuit(*nll);
    minimizer->setEvalErrorWall(1);
    minimizer->optimizeConst(1);
    minimizer->setStrategy(2);
    minimizer->migrad();
    fitResults = minimizer->save();

    if ((fitResults->minNll() - minNLL) < -20. || (fitResults->minNll() - minNLL) > 10000.)
    {
      cout << "[M=" << mass << "] WARNING : Delta(NLL) = " << fitResults->minNll() - minNLL << endl;
      DNll = 0.;
      fill = false;
    }
    else
    {
      DNll = fitResults->minNll() - minNLL;
      //std::cout << "nll: "  << DNll << std::endl;
      fill = true;
    }

    bool converged = fitResults->status() == 0 && fitResults->covQual() == 3;
    if (! converged)
    {
      std::cout << "[M=" << mass << "] WARNING: Fit has not converged!" << std::endl;
    }

    if (fill)
    {
      const double e = exp(-1. * DNll);
      results.likscan->SetBinContent(results.likscan->FindBin(x), DNll);
      results.pdfscan->SetBinContent(results.pdfscan->FindBin(x), e);
      // Positive part: exp shoud be lesser than 5e-6 AND x greater than MAX_LIKELIHOOD
      done = (e <= minExpNllValue && x >= MAX_LIKELIHOOD);
    }

    SAFE_DELETE(fitResults);
    SAFE_DELETE(minimizer);

    if (done)
    {
      std::cout << "Iterations stopped because likelihood is big enough" << std::endl;
      break;
    }
  }

  delete nll;

  // No systematic errors?
  if (systError < 1e-10) {

    for (double x = xLow; x < xHigh; x += steps) {
      double weight = results.pdfscan->GetBinContent(results.pdfscan->FindBin(x));
      results.pdfscan_wsyst->SetBinContent(results.pdfscan_wsyst->FindBin(x), weight);
      if (x >= 0) {
        results.pdfscan_wsyst_cut->SetBinContent(results.pdfscan_wsyst_cut->FindBin(x), weight);
      }
    }

  } else {
    // Perform smearing of pdfscan histogram

    // A seed of 0 means a random one
    TRandom3 *myRandom = new TRandom3(0);

    for (double x = xLow; x < xHigh; x += steps)
    {
      for (int i = 0; i < 50000; i++)
      {
        double xvalue = myRandom->Gaus(x, systError);
        results.pdfscan_wsyst->Fill(xvalue, results.pdfscan->GetBinContent(results.pdfscan->FindBin(x)));
      }
    }

    // Fill histogram with cut on x > 0. Don't do that when smearing, it'll cause weird effects around x = 0
    for (double x = xLow; x < xHigh; x += steps) {
      if (x >= 0) {
        double weight = results.pdfscan->GetBinContent(results.pdfscan->FindBin(x));
        results.pdfscan_wsyst_cut->SetBinContent(results.pdfscan_wsyst_cut->FindBin(x), weight);
      }
    }

    SAFE_DELETE(myRandom);

  }

  //calculate different 95% limits
  double areascan = 0.95;

  results.pdfscan->ComputeIntegral();
  results.pdfscan->Scale(1. / results.pdfscan->Integral());
  results.pdfscan->GetQuantiles(1, &results.scan_limit, &areascan);

  results.pdfscan_wsyst->ComputeIntegral();
  results.pdfscan_wsyst->Scale(1. / results.pdfscan_wsyst->Integral());
  results.pdfscan_wsyst->GetQuantiles(1, &results.scan_wsyst_limit, &areascan);

  results.pdfscan_wsyst_cut->ComputeIntegral();
  results.pdfscan_wsyst_cut->Scale(1. / results.pdfscan_wsyst_cut->Integral());
  results.pdfscan_wsyst_cut->GetQuantiles(1, &results.scan_wsyst_cut_limit, &areascan);
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
}

std::string buildCutFormula(RooAbsCategoryLValue& categories) {

  std::stringstream ss;
  RooSuperCategory* foo = dynamic_cast<RooSuperCategory*>(&categories);
  if (foo) {
    const RooArgSet& parentCategories = foo->inputCatList();
    TIterator *it2 = parentCategories.createIterator();
    RooAbsCategory* parentCategory = nullptr;
    bool first = true;
    while ((parentCategory = static_cast<RooAbsCategory*>(it2->Next()))) {
      ss << (first ? "" : " && ") << parentCategory->GetName() << "==" << parentCategory->getIndex();
      first = false;
    }
  } else {
    ss << categories.GetName() << "==" << categories.getIndex();
  }

  return ss.str();
}

void drawHistograms(RooAbsCategoryLValue& categories, RooRealVar& observable, int nBins, RooDataSet& dataset, RooSimultaneous& simPdfs, std::map<std::string, std::shared_ptr<BaseFunction>>& backgroundPdfs, int btag, bool savePlots, const std::string& prefix, const std::string& suffix, bool drawSignal, bool log, TFile* outputFile, bool drawOnlyData = false) {

  std::vector<std::shared_ptr<RooPlot>> plots;

  int n = categories.numTypes();
  int x = std::min(n, 2), y = (int) ceil((float) n / (float) x);
  std::cout << n << " categories. Dividing into " << x << "; " << y << std::endl;

  const int padWidth = 900;
  const int padHeight = 900;
  const float LUMI = 5.1;

  const int canvasWidth = padWidth * x;
  const int canvasHeight = padHeight * y;

  TCanvas *canvas = new TCanvas("canvas", "mTT fit", canvasWidth, canvasHeight);
  canvas->Divide(x, y);
  int currentPad = 1;

  Roo1DTable* table = nullptr;
  if (outputFile && !log) {
    table = dataset.table(categories);
  }

  TIterator* it = categories.typeIterator();
  RooCatType* type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    canvas->cd(currentPad++);

    RooPlot* plot = observable.frame(nBins);
    plot->SetTitle(""); //FIXME

    std::string category = type->GetName();
    std::string cleanedCategory = TString(category).ReplaceAll(";", "_").Data(); // Root does not like ';' in names

    categories = category.c_str();

    std::string cut = buildCutFormula(categories);

    RooDataSet* subData = static_cast<RooDataSet*>(dataset.reduce(cut.c_str()));
    subData->plotOn(plot);

    if (! drawOnlyData) {
      simPdfs.plotOn(plot, Slice(categories), ProjWData(*subData), Components(backgroundPdfs[category]->getPdf()), LineStyle(kDashed), LineColor(kRed), LineWidth(2));

      if (drawSignal)
        simPdfs.plotOn(plot, Slice(categories), ProjWData(*subData), LineColor(kBlue), LineWidth(2));
    }

    delete subData;

    std::string leptonName = TString(category).Contains("muon", TString::kIgnoreCase) ? "muons" : "electrons";
    float binningSize = (observable.getBinning().highBound() - observable.getBinning().lowBound()) / (float) nBins;

    plot->SetXTitle(TString::Format("#font[132]{#font[12]{m_{TT}} (GeV/#font[12]{c}^{2}), %s}", leptonName.c_str()));
    plot->SetYTitle(TString::Format("#font[132]{Events/(%0.2f GeV/#font[12]{c}^{2})}", binningSize));
    plot->SetTitleOffset(1.42, "Y");

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
    boost::regex regex("([0-9]+)-btag");
    boost::smatch regexResults;

    if (boost::regex_search(category, regexResults, regex)) {
      std::string result(regexResults[1].first, regexResults[2].second);
      btag = atoi(result.c_str());
    }

    TString btagLabel = "";
    if (btag == 2)
      btagLabel = "#geq 2 b-tags";
    else
      btagLabel = TString::Format("%d b-tag", btag);

    std::string leptonShortcutName = TString(category).Contains("muon", TString::kIgnoreCase) ? "#mu" : "e";

    TString legendLabel = TString::Format("#font[42]{%s, #geq 4 jets, %s}", leptonShortcutName.c_str(), btagLabel.Data());

    t.DrawLatex(0.53, 0.88, "#font[42]{CMS preliminary}");
    t.DrawLatex(0.53, 0.84, TString::Format("#font[42]{%0.2f fb^{-1} at #sqrt{s}=7 TeV}", LUMI));
    t.DrawLatex(0.53, 0.80, legendLabel);

    if (outputFile && !log) {
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
  }

  if (savePlots) {
    if (! log) {
      canvas->Print((OUTPUT_PATH + prefix + "_fitRes_" + suffix + ".pdf").c_str());
      canvas->Print((OUTPUT_PATH + prefix + "_fitRes_" + suffix + ".png").c_str());
    } else {
      canvas->Print((OUTPUT_PATH + prefix + "_fitRes_" + suffix + "_log.pdf").c_str());
      canvas->Print((OUTPUT_PATH + prefix + "_fitRes_" + suffix + "_log.png").c_str());
    }
  }
  
  delete table;
  delete canvas;
}

std::map<std::string, float> computeChi2(RooRealVar& observable, const RooSimultaneous& simPdf, RooAbsCategoryLValue& categories, RooAbsData& dataset, RooWorkspace& workspace) {

  std::map<std::string, float> results;

  const float resolution = 5.;
  const int nBinsForChi2 = (observable.getMax() - observable.getMin() + 0.5) / resolution;
  observable.setBins(nBinsForChi2);
  std::cout << "Binning dataset with " << nBinsForChi2 << " bins for chi2 computation (" << observable.getMin() << " -> " << observable.getMax() << " ; resolution: " << resolution << " GeV)" << std::endl;

  RooArgSet* floatingParameters = static_cast<RooArgSet*>(simPdf.getParameters(RooArgSet(observable))->selectByAttrib("Constant", false));
  int numberOfFloatingParams = floatingParameters->getSize() - 1;
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

    const RooArgList& coefs = pdf->coefList();

    std::cout << "\tNumber of signal events: " << ((RooAbsReal&) coefs[0]).getVal() << " out of " << table->get(type->GetName()) << std::endl;
    std::cout << "\tNumber of background events: " << ((RooAbsReal&) coefs[1]).getVal() << " out of " << table->get(type->GetName()) << std::endl;

    categories = type->GetName();
    std::string cut = buildCutFormula(categories);

    RooAbsData* reducedDataset = dataset.reduce(Cut(cut.c_str()));
    RooDataHist* binnedDataset = new RooDataHist("binnedDataset", "binned dataset", RooArgSet(observable), *reducedDataset);

    float chi2 = RooChi2Var("chi2", "chi2", *pdf, *binnedDataset, DataError(RooAbsData::Poisson)).getVal();
    float chi2NDF = (chi2 / (observable.getBins() - numberOfFloatingParams / numberOfCategories));
    std::cout << "\tChi2/NDF: " << chi2NDF << std::endl;

    results[type->GetName()] = chi2NDF;

    combinedChi2 += chi2;

    delete binnedDataset;
    delete reducedDataset;
  }

  std::cout << std::endl;
  delete table;

  combinedChi2 /= (numberOfCategories * observable.getBins() - numberOfFloatingParams);
  std::cout << "Combined Chi2: " << combinedChi2 << std::endl;
  results["combined"] = combinedChi2;

  return results;
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

void parseConfigFile(const std::string& filename, RooAbsCategoryLValue& categories, RooWorkspace& workspace) {

  std::fstream configFile((BASE_PATH + "/" + filename).c_str(), std::ios::in);
  std::string line;
  while (std::getline(configFile, line)) {
    if (line[0] == '#' || line.length() == 0)
      continue;

    workspace.factory(line.c_str());
  }

  workspace.Print("v");
}

void fitMtt(int massZprime, bool fit, string fitConfigurationFile, bool doLikScan, bool writeRootFile, bool writeTxtFile, bool saveFigures, bool doLimitCurve, int nToyExp, bool doLikScanInToys, int index, bool doDiscCurve, string syst_str, string systCBsign, string systCB, bool bkgOnly, bool muonsOnly, int btag)
{

  if ((syst_str != "nominal") && (syst_str != "JECup") && (syst_str != "JECdown"))
  {
    cout << "ERROR : Badly defined JEC configuration" << endl;
    return;
  }

  if (fitConfigurationFile == "auto" || fitConfigurationFile.empty())
  {
    // Use pdf function name from parameters.json
    fitConfigurationFile = "fit_pdf_faltb.json";
  }

  std::cout << "Loading fit configuration from '" << fitConfigurationFile << "'" << std::endl;

  const int maxBTag = 3;

  const bool combine = (btag == maxBTag);
  //btag = 2;

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
  Float_t minmTT = 500;
  Float_t maxmTT = 2000;
  Int_t nBins = 29;

  RooRealVar Mtt_KF_reco("Mtt_KF_reco", "Mtt_KF_reco", minmTT, maxmTT, "GeV/c^2");

  RooCategory whichLepton("whichLepton", "whichLepton");
  whichLepton.defineType("muon", 13);
  if (! muonsOnly) {
    whichLepton.defineType("electron", 11);
  }

  RooCategory btagCategory("btag", "btag");
  btagCategory.defineType("0-btag");
  btagCategory.defineType("1-btag");
  btagCategory.defineType("2-btag");

  RooSuperCategory superCategory("superCat", "superCat", RooArgList(whichLepton, btagCategory));

  RooAbsCategoryLValue& mainCategory = combine ? static_cast<RooAbsCategoryLValue&>(superCategory) : static_cast<RooAbsCategoryLValue&>(whichLepton);
  
  // Create main workspace for global pdf
  RooWorkspace mainWorkspace("mainWorkspace", "main workspace");

  std::string pdfSignalName;

  std::map<std::string, std::shared_ptr<BaseFunction>> backgroundPdfs = getCategoriesPdf(BASE_PATH + "/fit_configuration", fitConfigurationFile, Mtt_KF_reco, massZprime, "background", mainCategory, &pdfSignalName);

  for (auto& pdf: backgroundPdfs) {
    std::cout << "Background pdf: " << pdf.first << " ";
    pdf.second->getPdf().Print();
    mainWorkspace.import(pdf.second->getPdf());
  }

  TString prefix = TString::Format("data_2011_%s_%d", syst_str.c_str(), massZprime);
  //TString suffix = TString::Format("%s_%s", pdfSignalName.c_str(), bkgfit_str.c_str());
  TString suffix = TString::Format("%s", pdfSignalName.c_str());
  TString indexJob = TString::Format("job%d", index);


  // Updated 2012-3-20
  // See https://hypernews.cern.ch/HyperNews/CMS/get/physics-announcements/1531.html
  double lumi_mu            = 4678. * 1.066; // Original: 4678; +6.6%
  double lumi_e             = 4682. * 1.066; // Original: 4682; +6.6%
  double s_lumi_mu_percent  = 2.2 / 100.; // 2.2%

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

  std::map<int, double> s_sys_JEC; // JEC systematic error for each b-tag
  std::map<int, double> s_sys_PDF; // Signal PDF systematic error for each b-tag
  std::map<int, double> s_sys_PDF_CB; // Background PDF systematic error for each b-tag

  double sigma_ref = 0;; // Reference cross-section, obtained when fitting once
  loadSigmaRef(massZprime, btag, sigma_ref);

  if (! combine) {
    loadEfficiencies(massZprime, syst_str, btag, sel_eff_mu[btag], sel_eff_e[btag], hlt_eff_mu[btag], hlt_eff_e[btag], s_sel_eff_mu[btag], s_sel_eff_e[btag], s_hlt_eff_mu[btag], s_hlt_eff_e[btag]);
    loadSystematics(massZprime, btag, s_sys_JEC[btag], s_sys_PDF[btag], s_sys_PDF_CB[btag]);
  } else {
    // Load for 0, 1 and 2 btag
    for (int i = 0; i < maxBTag; i++) {
      loadEfficiencies(massZprime, syst_str, i, sel_eff_mu[i], sel_eff_e[i], hlt_eff_mu[i], hlt_eff_e[i], s_sel_eff_mu[i], s_sel_eff_e[i], s_hlt_eff_mu[i], s_hlt_eff_e[i]);
      loadSystematics(massZprime, i, s_sys_JEC[i], s_sys_PDF[i], s_sys_PDF_CB[i]);
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

  std::map<std::string, std::pair<std::shared_ptr<TFile>, RooAbsPdf*>> signalPdfs;

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

    TString workspaceFile = TString::Format("%s/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), syst_str.c_str(), massZprime, pdfSignalName.c_str(), categoryBTag);

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
    renameAndSetPdfParametersConst(RooArgSet(Mtt_KF_reco), *pdf, type->GetName());
    signalPdfs[category] = std::make_pair(file, pdf);
    mainWorkspace.import(*pdf);
  }

  std::cout << "Done." << std::endl;

  if (systCB != "none") // Signal systematics mode
  {
    RooAbsArg* arg = mainWorkspace.arg(systCB.c_str());
    if (! arg)
    {
      std::cerr << "ERROR: parameter " << systCB << " not found inside workspace!" << std::endl;
      exit(1);
    }

    RooRealVar* var = static_cast<RooRealVar*>(arg);
    double oldValue = var->getVal();
    if (systCBsign == "up")
    {
      var->setVal(var->getVal() + var->getErrorHi());
    }
    else
    {
      // Note: getErrorLo() returns something negative.
      var->setVal(var->getVal() + var->getErrorLo());
    }

    std::cout << systCB << " value changed from " << oldValue << " to " << var->getVal() << std::endl;
  }

  // Compute total systematic error
  std::map<int, double> eff_mu;
  std::map<int, double> eff_e;

  double combined_efficiency = 0;
  
  if (! combine) {

    eff_mu[btag] = computeEfficiency(sel_eff_mu[btag], hlt_eff_mu[btag]);
    eff_e[btag]  = computeEfficiency(sel_eff_e[btag], hlt_eff_e[btag]);
    combined_efficiency = eff_mu[btag];

  } else {

    for (int i = 0; i < maxBTag; i++) {
      eff_mu[i] = computeEfficiency(sel_eff_mu[i], hlt_eff_mu[i]);
      eff_e[i]  = computeEfficiency(sel_eff_e[i], hlt_eff_e[i]);
    }

    combined_efficiency = eff_mu[2]; // Our parameters is nSig_mu for 2 btag. Use its efficiency for sigma computation
  }

  std::cout << "Selection efficiency: " << combined_efficiency * 100 << " %" << std::endl;

  std::map<int, double> s_eff_mu_percent;
  std::map<int, double> s_eff_mu_pb;

  std::map<int, double> s_eff_e_percent;
  std::map<int, double> s_eff_e_pb;

  std::map<int, double> s_yield_percent;
  std::map<int, double> s_yield_pb;

  double s_lumi_mu_pb                   = 0.;
  double b_tagging_corr_error_pb        = 0.;

  double err_sys_pb                     = 0.;
  double err_sys_events                 = 0.;
  double err_sys_percent                = 0.;

  s_lumi_mu_pb = sigma_ref * s_lumi_mu_percent;
  b_tagging_corr_error_pb = sigma_ref * b_tagging_corr_error_relative;

  if (! combine) {
    s_eff_mu_percent[btag] = sqrt(s_sel_eff_mu[btag] * s_sel_eff_mu[btag]
        + s_hlt_eff_mu[btag] * s_hlt_eff_mu[btag]
        + trigger_corr_muons_error_relative * trigger_corr_muons_error_relative
        + muID_correction_error_relative * muID_correction_error_relative
        + muIso_correction_error_relative * muIso_correction_error_relative);

    s_eff_e_percent[btag] = sqrt(s_hlt_eff_e[btag] * s_hlt_eff_e[btag]
        + trigger_corr_ele_error_relative * trigger_corr_ele_error_relative
        + eleID_correction_error_relative * eleID_correction_error_relative
        + eleIso_correction_error_relative * eleIso_correction_error_relative); //dont put the sel eff cause it's already counted in s_eff_mu ?

    s_yield_percent[btag] = sqrt(s_sys_JEC[btag] * s_sys_JEC[btag] + s_sys_PDF[btag] * s_sys_PDF[btag] + s_sys_PDF_CB[btag] * s_sys_PDF_CB[btag]);

    s_yield_pb[btag] = sigma_ref * s_yield_percent[btag];
    s_eff_mu_pb[btag] = sigma_ref * s_eff_mu_percent[btag];
    s_eff_e_pb[btag] = sigma_ref * s_eff_e_percent[btag];

    if (ONLY_LUMI_SYST) {
      std::cout << "WARNING: Using only luminosity error for systematics" << std::endl;
      err_sys_percent = s_lumi_mu_percent;
    } else {
      err_sys_percent = sqrt(
          s_yield_percent[btag] * s_yield_percent[btag] +
          s_eff_mu_percent[btag] * s_eff_mu_percent[btag] +
          s_lumi_mu_percent * s_lumi_mu_percent +
          b_tagging_corr_error_relative * b_tagging_corr_error_relative +
          s_eff_e_percent[btag] * s_eff_e_percent[btag]);
    }

    err_sys_pb = sigma_ref * err_sys_percent;
    err_sys_events = eff_mu[btag] * lumi_mu * br_semil * err_sys_pb;

  } else {

    // FIXME: Is this even right?!

    if (ONLY_LUMI_SYST)
      std::cout << "WARNING: Using only luminosity error for systematics" << std::endl;

    for (int i = 0; i < maxBTag; i++) {
      s_eff_mu_percent[i] = sqrt(s_sel_eff_mu[i] * s_sel_eff_mu[i]
          + s_hlt_eff_mu[i] * s_hlt_eff_mu[i]
          + trigger_corr_muons_error_relative * trigger_corr_muons_error_relative
          + muID_correction_error_relative * muID_correction_error_relative
          + muIso_correction_error_relative * muIso_correction_error_relative);

      s_eff_e_percent[i] = sqrt(s_hlt_eff_e[i] * s_hlt_eff_e[i]
          + trigger_corr_ele_error_relative * trigger_corr_ele_error_relative
          + eleID_correction_error_relative * eleID_correction_error_relative
          + eleIso_correction_error_relative * eleIso_correction_error_relative); //dont put the sel eff cause it's already counted in s_eff_mu ?

      s_yield_percent[i] = sqrt(s_sys_JEC[i] * s_sys_JEC[i] + s_sys_PDF[i] * s_sys_PDF[i] + s_sys_PDF_CB[i] * s_sys_PDF_CB[i]);

      s_yield_pb[i] = sigma_ref * s_yield_percent[i];
      s_eff_mu_pb[i] = sigma_ref * s_eff_mu_percent[i];
      s_eff_e_pb[i] = sigma_ref * s_eff_e_percent[i];

      if (ONLY_LUMI_SYST) {
        err_sys_percent += s_lumi_mu_percent;
      } else {
        err_sys_percent += sqrt(
            s_yield_percent[i] * s_yield_percent[i] +
            s_eff_mu_percent[i] * s_eff_mu_percent[i] +
            s_lumi_mu_percent * s_lumi_mu_percent +
            b_tagging_corr_error_relative * b_tagging_corr_error_relative +
            s_eff_e_percent[i] * s_eff_e_percent[i]);
      }

    }

    err_sys_pb     = sigma_ref * err_sys_percent;
    //err_sys_events = eff_mu[i] * lumi_mu * br_semil * err_sys_pb;
  }

#ifdef NO_SYST
  err_sys_percent = err_sys_pb = err_sys_events = 0.;
  std::cout << "WARNING: Systematics are set to 0. If it's not wanted, please undef NO_SYST." << std::endl;
#endif

  // Output systematics error
  if (! combine) {
    std::cout << std::endl;

    double err_base = eff_mu[btag] * lumi_mu * br_semil;
    std::cout << "Reference cross-section: " << sigma_ref << " pb" << std::endl;
    std::cout << "Efficiencies: " << std::endl;
    std::cout << " - Eff: " << combined_efficiency * 100 << " %" << std::endl;
    std::cout << " - Lumi mu: " << lumi_mu << " /pb" << std::endl;
    std::cout << std::endl;
    std::cout << "Total syst errors: " <<
#ifdef NO_SYST
      Bash::set_color(Bash::Color::RED) <<
#endif
      err_sys_events << " events ; " << err_sys_pb << " pb ; " << err_sys_percent << " %" <<
#ifdef NO_SYST
      Bash::set_color() <<
#endif
      std::endl;
    std::cout << " - Eff mu: " << err_base * s_eff_mu_pb[btag] << " events ; " << s_eff_mu_pb[btag] << " pb ; " << s_eff_mu_percent[btag] << " %" << std::endl;
    std::cout << " - Eff e: " << err_base * s_eff_e_pb[btag] << " events ; " << s_eff_e_pb[btag] << " pb ; " << s_eff_e_percent[btag] << " %" << std::endl;
    std::cout << " - Lumi mu: " << err_base * s_lumi_mu_pb << " events ; " << s_lumi_mu_pb << " pb ; " << s_lumi_mu_percent << " %" << std::endl;
    std::cout << " - B-tagging: " << err_base * b_tagging_corr_error_pb << " events ; " << b_tagging_corr_error_pb << " pb ; " << b_tagging_corr_error_relative << " %" << std::endl;
    std::cout << " - Yield: " << std::endl;
    std::cout << "    - JEC: " << err_base * sigma_ref * s_sys_JEC[btag] << " events ; " << sigma_ref * s_sys_JEC[btag] << " pb ; " << s_sys_JEC[btag]  * 100 << " %" << std::endl;
    std::cout << "    - Bkg PDF: " << err_base * sigma_ref * s_sys_PDF[btag] << " events ; " << sigma_ref * s_sys_PDF[btag] << " pb ; " << s_sys_PDF[btag] * 100 << " %" << std::endl;
    std::cout << "    - Signal PDF: " << err_base * sigma_ref * s_sys_PDF_CB[btag] << " events ; " << sigma_ref * s_sys_PDF_CB[btag] << " pb ; " << s_sys_PDF_CB[btag] * 100 << " %" << std::endl;
    std::cout << std::endl;
  } else {
    std::cout << std::endl << Bash::set_color(Bash::Color::RED) << "Systematics are not yet implemented for combined analysis !!" << Bash::set_color() << std::endl << std::endl;
    err_sys_events = err_sys_pb = 0;
  }

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

  //Background::Function* bkgFct = Background::createFunction(bkgfit_str, Mtt_KF_reco);

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
  std::string configFile = (combine) ? "combined_btag_global_pdf.script" : "individual_btag_global_pdf.script";
  parseConfigFile(configFile, mainCategory, mainWorkspace);

  RooRealVar& nSig = *mainWorkspace.var("nSig");

  // Mtt_KF_reco global PDFs
  /*RooRealVar nSig_mu("nSig_mu", "number of sig events", 0., -1000., 2000.);
  RooFormulaVar nSig_e("nSig_e", "number of sig events", "nSig_mu*effRatio*lumiRatio", RooArgList(nSig_mu, effRatio, lumiRatio));
  RooRealVar nBkg_mu("nBkg_mu", "number of bkg events", 2000., 0., 50000);
  RooRealVar nBkg_e("nBkg_e", "number of bkg events", 2000., 0., 50000);*/

  //RooAddPdf globalPdf_mu("global_PDF_mu", "sigPdf+BkgPdf", RooArgList(*sigPdf_mu, backgroundPdfs["muon"]->getPdf()), RooArgList(nSig_mu, nBkg_mu));
  //RooAddPdf globalPdf_e("global_PDF_e", "sigPdf+BkgPdf", RooArgList(*sigPdf_e, backgroundPdfs["electron"]->getPdf()), RooArgList(nSig_e, nBkg_e));

  RooSimultaneous simPdf("simPdf", "simultaneous pdf", mainCategory) ;

  // If we are doing the limit curve, fit with background only
  //if (doLimitCurve) {
  //  simPdf.addPdf(*bkgPdf_mu, "muon");
  //  simPdf.addPdf(*bkgPdf_e, "electron");
  //} else {
  if (bkgOnly)
    std::cout << "Warning: fitting with background model only!" << std::endl;

  // NO shared_ptr<> here.
  std::map<std::string, RooAbsPdf*> globalPdfs;

  it = mainCategory.typeIterator();
  type = nullptr;
  while ((type = static_cast<RooCatType*>(it->Next()))) {
    std::string name = type->GetName();

    std::string cleanedCategory = TString(type->GetName())
      .ReplaceAll(";", "_").ReplaceAll("{", "").ReplaceAll("}", "").ReplaceAll("-", "").Data(); // For workspace names

    //if (muonsOnly && !TString(name).Contains("muon", TString::kIgnoreCase))
    //  continue;

    std::string workspaceName = "global_pdf_" + cleanedCategory;
    std::cout << "Looking for " << workspaceName << " inside workspace" << std::endl;
    globalPdfs[name] = mainWorkspace.pdf(workspaceName.c_str());

    const RooAbsPdf& pdf = (bkgOnly) ? backgroundPdfs[name]->getPdf() : *globalPdfs[name];

    std::cout << "Adding pdf ";
    pdf.Print();
    std::cout << " for category " << name << std::endl;
      
    simPdf.addPdf(pdf, name.c_str());
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
      TString datafilename_0btag = TString::Format("%s/data/ds_data_2011-%s_0_btag.txt", BASE_PATH.c_str(), syst_str.c_str());
      TString datafilename_1btag = TString::Format("%s/data/ds_data_2011-%s_1_btag.txt", BASE_PATH.c_str(), syst_str.c_str());
      TString datafilename_2btag = TString::Format("%s/data/ds_data_2011-%s_2_btag.txt", BASE_PATH.c_str(), syst_str.c_str());

      std::cout << "Combining b-tag from " << datafilename_0btag << ", " << datafilename_1btag << " and " << datafilename_2btag << std::endl;

      RooDataSet* dataset_0btag = RooDataSet::read(datafilename_0btag, RooArgList(Mtt_KF_reco, whichLepton));

      RooDataSet* dataset_1btag = RooDataSet::read(datafilename_1btag, RooArgList(Mtt_KF_reco, whichLepton));

      RooDataSet* dataset_2btag = RooDataSet::read(datafilename_2btag, RooArgList(Mtt_KF_reco, whichLepton));

      dataOrig = new RooDataSet("combData", "combined data", RooArgSet(Mtt_KF_reco, whichLepton), Index(btagCategory), Import("0-btag", *dataset_0btag), Import("1-btag", *dataset_1btag), Import("2-btag", *dataset_2btag));

      dataOrig->table(superCategory)->Print("v");

    } else {
      // Get dataset from external file
      TString datafilename = TString::Format("%s/data/ds_data_2011-%s_%d_btag.txt", BASE_PATH.c_str(), syst_str.c_str(), btag);
      cout << "Using data from " << datafilename << endl;
      dataOrig = RooDataSet::read(datafilename.Data(), RooArgList(Mtt_KF_reco, whichLepton));
    }

    // Reduce data set
    //TODO: This step is not needed, remove it
    TString cutstr = TString::Format("Mtt_KF_reco>%d&&Mtt_KF_reco<%d", (int) minmTT, (int) maxmTT);
    RooDataSet* RedData = (RooDataSet*) dataOrig->reduce(cutstr.Data());
    std::cout << "Dataset entries: " << RedData->numEntries() << std::endl;

    Double_t minmTTFit = minmTT + 0.0;
    Double_t maxmTTFit = maxmTT - 0.0;
    Mtt_KF_reco.setRange(minmTTFit, maxmTTFit);
    // RooFit::Optimize(1) is needed in RooFit 3.50, otherwise the fit
    // does NOT converge. It disables variable caching introduced by RooFit 3.50
    // Optimize(0) does also works, but it disable caching completely.
    // It seems fit does NOT converge when using NumCPU != 1. Awesome!
    RooFitResult *fitResult = simPdf.fitTo(*RedData, Save(), Optimize(1));
    fitResult->Print("v");
    
    TFile* outputFile = nullptr;
    if (writeRootFile) {
      outputFile = TFile::Open(OUTPUT_PATH + prefix + "_fitRes_" + suffix + ".root", "RECREATE");
    }

    drawHistograms(mainCategory, Mtt_KF_reco, nBins, *RedData, simPdf, backgroundPdfs, btag, saveFigures, std::string(prefix), std::string(suffix), !bkgOnly, false, outputFile);
    drawHistograms(mainCategory, Mtt_KF_reco, nBins, *RedData, simPdf, backgroundPdfs, btag, saveFigures, std::string(prefix), std::string(suffix), !bkgOnly, true, outputFile);

    std::map<std::string, float> chi2 = computeChi2(Mtt_KF_reco, simPdf, mainCategory, *RedData, mainWorkspace);

    if (outputFile) {
      outputFile->cd();
      fitResult->Write("fitResult");
    }

    if (outputFile) {
      outputFile->Close();
      delete outputFile;
    }

    //FIXME
    double sigmaZ = nSig.getVal() / (combined_efficiency * lumi_mu * br_semil);

    double errorqstat = nSig.getError() * nSig.getError() / (combined_efficiency * lumi_mu * br_semil * combined_efficiency * lumi_mu * br_semil);
    double errorqtot_pb = errorqstat + s_yield_pb[btag] * s_yield_pb[btag] + b_tagging_corr_error_pb * b_tagging_corr_error_pb + s_eff_mu_pb[btag] * s_eff_mu_pb[btag] + s_lumi_mu_pb * s_lumi_mu_pb;
    double Limit_Z_obs_pb = sigma_ref + 2. * sqrt(errorqtot_pb);

    if (!doLikScan) {
      cout << "The Zprime cross section is " << sigmaZ << " +- " << sqrt(errorqtot_pb) << " pb" << endl;
      cout << "The 95% C.L. upper limit on the Zprime cross section is " << Limit_Z_obs_pb << " pb" << endl;

      if (DO_SYST_COMPUTATION) {

        if (systCB != "none") {
          // Save the systematics
          saveSystematicParameter(massZprime, btag, "signal", systCB, systCBsign, nSig.getVal(), sigmaZ, chi2["combined"], fitResult);
        } else if (syst_str != "nominal") {
          saveSystematicParameter(massZprime, btag, "jec", syst_str, "", nSig.getVal(), sigmaZ, chi2["combined"], fitResult);
        } else {
          saveSystematicParameter(massZprime, btag, "background", fitConfigurationFile, "", nSig.getVal(), sigmaZ, chi2["combined"], fitResult);
        }
      }

      if (SAVE_SIGMA) {
        saveSigma(massZprime, btag, sigmaZ, chi2["combined"], nSig.getVal(), fitResult);
      }
    } else {

      LikelihoodResults results;
      //FIXME. Steps is 2
      doLikelihoodScan(*RedData, simPdf, nSig, massZprime, fitResult->minNll(), 20, err_sys_events, results);

      TFile likelihoodFile(OUTPUT_PATH + prefix + "_likscan_" + suffix + ".root", "recreate");
      likelihoodFile.cd();
      results.likscan->Write();
      results.pdfscan->Write();
      results.pdfscan_wsyst->Write();
      results.pdfscan_wsyst_cut->Write();
      likelihoodFile.Close();

      TCanvas * tmpCanvas = new TCanvas("tmpCanvas", "", 600, 600);

      results.likscan->Draw();
      tmpCanvas->Print(OUTPUT_PATH + prefix + "_likscan_" + suffix + ".pdf");

      tmpCanvas->Clear();
      results.pdfscan->Draw();
      tmpCanvas->Print(OUTPUT_PATH + prefix + "_pdfscan_" + suffix + ".pdf");

      tmpCanvas->Clear();
      results.pdfscan_wsyst->Draw();
      tmpCanvas->Print(OUTPUT_PATH + prefix + "_pdfscan_wsyst_" + suffix + ".pdf");

      tmpCanvas->Clear();
      results.pdfscan_wsyst_cut->Draw();
      tmpCanvas->Print(OUTPUT_PATH + prefix + "_pdfscan_wsyst_cut_" + suffix + ".pdf");

      delete tmpCanvas;

      cout << "The Zprime cross section is " << sigmaZ << " +- " << sqrt(errorqtot_pb) << " pb" << endl;
      cout << "The 95% C.L. upper limit on the Zprime cross section is " << Limit_Z_obs_pb << " pb" << endl;
      cout << "95% prob. limit from scan " << results.scan_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
      cout << "95% prob. limit from scan with systematics " << results.scan_wsyst_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
      cout << "95% prob. limit from scan with systematics and prior sigma>0. " << results.scan_wsyst_cut_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
      ofstream outlikscan(OUTPUT_PATH + prefix + "_likscan.txt");
      outlikscan << "The Zprime cross section is " << sigmaZ << " +- " << sqrt(errorqtot_pb) << " pb" << endl;
      outlikscan << "The 95% C.L. upper limit on the Zprime cross section is " << Limit_Z_obs_pb << " pb" << endl;
      outlikscan << "Systematics error used for gaussian smearing: " << err_sys_events << std::endl;
      outlikscan << "95% prob. limit from scan " << results.scan_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
      outlikscan << "95% prob. limit from scan with systematics " << results.scan_wsyst_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
      outlikscan << "95% prob. limit from scan with systematics and prior sigma>0. " << results.scan_wsyst_cut_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
      outlikscan.close();

      saveLikelihoodResults(massZprime, (combine) ? 3 : btag, results, combined_efficiency * lumi_mu * br_semil);

      results.release();
    }

    /*
    if (writeTxtFile)
    {
      ofstream outParFileNominalFit(OUTPUT_PATH + prefix + "_fitRes_" + suffix + ".txt");
      RooArgSet outFitPar(nSig_mu, nSig_e, nBkg_mu, nBkg_e, nEventsRed_mu, nEventsRed_e);
      outFitPar.writeToStream(outParFileNominalFit, kFALSE);
      outParFileNominalFit << "Chi2/Ndf for muons = "     << chi2mu << "   |   Chi2/Ndf for electrons = " << chi2e << endl;
      outParFileNominalFit << "Chi2/Ndf total = " << chi2tot << endl;
      outParFileNominalFit << "The minimum of the likelihood for the fit is " << fitResult->minNll() << endl;
      outParFileNominalFit << "The Zprime cross section is " << sigmaZ << " +- " << sqrt(errorqtot) << endl;
      outParFileNominalFit << "The 95% C.L. upper limit on the Zprime cross section is " << Limit_Z_obs << endl;
      outParFileNominalFit << endl;
      outParFileNominalFit.close();
    }
    */

    SAFE_DELETE(fitResult);
    SAFE_DELETE(RedData);
  }

  if (doLimitCurve)
  {

    std::cout << "Starting toys analysis..." << std::endl;

    // 0 means a random seed
    RooRandom::randomGenerator()->SetSeed(0);

    // Create PDF for toys generation
    std::map<std::string, std::shared_ptr<BaseFunction>> backgroundPdfsForToys = getCategoriesPdf(BASE_PATH, "fit_pdf.json", Mtt_KF_reco, massZprime, "background", mainCategory, nullptr, "toy");

    std::map<std::string, std::shared_ptr<RooAbsPdf>> globalPdfsForToys;
    std::map<int, std::shared_ptr<RooSimultaneous>> simPdfsForGeneration;
    if (combine) {
      for (int i = 0; i < maxBTag; i++) {
        TString name = TString::Format("simPdfToyForGeneration_%d_btag", i);
        simPdfsForGeneration[i] = std::shared_ptr<RooSimultaneous>(new RooSimultaneous(name, "simultaneous pdf for toys generation", whichLepton));
      }
    } else {
      simPdfsForGeneration[btag] = std::shared_ptr<RooSimultaneous>(new RooSimultaneous("simPdfToyForGeneration", "simultaneous pdf for toys generation", whichLepton));
    }
    
    // Create PDF for fitting
    RooSimultaneous simPdfToyFit("simPdfToyFit", "simultaneous pdf for toys fitting", mainCategory);

    int nEventsToy = 0;

    TFile *fitFile = TFile::Open(TString::Format("%s/%s_fitRes_%s.root", BASE_PATH.c_str(), prefix.Data(), suffix.Data()));
    if (! fitFile) {
      fitFile = TFile::Open(TString::Format("%s/%s_fitRes_%s.root", OUTPUT_PATH.c_str(), prefix.Data(), suffix.Data()));
    }

    RooFitResult *fitResult = static_cast<RooFitResult*>(fitFile->Get("fitResult"));
    double fitMinNll = fitResult->minNll();
    SAFE_DELETE(fitResult);

    std::map<std::string, RooRealVar*> nEvents;
    it = mainCategory.typeIterator();
    type = nullptr;
    while ((type = static_cast<RooCatType*>(it->Next()))) {
      std::string category = type->GetName();
      std::string cleanedCategory = TString(category).ReplaceAll(";", "_").Data(); // Root does not like ';' in names
      TString objectName = TString::Format("nEvents_%s", cleanedCategory.c_str());

      RooRealVar* events = (RooRealVar*) fitFile->Get(objectName);
      events->SetName(objectName);
      nEventsToy += (int) events->getVal();

      nEvents[category] = events;
    }

    SAFE_DELETE(fitFile);

    TFile* toyResFile = NULL;
    //if (writeRootFile) {
    toyResFile = new TFile(OUTPUT_PATH + prefix + "_toylimit_" + suffix + "_" + indexJob + ".root", "RECREATE");
    toyResFile->cd();
    //}

    // Book historams for toy MC results
    TH1F* hSigFrac = new TH1F("hSigFrac_mu", "signal pull (mu) in Toy MC", 20, -4, 4);
    TH1F* hResidual = new TH1F("hResidual_mu", "residual in Toy MC", 100, -500, 500);
    TH1F* hSigma = new TH1F("hSigma", "signal pull (sigma) in Toy MC", 15, -3, 3);
    TH1F* hLoErr = new TH1F("hLoErr_mu", "lo err (mu) in Toy MC", 40, 0, 200);
    TH1F* hHiErr = new TH1F("hHiErr_mu", "hi err (mu) in Toy MC", 40, 0, 200);
    TH1F* hLimit = new TH1F("hLimit_mu", "signal limit (mu) in Toy MC", 33, -32, 100);
    TH1F* hLimit_Z = new TH1F("hLimit_Z", "limit on the Z cross section", 990, -32, 100);
    TH1F* hMinNll = new TH1F("hMinNll", "minimum NLL for Toy MC", 20, fitMinNll - 0.3 * fitMinNll, fitMinNll + 0.3 * fitMinNll);

    // Set background variables values
    std::cout << "Parameters: " << std::endl;

    it = mainCategory.typeIterator();
    type = nullptr;
    while ((type = static_cast<RooCatType*>(it->Next()))) {

      std::string category = type->GetName();

      std::cout << "Category " << category << ": " << std::endl;

      std::map<std::string, std::shared_ptr<RooRealVar>>& parameters = backgroundPdfs[category]->getParameters();
      for (auto& parameter: parameters) {
        std::string parameterName = category + "_" + parameter.first;
        RooRealVar* var = mainWorkspace.var(parameterName.c_str());

        backgroundPdfsForToys[category]->getParameters()[parameter.first]->setVal(var->getVal());
        backgroundPdfsForToys[category]->getParameters()[parameter.first]->setAsymError(var->getErrorLo(), var->getErrorHi());
        backgroundPdfsForToys[category]->getParameters()[parameter.first]->setConstant(true);
        std::cout << parameterName << " = " << var->getVal() << std::endl;
      }

      // Create extended pdf
      globalPdfsForToys[category] = std::shared_ptr<RooAbsPdf>(new RooExtendPdf(("bkgExtended_" + category).c_str(), "Background extended pdf for toys generation", backgroundPdfsForToys[category]->getPdf(), *nEvents[category]));

      if (combine) {
        std::string leptonName = TString(category).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron";
        static boost::regex regex("([0-9]+)-btag");
        boost::smatch regexResults;

        int extractedBTag = 0;
        if (boost::regex_search(category, regexResults, regex)) {
          std::string result(regexResults[1].first, regexResults[2].second);
          extractedBTag = atoi(result.c_str());
        }

        simPdfsForGeneration[extractedBTag]->addPdf(*globalPdfsForToys[category], leptonName.c_str());
      } else {
        simPdfsForGeneration[btag]->addPdf(*globalPdfsForToys[category], category.c_str());
      }

      simPdfToyFit.addPdf(*globalPdfs[category], category.c_str());
    }

    std::cout << "Generating " << nEventsToy << " toys events" << std::endl;

    // Create parameters for toys generation
    std::map<int, std::shared_ptr<RooAbsPdf::GenSpec>> genSpecs;

    if (! combine) {
      genSpecs[btag] = std::shared_ptr<RooAbsPdf::GenSpec>(simPdfsForGeneration[btag]->prepareMultiGen(RooArgSet(Mtt_KF_reco, whichLepton), NumEvents(nEventsToy), Extended(true)/*, Verbose(true)*/));
    } else {
      for (int i = 0; i < maxBTag; i++) {
        genSpecs[i] = std::shared_ptr<RooAbsPdf::GenSpec>(simPdfsForGeneration[i]->prepareMultiGen(RooArgSet(Mtt_KF_reco, whichLepton, btagCategory), /*NumEvents(nEventsToy),*/ Extended(true)/*, Verbose(true)*/));
      }
    }

    double minmTTFit = minmTT + 0.0;
    double maxmTTFit = maxmTT - 0.0;
    Mtt_KF_reco.setRange(minmTTFit, maxmTTFit);

    // NLL for fitting
    RooAbsReal* nll = NULL;

    // Number of experiments with worse fits
    Int_t nWorseFit = 0;
    // Loop over experiments
    for (int i = 0; i < nToyExp; i++)
    {
      // Generating PDF
      std::cout << "Generating toy experiment number " << i << " of " << nToyExp << std::endl;

      std::cout << "Generating distribution..." << std::endl;
      RooDataSet* toyData = nullptr;
      std::map<int, RooDataSet*> toyDataForEachBTag;
      if (! combine) {
        toyData = simPdfsForGeneration[btag]->generate(*genSpecs[btag]);
      } else {
        toyData = new RooDataSet("toy_dataset", "toy dataset", RooArgSet(Mtt_KF_reco, whichLepton, btagCategory));
        for (int j = 0; j < maxBTag; j++) {
          btagCategory.setIndex(j);
          toyDataForEachBTag[i] =  simPdfsForGeneration[j]->generate(*genSpecs[j]);
          toyData->append(*toyDataForEachBTag[i]);
        }
      }
      std::cout << "done." << std::endl;

      TFile* myFile = TFile::Open(OUTPUT_PATH + prefix + "_toylimit_" + suffix + "_" + indexJob + "_" + (Long_t) i + ".root", "RECREATE");
      drawHistograms(mainCategory, Mtt_KF_reco, nBins, *toyData, simPdfToyFit, backgroundPdfsForToys, btag, false, std::string(prefix), std::string(suffix), false, false, myFile, true);
      myFile->Close();
      delete myFile;

      std::cout << "Fitting distribution ..." << std::endl;
      if (nll == NULL)
      {
        // Only create the nll the first time
        nll = simPdfToyFit.createNLL(*toyData, RooFit::Optimize(1));
      }
      else
      {
        nll->setData(*toyData);
      }

      //TODO: Maybe we need to reset all pdf parameters to their default values
      // Be sure that nSig_mu is not fixed anymore, and reset to 0
      nSig.setVal(0);
      nSig.setConstant(false);

      // Fit
      RooMinuit* minimizer = new RooMinuit(*nll);
      minimizer->setEvalErrorWall(1);
      minimizer->optimizeConst(1);
      minimizer->migrad();

      // Only compute errors for nSig
      minimizer->minos(RooArgSet(nSig));

      std::cout << "done." << std::endl;

      RooFitResult* toyFitRes = minimizer->save();
      toyFitRes->Print("v");

      double nSigVal   = nSig.getVal();
      double nSigErrHi = nSig.getAsymErrorHi();
      double nSigErrLo = nSig.getAsymErrorLo();

      hResidual->Fill(nSigVal);
      hHiErr->Fill(nSigErrHi);
      hLoErr->Fill(-1. * nSigErrLo);

      double nSigErr;

      if (nSigVal < 0.)
        nSigErr = nSigErrHi;
      else
        nSigErr = -1. * nSigErrLo;

      if (nSigErr > 0.)      
        hSigFrac->Fill(nSigVal / nSigErr);
      else
        hSigFrac->Fill(0);
      

      double Limit = (nSigVal + 2. * nSigErrHi);
      hLimit->Fill(Limit);

      if (nSigVal == 0)
        nSigVal = 0.00000001;

      double sigmaZl = nSigVal / (combined_efficiency * lumi_mu * br_semil);
      double errorqstatl_pb = nSigErrHi * nSigErrHi / (combined_efficiency * lumi_mu * br_semil * combined_efficiency * lumi_mu * br_semil);
      double errorqtotl_pb = errorqstatl_pb + s_yield_pb[btag] * s_yield_pb[btag] + b_tagging_corr_error_pb * b_tagging_corr_error_pb + s_eff_mu_pb[btag] * s_eff_mu_pb[btag] + s_lumi_mu_pb * s_lumi_mu_pb;
      double Limit_Z = sigmaZl + 2. * sqrt(errorqtotl_pb);

      hSigma->Fill(sigmaZl / sqrt(errorqstatl_pb));
      std::cout << "sigmaZ = " << sigmaZl << " +- " << sqrt(errorqtotl_pb) << " pb" << std::endl;

      if (!doLikScanInToys)
      {
        // fill "normally" with mean value+2*sigma
        hLimit_Z->Fill(Limit_Z);
      }
      else
      {

        LikelihoodResults results;
        //FIXME: It's 10 steps
        doLikelihoodScan(*toyData, simPdfToyFit, nSig, massZprime, toyFitRes->minNll(), 20, err_sys_events, results);

        TString dirName = TString::Format("likscans_%s_toy_%d", indexJob.Data(), i);
        toyResFile->mkdir(dirName);
        toyResFile->cd(dirName);

        results.likscan->Write();
        results.pdfscan->Write();
        results.pdfscan_wsyst->Write();
        results.pdfscan_wsyst_cut->Write();

        toyResFile->cd();

        //cout << "i fill with: " << results.scan_wsyst_cut_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
        hLimit_Z->Fill(results.scan_wsyst_cut_limit / (combined_efficiency * lumi_mu * br_semil));

        results.release();
      }

      hMinNll->Fill(toyFitRes->minNll());

      if (toyFitRes->minNll() > fitMinNll)
        nWorseFit++;

      delete toyData;

      for (auto& j: toyDataForEachBTag)
        delete j.second;

      delete toyFitRes;
      delete minimizer;

    }

    Float_t prob = ((Float_t)nWorseFit) / ((Float_t)nToyExp);
    std::cout << "Probability of having a worse fits = " << prob << std::endl;

    // Save results to file
    toyResFile->cd();
    hSigFrac->Write();
    hResidual->Write();
    hSigma->Write();
    hLoErr->Write();
    hHiErr->Write();
    hMinNll->Write();
    hLimit_Z->Write();

    toyResFile->Close();
    delete toyResFile;
  }
}

