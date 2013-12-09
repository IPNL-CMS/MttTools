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

void saveTemporaryResult(const std::string& outputFile, double result) {

  Json::Value root;
  root["result"] = result;

  Json::StyledWriter writer;
  const std::string json = writer.write(root);

  std::ofstream file(outputFile.c_str());
  file << json;
  file.close();
}

void fitMtt(std::map<int, TChain*> chains, int massZprime, bool fit, string bkgfit_str, bool doLikScan, bool writeRootFile, bool saveFigures, bool doLimitCurve, int nToyExp, bool doLikScanInToys, int index, string syst_str, string systCBsign, string systCB, bool bkgOnly, bool muonsOnly, int btag, bool useSharedMemory, key_t shm_key, const std::string& customWorkspaceFile, bool fixBackground, bool saveWorkspace,
    // Background systematics
    const std::string& temporaryResultFile,
    bool doBackgroundSyst, const std::string& backgroundParameterName, double backgroundParameterValue
    );

std::string BASE_PATH;
std::string OUTPUT_PATH;
std::string EFF_FILE;
bool DO_SYST_COMPUTATION = false;
bool SAVE_SIGMA = false;
bool ONLY_LUMI_SYST = false;
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

    TCLAP::ValueArg<int> massArg("m", "mass", "Zprime mass", true, 750, "integer");
    TCLAP::SwitchArg fitArg("", "no-fit", "Don't do the reference fit", true);
    TCLAP::SwitchArg doLikScanArg("", "scan", "Do the likelihood scan", false);
    TCLAP::SwitchArg writeRootArg("", "no-root-files", "Don't write root files", true);
    TCLAP::SwitchArg saveFiguresArg("", "no-figs", "Don't save figures", true);
    TCLAP::SwitchArg doLimitCurveArg("", "limit-curve", "Do the limit curve", false);
    TCLAP::ValueArg<int> nToyArg("", "toys", "Number of toys exp.", false, 1, "integer");
    TCLAP::SwitchArg doLikInToyArg("", "no-scan-in-toys", "Don't do the likelihood scan in toys", true);
    TCLAP::ValueArg<int> indexArg("", "index", "Index", false, 1, "integer");

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
    cmd.add(saveFiguresArg);
    cmd.add(doLimitCurveArg);
    cmd.add(nToyArg);
    cmd.add(doLikInToyArg);
    cmd.add(indexArg);
    cmd.add(systArg);
    cmd.add(systSignArg);
    cmd.add(systCBArg);

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    std::string default_base_path = "./analysis/" + getAnalysisUUID() + "/";
    TCLAP::ValueArg<std::string> pathArg("", "path", "Folder where to load files", false, default_base_path, "string", cmd);
    TCLAP::ValueArg<std::string> outputPathArg("", "output-path", "Folder where output files are stored", false, default_base_path, "string", cmd);
    TCLAP::ValueArg<std::string> effArg("", "eff-file", "File where efficiences are stored (JSON format)", false, "efficiencies.json", "string", cmd);
    TCLAP::SwitchArg doSystComputationArg("", "syst-computation", "Set this flag if you are running me to compute systematic", cmd);
    TCLAP::SwitchArg saveSigmaArg("", "save-sigma-ref", "Set this flag if you want to save the computed sigma as the ref one", cmd);
    TCLAP::SwitchArg bkgOnlyArg("", "bkg-only", "Fit using background model only", cmd);
    TCLAP::SwitchArg onlyMuonArg("", "muons-only", "Compute limits using only semi-mu data", cmd);
    TCLAP::SwitchArg onlyLumiSystArg("", "only-lumi-syst", "Only use luminosity error for systematics", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);
    TCLAP::SwitchArg verboseArg("v", "verbose", "Verbose mode", cmd);
    TCLAP::SwitchArg batchArg("", "batch", "Run in batch mode", cmd);
    TCLAP::ValueArg<std::string> fitConfigFileArg("", "config-file", "Configuration file name containing fit parameters", false, "fit_pdf_falt.json", "string", cmd);
    
    // Shared memory management
    TCLAP::SwitchArg useSharedMemoryArg("", "shared-memory", "Save fit results into a shared memory area. The shm key must be specified", cmd);
    TCLAP::ValueArg<key_t> sharedMemoryKeyArg("", "shm-key", "The shm key used to create the shm", false, 0, "integer", cmd);

    // Workspace
    TCLAP::ValueArg<std::string> workspaceArg("", "workspace", "Use a custom workspace containing the signal pdf to use", false, "", "string", cmd);

    // Fix background
    TCLAP::SwitchArg fixBackgroundArg("", "fix-background", "Fit once with background only, then fix the background and refit background + signal", cmd);

    TCLAP::SwitchArg saveWorkspaceArg("", "save-workspace", "Save the workspace for redoing plot after", cmd);

    // New background syst calculation
    TCLAP::ValueArg<std::string> backgroundParameterNameArg("", "bkg-parameter-name", "The background parameter name to vary", false, "", "string", cmd);
    TCLAP::ValueArg<double> backgroundParameterValueArg("", "bkg-parameter-value", "The background parameter value", false, 0., "double", cmd);
    TCLAP::ValueArg<std::string> temporaryResultFileArg("", "temp-output-file", "A temporary filename to store result", false, "", "string", cmd);

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

    bool doBackgroundSyst = backgroundParameterNameArg.isSet() && backgroundParameterValueArg.isSet() && temporaryResultFileArg.isSet();

    fitMtt(chains, massArg.getValue(), fitArg.getValue(), fitConfigFileArg.getValue(), doLikScanArg.getValue(), writeRootArg.getValue(), saveFiguresArg.getValue(), doLimitCurveArg.getValue(), nToyArg.getValue(), doLikInToyArg.getValue(), indexArg.getValue(), systArg.getValue(), systSignArg.getValue(), systCBArg.getValue(), bkgOnlyArg.getValue(), onlyMuonArg.getValue(), btagArg.getValue(), useSharedMemoryArg.getValue(), sharedMemoryKeyArg.getValue(), workspaceArg.getValue(), fixBackgroundArg.getValue(), saveWorkspaceArg.getValue(),
        // Background systematics
        temporaryResultFileArg.getValue(),
        doBackgroundSyst, backgroundParameterNameArg.getValue(), backgroundParameterValueArg.getValue()
        );

    for (auto& chain: chains)
      delete chain.second;

  }
  catch (TCLAP::ArgException& e)
  {
    std::cerr << e.error() << std::endl;
  }
}

void loadEfficiencies(int mass, const std::string& jecType, int btag, double& eff_mu, double& eff_e, double& error_eff_mu, double& error_eff_e)
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

  if (! effNode.isObject())
  {
    std::cerr << "ERROR: malformated JSON file. Exiting." << std::endl;
    exit(1);
  }

  eff_mu = effNode["eff_mu"].asDouble();
  eff_e = effNode["eff_e"].asDouble();
  error_eff_mu = effNode["error_eff_mu"].asDouble();
  error_eff_e = effNode["error_eff_e"].asDouble();

  std::cout << "Efficiencies loaded successfully for " << mass << " and systematic " << jecType << std::endl;
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

void saveSigma(int mass, int btag, double sigma, double chi, double events, RooFitResult* fitRes)
{

  FILE* lock = fopen((BASE_PATH + "sigma_reference.lock").c_str(), "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  std::ifstream file((BASE_PATH + "sigma_reference.json").c_str());
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  const std::string uuid = getAnalysisUUID(BASE_PATH);
  root[uuid][strMass][btagStr]["sigma"] = sigma;
  root[uuid][strMass][btagStr]["chi2"] = chi;
  root[uuid][strMass][btagStr]["events"] = events;
  root[uuid][strMass][btagStr]["fit_covQual"] = fitRes->covQual();
  root[uuid][strMass][btagStr]["fit_status"] = fitRes->status();


  std::ofstream ofile;
  ofile.open((BASE_PATH + "sigma_reference.json").c_str(), std::ios::out | std::ios::trunc);
  Json::StyledWriter writer;
  ofile << writer.write(root);
  ofile.close();

  fclose(lock);
}

void saveSystematicParameter(int mass, int btag, const std::string& type, const std::string& param, const std::string& subparam, double nEvents, double sigma, double chi2, RooFitResult* fitRes)
{

  FILE* lock = fopen((BASE_PATH + "systematics_parameters.lock").c_str(), "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  std::ifstream file((BASE_PATH + "systematics_parameters.json").c_str());
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  const std::string uuid = getAnalysisUUID(BASE_PATH);

  if (subparam.length() > 0)
  {
    root[uuid][strMass][btagStr][type][param]["sigma"][subparam] = sigma;
    root[uuid][strMass][btagStr][type][param]["events"][subparam] = nEvents;
    root[uuid][strMass][btagStr][type][param]["chi2"][subparam] = chi2;
    root[uuid][strMass][btagStr][type][param]["fit_covQual"][subparam] = fitRes->covQual();
    root[uuid][strMass][btagStr][type][param]["fit_status"][subparam] = fitRes->status();
  }
  else
  {
    root[uuid][strMass][btagStr][type][param]["sigma"] = sigma;
    root[uuid][strMass][btagStr][type][param]["events"] = nEvents;
    root[uuid][strMass][btagStr][type][param]["chi2"] = chi2;
    root[uuid][strMass][btagStr][type][param]["fit_covQual"] = fitRes->covQual();
    root[uuid][strMass][btagStr][type][param]["fit_status"] = fitRes->status();
  }

  //TODO: reswitch to ofstream
  FILE* fd = fopen((BASE_PATH + "systematics_parameters.json").c_str(), "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

#define MIN_LIKELIHOOD -200.
#define MAX_LIKELIHOOD 800.

/**
 * Compute the likelihood range
 */
void computeLikeloodRange(const RooRealVar& v, double& xLow, double& xHigh)
{

  double meanValue = v.getVal();
  xLow = meanValue - 300;
  xHigh = meanValue + 900;
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
  TH1* pdfscan_cut;
  TH1* pdfscan_wsyst;
  TH1* pdfscan_wsyst_cut;

  double scan_limit;
  double scan_cut_limit;
  double scan_wsyst_limit;
  double scan_wsyst_cut_limit;

  LikelihoodResults():
    likscan(NULL), pdfscan(NULL), pdfscan_cut(NULL), pdfscan_wsyst(NULL), pdfscan_wsyst_cut(NULL),
    scan_limit(0.), scan_cut_limit(0.), scan_wsyst_limit(0.), scan_wsyst_cut_limit(0.)
  {}

  void release()
  {
    SAFE_DELETE(likscan);
    SAFE_DELETE(pdfscan);
    SAFE_DELETE(pdfscan_cut);
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
  FILE* lock = fopen((BASE_PATH + "likelihood_scan.lock").c_str(), "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  std::ifstream file((BASE_PATH + "likelihood_scan.json").c_str());
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  const std::string uuid = getAnalysisUUID(BASE_PATH);
  root[uuid][strMass][btagStr]["scan_limit"] = results.scan_limit / denominator;
  root[uuid][strMass][btagStr]["scan_wsyst_limit"] = results.scan_wsyst_limit / denominator;
  root[uuid][strMass][btagStr]["scan_wsyst_cut_limit"] = results.scan_wsyst_cut_limit / denominator;

  std::ofstream ofile;
  ofile.open((BASE_PATH + "likelihood_scan.json").c_str(), std::ios::out | std::ios::trunc);
  Json::StyledWriter writer;
  ofile << writer.write(root);
  ofile.close();

  fclose(lock);
}

void doLikelihoodScan(RooAbsData& dataset, RooAbsPdf& pdf, RooRealVar& observable, const double mass, const double minNLL, const int steps, const double systError, LikelihoodResults& results)
{

  RooFitResult* fitResults = NULL;

  double xLow = 0, xHigh = 0;
  computeLikeloodRange(observable, xLow, xHigh);
  int nLHBins = (int)(xHigh - xLow) / steps;

  std::cout << "Likehood from " << xLow << " to " << xHigh << "; bins: " << nLHBins << std::endl;

  results.likscan = new TH1D("likscan", "likscan;Number of signal events;-#Delta Nll", nLHBins, xLow, xHigh);
  results.pdfscan = new TH1D("pdfscan", "pdfscan;Number of signal events;PDF", nLHBins, xLow, xHigh);
  results.pdfscan->GetYaxis()->SetTitleOffset(1.5);

  // Take into account gaussian smearing
  double xSmearedHigh = xHigh + 3 * systError;
  int nLHSmearedBins = (int)(xSmearedHigh - xLow) / steps;
  int nLHCuttedBins = (int) xSmearedHigh / steps;
  std::cout << "Likehood smearing from " << xLow << " to " << xSmearedHigh << "; bins: " << nLHSmearedBins << std::endl;
  results.pdfscan_wsyst = new TH1D("pdfscan_wsyst", "pdfscan_wsyst;Number of signal events;PDF", nLHSmearedBins, xLow, xSmearedHigh);
  results.pdfscan_wsyst->GetYaxis()->SetTitleOffset(1.5);
  results.pdfscan_cut = new TH1D("pdfscan_cut", "pdfscan_cut;Number of signal events;PDF", nLHCuttedBins, 0, xSmearedHigh);
  results.pdfscan_cut->GetYaxis()->SetTitleOffset(1.5);
  results.pdfscan_wsyst_cut = new TH1D("pdfscan_wsyst_cut", "pdfscan_wsyst_cut;Number of signal events;PDF", nLHCuttedBins, 0, xSmearedHigh);
  results.pdfscan_wsyst_cut->GetYaxis()->SetTitleOffset(1.5);

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
  RooAbsReal* nll = pdf.createNLL(dataset, RooFit::Optimize(0)/*, RooFit::NumCPU(2)*/);

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
    minimizer->setStrategy(1);
    minimizer->setEvalErrorWall(0);
    minimizer->optimizeConst(0);
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
    minimizer->setStrategy(1);
    minimizer->setEvalErrorWall(0);
    minimizer->optimizeConst(0);
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
        double weight = results.pdfscan_wsyst->GetBinContent(results.pdfscan_wsyst->FindBin(x));
        results.pdfscan_wsyst_cut->SetBinContent(results.pdfscan_wsyst_cut->FindBin(x), weight);

        weight = results.pdfscan->GetBinContent(results.pdfscan->FindBin(x));
        results.pdfscan_cut->SetBinContent(results.pdfscan_cut->FindBin(x), weight);
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

  results.pdfscan_cut->ComputeIntegral();
  results.pdfscan_cut->Scale(1. / results.pdfscan_cut->Integral());
  results.pdfscan_cut->GetQuantiles(1, &results.scan_cut_limit, &areascan);

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

    double x1 = 0.18, y1 = 0.14, x2 = 0.57, y2 = 0.52;

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

    if (! bkgOnly) {
      const RooArgList& coefs = pdf->coefList();

      std::cout << "\tNumber of signal events: " << ((RooAbsReal&) coefs[0]).getVal() << " out of " << table->get(type->GetName()) << std::endl;
      std::cout << "\tNumber of background events: " << ((RooAbsReal&) coefs[1]).getVal() << " out of " << table->get(type->GetName()) << std::endl;
    }

    categories = type->GetName();
    std::string cut = buildCutFormula(categories);

    std::shared_ptr<RooAbsData> reducedDataset(dataset.reduce(Cut(cut.c_str())));
    std::shared_ptr<RooDataHist> binnedDataset(new RooDataHist("binnedDataset", "binned dataset", RooArgSet(observable), *reducedDataset));

    float chi2 = RooChi2Var("chi2", "chi2", *pdf, *binnedDataset, DataError(RooAbsData::SumW2)).getVal();
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

void setPdfParametersRange(const RooArgSet& observables, const RooAbsPdf& pdf, double sigma)
{
  RooArgSet* params = pdf.getParameters(observables);
  TIterator* iter = params->createIterator();
  RooRealVar* var = NULL;
  while ((var = static_cast<RooRealVar*>(iter->Next())))
  {
    std::cout << "Setting range of " << var->GetName() << " to +/- " << sigma << " sigma" << std::endl;
    var->setRange(var->getVal() - sigma * var->getError(), var->getVal() + sigma * var->getError());
  }

  delete iter;
  delete params;
}

void mixEfficiencies(
    const std::map<int, double>& effs1,
    const std::map<int, double>& effs2,
    const std::map<int, double>& err_effs1,
    const std::map<int, double>& err_effs2,    
    const std::string& prefix1,
    const std::string& prefix2,
    RooWorkspace& workspace) {

  for (const auto& eff1: effs1) {
    for (const auto& eff2: effs2) {
      if (prefix1 == prefix2 && eff1.first == eff2.first)
        continue;

      std::stringstream ss;
      ss << "eff_ratio_" << prefix1 << "_" << eff1.first << "b_" << prefix2 << "_" << eff2.first << "b[" << eff1.second / eff2.second << "]";

      //std::cout << ss.str() << std::endl;

      workspace.factory(ss.str().c_str());
    }
  }

  if (prefix1 != prefix2) {
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

  //RooHistPdf* hist_pdf = new RooHistPdf(std::string("signal_" + categoryName).c_str(), std::string("signal_" + categoryName).c_str(), RooArgSet(observable), *binnedInterpolatedDataset);

  //interpolation->SetName(std::string("signal_" + categoryName).c_str());
  //renameAndSetPdfParametersConst(RooArgSet(observable), *interpolation, (suffix.length() == 0) ? cleanedCategory : suffix);
  renameAndSetPdfParametersConst(RooArgSet(observable), *keys_pdf, goodSuffix);

  observable.setBins(oldBinning);

  return keys_pdf;
}

void fitMtt(std::map<int, TChain*> eventChain, int massZprime, bool fit, string fitConfigurationFile, bool doLikScan, bool writeRootFile, bool saveFigures, bool doLimitCurve, int nToyExp, bool doLikScanInToys, int index, string syst_str, string systCBsign, string systCB, bool bkgOnly, bool muonsOnly, int btag, bool useSharedMemory, key_t shm_key, const std::string& customWorkspaceFile, bool fixBackground, bool saveWorkspace,
    // Background systematics
    const std::string& temporaryResultFile,
    bool doBackgroundSyst, const std::string& backgroundParameterName, double backgroundParameterValue
    )
{

  if ((syst_str != "nominal") && (syst_str != "JECup") && (syst_str != "JECdown"))
  {
    cout << "ERROR : Badly defined JEC configuration" << endl;
    return;
  }

  if (fitConfigurationFile == "auto" || fitConfigurationFile.empty())
    fitConfigurationFile = "fit_pdf_falt.json";

  if (! fixBackground) {
    fixBackground = analysisFixedBackground();
  } else {
    std::cout << Bash::set_color(Bash::Color::RED) << "Fitting with fixed background (--fix-background flag overrides analysis configuration)" << Bash::set_color() << std::endl;
  }

  std::cout << "Loading fit configuration from '" << fitConfigurationFile << "'" << std::endl;

  // Systematics?
  bool useSystematics = analysisUseSystematics(BASE_PATH);

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
    RooMsgService::instance().addStream(RooFit::ERROR); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
    RooMsgService::instance().setSilentMode(true);
  } else {
    RooMsgService::instance().addStream(RooFit::INFO); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
  }

  //fit region
  Float_t minmTT = 325;
  Float_t maxmTT = 1000;

  RooRealVar mtt("mtt", "M_{t#bar{t}}", minmTT, maxmTT, "GeV");
  RooRealVar weight("weight", "weight", 0, 100000);

  Double_t minmTTFit = minmTT + 0.0;
  Double_t maxmTTFit = maxmTT - 0.0;
  
  mtt.setRange(minmTTFit, maxmTTFit);

  mtt.setRange("R1", minmTTFit, massZprime - 50);
  mtt.setRange("R2", massZprime + 50, maxmTTFit);

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

  std::map<std::string, std::shared_ptr<BaseFunction>> backgroundPdfs = getCategoriesPdf(BASE_PATH + "/fit_configuration", fitConfigurationFile, mtt, NULL, massZprime, "background", mainCategory, NULL);

  std::map<std::string, RooAbsPdf*> backgroundPdfsFromWorkspace;

  for (auto& pdf: backgroundPdfs) {
    std::cout << "Background pdf: " << pdf.first << " ";
    pdf.second->getPdf().Print();
    mainWorkspace.import(pdf.second->getPdf());

    backgroundPdfsFromWorkspace[pdf.first] = mainWorkspace.pdf(pdf.second->getPdf().GetName());
  }

  TString prefix = TString::Format("data_2012_%s_%d", syst_str.c_str(), massZprime);
  //TString suffix = TString::Format("%s_%s", pdfSignalName.c_str(), bkgfit_str.c_str());
  TString suffix = TString::Format("%s", analysisName.c_str());
  TString indexJob = TString::Format("job%d", index);


  // Luminosity updated for full 2012 run, with PIXEL calculation, and Winter13 rereco
  double lumi_mu = 19700;
  double lumi_e  = 19667;

  double s_lumi_mu_percent  = 2.6 / 100.; // 2.6%

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

  std::map<int, double> s_sel_eff_mu; // = Δ(eff_mu)
  std::map<int, double> s_sel_eff_e; // = Δ(eff_e)

  std::map<int, double> s_hlt_eff_mu; // ?
  std::map<int, double> s_hlt_eff_e; // ?

  double s_sys_JEC; // JEC systematic error
  double s_sys_PDF; // Signal PDF systematic error
  double s_sys_PDF_CB; // Background PDF systematic

  double sigma_ref = 0;; // Reference cross-section, obtained when fitting once
  loadSigmaRef(massZprime, btag, sigma_ref);

  loadSystematics(massZprime, btag, s_sys_JEC, s_sys_PDF, s_sys_PDF_CB);

  if (! combine) {
    loadEfficiencies(massZprime, syst_str, btag, sel_eff_mu[btag], sel_eff_e[btag], s_sel_eff_mu[btag], s_sel_eff_e[btag]);
  } else {
    // Load for 0, 1 and 2 btag
    for (int i = minBTag; i <= maxBTag; i++) {
      loadEfficiencies(massZprime, syst_str, i, sel_eff_mu[i], sel_eff_e[i], s_sel_eff_mu[i], s_sel_eff_e[i]);
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

  std::map<
    std::string, // Category name
    SignalPdf 
  > signalPdfsFromWorkspace;

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

    if (massZprime == 500 || massZprime == 700 || massZprime == 1000 || massZprime == 1250 || massZprime == 1500 || massZprime == 2000) {

      TString workspaceFile = TString::Format("%s/frit/%s-Higgs%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), syst_str.c_str(), massZprime, analysisName.c_str(), categoryBTag);
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

      bool higgsMode = true;

      RooAbsPdf* positive_pdf = workspace->pdf(TString::Format(pdfName.c_str(), "positive"));
      if (! positive_pdf) {
        higgsMode = false;
      }

      if (higgsMode) {
        positive_pdf->SetName(std::string("positive_signal_" + std::string(type->GetName())).c_str());
        renameAndSetPdfParametersConst(RooArgSet(mtt), *positive_pdf, type->GetName());
        mainWorkspace.import(*positive_pdf);

        RooAbsPdf* negative_pdf = workspace->pdf(TString::Format(pdfName.c_str(), "negative"));
        if (! positive_pdf) {
          std::cerr << "ERROR: Signal pdf " << pdfName << " not found inside workspace." << std::endl;
          exit(1);
        }

        negative_pdf->SetName(std::string("negative_signal_" + std::string(type->GetName())).c_str());
        renameAndSetPdfParametersConst(RooArgSet(mtt), *negative_pdf, type->GetName());
        mainWorkspace.import(*negative_pdf);

        RooRealVar* positive_integral = workspace->var(TString::Format(integralName.c_str(), "positive"));
        std::cout << positive_integral->getVal() << std::endl;
        mainWorkspace.import(*positive_integral);

        RooRealVar* negative_integral = workspace->var(TString::Format(integralName.c_str(), "negative"));
        std::cout << negative_integral->getVal() << std::endl;
        mainWorkspace.import(*negative_integral);

        SignalPdf pdf = {
          .positive_pdf = mainWorkspace.pdf(positive_pdf->GetName()),
          .positive_integral = mainWorkspace.var(positive_integral->GetName()),
          .negative_pdf = mainWorkspace.pdf(negative_pdf->GetName()),
          .negative_integral = mainWorkspace.var(negative_integral->GetName())
        };

        signalPdfsFromWorkspace[category] = pdf;
      } else {
        std::string zprimePdfName = "signal_" + std::string((TString(type->GetName()).Contains("muon", TString::kIgnoreCase) ? "muon" : "electron"));

        RooAbsPdf* pdf = workspace->pdf(zprimePdfName.c_str());
        pdf->SetName(std::string("signal_" + std::string(type->GetName())).c_str());
        renameAndSetPdfParametersConst(RooArgSet(mtt), *pdf, type->GetName());
        mainWorkspace.import(*pdf);
      }
    } else {

      RooAbsPdf* interpolation = getInterpolatedPdf(mtt, massZprime, "nominal", categoryBTag, type->GetName());
      mainWorkspace.import(*interpolation);

      //signalPdfsFromWorkspace[category] = mainWorkspace.pdf(interpolation->GetName());
    }
  }

  std::cout << "Done." << std::endl;

  mainWorkspace.Print("v");

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

  // Compute total systematic error
  //std::map<int, double> eff_mu;
  //std::map<int, double> eff_e;

  double mu_efficiency = 0;
  double mu_efficiency_error = 0;
  
  // Scale factors are already included in event weight
  if (! combine) {

    mu_efficiency = sel_eff_mu[btag];
    mu_efficiency_error = s_sel_eff_mu[btag];

  } else {

    mu_efficiency = sel_eff_mu[2]; // Our parameters is nSig_mu for 2 btag. Use its efficiency for sigma computation
    mu_efficiency_error = s_sel_eff_mu[2];
  }

  std::cout << std::endl;
  std::cout << "Analysis efficiencies" << std::endl;
  std::cout << "-------------------------------------" << std::endl << std::endl;

  std::cout << "Luminosity: " << Bash::set_color(Bash::Color::MAGENTA) << lumi_mu << " /pb" << Bash::set_color() << std::endl;
  std::cout << "Reference cross-section: " << Bash::set_color(Bash::Color::MAGENTA) << sigma_ref << " pb" << Bash::set_color() << std::endl;
  std::cout << "Total efficiency (selection * SFs): " << std::endl;

  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << "\tMuonic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << sel_eff_mu[i] * 100 << " % +/- " << s_sel_eff_mu[i] * 100 << " %" << Bash::set_color() << std::endl;
    std::cout << "\tElectronic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << sel_eff_e[i] * 100 << " % +/- " << s_sel_eff_e[i] * 100 << " %" << Bash::set_color() << std::endl;
  }
  std::cout << "Expected number of signal events: " << std::endl;

  for (int i = minBTag; i <= maxBTag; i++) {
    std::cout << "\tMuonic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << lumi_mu * sel_eff_mu[i] << Bash::set_color() << std::endl;
    std::cout << "\tElectronic channel, " << i << " b-tag: " << Bash::set_color(Bash::Color::MAGENTA) << lumi_e * sel_eff_e[i] << Bash::set_color() << std::endl;
  }

  if (combine) {
    mixEfficiencies(sel_eff_e, sel_eff_mu, s_sel_eff_e, s_sel_eff_mu, "e", "mu", mainWorkspace);
    mixEfficiencies(sel_eff_e, sel_eff_e, s_sel_eff_e, s_sel_eff_e, "e", "e", mainWorkspace);
    mixEfficiencies(sel_eff_mu, sel_eff_mu, s_sel_eff_mu, s_sel_eff_mu, "mu", "mu", mainWorkspace);
  } else {
    std::stringstream ss;
    ss << "eff_ratio_e_mu[" << sel_eff_e[btag] / sel_eff_mu[btag]<< "]";
    mainWorkspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "eff_e[" << sel_eff_e[btag] << "]";
    mainWorkspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "eff_mu[" << sel_eff_mu[btag] << "]";
    mainWorkspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "error_eff_e[" << s_sel_eff_e[btag] << "]";
    mainWorkspace.factory(ss.str().c_str());

    ss.str(std::string());
    ss << "error_eff_mu[" << s_sel_eff_mu[btag] << "]";
    mainWorkspace.factory(ss.str().c_str());
  }

  RooRealVar lumiRatio("lumiRatio", "luminosity ratio", lumi_e / lumi_mu, "");
  mainWorkspace.import(lumiRatio);

  // Import luminosity inside workspace
  mainWorkspace.factory(TString::Format("lumi_e[%.5f]", lumi_e));
  mainWorkspace.factory(TString::Format("lumi_mu[%.5f]", lumi_mu));
  mainWorkspace.factory(TString::Format("error_lumi_e[%.5f]", lumi_e * s_lumi_mu_percent));
  mainWorkspace.factory(TString::Format("error_lumi_mu[%.5f]", lumi_mu * s_lumi_mu_percent));

  if (! combine) {
    std::cout << "Efficiencies ratio x lumis ratio: " << sel_eff_e[btag] / sel_eff_mu[btag] * lumi_e / lumi_mu << std::endl;
  }

  // Read config file for global pdf
  std::string configFile = (combine)
    ? TString::Format("combined_%s_global_pdf.script", scriptNamePrefix.c_str()).Data()
    : "individual_btag_global_pdf.script";

  // Define cross section
  
  RooConstVar sigma("sigma", "theoretical cross-section", 1);
  mainWorkspace.import(sigma);

  parseConfigFile(configFile, /*mainCategory,*/ mainWorkspace);

  // mu is the signal strength, sigma_signal / sigma_theoretical
  RooRealVar& mu = *mainWorkspace.var("mu");

  //FIXME
  RooRealVar& nSig_mu = *mainWorkspace.var("nSig_mu");
  RooRealVar& nSig_e = *mainWorkspace.var("nSig_e");

  //FIXME: Correctly handle 1+2 btag case
  // Retrieve contraints
  RooRealVar* eff_e_constrained = mainWorkspace.var("eff_e_constrained");
  if (eff_e_constrained) {
    eff_e_constrained->setMax(1);
    eff_e_constrained->setMin(0);
    eff_e_constrained->setVal(sel_eff_e[btag]);
    //eff_e_constrained->setConstant(true);
  }

  RooRealVar* eff_mu_constrained = mainWorkspace.var("eff_mu_constrained");
  if (eff_mu_constrained) {
    eff_mu_constrained->setMax(1);
    eff_mu_constrained->setMin(0);
    eff_mu_constrained->setVal(sel_eff_mu[btag]);
    //eff_mu_constrained->setConstant(true);
  }

  RooRealVar* lumi_e_constrained = mainWorkspace.var("lumi_e_constrained");
  if (lumi_e_constrained) {
    lumi_e_constrained->setMax(lumi_e + 10 * s_lumi_mu_percent * lumi_e);
    lumi_e_constrained->setMin(lumi_e - 10 * s_lumi_mu_percent * lumi_e);
    lumi_e_constrained->setVal(lumi_e);
    //lumi_e_constrained->setConstant(true);
  }

  RooRealVar* lumi_mu_constrained = mainWorkspace.var("lumi_mu_constrained");
  if (lumi_mu_constrained) {
    lumi_mu_constrained->setMax(lumi_mu + 10 * lumi_mu * s_lumi_mu_percent);
    lumi_mu_constrained->setMin(lumi_mu - 10 * lumi_mu * s_lumi_mu_percent);
    lumi_mu_constrained->setVal(lumi_mu);
    //lumi_mu_constrained->setConstant(true);
  }

  // mtt global PDFs
  /*RooRealVar nSig_mu("nSig_mu", "number of sig events", 0., -1000., 2000.);
  RooFormulaVar nSig_e("nSig_e", "number of sig events", "nSig_mu*effRatio*lumiRatio", RooArgList(nSig_mu, effRatio, lumiRatio));
  RooRealVar nBkg_mu("nBkg_mu", "number of bkg events", 2000., 0., 50000);
  RooRealVar nBkg_e("nBkg_e", "number of bkg events", 2000., 0., 50000);*/

  //RooAddPdf globalPdf_mu("global_PDF_mu", "sigPdf+BkgPdf", RooArgList(*sigPdf_mu, backgroundPdfs["muon"]->getPdf()), RooArgList(nSig_mu, nBkg_mu));
  //RooAddPdf globalPdf_e("global_PDF_e", "sigPdf+BkgPdf", RooArgList(*sigPdf_e, backgroundPdfs["electron"]->getPdf()), RooArgList(nSig_e, nBkg_e));

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

    //if (muonsOnly && !TString(name).Contains("muon", TString::kIgnoreCase))
    //  continue;

    std::string workspaceName = "global_pdf_" + cleanedCategory;
    std::cout << "Looking for " << workspaceName << " inside workspace" << std::endl;
    globalPdfs[name] = mainWorkspace.pdf(workspaceName.c_str());

    //const RooAbsPdf& pdf = (bkgOnly) ? *backgroundPdfsFromWorkspace[name] : *globalPdfs[name];
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
    std::cout << "Fitting..." << std::endl;

    TFile* outputFile = nullptr;
    if (writeRootFile) {
      outputFile = TFile::Open(OUTPUT_PATH + prefix + "_fitRes_" + suffix + ".root", "RECREATE");
    }


    RooFitResult *fitResult = nullptr;

    //Set range parameters

    /*
    // Do a baground only fit, and limit parameters to +/- 10 sigmas
    fitResult = simPdfBackgroundOnly.fitTo(*datasetToFit, Save(), Optimize(0));
    fitResult->Print("v");
    delete fitResult;

    it = mainCategory.typeIterator();
    type = nullptr;
    while ((type = static_cast<RooCatType*>(it->Next()))) {
      TString name = TString::Format("background_%s", type->GetName());
      RooAbsPdf* pdf = mainWorkspace.pdf(name);
      setPdfParametersRange(RooArgSet(mtt), *pdf, 10);
    }
    */
   
    std::cout << "Background only ..." << std::endl;

    //simPdfBackgroundOnly.fitTo(*datasetToFit, [>, Optimize(0),<] Strategy(1), Minimizer(
          //"Minuit", "migrad")
        //);
    fitResult = simPdfBackgroundOnly.fitTo(*datasetToFit, Save(),/*, Optimize(0),*/ Strategy(1), Minimizer(
          "Minuit", "migrad")
        );

    fitResult->Print("v");
    delete fitResult;

    drawHistograms(mainCategory, mtt, *dataOrig, simPdfBackgroundOnly, backgroundPdfsFromWorkspace, btag, saveFigures, std::string(prefix), std::string(suffix) + "_bkg_only", false, true, nullptr);

    it = mainCategory.typeIterator();
    type = nullptr;
    //while ((type = static_cast<RooCatType*>(it->Next()))) {
      //setPdfParametersRange(RooArgSet(mtt), *simPdfBackgroundOnly.getPdf(type->GetName()), 0);
    //}

    if (bkgOnly) {
      computeChi2(mtt, simPdfBackgroundOnly, mainCategory, *dataOrig, mainWorkspace, true);
    } else {
      std::cout << "Background (floating) + signal ..." << std::endl;

      // And refit
      //simPdf.fitTo(*datasetToFit, Strategy(1), Minimizer("Minuit2", "migrad"),
          //ExternalConstraints(RooArgSet(*mainWorkspace.pdf("eff_e_constraint"), *mainWorkspace.pdf("eff_mu_constraint"), *mainWorkspace.pdf("lumi_e_constraint"), *mainWorkspace.pdf("lumi_mu_constraint")))
          //);
      fitResult = simPdf.fitTo(*datasetToFit, Save(), Strategy(1), Minimizer("Minuit2", "migrad"),
          ExternalConstraints(RooArgSet(*mainWorkspace.pdf("eff_e_constraint"), *mainWorkspace.pdf("eff_mu_constraint"), *mainWorkspace.pdf("lumi_e_constraint"), *mainWorkspace.pdf("lumi_mu_constraint")))
          );
      //fitResult = simPdf.fitTo(*datasetToFit, Save(), Strategy(1), Minimizer("Minuit2", "migrad"), Minos(mu));
      fitResult->Print("v");

      std::cout << "Done." << std::endl;

      if (saveWorkspace) {
        // Save fitted pdf and datasets in order to redo some plots
        RooWorkspace higgsWorkspace("w");

        fitResult->SetName("fit_results");
        higgsWorkspace.import(*fitResult);
        higgsWorkspace.import(*dataOrig);

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

          TString name = TString::Format("n_signal_%s", category.c_str());

          // Get number of fitted signal events
          RooAddPdf* globalPdf = dynamic_cast<RooAddPdf*>(mainWorkspace.pdf(workspaceName.c_str()));
          higgsWorkspace.import(
              *globalPdf,
              RooFit::RecycleConflictNodes(),
              RooFit::RenameVariable(workspaceName.c_str(), TString::Format("global_pdf_%s", workspace_suffix.Data()))
              );

          RooAbsPdf* pdf = NULL;
          /*
          // Import background
          name = TString::Format("background_%s", category.c_str());
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
           */


          if (massZprime == 500 || massZprime == 750 || massZprime == 1000 || massZprime == 1250 || massZprime == 1500 || massZprime == 2000) {

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

            // Import PDF up
            //workspaceFile = TString::Format("%s/frit/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), "pdfUp", massZprime, analysisName.c_str(), extractedBTag);
            //f.reset(TFile::Open(workspaceFile.Data()));
            //pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(name);
            //pdf->SetName(TString::Format("signal_%s_pdfUp", workspace_suffix.Data()));
            //higgsWorkspace.import(*pdf);

            //// Import PDF down
            //workspaceFile = TString::Format("%s/frit/%s-Zprime%d_%s_%d_btag_workspace.root", BASE_PATH.c_str(), "pdfDown", massZprime, analysisName.c_str(), extractedBTag);
            //f.reset(TFile::Open(workspaceFile.Data()));
            //pdf = static_cast<RooWorkspace*>(f->Get("w"))->pdf(name);
            //pdf->SetName(TString::Format("signal_%s_pdfDown", workspace_suffix.Data()));
            //higgsWorkspace.import(*pdf);


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

            name = TString::Format("signal_%s_jerUp", workspace_suffix.Data());
            pdf = getInterpolatedPdf(mtt, massZprime, "JERup", extractedBTag, category, name.Data());
            pdf->SetName(name);
            higgsWorkspace.import(*pdf);

            name = TString::Format("signal_%s_jerDown", workspace_suffix.Data());
            pdf = getInterpolatedPdf(mtt, massZprime, "JERdown", extractedBTag, category, name.Data());
            pdf->SetName(name);
            higgsWorkspace.import(*pdf);

            name = TString::Format("signal_%s_puUp", workspace_suffix.Data());
            pdf = getInterpolatedPdf(mtt, massZprime, "puUp", extractedBTag, category, name.Data());
            pdf->SetName(name);
            higgsWorkspace.import(*pdf);

            name = TString::Format("signal_%s_puDown", workspace_suffix.Data());
            pdf = getInterpolatedPdf(mtt, massZprime, "puDown", extractedBTag, category, name.Data());
            pdf->SetName(name);
            higgsWorkspace.import(*pdf);

            name = TString::Format("signal_%s_pdfUp", workspace_suffix.Data());
            pdf = getInterpolatedPdf(mtt, massZprime, "pdfUp", extractedBTag, category, name.Data());
            pdf->SetName(name);
            higgsWorkspace.import(*pdf);

            name = TString::Format("signal_%s_pdfDown", workspace_suffix.Data());
            pdf = getInterpolatedPdf(mtt, massZprime, "pdfDown", extractedBTag, category, name.Data());
            pdf->SetName(name);
            higgsWorkspace.import(*pdf);
          }
        }

        TString outputFileName = TString::Format("zprime_%d_workspace_after_fit.root", massZprime);
        //if (outputFile.length() > 0)
        //outputFileName = outputFile;

        higgsWorkspace.writeToFile(outputFileName);
      }

      drawHistograms(mainCategory, mtt, *dataOrig, simPdf, backgroundPdfsFromWorkspace, btag, saveFigures, std::string(prefix), std::string(suffix), !bkgOnly, true, outputFile);

      if (SAVE_SIGMA) {
        // Save background parameters for background systematics
        TString backgroundPdfWorkspaceFile = TString::Format("%s/background_parameters_%d_%d_btag.root", BASE_PATH.c_str(), massZprime, btag);
        RooWorkspace backgroundPdfWorkspace("w", "Background PDF workspace");
        it = mainCategory.typeIterator();
        type = nullptr;
        while ((type = static_cast<RooCatType*>(it->Next()))) {
          backgroundPdfWorkspace.import(*backgroundPdfsFromWorkspace[type->GetName()]->getParameters(RooArgSet(mtt)));
        }

        backgroundPdfWorkspace.writeToFile(backgroundPdfWorkspaceFile);
      }

      if (doBackgroundSyst) {
        delete fitResult;

        // Fix background parameters
        it = mainCategory.typeIterator();
        type = nullptr;
        while ((type = static_cast<RooCatType*>(it->Next()))) {
          setPdfParametersConst(mtt, *backgroundPdfsFromWorkspace[type->GetName()], true);
        }

        RooRealVar* var = static_cast<RooRealVar*>(mainWorkspace.var(backgroundParameterName.c_str()));
        if (! var) {
          std::cout << "ERROR: Someone asked me to set the value of '" << backgroundParameterName << "', but I can't find it..." << std::endl;
        } else {
          var->setVal(backgroundParameterValue);
          std::cout << "Value of '" << backgroundParameterName << "' set to " << backgroundParameterValue << std::endl;
        }

        mu.setVal(0);

        std::cout << "Background (fixed) + signal ..." << std::endl;

        simPdf.fitTo(*datasetToFit, Optimize(0), Strategy(1));
        simPdf.fitTo(*datasetToFit, Optimize(0), Strategy(1));
        fitResult = simPdf.fitTo(*datasetToFit, Save(), Optimize(0), Strategy(1));
        fitResult->Print("v");

        std::cout << "Done." << std::endl;
      }

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

      std::cout << "mu: " << mu.getVal() << std::endl;
      double sigma_signal = mu.getVal() * sigma.getVal();
      double sigma_error = mu.getError() / mu.getVal() * sigma_signal;

      double Limit_Z_obs_pb = sigma_signal + 2. * sigma_error;

      if (!doLikScan) {
        cout << "The signal cross section is " << sigma_signal << " +/- " << sigma_error << " pb" << endl;
        cout << "The estimated 95% C.L. upper limit on the signal cross section is " << Limit_Z_obs_pb << " pb" << endl;

        //if (DO_SYST_COMPUTATION) {

        //if (systCB != "none") {
        //// Save the systematics
        //saveSystematicParameter(massZprime, btag, "signal", systCB, systCBsign, nSig.getVal(), sigmaZ, chi2["combined"], fitResult);
        //} else if (syst_str != "nominal") {
        //saveSystematicParameter(massZprime, btag, "jec", syst_str, "", nSig.getVal(), sigmaZ, chi2["combined"], fitResult);
        //} else {
        //saveSystematicParameter(massZprime, btag, "background", fitConfigurationFile, "", nSig.getVal(), sigmaZ, chi2["combined"], fitResult);
        //}
        //}

        if (SAVE_SIGMA) {
          saveSigma(massZprime, btag, sigma_signal, chi2["combined"], 0, fitResult);
        }

        //if (useSharedMemory) {

        //SHMFitResults results = {nSig.getVal(), nSig.getError(), sigmaZ, sqrt(errorqtot_pb), chi2["combined"], fitResult->covQual(), fitResult->status()};

        //int shmid;
        //if ((shmid = shmget(shm_key, sizeof(SHMFitResults) * 8, 0666)) < 0) {
        //perror("shmget");
        //exit(1);
        //}

        //void* shm = NULL;
        //if ((shm = shmat(shmid, NULL, 0)) == (void *) -1) {
        //perror("Can't map shared memory to local memory");
        //exit(1);
        //}

        //memcpy(shm, (void*) &results, sizeof(SHMFitResults));

        //shmdt(shm);
        //}

        if (temporaryResultFile.length() > 0) {
          saveTemporaryResult(temporaryResultFile, sigma_signal);
        }

      } else {

        LikelihoodResults results;
        //FIXME. Steps is 2
        //doLikelihoodScan(*datasetToFit, simPdf, nSig, massZprime, fitResult->minNll(), 2, systematics_error_events, results);

        TFile likelihoodFile(OUTPUT_PATH + prefix + "_likscan_" + suffix + ".root", "recreate");
        likelihoodFile.cd();
        results.likscan->Write();
        results.pdfscan->Write();
        results.pdfscan_cut->Write();
        results.pdfscan_wsyst->Write();
        results.pdfscan_wsyst_cut->Write();
        likelihoodFile.Close();

        TCanvas * tmpCanvas = new TCanvas("tmpCanvas", "", 600, 600);
        tmpCanvas->SetLeftMargin(0.18);

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

        //cout << "The signal cross section is " << sigma_signal << " +/- " << sigma_error << " pb" << endl;
        //cout << "The estimated 95% C.L. upper limit on the Zprime cross section is " << Limit_Z_obs_pb << " pb" << endl;
        //cout << "95% prob. limit from scan " << results.scan_limit / (total_efficiency * lumi_mu * br_semil) << endl;
        //cout << "95% prob. limit from scan with systematics " << results.scan_wsyst_limit / (total_efficiency * lumi_mu * br_semil) << endl;
        //cout << "95% prob. limit from scan without systematics and prior sigma>0. " << results.scan_cut_limit / (total_efficiency * lumi_mu * br_semil) << endl;
        //cout << "95% prob. limit from scan with systematics and prior sigma>0. " << results.scan_wsyst_cut_limit / (total_efficiency * lumi_mu * br_semil) << endl;

        /*
           ofstream outlikscan(OUTPUT_PATH + prefix + "_likscan.txt");
           outlikscan << "The Zprime cross section is " << sigmaZ << " +- " << sqrt(errorqtot_pb) << " pb" << endl;
           outlikscan << "The 95% C.L. upper limit on the Zprime cross section is " << Limit_Z_obs_pb << " pb" << endl;
           outlikscan << "Systematics error used for gaussian smearing: " << err_sys_events << std::endl;
           outlikscan << "95% prob. limit from scan " << results.scan_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
           outlikscan << "95% prob. limit from scan with systematics " << results.scan_wsyst_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
           outlikscan << "95% prob. limit from scan with systematics and prior sigma>0. " << results.scan_wsyst_cut_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
           outlikscan.close();
           */

        //saveLikelihoodResults(massZprime, btag, results, total_efficiency * lumi_mu * br_semil);

        results.release();
      }
      SAFE_DELETE(fitResult);
    }
  }

  if (doLimitCurve)
  {

    std::cout << "Starting toys analysis..." << std::endl;

    // 0 means a random seed
    RooRandom::randomGenerator()->SetSeed(0);

    // Create PDF for toys generation
    std::map<std::string, std::shared_ptr<BaseFunction>> backgroundPdfsForToys = getCategoriesPdf(BASE_PATH + "/fit_configuration", fitConfigurationFile, mtt, NULL, massZprime, "background", mainCategory, nullptr, "toy");

    std::map<std::string, std::shared_ptr<RooAbsPdf>> globalPdfsForToys;
    std::map<int, std::shared_ptr<RooSimultaneous>> simPdfsForGeneration;
    if (combine) {
      for (int i = minBTag; i <= maxBTag; i++) {
        TString name = TString::Format("simPdfToyForGeneration_%d_btag", i);
        simPdfsForGeneration[i] = std::shared_ptr<RooSimultaneous>(new RooSimultaneous(name, "simultaneous pdf for toys generation", lepton_type));
      }
    } else {
      simPdfsForGeneration[btag] = std::shared_ptr<RooSimultaneous>(new RooSimultaneous("simPdfToyForGeneration", "simultaneous pdf for toys generation", lepton_type));
    }
    
    // Create PDF for fitting
    RooSimultaneous simPdfToyFit("simPdfToyFit", "simultaneous pdf for toys fitting", mainCategory);
    RooSimultaneous simPdfToyFitBackgroundOnly("simPdfToyFitBackgroundOnly", "simultaneous pdf for toys fitting (background only)", mainCategory);

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
    TH1F* hLimit_Z = new TH1F("hLimit_Z", "limit on the Z cross section", 500, 0, 7);
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
      simPdfToyFitBackgroundOnly.addPdf(*backgroundPdfsFromWorkspace[category], category.c_str());
    }

    std::cout << "Generating " << nEventsToy << " toys events" << std::endl;

    // Create parameters for toys generation
    std::map<int, std::shared_ptr<RooAbsPdf::GenSpec>> genSpecs;

    if (! combine) {
      genSpecs[btag] = std::shared_ptr<RooAbsPdf::GenSpec>(simPdfsForGeneration[btag]->prepareMultiGen(RooArgSet(mtt, lepton_type), NumEvents(nEventsToy), Extended(true)/*, Verbose(true)*/));
    } else {
      for (int i = minBTag; i <= maxBTag; i++) {
        genSpecs[i] = std::shared_ptr<RooAbsPdf::GenSpec>(simPdfsForGeneration[i]->prepareMultiGen(RooArgSet(mtt, lepton_type, btagCategory), /*NumEvents(nEventsToy),*/ Extended(true)/*, Verbose(true)*/));
      }
    }

    // NLL for fitting
    RooAbsReal* nll = NULL;
    RooAbsReal* nll_background = NULL;

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
        toyData = new RooDataSet("toy_dataset", "toy dataset", RooArgSet(mtt, lepton_type, btagCategory));
        int btagIndex = 0;
        for (int j = minBTag; j <= maxBTag; j++) {
          btagCategory.setIndex(btagIndex++);
          toyDataForEachBTag[j] =  simPdfsForGeneration[j]->generate(*genSpecs[j]);
          toyData->append(*toyDataForEachBTag[j]);
        }

        toyData->table(superCategory)->Print("v");
      }
      std::shared_ptr<RooDataHist> binnedDatasetForToys = std::shared_ptr<RooDataHist>(toyData->binnedClone());
      std::cout << "done." << std::endl;
      std::cout << "Dataset entries: " << toyData->numEntries() << std::endl;
      binnedDatasetForToys->Print();

      /*
      TFile* myFile = TFile::Open(OUTPUT_PATH + prefix + "_toylimit_" + suffix + "_" + indexJob + "_" + (Long_t) i + ".root", "RECREATE");
      drawHistograms(mainCategory, mtt, nBins, *toyData, simPdfToyFit, backgroundPdfsForToys, btag, false, std::string(prefix), std::string(suffix), false, false, myFile, true);
      myFile->Close();
      delete myFile;
      */
      if (fixBackground) {
        // Release background parameter for the fit
        it = mainCategory.typeIterator();
        type = nullptr;
        while ((type = static_cast<RooCatType*>(it->Next()))) {
          setPdfParametersConst(mtt, *backgroundPdfsFromWorkspace[type->GetName()], false);
        }
      }


      if (nll == NULL)
      {
        // Only create the nll the first time
        nll = simPdfToyFit.createNLL(*binnedDatasetForToys, RooFit::Optimize(0));
        nll_background = simPdfToyFitBackgroundOnly.createNLL(*binnedDatasetForToys, RooFit::Optimize(0));
      }
      else
      {
        nll->setData(*binnedDatasetForToys);
        nll_background->setData(*binnedDatasetForToys);
      }

      // Be sure that nSig_mu is not fixed anymore, and reset to 0
      mu.setVal(0);
      mu.setConstant(false);

      // Background fit only
      std::cout << "Fitting distribution (background only) ..." << std::endl;
      RooMinuit* minimizer = new RooMinuit(*nll_background);
      minimizer->setStrategy(1);
      minimizer->setEvalErrorWall(0);
      minimizer->optimizeConst(0);
      minimizer->migrad();
      minimizer->migrad();
      minimizer->migrad();

      RooFitResult* tmp = minimizer->save();
      tmp->Print("v");
      delete tmp;

      std::cout << "done.";
      delete minimizer;

      // Fit
      std::cout << "Fitting distribution (background floating + signal) ..." << std::endl;
      minimizer = new RooMinuit(*nll);
      minimizer->setStrategy(1);
      minimizer->setEvalErrorWall(0);
      minimizer->optimizeConst(0);
      minimizer->migrad();
      minimizer->migrad();
      minimizer->migrad();
      std::cout << "done.";

      std::cout << " Minos:" << std::endl;

      // Only compute errors for nSig
      minimizer->minos(RooArgSet(mu));

      std::cout << "done." << std::endl;

      RooFitResult* toyFitRes = minimizer->save();
      toyFitRes->Print("v");

      if (fixBackground) {
        // Set background parameter constant for the likelihood scan
        it = mainCategory.typeIterator();
        type = nullptr;
        while ((type = static_cast<RooCatType*>(it->Next()))) {
          setPdfParametersConst(mtt, *backgroundPdfsFromWorkspace[type->GetName()], true);
        }
      }

      //double nSigVal   = nSig.getVal();
      //double nSigErrHi = nSig.getAsymErrorHi();
      //double nSigErrLo = nSig.getAsymErrorLo();

      //hResidual->Fill(nSigVal);
      //hHiErr->Fill(nSigErrHi);
      //hLoErr->Fill(-1. * nSigErrLo);

      double nSigErr;

      //if (nSigVal < 0.)
        //nSigErr = nSigErrHi;
      //else
        //nSigErr = -1. * nSigErrLo;

      //if (nSigErr > 0.)      
        //hSigFrac->Fill(nSigVal / nSigErr);
      //else
        //hSigFrac->Fill(0);
      

      //double Limit = (nSigVal + 2. * nSigErrHi);
      //hLimit->Fill(Limit);

      //if (nSigVal == 0)
        //nSigVal = 0.00000001;

      //double sigmaZl = nSigVal / (total_efficiency * lumi_mu * br_semil);
      //double errorqstatl_pb = nSigErrHi * nSigErrHi / (total_efficiency * lumi_mu * br_semil * total_efficiency * lumi_mu * br_semil);
      //double errorqtotl_pb = errorqstatl_pb + s_yield_pb[btag] * s_yield_pb[btag] + b_tagging_corr_error_pb * b_tagging_corr_error_pb + s_eff_mu_pb[btag] * s_eff_mu_pb[btag] + s_lumi_mu_pb * s_lumi_mu_pb;
      //double Limit_Z = sigmaZl + 2. * sqrt(errorqtotl_pb);

      //hSigma->Fill(sigmaZl / sqrt(errorqstatl_pb));
      //std::cout << "sigmaZ = " << sigmaZl << " +- " << sqrt(errorqtotl_pb) << " pb" << std::endl;

      //if (!doLikScanInToys)
      //{
        //// fill "normally" with mean value+2*sigma
        //hLimit_Z->Fill(Limit_Z);
      //}
      //else
      //{

        //LikelihoodResults results;
        ////FIXME: It's 10 steps
        //doLikelihoodScan(*binnedDatasetForToys, simPdfToyFit, nSig, massZprime, toyFitRes->minNll(), 10, systematics_error_events, results);

        //TString dirName = TString::Format("likscans_%s_toy_%d", indexJob.Data(), i);
        //toyResFile->mkdir(dirName);
        //toyResFile->cd(dirName);

        //results.likscan->Write();
        //results.pdfscan->Write();
        //results.pdfscan_wsyst->Write();
        //results.pdfscan_wsyst_cut->Write();

        //toyResFile->cd();

        ////cout << "i fill with: " << results.scan_wsyst_cut_limit / (combined_efficiency * lumi_mu * br_semil) << endl;
        //hLimit_Z->Fill(results.scan_wsyst_cut_limit / (total_efficiency * lumi_mu * br_semil));

        //results.release();
      //}

      //hMinNll->Fill(toyFitRes->minNll());

      //if (toyFitRes->minNll() > fitMinNll)
        //nWorseFit++;

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

