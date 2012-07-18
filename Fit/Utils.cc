#include "Utils.h"

double b_tagging_correction = 0.90; // (0.95^2)
double b_tagging_corr_error_relative = 0.08; //0.30; 2010

double trigger_correction_muons = 1.; // PRELIMINRAY -- was 0.9946
double trigger_corr_muons_error_relative = 0.0001;
double trigger_correction_ele = 1.; //PRELIMINARY!!
double trigger_corr_ele_error_relative = 0.0001;

double muID_correction = 0.997;
double muID_correction_error_relative = 0.0012;
double muIso_correction = 0.998;
double muIso_correction_error_relative = 0.0002;

double eleID_correction = 1.002;
double eleID_correction_error_relative = 0.0004;
double eleIso_correction = 0.9999;
double eleIso_correction_error_relative = 0.0004;

bool fileExists(const std::string& filename) {
  struct stat buf;
  if (stat(filename.c_str(), &buf) != -1)
    return true;

  return false;
}

void getJsonRoot(const std::string& filename, Json::Value& root) {
  std::ifstream file(filename.c_str());
  Json::Reader reader;
  if (! reader.parse(file, root)) {
    std::cerr << "ERROR: Failed to parse " << filename << "." << std::endl;
    exit(1);
  }
  file.close();
}

std::string formatPath(const std::string& base, const std::string& filename) {
  std::stringstream ss;
  ss << base << "/" << filename;

  return ss.str();
}

std::string getSignalPdfName(const std::string& base) {
  Json::Value root;
  getJsonRoot(formatPath(base, "parameters.json"), root);

  //FIXME: Check
  return root["parameters"]["pdf"]["signal"].asString();
}

std::string getSignalPdfName() {
  return getSignalPdfName(".");
}

std::string getFitBackgroundPdfName(const std::string& base) {
  Json::Value root;
  getJsonRoot(formatPath(base, "parameters.json"), root);

  //FIXME: Check
  return root["parameters"]["pdf"]["background"]["fit"].asString();
}

std::string getFitBackgroundPdfName() {
  return getFitBackgroundPdfName(".");
}

std::string getFritBackgroundPdfName(const std::string& base) {
  Json::Value root;
  getJsonRoot(formatPath(base, "parameters.json"), root);

  //FIXME: Check
  return root["parameters"]["pdf"]["background"]["frit"].asString();
}

std::string getFritBackgroundPdfName() {
  return getFritBackgroundPdfName(".");
}

char** getSystCLParameters(const std::string& mass, bool muonsOnly, int btag, ...) {
 // fitMtt", "-m", ss.str().c_str(), "--syst", (*param).c_str(), "--syst-computation", "--no-figs", "--no-text-files", "--no-root-files", "--muons-only  
 
  va_list paramList;
  va_start(paramList, btag);

  std::vector<const char*> params = { "fitMtt", "-m", mass.c_str(), "--syst-computation", "--no-figs", "--no-text-files", "--no-root-files", "--b-tag" };

  std::stringstream ss;
  ss << btag;

  params.push_back(strdup(ss.str().c_str()));

  if (muonsOnly)
    params.push_back("--muons-only");

  char * arg = NULL;
  while ((arg = va_arg(paramList, char*))) {
    params.push_back(arg);
  }
  params.push_back(NULL);

  char** array = new char*[params.size()];
  for (size_t i = 0; i < params.size(); i++) {
    array[i] = const_cast<char*>(params[i]);
  }

  return array;
}

double computeEfficiency(double selEfficiency, double hltEfficiency) {
  return trigger_correction_muons * muID_correction  * muIso_correction  * b_tagging_correction * selEfficiency * hltEfficiency;
}
