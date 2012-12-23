#pragma once

#include <sys/stat.h>

#include <sstream>
#include <fstream>
#include <json/json.h>

#include <string.h>
#include <cstdarg>

struct SHMFitResults {
  double nSignalEvents;
  double nSignalEvents_error;

  double sigma;
  double sigma_error;

  double chiSquare;
  
  int fit_coQual;
  int fit_status;
};

extern double b_tagging_scale_factor;
extern double b_tagging_scale_factor_error;

extern double trigger_scale_factor_muons;
extern double trigger_scale_factor_muons_error;

extern double trigger_scale_factor_electrons;
extern double trigger_scale_factor_electrons_error;

extern double trigger_correction_muons;
extern double trigger_corr_muons_error_relative;
extern double trigger_correction_ele;
extern double trigger_corr_ele_error_relative;

extern double muID_correction;
extern double muID_correction_error_relative;
extern double muIso_correction;
extern double muIso_correction_error_relative;
  //
extern double eleID_correction;
extern double eleID_correction_error_relative;
extern double eleIso_correction;
extern double eleIso_correction_error_relative;

extern bool fileExists(const std::string& filename);
extern void getJsonRoot(const std::string& filename, Json::Value& root, bool exitOnError = true);
extern std::string formatPath(const std::string& base, const std::string& filename);
extern char** getSystCLParameters(const std::string& mass, const std::string& file, bool singleFile, bool muonsOnly, int btag, ...);

extern int         getAnalysisId  (const std::string& base = ".");
extern std::string getAnalysisUUID(const std::string& base = ".");
extern std::string getAnalysisName(const std::string& base = ".");

extern bool        analysisUseSystematics (const std::string& base = ".");
extern bool        analysisUseSignalSyst  (const std::string& base = ".");
extern bool        analysisUseBkgSyst     (const std::string& base = ".");
extern bool        analysisUseJECSyst     (const std::string& base = ".");

extern "C" double computeEfficiencyMuons(double selEfficiency, double hltEfficiency);
extern "C" double computeEfficiencyElectrons(double selEfficiency, double hltEfficiency);
