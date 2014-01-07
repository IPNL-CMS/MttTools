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

enum AnalysisType {
  HIGGS,
  ZPRIME
};

extern double b_tagging_scale_factor;
extern double b_tagging_scale_factor_error;

extern double b_tagging_efficiency;
extern double b_tagging_efficiency_error;

extern double trigger_scale_factor_muons;
extern double trigger_scale_factor_muons_error;

extern double trigger_scale_factor_electrons;
extern double trigger_scale_factor_electrons_error;

extern double muonID_scale_factor;
extern double muonID_scale_factor_error;
extern double muonIso_scale_factor;
extern double muonIso_scale_factor_error;

extern double electron_scale_factor;
extern double electron_scale_factor_error;

extern bool fileExists(const std::string& filename);
extern void getJsonRoot(const std::string& filename, Json::Value& root, bool exitOnError = true);
extern std::string formatPath(const std::string& base, const std::string& filename);
extern char** getSystCLParameters(const std::string& mass, const std::string& file, bool singleFile, bool fixBackground, bool muonsOnly, int btag, ...);

extern int         getAnalysisId  (const std::string& base = ".");
extern std::string getAnalysisUUID(const std::string& base = ".");
extern std::string getAnalysisName(const std::string& base = ".");

extern AnalysisType getAnalysisType(const std::string& base = ".");
extern const char* getAnalysisPrefix(const std::string& basee = ".");
