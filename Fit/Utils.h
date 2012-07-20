#pragma once

#include <sys/stat.h>

#include <sstream>
#include <fstream>
#include <json/json.h>

#include <string.h>
#include <cstdarg>

extern double b_tagging_correction;
extern double b_tagging_corr_error_relative;

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
extern std::string getSignalPdfName(const std::string& base);
extern std::string getSignalPdfName();
extern std::string getFitBackgroundPdfName(const std::string& base);
extern std::string getFitBackgroundPdfName();
extern std::string getFritBackgroundPdfName(const std::string& base);
extern std::string getFritBackgroundPdfName();
extern char** getSystCLParameters(const std::string& mass, bool muonsOnly, int btag, ...);

extern "C" double computeEfficiency(double selEfficiency, double hltEfficiency);
