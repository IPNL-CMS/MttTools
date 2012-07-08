#pragma once

#include <sstream>
#include <fstream>
#include <json/json.h>

#include <string.h>
#include <cstdarg>

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
