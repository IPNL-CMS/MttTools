#pragma once

#include <memory>

#include <RooAbsCategory.h>
#include <RooCatType.h>

#include "BaseFunction.h"

#include "SignalFunctions.h"
#include "BackgroundFunctions.h"

#include "Utils.h"

std::shared_ptr<BaseFunction> getPdf(const std::string& name, const std::string& pdfName) {

  if (name == "crystalball")
    return std::shared_ptr<BaseFunction>(new CrystalBall(pdfName));

  if (name == "faltB")
    return std::shared_ptr<BaseFunction>(new FAltB(pdfName));

  if (name == "falt")
    return std::shared_ptr<BaseFunction>(new FAlt(pdfName));

  if (name == "exp")
    return std::shared_ptr<BaseFunction>(new Exp(pdfName));

  return nullptr;
}

std::shared_ptr<BaseFunction> createPdf(const Json::Value& value, const std::string& inputPath, RooRealVar& observable, int mass, const std::string& type, const std::string& prefix) {
  if (! value.isMember("name") || ! value.isMember("parameters"))
    assert(false);

  std::string name = value["name"].asString();
  std::string parameters = value["parameters"].asString();

  std::shared_ptr<BaseFunction> pdf = getPdf(name, type + "_" + prefix);
  if (! pdf.get()) {
    std::cout << "Error: PDF '" << name << "' not found." << std::endl;
    assert(false);
    return nullptr;
  }

  Json::Value parametersRoot;
  getJsonRoot(inputPath + "/pdf_parameters.json", parametersRoot);

  if (! parametersRoot.isMember(parameters)) {
    std::cout << "Error: Parameters set '" << parameters << "' not found." << std::endl;
    assert(false);
    return nullptr;
  }

  pdf->setParameters(jsonToRealVars(observable, mass, prefix, parametersRoot[parameters]));
  assert(pdf->createPdf(observable));

  return pdf;
}

std::shared_ptr<BaseFunction> jsonToPdf(const std::string& inputPath, RooRealVar& observable, int mass, const std::string& prefix, const std::string& type, const Json::Value& root) {

  if (! root.isMember(type)) {
    std::cout << "Error: '" << type << "' node not found" << std::endl;
    assert(false);
  }

  return createPdf(root[type], inputPath, observable, mass, type, prefix);
}

std::map<std::string, std::shared_ptr<BaseFunction>> getCategoriesPdf(const std::string& inputPath, const std::string& inputFile, RooRealVar& observable, int mass, const std::string& type, const RooAbsCategory& categories, const std::string& prefix = "") {

  std::map<std::string, std::shared_ptr<BaseFunction>> results;

  Json::Value root;
  getJsonRoot(inputPath + "/" + inputFile, root);

  TIterator* it = categories.typeIterator();
  const RooCatType* catType = nullptr;

  while ((catType = static_cast<RooCatType*>(it->Next()))) {
    if (! catType)
      continue;

    std::string name = catType->GetName();

    if (! root.isMember(name)) {
      std::cout << "Error: Looking for Pdf '" << name << "' inside 'pdf.json', but nothing found." << std::endl;
      assert(false);
    }

    std::string pdfsPrefix = prefix;
    if (pdfsPrefix.length() != 0)
      pdfsPrefix += "_" + name;
    else
      pdfsPrefix += name;

    auto pdf = jsonToPdf(inputPath, observable, mass, pdfsPrefix, type, root[name]);
    assert(pdf.get());

    results[name] = pdf;
  }

  return results;
}
