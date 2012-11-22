#pragma once

#include <memory>

#include <RooAbsCategory.h>
#include <RooCatType.h>
#include <RooAbsCategoryLValue.h>

#include "BaseFunction.h"

#include "SignalFunctions.h"
#include "BackgroundFunctions.h"

#include "Utils.h"

std::shared_ptr<BaseFunction> getPdf(const std::string& name, const std::string& pdfName) {

  if (name == "crystalball")
    return std::shared_ptr<BaseFunction>(new CrystalBall(pdfName));

  if (name == "keyspdf")
    return std::shared_ptr<BaseFunction>(new KeysPdf(pdfName));

  if (name == "faltB")
    return std::shared_ptr<BaseFunction>(new FAltB(pdfName));

  if (name == "falt")
    return std::shared_ptr<BaseFunction>(new FAlt(pdfName));

  if (name == "exp")
    return std::shared_ptr<BaseFunction>(new Exp(pdfName));

  return nullptr;
}

/**
 * Build a cut formula for all categories
 *
 * @return A cut formula like "whichLepton == 11 && btag == 2"
 */
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

std::shared_ptr<BaseFunction> createPdf(const Json::Value& value, const std::string& inputPath, RooRealVar& observable, RooDataSet* dataset, int mass, const std::string& type, const std::string& prefix) {
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

  if (parameters != "none") {
    Json::Value parametersRoot;
    getJsonRoot(inputPath + "/pdf_parameters.json", parametersRoot);

    if (! parametersRoot.isMember(parameters)) {
      std::cout << "Error: Parameters set '" << parameters << "' not found." << std::endl;
      assert(false);
      return nullptr;
    }

    pdf->setParameters(jsonToRealVars(observable, mass, prefix, parametersRoot[parameters]));
  }

  if (dataset) {
    pdf->setDataset(dataset);
  }

  assert(pdf->createPdf(observable));

  return pdf;
}

std::shared_ptr<BaseFunction> jsonToPdf(const std::string& inputPath, RooRealVar& observable, RooDataSet* dataset, int mass, const std::string& prefix, const std::string& type, const Json::Value& root) {

  if (! root.isMember(type)) {
    std::cout << "Error: '" << type << "' node not found" << std::endl;
    assert(false);
  }

  return createPdf(root[type], inputPath, observable, dataset, mass, type, prefix);
}

/**
 * Parse a config file and create PDFs, one for each category of the analysis
 *
 * @param inputPath Where is the inputFile located
 * @param inputFile The config file (JSON format) containing informations about which pdf use for a given category
 * @param observable The main observable of the analysis
 * @param dataset The dataset of the analysis. Can be NULL if not needed by any PDF.
 * @param mass The current Zprime mass
 * @param type The PDF type you want to retrieve (signal or background)
 * @param categories List of possible category for this analysis
 * @param[out] id If not NULL, the id of the inputFile is stored
 * @param prefix A string to use as a prefix for each variable name created
 *
 * @return A map containing the PDF for each categories
 */
std::map<std::string, std::shared_ptr<BaseFunction>> getCategoriesPdf(const std::string& inputPath, const std::string& inputFile, RooRealVar& observable, RooDataSet* dataset, int mass, const std::string& type, RooAbsCategoryLValue& categories, std::string* id, const std::string& prefix = "") {

  std::map<std::string, std::shared_ptr<BaseFunction>> results;

  Json::Value root;
  getJsonRoot(inputPath + "/" + inputFile, root);

  if (id) {
    *id = root["id"].asString();
  }

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

    RooDataSet* reducedDataset = nullptr;
    if (dataset) {
      // Reduce dataset to keep only data for this category
      categories = name.c_str();
      std::string cutFormula = buildCutFormula(categories);
      reducedDataset = static_cast<RooDataSet*>(dataset->reduce(cutFormula.c_str()));
    }

    auto pdf = jsonToPdf(inputPath, observable, reducedDataset, mass, pdfsPrefix, type, root[name]);
    assert(pdf.get());

    results[name] = pdf;
  }

  return results;
}
