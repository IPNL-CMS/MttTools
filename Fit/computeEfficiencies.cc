#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

#include <json/json.h>
#include <tclap/CmdLine.h>

#include <TString.h>

#include "Utils.h"

#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TCanvas.h>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

std::string base_path = "";

namespace fs = boost::filesystem;

void loadSelection(int mass, const std::string& syst, int btag, float& nSelectionMu, float& errNSelectionMu, float& nSelectionE, float& errNSelectionE) {

  Json::Reader reader;
  Json::Value root;
  std::ifstream file((base_path + "/frit_efficiencies.json").c_str());
  bool success = reader.parse(file, root);
  file.close();

  if (! success) {
    std::cerr << "ERROR: Can't parse " << "frit_efficiencies.json" << ". Exiting" << std::endl;
    exit(1);
  }

  root = root[getAnalysisUUID()];

  std::stringstream ss;
  ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  ss.str(std::string());
  ss << mass;
  std::string massStr = ss.str();

  if (! root.isMember(massStr)) {
    std::cerr << "ERROR: Mass not found inside efficiencies files" << std::endl;
    exit(1);
  }

  Json::Value massNode = root[massStr][btagStr];

  if (! massNode.isMember(syst)) {
    std::cerr << "ERROR: '" << syst << "' not found for m=" << mass << " in JSON file. Setting to 0." << std::endl;

    nSelectionMu = 0;
    errNSelectionMu = 0;
    nSelectionE = 0;
    errNSelectionE = 0;

    return;
  }

  Json::Value systNode = massNode[syst];

  nSelectionMu = systNode["muon"]["integral"].asDouble();
  errNSelectionMu = systNode["muon"]["error"].asDouble();

  nSelectionE = systNode["electron"]["integral"].asDouble();
  errNSelectionE = systNode["electron"]["error"].asDouble();
}

class Efficiency {
  public:
    int mass;

    float effTrig_mu;
    float effTrig_e;
    float error_effTrig_mu;
    float error_effTrig_e;

    float selectionEff_mu;
    float selectionEff_e;
    float error_selectionEff_mu;
    float error_selectionEff_e;

    float hasTriggerEfficiencies;

    void copy(const Efficiency& from) {
      effTrig_mu = from.effTrig_mu;
      effTrig_e  = from.effTrig_e;
      error_effTrig_mu = from.error_effTrig_mu;
      error_effTrig_e  = from.error_effTrig_e;

      selectionEff_mu  = from.selectionEff_mu;
      selectionEff_e   = from.selectionEff_e;
      error_selectionEff_mu = from.error_selectionEff_mu;
      error_selectionEff_e  = from.error_selectionEff_e;
    }

    Efficiency(int m) {
      effTrig_mu = 1.; effTrig_e = 1.;
      selectionEff_mu = 0; selectionEff_e = 0;

      error_effTrig_mu = error_effTrig_e = error_selectionEff_e = error_selectionEff_mu = 0;

      mass = m;

      hasTriggerEfficiencies = false;
    }

    Efficiency(const Efficiency& from) {
      operator=(from);
    }

    Efficiency& operator=(const Efficiency& from) {
      mass = from.mass;

      copy(from);

      return *this;
    }

    void compute(AnalysisType type, float n_gen, float n_sel_mu, float n_sel_e, float err_n_sel_mu, float err_n_sel_e) {

      selectionEff_mu = n_sel_mu / n_gen;
      selectionEff_e  = n_sel_e / n_gen;

      if (type == HIGGS) {
        // We don't perform a fit for the signal shape determination.
        // Use standard error on efficiency
        error_selectionEff_mu = std::sqrt( selectionEff_mu * (1 - selectionEff_e) / n_gen );
        error_selectionEff_e = std::sqrt( selectionEff_e * (1 - selectionEff_e) / n_gen );
      } else {
        error_selectionEff_mu = err_n_sel_mu / n_gen;
        error_selectionEff_e = err_n_sel_e / n_gen;
      }
    }

    void loadTriggerEfficiencies(const std::string& btag, fs::path filename) {

      Json::Value root;
      getJsonRoot(filename.string(), root, false);

      std::stringstream ss;
      ss << mass;
      std::string massStr = ss.str();

      const std::string btagKey = (boost::format("%sbtag") % btag).str();

      if (! root.isMember(btagKey)) {
        std::cout << "Warning: btag key '" << btagKey << "' not found inside trigger efficiencies JSON file" << std::endl;
        return;
      }

      Json::Value btagNode = root[btagKey];

      if (! btagNode.isMember(massStr)) {
        std::cout << "Warning: mass '" << mass << "' not found inside trigger efficiencies JSON file" << std::endl;
        return;
      }

      Json::Value massNode = btagNode[massStr];
      effTrig_mu = massNode["efficiency_mu"].asDouble();
      effTrig_e = massNode["efficiency_e"].asDouble();

      error_effTrig_mu = massNode["error_efficiency_mu"].asDouble();
      error_effTrig_e = massNode["error_efficiency_e"].asDouble();

      hasTriggerEfficiencies = true;
    }

    Json::Value getAsJSON() {

      Json::Value array(Json::objectValue);
      array["eff_mu"] = selectionEff_mu;
      array["eff_e"] = selectionEff_e;
      array["error_eff_mu"] = error_selectionEff_mu;
      array["error_eff_e"] = error_selectionEff_e;

      array["trigger_eff_mu"] = effTrig_mu,
      array["trigger_eff_e"] = effTrig_e;
      array["trigger_error_eff_mu"] = error_effTrig_mu;
      array["trigger_error_eff_e"] = error_effTrig_e;

      return array;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Efficiency& eff);
};

std::ostream& operator<<(std::ostream& stream, const Efficiency& eff) {
  stream << "Selection efficiency for M = " << eff.mass<< std::endl;
  stream << '\t' << " - semi-mu channel: " << eff.selectionEff_mu * 100 << "% +/- " << eff.error_selectionEff_mu * 100 << "%" << std::endl;
  stream << '\t' << " - semi-e channel:  " << eff.selectionEff_e * 100 << "% +/- " << eff.error_selectionEff_e * 100 << "%" << std::endl;

  if (eff.hasTriggerEfficiencies) {
    stream << "Trigger efficiency for M = " << eff.mass<< std::endl;
    stream << '\t' << " - semi-mu channel: " << eff.effTrig_mu * 100 << "% +/- " << eff.error_effTrig_mu * 100 << "%" << std::endl;
    stream << '\t' << " - semi-e channel:  " << eff.effTrig_e * 100 << "% +/- " << eff.error_effTrig_e * 100 << "%" << std::endl;
  }

  stream << std::endl;

  return stream;
}

float loadNumberOfGeneratedEvents(fs::path filename, int mass) {
  Json::Value root;
  getJsonRoot(filename.string(), root, false);

  std::stringstream ss;
  ss << mass;
  std::string massStr = ss.str();

  if (! root.isMember(massStr)) {
    std::cout << "ERROR: mass '" << mass << "' not found inside generated events JSON file. Exiting" << std::endl;
    exit(1);
  }

  return root[massStr].asDouble();
}

int main(int argc, char** argv) {

  try {
    TCLAP::CmdLine cmd("Compute efficiencies", ' ', "0.1");

    TCLAP::ValueArg<int> massArg("m", "mass", "Signal mass", true, 500, "int", cmd);

    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);

    TCLAP::ValueArg<std::string> triggerEffFileArg("", "trigger-eff", "JSON file containing trigger efficiencies", false, "data/zprime_trigger_efficiencies.json", "filename", cmd);
    TCLAP::ValueArg<std::string> nGenFileArg("", "gen-file", "JSON file containing number of generated events for each signal point", true, "data/", "filename", cmd);

    TCLAP::ValueArg<std::string> jecArg("", "jec", "JEC", false, "nominal", "nominal/JECup/JECdown", cmd);
    TCLAP::ValueArg<std::string> jerArg("", "jer", "JER", false, "nominal", "nominal/JECup/JECdown", cmd);
    TCLAP::ValueArg<std::string> puArg("", "pu", "PU", false, "nominal", "nominal/JECup/JECdown", cmd);
    TCLAP::ValueArg<std::string> pdfArg("", "pdf", "PDF", false, "nominal", "nominal/JECup/JECdown", cmd);

    // TCLAP::SwitchArg ignoreInterpolatedArg("", "ignore-interpolated", "Ignore interpolated mass for extrapolation", cmd, false);

    cmd.parse(argc, argv);

    int mass = massArg.getValue();
    int btag = btagArg.getValue();

    std::string jec = jecArg.getValue();
    std::string jer = jerArg.getValue();
    std::string pu = puArg.getValue();
    std::string pdf = pdfArg.getValue();

    AnalysisType analysisType = getAnalysisType();

    bool computeTriggerEfficiencies = analysisType == ZPRIME;
    fs::path triggerFilename(triggerEffFileArg.getValue());

    if (computeTriggerEfficiencies && !fs::exists(triggerFilename)) {
      std::cout << "Warning: trigger efficiencies JSON file '" << triggerFilename << "' does not exist. No trigger efficiency will be computed" << std::endl;
      computeTriggerEfficiencies = false;
    }

    fs::path genFilename(nGenFileArg.getValue());
    if (! fs::exists(genFilename)) {
      std::cout << "ERROR: JSON file containing the number of generated events does not exists." << std::endl;

      return 1;
    }

    std::string syst = "nominal";
    if (jec != "nominal" || jer != "nominal" || pu != "nominal" || pdf != "nominal") {
      if (jec != "nominal") {
        if (jec == "up")
          syst = "JECup";
        else
          syst = "JECdown";
      } else if (jer != "nominal") {
        if (jer == "up")
          syst = "JERup";
        else
          syst = "JERdown";
      } else if (pu != "nominal") {
        if (pu == "up")
          syst = "puUp";
        else
          syst = "puDown";
      } else if (pdf != "nominal") {
        if (pdf == "up")
          syst = "pdfUp";
        else
          syst = "pdfDown";
      }
    }

    std::stringstream stream;
    stream << btag;

    std::string btagStr = stream.str();

    base_path = "./analysis/" + getAnalysisUUID();

    float Nsel_mu;
    float ErrNsel_mu;
    float Nsel_e;
    float ErrNsel_e;
    loadSelection(mass, syst, btag, Nsel_mu, ErrNsel_mu, Nsel_e, ErrNsel_e);

    float Ngen = loadNumberOfGeneratedEvents(genFilename, mass);

    std::cout << "Number of selected events for M=" << mass << std::endl;
    std::cout << "\t- semi-mu channel: " << Nsel_mu << " +/- " << ErrNsel_mu << std::endl;
    std::cout << "\t- semi-e channel: " << Nsel_e << " +/- " << ErrNsel_e << std::endl;
    std::cout << std::endl;

    Efficiency eff(mass);
    eff.compute(analysisType, Ngen, Nsel_mu, Nsel_e, ErrNsel_mu, ErrNsel_e);

    if (computeTriggerEfficiencies) {
      eff.loadTriggerEfficiencies(btagStr, triggerFilename);
    }

    Json::Value root;
    getJsonRoot(base_path + "/efficiencies.json", root, false);

    {
      std::stringstream ss;
      ss << mass;
      std::string massStr = ss.str();

      root[getAnalysisUUID()][massStr][btagStr][syst] = eff.getAsJSON();
      std::cout << eff << std::endl;
      std::cout << std::endl;
    }

    Json::StyledWriter writer;
    std::ofstream output(base_path + "/efficiencies.json");
    output << writer.write(root);
    output.close();
    std::cout << "Efficiencies saved as 'efficiences.json'" << std::endl;

  } catch (TCLAP::ArgException& e) {

  }

}
