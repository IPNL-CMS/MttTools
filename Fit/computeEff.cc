#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <json/json.h>
#include <tclap/CmdLine.h>

#include <TString.h>

#include "Utils.h"

#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TCanvas.h>

std::string base_path = "";

void loadSelection(const std::string& syst, int btag, const std::vector<int>& masses, std::vector<float> &nSelectionMu, std::vector<float>& errNSelectionMu, std::vector<float> &nSelectionE, std::vector<float>& errNSelectionE) {

  Json::Reader reader;
  Json::Value root;
  std::ifstream file((base_path + "/frit_efficiencies.json").c_str());
  bool success = reader.parse(file, root);
  file.close();

  if (! success) {
    std::cerr << "ERROR: Can't parse " << "frit_efficiencies.json" << ". Exiting" << std::endl;
    exit(1);
  }

  nSelectionMu.resize(masses.size());
  nSelectionE.resize(masses.size());
  errNSelectionMu.resize(masses.size());
  errNSelectionE.resize(masses.size());

  root = root[getAnalysisUUID()];

  std::stringstream ss;
  ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  for (auto& member: root.getMemberNames()) {
    int index = -1;
    for (uint32_t i = 0; i < masses.size(); i++) {
      if (masses[i] == std::stoi(member)) {
        index = i;
        break;
      }
    }

    if (index == -1) {
      std::cerr << "ERROR: mass " << member << " is not supported" << std::endl;
      continue;
    }

    Json::Value massNode = root[member][btagStr];

    if (! massNode.isMember(syst)) {
      std::cerr << "ERROR: '" << syst << "' not found for m=" << member << " in JSON file. Setting to 0." << std::endl;
      
      nSelectionMu[index] = 0;
      errNSelectionMu[index] = 0;
      nSelectionE[index] = 0;
      errNSelectionE[index] = 0;

      continue;
    }

    Json::Value systNode = massNode[syst];

    nSelectionMu[index] = systNode["muon"]["integral"].asDouble();
    errNSelectionMu[index] = systNode["muon"]["error"].asDouble();

    nSelectionE[index] = systNode["electron"]["integral"].asDouble();
    errNSelectionE[index] = systNode["electron"]["error"].asDouble();
  }
}

class Efficiencies {
  public: 
    bool   isInterpolated;
    int    mass;

    float effTrig_mu;
    float effTrig_e;
    float error_effTrig_mu;
    float error_effTrig_e;

    float selectionEff_mu;
    float selectionEff_e;
    float error_selectionEff_mu;
    float error_selectionEff_e;

    void copy(const Efficiencies& from) {
      effTrig_mu = from.effTrig_mu;
      effTrig_e  = from.effTrig_e;
      error_effTrig_mu = from.error_effTrig_mu;
      error_effTrig_e  = from.error_effTrig_e;

      selectionEff_mu  = from.selectionEff_mu;
      selectionEff_e   = from.selectionEff_e;
      error_selectionEff_mu = from.error_selectionEff_mu;
      error_selectionEff_e  = from.error_selectionEff_e;
    }

    Efficiencies(int m) {
      effTrig_mu = 1.; effTrig_e = 1.;
      selectionEff_mu = 0; selectionEff_e = 0;

      mass = m;
      isInterpolated = (mass != 500 && mass != 700);
    }

    Efficiencies() {

    }

    Efficiencies(const Efficiencies& from) {
      operator=(from);
    }

    Efficiencies& operator=(const Efficiencies& from) {
      isInterpolated = from.isInterpolated;
      mass = from.mass;

      copy(from);

      return *this;
    }

    Json::Value getAsJSON() {
      
      Json::Value array(Json::objectValue);
      array["eff_mu"] = selectionEff_mu;
      array["eff_e"] = selectionEff_e;
      array["error_eff_mu"] = error_selectionEff_mu;
      array["error_eff_e"] = error_selectionEff_e;

      return array;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Efficiencies& eff);
};

std::ostream& operator<<(std::ostream& stream, const Efficiencies& eff) {
  stream << "Selection efficiency for M = " << eff.mass<< std::endl;
  stream << '\t' << " - semi-mu channel: " << eff.selectionEff_mu * 100 << "% +/- " << eff.error_selectionEff_mu * 100 << "%" << std::endl;
  stream << '\t' << " - semi-e channel:  " << eff.selectionEff_e * 100 << "% +/- " << eff.error_selectionEff_e * 100 << "%" << std::endl;
  stream << std::endl;

  return stream;
}

int main(int argc, char** argv) {

  try {
    TCLAP::CmdLine cmd("Compute efficiencies", ' ', "0.1");

    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);

    TCLAP::ValueArg<std::string> jecArg("", "jec", "JEC", false, "nominal", "nominal/JECup/JECdown", cmd);
    TCLAP::ValueArg<std::string> jerArg("", "jer", "JER", false, "nominal", "nominal/JECup/JECdown", cmd);
    TCLAP::ValueArg<std::string> puArg("", "pu", "PU", false, "nominal", "nominal/JECup/JECdown", cmd);
    TCLAP::ValueArg<std::string> pdfArg("", "pdf", "PDF", false, "nominal", "nominal/JECup/JECdown", cmd);

    // TCLAP::SwitchArg ignoreInterpolatedArg("", "ignore-interpolated", "Ignore interpolated mass for extrapolation", cmd, false);

    cmd.parse(argc, argv);

    int btag = btagArg.getValue();

    std::string jec = jecArg.getValue();
    std::string jer = jerArg.getValue();
    std::string pu = puArg.getValue();
    std::string pdf = pdfArg.getValue();

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
    bool ignoreInterpolated = !analysisUseInterpolation();

    std::map<int, Efficiencies> efficiencies;
    //for (int i = 500; i <= 2000; i += 125) {
      //efficiencies[i] = Efficiencies(i);
    //}
    
    efficiencies = {
      {500, Efficiencies(500)},
      {700, Efficiencies(700)}
    };

    std::vector<int> M = {500, 700};

    //--- selection efficiencies
    std::vector<float> Ngen = {
      259982,
      379963
    };

    std::vector<float> Nsel_mu;
    std::vector<float> ErrNsel_mu;
    std::vector<float> Nsel_e;
    std::vector<float> ErrNsel_e;
    loadSelection(syst, btag, M, Nsel_mu, ErrNsel_mu, Nsel_e, ErrNsel_e);

    for (uint32_t i = 0; i <M.size(); i++) {
      std::cout << "Number of selected events for M=" << M[i] << std::endl;
      std::cout << "\t- semi-mu channel: " << Nsel_mu[i] << " +/- " << ErrNsel_mu[i] << std::endl;
      std::cout << "\t- semi-e channel: " << Nsel_e[i] << " +/- " << ErrNsel_e[i] << std::endl;
      std::cout << std::endl;
    }

    TGraphErrors e_mu;
    TGraphErrors e_mu_low;
    TGraphErrors e_mu_high;

    TGraphErrors e_e;
    TGraphErrors e_e_low;
    TGraphErrors e_e_high;

    TF1 selectionEff_fit_mu("sel_eff_fit_mu", "pol3", 750, 2000);
    TF1* selectionEff_fit_mu_low = (TF1*) selectionEff_fit_mu.Clone("sel_eff_fit_mu_low");
    TF1* selectionEff_fit_mu_high = (TF1*) selectionEff_fit_mu.Clone("sel_eff_fit_mu_high");

    TF1 selectionEff_fit_e("sel_eff_fit_e", "pol3", 750, 2000);
    TF1* selectionEff_fit_e_low = (TF1*) selectionEff_fit_e.Clone("sel_eff_fit_e_low");
    TF1* selectionEff_fit_e_high = (TF1*) selectionEff_fit_e.Clone("sel_eff_fit_e_high");

    int index = 0;
    Efficiencies* lowMass_eff = nullptr;
    for (auto& i: efficiencies) {
      Efficiencies& eff = i.second;

      if (! i.second.isInterpolated) {
        lowMass_eff = &eff;

        eff.selectionEff_mu = Nsel_mu[index] / Ngen[index];
        eff.selectionEff_e  = Nsel_e[index] / Ngen[index];
        eff.error_selectionEff_mu = std::sqrt( eff.selectionEff_mu * (1 - eff.selectionEff_e) / Ngen[index] );
        eff.error_selectionEff_e = std::sqrt( eff.selectionEff_e * (1 - eff.selectionEff_e) / Ngen[index] );

        if (! ignoreInterpolated && i.first != 500) {
          e_mu.SetPoint(index, i.first, i.second.selectionEff_mu);
          e_mu.SetPointError(index, 0, i.second.error_selectionEff_mu);

          e_mu_low.SetPoint(index, i.first, i.second.selectionEff_mu - i.second.error_selectionEff_mu);
          e_mu_high.SetPoint(index, i.first, i.second.selectionEff_mu + i.second.error_selectionEff_mu);

          e_e.SetPoint(index, i.first, i.second.selectionEff_e);
          e_e.SetPointError(index, 0, i.second.error_selectionEff_e);

          e_e_low.SetPoint(index, i.first, i.second.selectionEff_e - i.second.error_selectionEff_e);
          e_e_high.SetPoint(index, i.first, i.second.selectionEff_e + i.second.error_selectionEff_e);
        }

        index++;
      } else if (ignoreInterpolated) {
        eff.copy(*lowMass_eff);
      }
    }

    if (! ignoreInterpolated) {
      e_mu.Fit(&selectionEff_fit_mu, "QR");
      e_mu_low.Fit(selectionEff_fit_mu_low, "QR");
      e_mu_high.Fit(selectionEff_fit_mu_high, "QR");

      e_e.Fit(&selectionEff_fit_e, "QR");
      e_e_low.Fit(selectionEff_fit_e_low, "QR");
      e_e_high.Fit(selectionEff_fit_e_high, "QR");
    }

    if (! ignoreInterpolated) {
      for (auto& i: efficiencies) {
        Efficiencies& eff = i.second;

        if (i.second.isInterpolated) {
          eff.selectionEff_mu = selectionEff_fit_mu.Eval(i.first);
          eff.selectionEff_e  = selectionEff_fit_e.Eval(i.first);

          eff.error_selectionEff_mu = fabs(selectionEff_fit_mu_high->Eval(i.first) - selectionEff_fit_mu_low->Eval(i.first)) / 2.;
          eff.error_selectionEff_e = fabs(selectionEff_fit_e_high->Eval(i.first) - selectionEff_fit_e_low->Eval(i.first)) / 2.; 
        }
      }
    }

    Json::Value root;
    getJsonRoot(base_path + "/efficiencies.json", root, false);

    for (auto& i: efficiencies) {

      if (ignoreInterpolated && i.second.isInterpolated)
        continue;

      std::stringstream ss;
      ss << i.first;
      std::string mass = ss.str();

      root[getAnalysisUUID()][mass][btagStr][syst] = i.second.getAsJSON();
      std::cout << i.second << std::endl;
      std::cout << std::endl;
    }

    Json::StyledWriter writer;
    std::ofstream output(base_path + "/efficiencies.json");
    output << writer.write(root);
    output.close();
    std::cout << "Efficiencies saved as 'efficiences.json'" << std::endl;

    if (jec == "nominal" && pu == "nominal" && jer == "nominal" && pdf == "nominal") {
      TString noteFilename = TString::Format("%s/efficiencies_table_%s_%d_btag.tex", base_path.c_str(), getAnalysisName().c_str(), btagArg.getValue());

      // table latex pour la note :
      std::ofstream latex(noteFilename);
      latex << "\\mtt";

      for (auto& i: efficiencies) {
        if (i.second.isInterpolated)
          continue;
        latex << " & " << i.first << " GeV";
      }

      latex << "\\\\" << std::endl << "\\hline" << std::endl;
      latex << std::setiosflags(std::ios::fixed) << std::setprecision(2);

      latex << "$\\epsilon(Z^{\\prime}), semi-mu$ (\\%)";
      for (auto& i: efficiencies) {
        if (i.second.isInterpolated)
          continue;
        latex << " & " << i.second.selectionEff_mu * 100 << " $\\pm$ " << std::setprecision(4) << i.second.error_selectionEff_mu * 100 << std::setprecision(2);
      }
      latex << "\\\\" << std::endl;

      latex << "$\\epsilon(Z^{\\prime}), semi-e$ (\\%)"; 
      for (auto& i: efficiencies) {
        if (i.second.isInterpolated)
          continue;
        latex << " & " << i.second.selectionEff_e * 100 << " $\\pm$ " << std::setprecision(4) << i.second.error_selectionEff_e * 100 << std::setprecision(2);
      }
      latex << "\\\\" << std::endl;

      latex.close();

      std::cout << "Latex table saved as '" << noteFilename << "'" << std::endl;
    }

  } catch (TCLAP::ArgException& e) {

  }

}
