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

void loadSelection(const std::string& jecType, int btag, const int (&masses)[4], float (&nSelectionMu)[4], float (&errNSelectionMu)[4], float (&nSelectionE)[4], float (&errNSelectionE)[4]) {

  Json::Reader reader;
  Json::Value root;
  std::ifstream file("frit_efficiencies.json");
  bool success = reader.parse(file, root);
  file.close();

  if (! success) {
    std::cerr << "ERROR: Can't parse " << "frit_efficiencies.json" << ". Exiting" << std::endl;
    exit(1);
  }

  root = root[getAnalysisUUID()];

  for (int i = 0; i < 4; i++) {
    const int& mass = masses[i];

    std::stringstream ss;
    ss << mass;
    std::string strMass = ss.str();

    ss.clear(); ss.str(std::string());
    ss << btag;
    std::string btagStr = ss.str();

    if (! root.isMember(strMass)) {
      std::cerr << "ERROR: mass '" << mass << "' not found in JSON file. Exiting." << std::endl;
      exit(1);
    }

    Json::Value massNode = root[strMass][btagStr];

    if (! massNode.isMember(jecType)) {
      std::cerr << "ERROR: '" << jecType << "' not found for m=" << mass << " in JSON file. Setting to 0." << std::endl;
      
      nSelectionMu[i] = 0;
      errNSelectionMu[i] = 0;
      nSelectionE[i] = 0;
      errNSelectionE[i] = 0;

      continue;
    }

    Json::Value jecNode = massNode[jecType];

    nSelectionMu[i] = jecNode["muon"]["events"].asDouble();
    errNSelectionMu[i] = jecNode["muon"]["error"].asDouble();

    nSelectionE[i] = jecNode["electron"]["events"].asDouble();
    errNSelectionE[i] = jecNode["electron"]["error"].asDouble();
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
      isInterpolated = (mass != 750 && mass != 1000 && mass != 1250 && mass != 1500);
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
      
      Json::Value array(Json::arrayValue);
      array.append(selectionEff_mu);
      array.append(selectionEff_e);
      array.append(effTrig_mu);
      array.append(effTrig_e);
      array.append(error_selectionEff_mu / selectionEff_mu);
      array.append(error_selectionEff_e / selectionEff_e);
      array.append(error_effTrig_mu / effTrig_mu);
      array.append(error_effTrig_e / effTrig_e);

      return array;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Efficiencies& eff);
};

std::ostream& operator<<(std::ostream& stream, const Efficiencies& eff) {
  stream << "M_Z' = " << eff.mass<< std::endl;
  stream << "Selection efficiency: "
    << '\t' << "muon = " << eff.selectionEff_mu << " +/- " << eff.error_selectionEff_mu 
    << '\t' << "electron = " << eff.selectionEff_e << " +/- " << eff.error_selectionEff_e;
  stream << std::endl;
  stream << "Trigger efficiency: "
    << '\t' << "muon = " << eff.effTrig_mu << " +/- " << eff.error_effTrig_mu 
    << '\t' << "electron = " << eff.effTrig_e << " +/- " << eff.error_effTrig_e;

  return stream;
}

int main(int argc, char** argv) {

  try {
    TCLAP::CmdLine cmd("Compute efficiencies", ' ', "0.1");

    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);
    TCLAP::ValueArg<std::string> jecArg("", "jec", "Type of JEC", false, "nominal", "nominal/JECup/JECdown", cmd);
    TCLAP::SwitchArg ignoreInterpolatedArg("", "ignore-interpolated", "Ignore interpolated mass for extrapolation", cmd, false);

    cmd.parse(argc, argv);

    int btag = btagArg.getValue();
    std::string jec = jecArg.getValue();

    std::stringstream stream;
    stream << btag;

    std::string btagStr = stream.str();

    std::map<int, Efficiencies> efficiencies;
    for (int i = 750; i <= 1500; i += 50) {
      efficiencies[i] = Efficiencies(i);
    }

    const int M[4] = {750, 1000, 1250, 1500};

    // HLT efficiencies
    // See https://docs.google.com/spreadsheet/ccc?key=0AsI4zLOlSqcUdHhiYmFKbDIxY3YwWlJHdE9NSVhMSnc
    // for details about effiencies

    TGraphErrors trig_e;
    TGraphErrors trig_e_low;
    TGraphErrors trig_e_high;

    TGraphErrors trig_mu;
    TGraphErrors trig_mu_low;
    TGraphErrors trig_mu_high;

    int index = 0;

    for (auto& i: efficiencies) {

      if (i.second.isInterpolated)
        continue;

      switch (i.first) {
        case 750:

          i.second.effTrig_mu = 0.915380675350528;
          i.second.effTrig_e  = 0.96626144360191;
          i.second.error_effTrig_mu = 0.002953712440538;
          i.second.error_effTrig_e  = 0.001977197511629;
          break;
        case 1000:

          i.second.effTrig_mu = 0.919285739933724;
          i.second.effTrig_e  = 0.96382348703396;
          i.second.error_effTrig_mu = 0.002772009485335;
          i.second.error_effTrig_e  = 0.001970065457726;
          break;
        case 1250:

          i.second.effTrig_mu = 0.911219301358804;
          i.second.effTrig_e  = 0.95970579322911;
          i.second.error_effTrig_mu = 0.002954293566453;
          i.second.error_effTrig_e  = 0.002040906650861;
          break;
        case 1500:

          i.second.effTrig_mu = 0.919850794664936;
          i.second.effTrig_e  = 0.95589371607159;
          i.second.error_effTrig_mu = 0.002964341074836;
          i.second.error_effTrig_e  = 0.00222098605711;
          break;
      }

      if (! ignoreInterpolatedArg.getValue()) {
        trig_mu.SetPoint(index, i.first, i.second.effTrig_mu);
        trig_mu.SetPointError(index, 0., i.second.error_effTrig_mu);

        trig_mu_low.SetPoint(index, i.first, i.second.effTrig_mu - i.second.error_effTrig_mu);
        trig_mu_low.SetPointError(index, 0., i.second.error_effTrig_mu); // Needed for fit
        trig_mu_high.SetPoint(index, i.first, i.second.effTrig_mu + i.second.error_effTrig_mu);
        trig_mu_high.SetPointError(index, 0., i.second.error_effTrig_mu); // Needed for fit

        trig_e.SetPoint(index, i.first, i.second.effTrig_e);
        trig_e.SetPointError(index, 0., i.second.error_effTrig_e);

        trig_e_low.SetPoint(index, i.first, i.second.effTrig_e - i.second.error_effTrig_e);
        trig_e_low.SetPointError(index, 0., i.second.error_effTrig_e); // Needed for fit
        trig_e_high.SetPoint(index, i.first, i.second.effTrig_e + i.second.error_effTrig_e);
        trig_e_high.SetPointError(index++, 0., i.second.error_effTrig_e); // Needed for fit
      }
    }

    TF1 triggerEff_fit_mu("sel_eff_fit_mu", "[0] + [1] * x + [2] * x * x + [3] * x * x * x", 750, 1500);
    TF1* triggerEff_fit_mu_low = (TF1*) triggerEff_fit_mu.Clone("triggerEff_fit_mu_low");
    TF1* triggerEff_fit_mu_high = (TF1*) triggerEff_fit_mu.Clone("triggerEff_fit_mu_high");

    TF1 triggerEff_fit_e("sel_eff_fit_e", "[0] + [1] * x + [2] * x * x + [3] * x * x * x", 750, 1500);
    TF1* triggerEff_fit_e_low = (TF1*) triggerEff_fit_e.Clone("triggerEff_fit_e_low");
    TF1* triggerEff_fit_e_high = (TF1*) triggerEff_fit_e.Clone("triggerEff_fit_e_high");

    if (! ignoreInterpolatedArg.getValue()) {
      trig_mu.Fit(&triggerEff_fit_mu, "QMR");
      trig_mu_low.Fit(triggerEff_fit_mu_low, "QMR");
      trig_mu_high.Fit(triggerEff_fit_mu_high, "QMR");

      trig_e.Fit(&triggerEff_fit_e, "QR");
      trig_e_low.Fit(triggerEff_fit_e_low, "QR");
      trig_e_high.Fit(triggerEff_fit_e_high, "QR");
    }

    TCanvas c("c", "c", 800, 800);
    
    /*
    {
      TMultiGraph *mg = new TMultiGraph();

      mg->Add(&trig_e, "lp");
      mg->Add(&trig_e_low, "lp");
      mg->Add(&trig_e_high, "lp");

      mg->Draw("a");

      c.Print("trig_eff.root");

      delete mg;
    }
    */

    //--- selection efficiencies
    const float N0[4] = {108827, 102411, 96994, 96194};

    float Nsel_mu[4];
    float ErrNsel_mu[4];
    float Nsel_e[4];
    float ErrNsel_e[4];
    loadSelection(jec, btag, M, Nsel_mu, ErrNsel_mu, Nsel_e, ErrNsel_e);

    TGraphErrors e_mu;
    TGraphErrors e_mu_low;
    TGraphErrors e_mu_high;

    TGraphErrors e_e;
    TGraphErrors e_e_low;
    TGraphErrors e_e_high;

    TF1 selectionEff_fit_mu("sel_eff_fit_mu", "[0] + [1] * x + [2] * x * x", 750, 1500);
    TF1* selectionEff_fit_mu_low = (TF1*) selectionEff_fit_mu.Clone("sel_eff_fit_mu_low");
    TF1* selectionEff_fit_mu_high = (TF1*) selectionEff_fit_mu.Clone("sel_eff_fit_mu_high");

    TF1 selectionEff_fit_e("sel_eff_fit_e", "[0] + [1] * x + [2] * x * x", 750, 1500);
    TF1* selectionEff_fit_e_low = (TF1*) selectionEff_fit_e.Clone("sel_eff_fit_e_low");
    TF1* selectionEff_fit_e_high = (TF1*) selectionEff_fit_e.Clone("sel_eff_fit_e_high");

    index = 0;
    Efficiencies* lowMass_eff = nullptr;
    for (auto& i: efficiencies) {
      Efficiencies& eff = i.second;

      if (! i.second.isInterpolated) {
        lowMass_eff = &eff;

        eff.selectionEff_mu = Nsel_mu[index] / N0[index];
        eff.selectionEff_e  = Nsel_e[index] / N0[index];
        eff.error_selectionEff_mu = ErrNsel_mu[index] / N0[index];
        eff.error_selectionEff_e = ErrNsel_e[index] / N0[index];

        if (! ignoreInterpolatedArg.getValue()) {
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
      } else if (ignoreInterpolatedArg.getValue()) {
        eff.copy(*lowMass_eff);
      }
    }

    if (! ignoreInterpolatedArg.getValue()) {
      e_mu.Fit(&selectionEff_fit_mu, "QR");
      e_mu_low.Fit(selectionEff_fit_mu_low, "QR");
      e_mu_high.Fit(selectionEff_fit_mu_high, "QR");

      e_e.Fit(&selectionEff_fit_e, "QR");
      e_e_low.Fit(selectionEff_fit_e_low, "QR");
      e_e_high.Fit(selectionEff_fit_e_high, "QR");
    }

    /*
    {
      TMultiGraph *mg = new TMultiGraph();

      mg->Add(&e_mu, "lp");
      mg->Add(&e_mu_low, "lp");
      mg->Add(&e_mu_high, "lp");

      mg->Draw("a");

      c.Print("sel_eff.root");

      delete mg;
    }
    */


    if (! ignoreInterpolatedArg.getValue()) {
      for (auto& i: efficiencies) {
        Efficiencies& eff = i.second;

        if (i.second.isInterpolated) {
          eff.selectionEff_mu = selectionEff_fit_mu.Eval(i.first);
          eff.selectionEff_e  = selectionEff_fit_e.Eval(i.first);

          eff.error_selectionEff_mu = fabs(selectionEff_fit_mu_high->Eval(i.first) - selectionEff_fit_mu_low->Eval(i.first)) / 2.;
          eff.error_selectionEff_e = fabs(selectionEff_fit_e_high->Eval(i.first) - selectionEff_fit_e_low->Eval(i.first)) / 2.; 

          eff.effTrig_mu = triggerEff_fit_mu.Eval(i.first);
          eff.effTrig_e = triggerEff_fit_e.Eval(i.first);

          eff.error_effTrig_mu = fabs(triggerEff_fit_mu_high->Eval(i.first) - triggerEff_fit_mu_low->Eval(i.first)) / 2.;
          eff.error_effTrig_e = fabs(triggerEff_fit_e_high->Eval(i.first) - triggerEff_fit_e_low->Eval(i.first)) / 2.;
        }
      }
    }

    Json::Value root;
    getJsonRoot("efficiencies.json", root, false);

    for (auto& i: efficiencies) {

      if (ignoreInterpolatedArg.getValue() && i.second.isInterpolated)
        continue;

      std::stringstream ss;
      ss << i.first;
      std::string mass = ss.str();

      root[getAnalysisUUID()][mass][btagStr][jec] = i.second.getAsJSON();
      std::cout << i.second << std::endl;
      std::cout << std::endl;
    }

    Json::StyledWriter writer;
    std::ofstream output("efficiencies.json");
    output << writer.write(root);
    output.close();
    std::cout << "Efficiencies saved as 'efficiences.json'" << std::endl;

    if (jec == "nominal") {
      TString noteFilename = TString::Format("efficiencies_table_%s_%d_btag.tex", getAnalysisName().c_str(), btagArg.getValue());

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
