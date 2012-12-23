#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <json/json.h>
#include <tclap/CmdLine.h>

#include <TString.h>

#include "Utils.h"

#include <TGraphErrors.h>
#include <TF1.h>

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
      array.append(0.009 / effTrig_mu);
      array.append(0.004 / effTrig_e);

      return array;
    }

    friend std::ostream& operator<<(std::ostream& stream, const Efficiencies& eff);
};

std::ostream& operator<<(std::ostream& stream, const Efficiencies& eff) {
  stream << "M_Z' = " << eff.mass
    << '\t' << "eff_muon = " << eff.selectionEff_mu << " +/- " << eff.error_selectionEff_mu 
    << '\t' << "eff_electron = " << eff.selectionEff_e << " +/- " << eff.error_selectionEff_e;

  return stream;
}

int main(int argc, char** argv) {

  try {
    TCLAP::CmdLine cmd("Compute efficiencies", ' ', "0.1");

    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);
    TCLAP::ValueArg<std::string> jecArg("", "jec", "Type of JEC", false, "nominal", "nominal/JECup/JECdown", cmd);

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


    //--- first compute trigger efficiencies
    // lumi : cf https://lyosvn.in2p3.fr/cms_top/wiki/Fall11_428
    const float lumimu_A1 = 211.599;
    const float lumimu_A2 = 929.748;
    const float lumimu_A3 = 368.037;
    const float lumimu_A4 = 412.359;
    const float lumimu_A5 = 246.527;
    const float lumimu_B6 = 1698. ;
    const float lumimu_B7 = 812.470;
    const float lumimu_tot = lumimu_A1+lumimu_A2+lumimu_A3+lumimu_A4+lumimu_A5+lumimu_B6+lumimu_B7;

    const float lumie_A1 = 216.240;
    const float lumie_A2 = 139.078;
    const float lumie_A3 = 790.670;
    const float lumie_A4 = 368.017;
    const float lumie_A5 = 658.886;
    const float lumie_B6 = 1697. ;
    const float lumie_B7 = 812.47 ;
    const float lumie_tot = lumie_A1+lumie_A2+lumie_A3+lumie_A4+lumie_A5+lumie_B6+lumie_B7;

    // trg efficiencies : cfg mail nicolas
    const float eff_IsoMu17[4] = {88.3, 84.2, 82.4, 79.3};
    const float eff_IsoMu17_DiCentralJet30[4] = {87.3, 83.5, 80.1 , 81.8}; 
    const float eff_IsoMu17_TriCentralJet30[4] = {84.0, 81.8, 78.5, 78.8};
    const float eff_IsoMu17_TriCentralPFJet30[4] = {84.1, 81.5, 78.2, 77.2};

    const float eff_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT[4] = {97.6, 97.5, 97.4, 97.7};
    const float eff_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT[4] = {93.7, 94.6, 95.0, 94.5};
    const float eff_Ele25_CaloIdVT_TrkIdT_TriCentralJet30[4] = {97.0, 96.6, 96.7, 96.8};
    const float eff_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30[4] = {96.6, 97.1, 96.8, 97.}; 
    const float eff_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30[4] = {96.7, 97.0, 96.5, 96.7};

    // compute integrated trigger efficiencies
    float eff_trg_mu[4];
    float eff_trg_e[4];
    std::cout << "---- HLT efficiencies ----" << std::endl;
    for (int i=0 ; i < 4 ; i++)
    {
      //eff_trg_mu[i]= ( (lumimu_A1 + lumimu_A2) * eff_IsoMu17[i] 
          //+ (lumimu_A3 + lumimu_A4) * eff_IsoMu17_DiCentralJet30[i] 
          //+ (lumimu_A5 + lumimu_B6) * eff_IsoMu17_TriCentralJet30[i] 
          //+ lumimu_B7 * eff_IsoMu17_TriCentralPFJet30[i] 
          //) / lumimu_tot / 100;
      eff_trg_mu[i] = 1.; //FIXME
      //eff_trg_e[i]= ( lumie_A1 * eff_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT[i] 
          //+ lumie_A2 * eff_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT[i] 
          //+ lumie_A3 * eff_Ele25_CaloIdVT_TrkIdT_TriCentralJet30[i] 
          //+ (lumie_A4 + lumie_A5 + lumie_B6) * eff_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30[i] 
          //+ lumie_B7 * eff_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30[i] 
          //) /lumie_tot / 100;
      eff_trg_e[i] = 1.; //FIXME
      std::cout << "M_Z' = " << M[i] 
        << '\t' << "eff_muon = " << eff_trg_mu[i]  
        << '\t' << "eff_electron = " << eff_trg_e[i] << std::endl;
    }
    std::cout << std::endl;

    for (auto& i: efficiencies) {
      i.second.effTrig_mu = 1.;
      i.second.effTrig_e  = 1.;
      i.second.error_effTrig_mu = 0;
      i.second.error_effTrig_e  = 0;
    }

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
    TF1* selectionEff_fit_mu_low = (TF1*) selectionEff_fit_mu.Clone();
    TF1* selectionEff_fit_mu_high = (TF1*) selectionEff_fit_mu.Clone();

    TF1 selectionEff_fit_e("sel_eff_fit_e", "[0] + [1] * x + [2] * x * x", 750, 1500);
    TF1* selectionEff_fit_e_low = (TF1*) selectionEff_fit_e.Clone();
    TF1* selectionEff_fit_e_high = (TF1*) selectionEff_fit_e.Clone();

    int index = 0;
    Efficiencies* lowMass_eff = nullptr;
    for (auto& i: efficiencies) {
      Efficiencies& eff = i.second;

      if (! i.second.isInterpolated) {
        lowMass_eff = &eff;

        eff.selectionEff_mu = Nsel_mu[index] / N0[index];
        eff.selectionEff_e  = Nsel_e[index] / N0[index];
        eff.error_selectionEff_mu = ErrNsel_mu[index] / N0[index];
        eff.error_selectionEff_e = ErrNsel_e[index] / N0[index];

        e_mu.SetPoint(index, i.first, i.second.selectionEff_mu);
        e_mu.SetPointError(index, 0, i.second.error_selectionEff_mu);

        e_mu_low.SetPoint(index, i.first, i.second.selectionEff_mu - i.second.error_selectionEff_mu);
        e_mu_high.SetPoint(index, i.first, i.second.selectionEff_mu + i.second.error_selectionEff_mu);

        e_e.SetPoint(index, i.first, i.second.selectionEff_e);
        e_e.SetPointError(index, 0, i.second.error_selectionEff_e);

        e_e_low.SetPoint(index, i.first, i.second.selectionEff_e - i.second.error_selectionEff_e);
        e_e_high.SetPoint(index, i.first, i.second.selectionEff_e + i.second.error_selectionEff_e);

        index++;
      } else {
        eff.copy(*lowMass_eff);
      }
    }

    e_mu.Fit(&selectionEff_fit_mu, "QR");
    e_mu_low.Fit(selectionEff_fit_mu_low, "QR");
    e_mu_high.Fit(selectionEff_fit_mu_high, "QR");

    e_e.Fit(&selectionEff_fit_e, "QR");
    e_e_low.Fit(selectionEff_fit_e_low, "QR");
    e_e_high.Fit(selectionEff_fit_e_high, "QR");

    for (auto& i: efficiencies) {
      Efficiencies& eff = i.second;

      if (i.second.isInterpolated) {
        eff.selectionEff_mu = selectionEff_fit_mu.Eval(i.first);
        eff.selectionEff_e  = selectionEff_fit_e.Eval(i.first);

        eff.error_selectionEff_mu = fabs(selectionEff_fit_mu_high->Eval(i.first) - selectionEff_fit_mu_low->Eval(i.first)) / 2.;
        eff.error_selectionEff_e = fabs(selectionEff_fit_e_high->Eval(i.first) - selectionEff_fit_e_low->Eval(i.first)) / 2.; 
      }
    }

    Json::Value root;
    getJsonRoot("efficiencies.json", root, false);

    for (auto& i: efficiencies) {
      std::stringstream ss;
      ss << i.first;
      std::string mass = ss.str();

      root[getAnalysisUUID()][mass][btagStr][jec] = i.second.getAsJSON();
      std::cout << i.second << std::endl;
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
        latex << " & " << i.first << " GeV";
      }

      latex << "\\\\" << std::endl << "\\hline" << std::endl;
      latex << std::setiosflags(std::ios::fixed) << std::setprecision(2);

      latex << "$\\epsilon(Z^{\\prime}), semi-mu$ (\\%)";
      for (auto& i: efficiencies) {
        latex << " & " << i.second.selectionEff_mu * 100 << " $\\pm$ " << std::setprecision(4) << i.second.error_selectionEff_mu * 100 << std::setprecision(2);
      }
      latex << "\\\\" << std::endl;

      latex << "$\\epsilon(Z^{\\prime}), semi-e$ (\\%)"; 
      for (auto& i: efficiencies) {
        latex << " & " << i.second.selectionEff_e * 100 << " $\\pm$ " << std::setprecision(4) << i.second.error_selectionEff_e * 100 << std::setprecision(2);
      }
      latex << "\\\\" << std::endl;

      latex.close();

      std::cout << "Latex table saved as '" << noteFilename << "'" << std::endl;
    }

  } catch (TCLAP::ArgException& e) {

  }

}
