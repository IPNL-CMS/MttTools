#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <json/json.h>
#include <tclap/CmdLine.h>

#include <TString.h>

#include "Utils.h"

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

    Json::Value massNode = root[getAnalysisUUID()][strMass][btagStr];

    if (! massNode.isMember(jecType)) {
      std::cerr << "ERROR: '" << jecType << "' not found for m=" << mass << " in JSON file. Setting to 0." << std::endl;
      
      nSelectionMu[i] = 0;
      errNSelectionMu[i] = 0;
      nSelectionE[i] = 0;
      errNSelectionE[i] = 0;

      continue;
    }

    Json::Value jecNode = massNode[jecType];

    nSelectionMu[i] = jecNode["mu"]["events"].asDouble();
    errNSelectionMu[i] = jecNode["mu"]["error"].asDouble();

    nSelectionE[i] = jecNode["e"]["events"].asDouble();
    errNSelectionE[i] = jecNode["e"]["error"].asDouble();
  }
}

int main(int argc, char** argv) {

  try {
    TCLAP::CmdLine cmd("Compute efficiencies", ' ', "0.1");

    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "int", cmd);

    cmd.parse(argc, argv);

    int btag = btagArg.getValue();

    std::stringstream stream;
    stream << btag;

    std::string btagStr = stream.str();

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

    //--- selection efficiencies
    const float N0[4] = {108827, 102411, 96994, 96194};

    float Nsel_mu_nominal[4];
    float ErrNsel_mu_nominal[4];
    float Nsel_e_nominal[4];
    float ErrNsel_e_nominal[4];
    loadSelection("nominal", btag, M, Nsel_mu_nominal, ErrNsel_mu_nominal, Nsel_e_nominal, ErrNsel_e_nominal);

    float Nsel_mu_JECup[4];
    float ErrNsel_mu_JECup[4];
    float Nsel_e_JECup[4];
    float ErrNsel_e_JECup[4];
    loadSelection("JECup", btag, M, Nsel_mu_JECup, ErrNsel_mu_JECup, Nsel_e_JECup, ErrNsel_e_JECup);

    float Nsel_mu_JECdown[4];
    float ErrNsel_mu_JECdown[4];
    float Nsel_e_JECdown[4];
    float ErrNsel_e_JECdown[4];
    loadSelection("JECdown", btag, M, Nsel_mu_JECdown, ErrNsel_mu_JECdown, Nsel_e_JECdown, ErrNsel_e_JECdown);

    float eff_mu[4];
    float eff_e[4];
    float s_eff_mu[4];
    float s_eff_e[4];
    float eff_mu_JECup[4];
    float eff_e_JECup[4];
    float s_eff_mu_JECup[4];
    float s_eff_e_JECup[4];
    float eff_mu_JECdown[4];
    float eff_e_JECdown[4];
    float s_eff_mu_JECdown[4];
    float s_eff_e_JECdown[4];
    std::cout << "---- Selection efficiencies (w/o HLT)----" << std::endl;
    for (int i=0 ; i<4 ; i++)
    {
      eff_mu[i] = Nsel_mu_nominal[i] / N0[i];
      s_eff_mu[i] = ErrNsel_mu_nominal[i] / N0[i];
      eff_e[i] = Nsel_e_nominal[i] / N0[i];
      s_eff_e[i] = ErrNsel_e_nominal[i] / N0[i];
      std::cout << "nominal M_Z' = " << M[i]
        << '\t' << "eff_muon = " << eff_mu[i] << " +/- " << s_eff_mu[i] 
        << '\t' << "eff_electron = " << eff_e[i] << " +/- " << s_eff_e[i] << std::endl;
      eff_mu_JECup[i] = Nsel_mu_JECup[i]/N0[i];
      s_eff_mu_JECup[i] = ErrNsel_mu_JECup[i]/N0[i];
      eff_e_JECup[i] = Nsel_e_JECup[i]/N0[i];
      s_eff_e_JECup[i] = ErrNsel_e_JECup[i]/N0[i];
      std::cout << "JECup   M_Z' = " << M[i]
        << '\t' << "eff_muon = " << eff_mu_JECup[i] << " +/- " << s_eff_mu_JECup[i] 
        << '\t' << "eff_electron = " << eff_e_JECup[i] << " +/- " << s_eff_e_JECup[i] << std::endl;
      eff_mu_JECdown[i] = Nsel_mu_JECdown[i]/N0[i];
      s_eff_mu_JECdown[i] = ErrNsel_mu_JECdown[i]/N0[i];
      eff_e_JECdown[i] = Nsel_e_JECdown[i]/N0[i];
      s_eff_e_JECdown[i] = ErrNsel_e_JECdown[i]/N0[i];
      std::cout << "JECdown M_Z' = " << M[i]
        << '\t' << "eff_muon = " << eff_mu_JECdown[i] << " +/- " << s_eff_mu_JECdown[i] 
        << '\t' << "eff_electron = " << eff_e_JECdown[i] << " +/- " << s_eff_e_JECdown[i] << std::endl;
    }
    std::cout << std::endl;

    // nominal

    Json::Value root;
    getJsonRoot("efficiencies.json", root, false);
    for (int i=0 ; i<4 ; i++) {

      std::stringstream ss;
      ss << M[i];

      std::string mass = ss.str();

      //root[mass]["nominal"] = 
      Json::Value array(Json::arrayValue);
      array.append(eff_mu[i]);
      array.append(eff_e[i]);
      array.append(eff_trg_mu[i]);
      array.append(eff_trg_e[i]);
      array.append(s_eff_mu[i] / eff_mu[i]);
      array.append(s_eff_e[i] / eff_e[i]);
      array.append(0.009 / eff_trg_mu[i]);
      array.append(0.004 / eff_trg_e[i]);

      root[getAnalysisUUID()][mass][btagStr]["nominal"] = array;
      array.clear();

      array.append(eff_mu_JECup[i]);
      array.append(eff_e_JECup[i]);
      array.append(eff_trg_mu[i]);
      array.append(eff_trg_e[i]);
      array.append(eff_mu_JECup[i] == 0 ? 0 : s_eff_mu_JECup[i] / eff_mu_JECup[i]);
      array.append(eff_e_JECup[i] == 0 ? 0 : s_eff_e_JECup[i] / eff_e_JECup[i]);
      array.append(0.009 / eff_trg_mu[i]);
      array.append(0.004 / eff_trg_e[i]);

      root[getAnalysisUUID()][mass][btagStr]["JECup"] = array;
      array.clear();

      array.append(eff_mu_JECdown[i]);
      array.append(eff_e_JECdown[i]);
      array.append(eff_trg_mu[i]);
      array.append(eff_trg_e[i]);
      array.append(eff_mu_JECdown[i] == 0 ? 0 : s_eff_mu_JECdown[i] / eff_mu_JECdown[i]);
      array.append(eff_e_JECdown[i] == 0 ? 0 : s_eff_e_JECdown[i] / eff_e_JECdown[i]);
      array.append(0.009 / eff_trg_mu[i]);
      array.append(0.004 / eff_trg_e[i]);

      root[getAnalysisUUID()][mass][btagStr]["JECdown"] = array;
    }

    Json::StyledWriter writer;
    std::ofstream output("efficiencies.json");
    output << writer.write(root);
    output.close();
    std::cout << "Efficiencies saved as 'efficiences.json'" << std::endl;

    TString noteFilename = TString::Format("efficiencies_table_%s_%d_btag.tex", getAnalysisName().c_str(), btagArg.getValue());

    // table latex pour la note :
    std::ofstream latex(noteFilename);
    latex << "\\mtt & 750 GeV & 1000 GeV & 1250 GeV & 1500 GeV \\\\ " << std::endl;
    latex << "\\hline" << std::endl;
    latex << std::setiosflags(std::ios::fixed) << std::setprecision(2) ;
    latex << "$\\epsilon(Z^{\\prime}), semi-mu$ (\\%)         &  " 
      << eff_mu[1]*100 << "$\\pm$" << s_eff_mu[0]*100 << " & "
      << eff_mu[2]*100 << "$\\pm$" << s_eff_mu[1]*100 << " & "
      << eff_mu[3]*100 << "$\\pm$" << s_eff_mu[2]*100 << " & "
      << eff_mu[4]*100 << "$\\pm$" << s_eff_mu[3]*100 << " \\\\" << std::endl;
    latex << "$\\epsilon(Z^{\\prime}), semi-e$ (\\%)          &  " 
      << eff_e[1]*100 << "$\\pm$" << s_eff_e[0]*100 << " & "
      << eff_e[2]*100 << "$\\pm$" << s_eff_e[1]*100 << " & "
      << eff_e[3]*100 << "$\\pm$" << s_eff_e[2]*100 << " & "
      << eff_e[4]*100 << "$\\pm$" << s_eff_e[3]*100 << " \\\\" << std::endl;
    latex.close();

    std::cout << "Latex table saved as '" << noteFilename << "'" << std::endl;

  } catch (TCLAP::ArgException& e) {

  }

}
