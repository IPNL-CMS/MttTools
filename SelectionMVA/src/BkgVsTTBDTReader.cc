#include <BkgVsTTBDTReader.h>

#include <Math/Vector3Dfwd.h>
#include <TMath.h>

namespace MVAReaderUtils {
  TLorentzVector* getP4(TClonesArray* array, int index) {
    return (TLorentzVector*) array->At(index);
  }
};

BkgVsTTBDTReader::BkgVsTTBDTReader(const std::vector<std::string>& inputFiles, bool isMC):
  MVAReader(inputFiles, isMC) {

    m_btagCalculator = std::make_shared<NBTagCalculator>(inputFiles, NOMINAL);
  }


void BkgVsTTBDTReader::initTrees() {
  m_mtt = loadChain("Mtt");
  m_jets = loadChain("jet_PF");

  setBranchAddress(m_mtt, "MET", MET);
  setBranchAddress(m_mtt, "aplanarity", aplanarity);
  setBranchAddress(m_mtt, "circularity", circularity);
  setBranchAddress(m_mtt, "sphericity", sphericity);
  setBranchAddress(m_mtt, "nBtaggedJets_CSVM", n_btagged_jets);

  setBranchAddress(m_mtt, "selectedLeptonP4", selectedLeptonP4_AfterChi2);
  setBranchAddress(m_mtt, "selectedFirstJetP4_AfterChi2", selectedFirstJetP4_AfterReco);
  setBranchAddress(m_mtt, "selectedSecondJetP4_AfterChi2", selectedSecondJetP4_AfterReco);
  setBranchAddress(m_mtt, "selectedNeutrinoP4_AfterChi2", selectedNeutrinoP4_AfterReco);
  setBranchAddress(m_mtt, "selectedHadronicBP4_AfterChi2", selectedHadronicBP4_AfterReco);
  setBranchAddress(m_mtt, "selectedLeptonicBP4_AfterChi2", selectedLeptonicBP4_AfterReco);  

  setBranchAddress(m_jets, "n_jets", n_jets);
  setBranchAddress(m_jets, "jet_4vector", jets_p4);
  setBranchAddress(m_jets, "jet_btag_CSV", jets_CSV_discriminant);
}

void BkgVsTTBDTReader::setupVariables() {

  m_reader->AddVariable("aplanarity", &aplanarity);
  m_reader->AddVariable("circularity", &circularity);
  m_reader->AddVariable("sphericity", &sphericity);
  //m_reader->AddVariable("mean_csv", &mean_csv);
  m_reader->AddVariable("n_btagged_jets", &n_btagged_jets_float);
  m_reader->AddVariable("st", &st);
  m_reader->AddVariable("theta_lepton", &theta_lepton);
  m_reader->AddVariable("met", &MET);
  //m_reader->AddVariable("neutrino_pt", &neutrino_pt);
  m_reader->AddVariable("lepton_pt", &lepton_pt);
  m_reader->AddVariable("lepton_eta", &lepton_eta);
  m_reader->AddVariable("leptonic_B_eta", &leptonic_B_eta);
  m_reader->AddVariable("hadronic_B_eta", &hadronic_B_eta);
  m_reader->AddVariable("hadronic_first_jet_pt + hadronic_second_jet_pt", &hadronic_first_jet_pt_plus_hadronic_second_jet_pt);
  m_reader->AddVariable("hadronic_first_jet_eta + hadronic_second_jet_eta", &hadronic_first_jet_eta_plus_hadronic_second_jet_eta);
  m_reader->AddVariable("leptonic_T_pt", &leptonic_T_pt);
  m_reader->AddVariable("leptonic_T_transverse_mass", &leptonic_T_transverse_mass);
  m_reader->AddVariable("hadronic_T_pt", &hadronic_T_pt);
  //m_reader->AddVariable("resonance_pt", &resonance_pt);
  //m_reader->AddVariable("resonance_eta", &resonance_eta);
  m_reader->AddVariable("neutrino_leptonic_B_delta_R", &neutrino_leptonic_B_delta_R);
  m_reader->AddVariable("neutrino_lepton_delta_R", &neutrino_lepton_delta_R);
  m_reader->AddVariable("neutrino_hadronic_B_delta_R", &neutrino_hadronic_B_delta_R);
  m_reader->AddVariable("neutrino_hadronic_first_jet_delta_R+neutrino_hadronic_second_jet_delta_R", &neutrino_hadronic_first_jet_delta_R_plus_neutrino_hadronic_second_jet_delta_R);
  //m_reader->AddVariable("neutrino_hadronic_second_jet_delta_R", &neutrino_hadronic_second_jet_delta_R);
  m_reader->AddVariable("lepton_leptonic_B_delta_R", &lepton_leptonic_B_delta_R);
  m_reader->AddVariable("lepton_hadronic_B_delta_R", &lepton_hadronic_B_delta_R);
  m_reader->AddVariable("lepton_hadronic_first_jet_delta_R+lepton_hadronic_second_jet_delta_R", &lepton_hadronic_first_jet_delta_R_plus_lepton_hadronic_second_jet_delta_R);
  //m_reader->AddVariable("lepton_hadronic_second_jet_delta_R", &lepton_hadronic_second_jet_delta_R);
  m_reader->AddVariable("leptonic_B_hadronic_B_delta_R", &leptonic_B_hadronic_B_delta_R);
  m_reader->AddVariable("leptonic_B_hadronic_first_jet_delta_R+leptonic_B_hadronic_second_jet_delta_R", &leptonic_B_hadronic_first_jet_delta_R_plus_leptonic_B_hadronic_second_jet_delta_R);
  //m_reader->AddVariable("leptonic_B_hadronic_second_jet_delta_R", &leptonic_B_hadronic_second_jet_delta_R);
  m_reader->AddVariable("hadronic_B_hadronic_first_jet_delta_R+hadronic_B_hadronic_second_jet_delta_R", &hadronic_B_hadronic_first_jet_delta_R_plus_hadronic_B_hadronic_second_jet_delta_R);
  //m_reader->AddVariable("hadronic_B_hadronic_second_jet_delta_R", &hadronic_B_hadronic_second_jet_delta_R);
  m_reader->AddVariable("hadronic_first_jet_hadronic_second_jet_delta_R", &hadronic_first_jet_hadronic_second_jet_delta_R);
  //m_reader->AddVariable("leptonic_W_hadronic_W_delta_R", &leptonic_W_hadronic_W_delta_R);
  //m_reader->AddVariable("leptonic_T_hadronic_T_delta_R", &leptonic_T_hadronic_T_delta_R);
  //m_reader->AddVariable("cos_theta_leading_top_resonance", &cos_theta_leading_top_resonance);
}

void BkgVsTTBDTReader::computeVariables(uint64_t entry) {

  //uint32_t n = 0;
  float HT = 0;
  //mean_csv = 0;
  for (uint32_t i = 0; i < (uint32_t) jets_p4->GetEntriesFast(); i++) {
    float pt = ((TLorentzVector*) (*jets_p4)[i])->Pt();
    //if (pt > 30 && jets_CSV_discriminant[i] > 0.244) {
      //n++;
      //mean_csv += jets_CSV_discriminant[i];
    //}

    HT += pt;
  }

  //if (n > 0)
    //mean_csv /= n;
 
  if (m_isMC)
    n_btagged_jets = m_btagCalculator->getNumberOfBTaggedJets(entry);

  n_btagged_jets_float = n_btagged_jets;

  LorentzVector selectedLeptonP4_AfterReco;
  TLorentzVector* leptonP4 = MVAReaderUtils::getP4(selectedLeptonP4_AfterChi2, 0);
  selectedLeptonP4_AfterReco = LorentzVector(leptonP4->Pt(), leptonP4->Eta(), leptonP4->Phi(), leptonP4->E());

  LorentzVector leptonic_W_p4 = selectedLeptonP4_AfterReco + *selectedNeutrinoP4_AfterReco;
  LorentzVector hadronic_W_p4 = *selectedFirstJetP4_AfterReco + *selectedHadronicBP4_AfterReco;
  LorentzVector leptonic_T_p4 = leptonic_W_p4 + *selectedLeptonicBP4_AfterReco;
  LorentzVector hadronic_T_p4 = hadronic_W_p4 + *selectedHadronicBP4_AfterReco;
  LorentzVector resonance_p4 = hadronic_T_p4 + leptonic_T_p4;

  st = HT + selectedLeptonP4_AfterReco.Pt();
  //neutrino_pt = selectedNeutrinoP4_AfterReco->Pt();
  lepton_pt = selectedLeptonP4_AfterReco.Pt();
  lepton_eta = selectedLeptonP4_AfterReco.Eta();
  leptonic_B_eta = selectedLeptonicBP4_AfterReco->Eta();
  hadronic_B_eta = selectedHadronicBP4_AfterReco->Eta();
  hadronic_first_jet_pt_plus_hadronic_second_jet_pt = selectedFirstJetP4_AfterReco->Pt() + selectedSecondJetP4_AfterReco->Pt();
  hadronic_first_jet_eta_plus_hadronic_second_jet_eta = selectedFirstJetP4_AfterReco->Eta() + selectedSecondJetP4_AfterReco->Eta();
  leptonic_T_pt = leptonic_T_p4.Pt();
  leptonic_T_transverse_mass = leptonic_T_p4.Mt();
  hadronic_T_pt = hadronic_T_p4.Pt();
  //resonance_pt = resonance_p4.Pt();
  //resonance_eta = resonance_p4.Eta();
  neutrino_leptonic_B_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedNeutrinoP4_AfterReco, *selectedLeptonicBP4_AfterReco);
  neutrino_lepton_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedNeutrinoP4_AfterReco, selectedLeptonP4_AfterReco);
  neutrino_hadronic_B_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedNeutrinoP4_AfterReco, *selectedHadronicBP4_AfterReco);
  neutrino_hadronic_first_jet_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedNeutrinoP4_AfterReco, *selectedFirstJetP4_AfterReco);
  neutrino_hadronic_second_jet_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedNeutrinoP4_AfterReco, *selectedSecondJetP4_AfterReco);
  neutrino_hadronic_first_jet_delta_R_plus_neutrino_hadronic_second_jet_delta_R = neutrino_hadronic_first_jet_delta_R + neutrino_hadronic_second_jet_delta_R;
  lepton_leptonic_B_delta_R = ROOT::Math::VectorUtil::DeltaR(selectedLeptonP4_AfterReco, *selectedLeptonicBP4_AfterReco);
  lepton_hadronic_B_delta_R = ROOT::Math::VectorUtil::DeltaR(selectedLeptonP4_AfterReco, *selectedHadronicBP4_AfterReco);
  lepton_hadronic_first_jet_delta_R = ROOT::Math::VectorUtil::DeltaR(selectedLeptonP4_AfterReco, *selectedFirstJetP4_AfterReco);
  lepton_hadronic_second_jet_delta_R = ROOT::Math::VectorUtil::DeltaR(selectedLeptonP4_AfterReco, *selectedSecondJetP4_AfterReco);
  lepton_hadronic_first_jet_delta_R_plus_lepton_hadronic_second_jet_delta_R = lepton_hadronic_first_jet_delta_R + lepton_hadronic_second_jet_delta_R;
  leptonic_B_hadronic_B_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedLeptonicBP4_AfterReco, *selectedHadronicBP4_AfterReco);
  leptonic_B_hadronic_first_jet_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedLeptonicBP4_AfterReco, *selectedFirstJetP4_AfterReco);
  leptonic_B_hadronic_second_jet_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedLeptonicBP4_AfterReco, *selectedSecondJetP4_AfterReco);
  leptonic_B_hadronic_first_jet_delta_R_plus_leptonic_B_hadronic_second_jet_delta_R = leptonic_B_hadronic_first_jet_delta_R + leptonic_B_hadronic_second_jet_delta_R;
  hadronic_B_hadronic_first_jet_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedHadronicBP4_AfterReco, *selectedFirstJetP4_AfterReco);
  hadronic_B_hadronic_second_jet_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedHadronicBP4_AfterReco, *selectedSecondJetP4_AfterReco);
  hadronic_B_hadronic_first_jet_delta_R_plus_hadronic_B_hadronic_second_jet_delta_R = hadronic_B_hadronic_first_jet_delta_R + hadronic_B_hadronic_second_jet_delta_R;
  hadronic_first_jet_hadronic_second_jet_delta_R = ROOT::Math::VectorUtil::DeltaR(*selectedFirstJetP4_AfterReco, *selectedSecondJetP4_AfterReco);
  //leptonic_W_hadronic_W_delta_R = ROOT::Math::VectorUtil::DeltaR(leptonic_W_p4, hadronic_W_p4);
  //leptonic_T_hadronic_T_delta_R = ROOT::Math::VectorUtil::DeltaR(leptonic_T_p4, hadronic_T_p4);

  //LorentzVector& leading_top = (leptonic_T_p4.Pt() > hadronic_T_p4.Pt())
    //? leptonic_T_p4
    //: hadronic_T_p4;
  //cos_theta_leading_top_resonance = ROOT::Math::VectorUtil::CosTheta(leading_top.Vect(), resonance_p4.Vect());

  ROOT::Math::Boost resonance_in_rest_frame_boost = ROOT::Math::Boost(resonance_p4.BoostToCM());
  LorentzVector leptonic_T_p4_in_resonance_rest_frame = resonance_in_rest_frame_boost(leptonic_T_p4);

  ROOT::Math::XYZVector leptonic_T_direction = leptonic_T_p4_in_resonance_rest_frame.Vect();

  // Find the boost to go in the top rest frame, and apply it to the lepton to have the lepton direction
  ROOT::Math::Boost leptonic_T_in_rest_frame_boost = ROOT::Math::Boost(leptonic_T_p4.BoostToCM());
  LorentzVector lepton_p4_in_resonance_rest_frame = leptonic_T_in_rest_frame_boost(selectedLeptonP4_AfterReco);

  ROOT::Math::XYZVector lepton_direction = lepton_p4_in_resonance_rest_frame.Vect();
  theta_lepton = ROOT::Math::VectorUtil::Angle(leptonic_T_direction, lepton_direction);
}

void BkgVsTTBDTReader::createPlots(TFile* f) {
  f->cd();

  h_aplanarity = new TH1F("aplanarity_bdt", "", 25, 0., 0.5);
  h_circularity = new TH1F("circularity_bdt", "", 25, 0., 1.);
  h_sphericity = new TH1F("sphericity_bdt", "", 25, 0., 1.);
  h_mean_csv = new TH1F("mean_csv_bdt", "", 25, 0., 1.);
  h_n_btagged_jets = new TH1F("n_btagged_jets_bdt", "", 6, -1, 5);
  h_st = new TH1F("st_bdt", "", 40, 0, 2000);
  h_theta_lepton = new TH1F("theta_lepton_bdt", "", 5, 0, M_PI);
  h_MET = new TH1F("MET_bdt", "", 20, 0., 400.);
  //h_neutrino_pt;
  h_lepton_pt = new TH1F("lepton_pt_bdt", "", 20, 0., 400.);
  h_lepton_eta = new TH1F("lepton_eta_bdt", "", 25, -2*M_PI, 2*M_PI);
  h_leptonic_B_eta = new TH1F("leptonic_B_eta_bdt", "", 25, -2*M_PI, 2*M_PI);
  h_hadronic_B_eta = new TH1F("hadronic_B_eta_bdt", "", 25, -2*M_PI, 2*M_PI);
  h_hadronic_first_jet_pt_plus_hadronic_second_jet_pt = new TH1F("hadronic_first_jet_pt_plus_hadronic_second_jet_pt_bdt", "", 15, 0., 600.);
  h_hadronic_first_jet_eta_plus_hadronic_second_jet_eta = new TH1F("hadronic_first_jet_eta_plus_hadronic_second_jet_eta_bdt", "", 25, -2*M_PI, 2*M_PI);
  h_leptonic_T_pt = new TH1F("leptonic_T_pt_bdt", "", 15, 0., 600.);
  h_leptonic_T_transverse_mass = new TH1F("leptonic_T_transverse_mass_bdt", "", 15, 0., 600.);
  h_hadronic_T_pt = new TH1F("hadronic_T_pt_bdt", "", 15, 0., 600.);
  //h_resonance_pt;
  //h_resonance_eta;
  h_neutrino_leptonic_B_delta_R = new TH1F("neutrino_leptonic_B_delta_R_bdt", "", 50, 0, 10);
  h_neutrino_lepton_delta_R = new TH1F("neutrino_lepton_delta_R _bdt", "", 50, 0, 10);
  h_neutrino_hadronic_B_delta_R = new TH1F("neutrino_hadronic_B_delta_R_bdt", "", 50, 0, 10);
  h_neutrino_hadronic_first_jet_delta_R_plus_neutrino_hadronic_second_jet_delta_R = new TH1F("neutrino_hadronic_first_jet_delta_R_plus_neutrino_hadronic_second_jet_delta_R_bdt", "", 100, 0, 20);
  h_neutrino_hadronic_first_jet_delta_R = new TH1F("neutrino_hadronic_first_jet_delta_R_bdt", "", 50, 0, 10);
  h_neutrino_hadronic_second_jet_delta_R = new TH1F("neutrino_hadronic_second_jet_delta_R_bdt", "", 50, 0, 10);
  h_lepton_leptonic_B_delta_R = new TH1F("lepton_leptonic_B_delta_R_bdt", "", 50, 0, 10);
  h_lepton_hadronic_B_delta_R = new TH1F("lepton_hadronic_B_delta_R_bdt", "", 50, 0, 10);
  h_lepton_hadronic_first_jet_delta_R_plus_lepton_hadronic_second_jet_delta_R = new TH1F("lepton_hadronic_first_jet_delta_R_plus_lepton_hadronic_second_jet_delta_R_bdt", "", 100, 0, 20);
  h_lepton_hadronic_first_jet_delta_R = new TH1F("lepton_hadronic_first_jet_delta_R_bdt", "", 50, 0, 10);
  h_lepton_hadronic_second_jet_delta_R = new TH1F("lepton_hadronic_second_jet_delta_R_bdt", "", 50, 0, 10);
  h_leptonic_B_hadronic_B_delta_R = new TH1F("leptonic_B_hadronic_B_delta_R_bdt", "", 50, 0, 10);
  h_leptonic_B_hadronic_first_jet_delta_R_plus_leptonic_B_hadronic_second_jet_delta_R = new TH1F("leptonic_B_hadronic_first_jet_delta_R_plus_leptonic_B_hadronic_second_jet_delta_R_bdt", "", 100, 0, 20);
  h_leptonic_B_hadronic_first_jet_delta_R = new TH1F("leptonic_B_hadronic_first_jet_delta_R_bdt", "", 50, 0, 10);
  h_leptonic_B_hadronic_second_jet_delta_R = new TH1F("leptonic_B_hadronic_second_jet_delta_R_bdt", "", 50, 0, 10);
  h_hadronic_B_hadronic_first_jet_delta_R_plus_hadronic_B_hadronic_second_jet_delta_R = new TH1F("hadronic_B_hadronic_first_jet_delta_R_plus_hadronic_B_hadronic_second_jet_delta_R_bdt", "", 100, 0, 20);
  h_hadronic_B_hadronic_first_jet_delta_R = new TH1F("hadronic_B_hadronic_first_jet_delta_R_bdt", "", 50, 0, 10);
  h_hadronic_B_hadronic_second_jet_delta_R = new TH1F("hadronic_B_hadronic_second_jet_delta_R_bdt", "", 50, 0, 10);
  h_hadronic_first_jet_hadronic_second_jet_delta_R = new TH1F("hadronic_first_jet_hadronic_second_jet_delta_R_bdt", "", 50, 0, 10);
  //h_leptonic_W_hadronic_W_delta_R;
  //h_leptonic_T_hadronic_T_delta_R;
  //h_cos_theta_leading_top_resonance;

  h_aplanarity->SetXTitle("Aplanarity");
  h_circularity->SetXTitle("Circularity");
  h_sphericity->SetXTitle("Sphericity");
  h_mean_csv->SetXTitle("Mean CSV");
  h_n_btagged_jets->SetXTitle("Number of b-tagged jets");
  h_st->SetXTitle("ST = HT + p_{T, Lepton} (GeV)");
  h_theta_lepton->SetXTitle("#theta (leptonic top, lepton)_{resonance rest frame}");
  h_MET->SetXTitle("MET (GeV)");
  //h_neutrino_pt->SetXTitle("");
  h_lepton_pt->SetXTitle("Lepton p_{T} (GeV)");
  h_lepton_eta->SetXTitle("Lepton #eta");
  h_leptonic_B_eta->SetXTitle("Leptonic B #eta");
  h_hadronic_B_eta->SetXTitle("Hadronic B #eta");
  h_hadronic_first_jet_pt_plus_hadronic_second_jet_pt->SetXTitle("p_{T}^{hadronic first jet} + p_{T}^{hadronic second jet} (GeV)");
  h_hadronic_first_jet_eta_plus_hadronic_second_jet_eta->SetXTitle("#eta^{hadronic first jet} + #eta^{hadronic second jet}");
  h_leptonic_T_pt->SetXTitle("Leptonic top p_{T} (GeV)");
  h_leptonic_T_transverse_mass->SetXTitle("Leptonic top m_{T} (GeV)");
  h_hadronic_T_pt->SetXTitle("Hadronic top p_{T} (GeV)");
  //h_resonance_pt->SetXTitle("");
  //h_resonance_eta->SetXTitle("");
  h_neutrino_leptonic_B_delta_R->SetXTitle("#Delta R (neutrino, leptonic B)");
  h_neutrino_lepton_delta_R->SetXTitle("#Delta R (neutrino, lepton)");
  h_neutrino_hadronic_B_delta_R->SetXTitle("#Delta R (neutrino, hadronic B)");
  h_neutrino_hadronic_first_jet_delta_R_plus_neutrino_hadronic_second_jet_delta_R->SetXTitle("#Delta R (#nu, hadronic 1^{st} jet) + #Delta R (#nu, hadronic 2^{nd} jet)");
  h_neutrino_hadronic_first_jet_delta_R->SetXTitle("#Delta R (#nu, hadronic 1^{st} jet)");
  h_neutrino_hadronic_second_jet_delta_R->SetXTitle("#Delta R (#nu, hadronic 2^{nd} jet)");
  h_lepton_leptonic_B_delta_R->SetXTitle("#Delta R (lepton, leptonic B)");
  h_lepton_hadronic_B_delta_R->SetXTitle("#Delta R (lepton, hadronic B)");
  h_lepton_hadronic_first_jet_delta_R_plus_lepton_hadronic_second_jet_delta_R->SetXTitle("#Delta R (lepton, hadronic 1^{st} jet) + #Delta R (lepton, hadronic 2^{nd} jet)");
  h_lepton_hadronic_first_jet_delta_R->SetXTitle("#Delta R (lepton, hadronic 1^{st} jet)");
  h_lepton_hadronic_second_jet_delta_R->SetXTitle("#Delta R (lepton, hadronic 2^{nd} jet)");
  h_leptonic_B_hadronic_B_delta_R->SetXTitle("#Delta R (leptonic B, hadronic B)");
  h_leptonic_B_hadronic_first_jet_delta_R_plus_leptonic_B_hadronic_second_jet_delta_R->SetXTitle("#Delta R (leptonic B, hadronic 1^{st} jet) + #Delta R (leptonic B, hadronic 2^{nd} jet)");
  h_leptonic_B_hadronic_first_jet_delta_R->SetXTitle("#Delta R (leptonic B, hadronic 1^{st} jet)");
  h_leptonic_B_hadronic_second_jet_delta_R->SetXTitle("#Delta R (leptonic B, hadronic 2^{nd} jet)");
  h_hadronic_B_hadronic_first_jet_delta_R_plus_hadronic_B_hadronic_second_jet_delta_R->SetXTitle("#Delta R (hadronic B, hadronic 1^{st} jet) + #Delta R (hadronic B, hadronic 2^{nd} jet)");
  h_hadronic_B_hadronic_first_jet_delta_R->SetXTitle("#Delta R (hadronic B, hadronic 1^{st} jet)");
  h_hadronic_B_hadronic_second_jet_delta_R->SetXTitle("#Delta R (hadronic B, hadronic 2^{nd} jet)");
  h_hadronic_first_jet_hadronic_second_jet_delta_R->SetXTitle("#Delta R (hadronic 1^{st} jet, hadronic 2^{nd} jet)");
  //h_leptonic_W_hadronic_W_delta_R->SetXTitle("");
  //h_leptonic_T_hadronic_T_delta_R->SetXTitle("");
  //h_cos_theta_leading_top_resonance->SetXTitle("");

}

void BkgVsTTBDTReader::fillPlots(float weight) {
  h_aplanarity->Fill(aplanarity, weight);
  h_circularity->Fill(circularity, weight);
  h_sphericity->Fill(sphericity, weight);
  h_mean_csv->Fill(mean_csv, weight);
  h_n_btagged_jets->Fill(n_btagged_jets, weight);
  h_st->Fill(st, weight);
  h_theta_lepton->Fill(theta_lepton, weight);
  h_MET->Fill(MET, weight);
  //h_neutrino_pt->Fill(, weight);
  h_lepton_pt->Fill(lepton_pt, weight);
  h_lepton_eta->Fill(lepton_eta, weight);
  h_leptonic_B_eta->Fill(leptonic_B_eta, weight);
  h_hadronic_B_eta->Fill(hadronic_B_eta, weight);
  h_hadronic_first_jet_pt_plus_hadronic_second_jet_pt->Fill(hadronic_first_jet_pt_plus_hadronic_second_jet_pt, weight);
  h_hadronic_first_jet_eta_plus_hadronic_second_jet_eta->Fill(hadronic_first_jet_eta_plus_hadronic_second_jet_eta, weight);
  h_leptonic_T_pt->Fill(leptonic_T_pt, weight);
  h_leptonic_T_transverse_mass->Fill(leptonic_T_transverse_mass, weight);
  h_hadronic_T_pt->Fill(hadronic_T_pt, weight);
  //h_resonance_pt->Fill(, weight);
  //h_resonance_eta->Fill(, weight);
  h_neutrino_leptonic_B_delta_R->Fill(neutrino_leptonic_B_delta_R, weight);
  h_neutrino_lepton_delta_R->Fill(neutrino_lepton_delta_R, weight);
  h_neutrino_hadronic_B_delta_R->Fill(neutrino_hadronic_B_delta_R, weight);
  h_neutrino_hadronic_first_jet_delta_R_plus_neutrino_hadronic_second_jet_delta_R->Fill(neutrino_hadronic_first_jet_delta_R_plus_neutrino_hadronic_second_jet_delta_R, weight);
  h_neutrino_hadronic_first_jet_delta_R->Fill(neutrino_hadronic_first_jet_delta_R, weight);
  h_neutrino_hadronic_second_jet_delta_R->Fill(neutrino_hadronic_second_jet_delta_R, weight);
  h_lepton_leptonic_B_delta_R->Fill(lepton_leptonic_B_delta_R, weight);
  h_lepton_hadronic_B_delta_R->Fill(lepton_hadronic_B_delta_R, weight);
  h_lepton_hadronic_first_jet_delta_R_plus_lepton_hadronic_second_jet_delta_R->Fill(lepton_hadronic_first_jet_delta_R_plus_lepton_hadronic_second_jet_delta_R, weight);
  h_lepton_hadronic_first_jet_delta_R->Fill(lepton_hadronic_first_jet_delta_R, weight);
  h_lepton_hadronic_second_jet_delta_R->Fill(lepton_hadronic_second_jet_delta_R, weight);
  h_leptonic_B_hadronic_B_delta_R->Fill(leptonic_B_hadronic_B_delta_R, weight);
  h_leptonic_B_hadronic_first_jet_delta_R_plus_leptonic_B_hadronic_second_jet_delta_R->Fill(leptonic_B_hadronic_first_jet_delta_R_plus_leptonic_B_hadronic_second_jet_delta_R, weight);
  h_leptonic_B_hadronic_first_jet_delta_R->Fill(leptonic_B_hadronic_first_jet_delta_R, weight);
  h_leptonic_B_hadronic_second_jet_delta_R->Fill(leptonic_B_hadronic_second_jet_delta_R, weight);
  h_hadronic_B_hadronic_first_jet_delta_R_plus_hadronic_B_hadronic_second_jet_delta_R->Fill(hadronic_B_hadronic_first_jet_delta_R_plus_hadronic_B_hadronic_second_jet_delta_R, weight);
  h_hadronic_B_hadronic_first_jet_delta_R->Fill(hadronic_B_hadronic_first_jet_delta_R, weight);
  h_hadronic_B_hadronic_second_jet_delta_R->Fill(hadronic_B_hadronic_second_jet_delta_R, weight);
  h_hadronic_first_jet_hadronic_second_jet_delta_R->Fill(hadronic_first_jet_hadronic_second_jet_delta_R, weight);
  //h_leptonic_W_hadronic_W_delta_R->Fill(, weight);
  //h_leptonic_T_hadronic_T_delta_R->Fill(, weight);
  //h_cos_theta_leading_top_resonance->Fill(, weight);
}

void BkgVsTTBDTReader::writePlots(TFile* f) {
  f->cd();

  h_aplanarity->Write();
  h_circularity->Write();
  h_sphericity->Write();
  h_mean_csv->Write();
  h_n_btagged_jets->Write();
  h_st->Write();
  h_theta_lepton->Write();
  h_MET->Write();
  //h_neutrino_pt->Write();
  h_lepton_pt->Write();
  h_lepton_eta->Write();
  h_leptonic_B_eta->Write();
  h_hadronic_B_eta->Write();
  h_hadronic_first_jet_pt_plus_hadronic_second_jet_pt->Write();
  h_hadronic_first_jet_eta_plus_hadronic_second_jet_eta->Write();
  h_leptonic_T_pt->Write();
  h_leptonic_T_transverse_mass->Write();
  h_hadronic_T_pt->Write();
  //h_resonance_pt->Write();
  //h_resonance_eta->Write();
  h_neutrino_leptonic_B_delta_R->Write();
  h_neutrino_lepton_delta_R->Write();
  h_neutrino_hadronic_B_delta_R->Write();
  h_neutrino_hadronic_first_jet_delta_R_plus_neutrino_hadronic_second_jet_delta_R->Write();
  h_neutrino_hadronic_first_jet_delta_R->Write();
  h_neutrino_hadronic_second_jet_delta_R->Write();
  h_lepton_leptonic_B_delta_R->Write();
  h_lepton_hadronic_B_delta_R->Write();
  h_lepton_hadronic_first_jet_delta_R_plus_lepton_hadronic_second_jet_delta_R->Write();
  h_lepton_hadronic_first_jet_delta_R->Write();
  h_lepton_hadronic_second_jet_delta_R->Write();
  h_leptonic_B_hadronic_B_delta_R->Write();
  h_leptonic_B_hadronic_first_jet_delta_R_plus_leptonic_B_hadronic_second_jet_delta_R->Write();
  h_leptonic_B_hadronic_first_jet_delta_R->Write();
  h_leptonic_B_hadronic_second_jet_delta_R->Write();
  h_hadronic_B_hadronic_first_jet_delta_R_plus_hadronic_B_hadronic_second_jet_delta_R->Write();
  h_hadronic_B_hadronic_first_jet_delta_R->Write();
  h_hadronic_B_hadronic_second_jet_delta_R->Write();
  h_hadronic_first_jet_hadronic_second_jet_delta_R->Write();
  //h_leptonic_W_hadronic_W_delta_R->Write();
  //h_leptonic_T_hadronic_T_delta_R->Write();
  //h_cos_theta_leading_top_resonance->Write();

}
