#pragma once

#include <MVAReader.h>

#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
typedef ROOT::Math::PtEtaPhiEVector LorentzVector;

class BkgVsTTBDTReader: public MVAReader {
  public:
    BkgVsTTBDTReader(const std::vector<std::string>& inputFiles);

  protected:
    virtual void initTrees() override;
    virtual void setupVariables() override;
    virtual void computeVariables() override;

  private:

    TChain* m_mtt;
    TChain* m_jets;

    // Tree variables
    TClonesArray* selectedLeptonP4_AfterChi2 = nullptr;
    LorentzVector* selectedFirstJetP4_AfterReco = nullptr;
    LorentzVector* selectedSecondJetP4_AfterReco = nullptr;
    LorentzVector* selectedHadronicBP4_AfterReco = nullptr;
    LorentzVector* selectedLeptonicBP4_AfterReco = nullptr;
    LorentzVector* selectedNeutrinoP4_AfterReco = nullptr;
    TClonesArray* jets_p4 = nullptr;
    uint32_t n_jets;
    float jets_CSV_discriminant[100];
    
    
    // BDT variables
    float aplanarity;
    float circularity;
    float sphericity;
    float mean_csv;
    int   n_btagged_jets;
    float n_btagged_jets_float;
    float st;
    float theta_lepton;
    float MET;
    //float neutrino_pt;
    float lepton_pt;
    float lepton_eta;
    float leptonic_B_eta;
    float hadronic_B_eta;
    float hadronic_first_jet_pt_plus_hadronic_second_jet_pt;
    float hadronic_first_jet_eta_plus_hadronic_second_jet_eta;
    float leptonic_T_pt;
    float leptonic_T_transverse_mass;
    float hadronic_T_pt;
    //float resonance_pt;
    //float resonance_eta;
    float neutrino_leptonic_B_delta_R;
    float neutrino_lepton_delta_R;
    float neutrino_hadronic_B_delta_R;
    float neutrino_hadronic_first_jet_delta_R_plus_neutrino_hadronic_second_jet_delta_R;
    float neutrino_hadronic_first_jet_delta_R;
    float neutrino_hadronic_second_jet_delta_R;
    float lepton_leptonic_B_delta_R;
    float lepton_hadronic_B_delta_R;
    float lepton_hadronic_first_jet_delta_R_plus_lepton_hadronic_second_jet_delta_R;
    float lepton_hadronic_first_jet_delta_R;
    float lepton_hadronic_second_jet_delta_R;
    float leptonic_B_hadronic_B_delta_R;
    float leptonic_B_hadronic_first_jet_delta_R_plus_leptonic_B_hadronic_second_jet_delta_R;
    float leptonic_B_hadronic_first_jet_delta_R;
    float leptonic_B_hadronic_second_jet_delta_R;
    float hadronic_B_hadronic_first_jet_delta_R_plus_hadronic_B_hadronic_second_jet_delta_R;
    float hadronic_B_hadronic_first_jet_delta_R;
    float hadronic_B_hadronic_second_jet_delta_R;
    float hadronic_first_jet_hadronic_second_jet_delta_R;
    //float leptonic_W_hadronic_W_delta_R;
    //float leptonic_T_hadronic_T_delta_R;
    //float cos_theta_leading_top_resonance;
};
