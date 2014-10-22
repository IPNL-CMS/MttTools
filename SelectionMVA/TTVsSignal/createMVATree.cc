#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <algorithm>
#include <string>
#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TRandom2.h>

#include <tclap/CmdLine.h>

#include <PUReweighter.h>
#include "TopTriggerEfficiencyProvider.h"
#include "ExtractorPostprocessing.h"

#include <BkgVsTTBDTReader.h>
#include <BDTCuts.h>

#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <Math/Vector3Dfwd.h>
typedef ROOT::Math::PtEtaPhiEVector LorentzVector;

PUProfile puProfile;

template<typename T>
void createBranch(TTree* t, const std::string name, T& ptr) {
  t->Branch(name.c_str(), &ptr, (name + "/F").c_str());
}

struct ObjectP4 {
  LorentzVector* p4;
  std::string name;
};

struct Object {
  float pt;
  float eta;
  float phi;
  float rapidity;
  float mass;
  float transverse_mass;

  ObjectP4* p4;

  void create(TTree* tree, const std::string& prefix, ObjectP4* p4) {
    createBranch(tree, prefix + "_pt", pt);
    createBranch(tree, prefix + "_eta", eta);
    createBranch(tree, prefix + "_phi", phi);
    createBranch(tree, prefix + "_rapidity", rapidity);
    createBranch(tree, prefix + "_mass", mass);
    createBranch(tree, prefix + "_transverse_mass", transverse_mass);

    this->p4 = p4;
  }

  void fill() {
    pt = p4->p4->Pt();
    eta = p4->p4->Eta();
    phi = p4->p4->Phi();
    rapidity = p4->p4->Rapidity();
    mass = p4->p4->M();
    transverse_mass = p4->p4->Mt();
  }
};

struct Relation {
  float delta_R;
  float delta_eta;
  float delta_phi;
  float delta_rapidity;

  ObjectP4* first_p4;
  ObjectP4* second_p4;

  void create(TTree* tree, const std::string& prefix1, const std::string& prefix2, ObjectP4* first_p4, ObjectP4* second_p4) {
  
    const std::string prefix = prefix1 + "_" + prefix2;
    createBranch(tree, prefix + "_delta_R", delta_R);
    createBranch(tree, prefix + "_delta_eta", delta_eta);
    createBranch(tree, prefix + "_delta_phi", delta_phi);
    createBranch(tree, prefix + "_delta_rapidity", delta_rapidity);

    this->first_p4 = first_p4;
    this->second_p4 = second_p4;
  }

  void fill() {
    delta_R = ROOT::Math::VectorUtil::DeltaR(*first_p4->p4, *second_p4->p4);
    delta_phi = ROOT::Math::VectorUtil::DeltaPhi(*first_p4->p4, *second_p4->p4);
    delta_eta = first_p4->p4->Eta() - second_p4->p4->Eta();
    delta_rapidity = first_p4->p4->Rapidity() - second_p4->p4->Rapidity();
  }
};

void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

TChain* loadChain(const std::vector<std::string>& inputFiles, const std::string& name) {

  TChain* chain = new TChain(name.c_str());

  for (const std::string& file: inputFiles) {
    chain->Add(file.c_str());
  }

  chain->SetCacheSize(30*1024*1024);
  chain->SetBranchStatus("*", 0);

  return chain;
}

void activateBranch(TBranch* branch) {
  branch->SetStatus(1);
  TObjArray* objArray = branch->GetListOfBranches();
  for (int i = 0; i < objArray->GetEntries(); i++) {
    TBranch* b = static_cast<TBranch*> (objArray->At(i));
    if (b) {
      activateBranch(b);
    }
  }
}

template<typename T>
void setBranchAddress(TChain* chain, const std::string& name, T& address) {

  TBranch* branch = chain->FindBranch(name.c_str());
  if (!branch)
    return;

  activateBranch(branch);
  chain->AddBranchToCache(name.c_str(), true);
  chain->SetBranchAddress(name.c_str(), &address);
}

void createTree(const std::vector<std::string>& inputFiles, const std::string& outputFile, int64_t maxEntries, float weight, bool isSignal, bool isSemiMu, int btag, bool mva, bool chi2, bool kf, bool hybrid, const std::string& bdtWeights) {
  TChain* mtt = loadChain(inputFiles, "Mtt");
  TChain* jets = loadChain(inputFiles, "jet_PF");
  TChain* events = loadChain(inputFiles, "event");
  TChain* vertices = loadChain(inputFiles, "Vertices");

  BkgVsTTBDTReader bkgVsTTBDTReader(inputFiles);
  bkgVsTTBDTReader.initMVA(bdtWeights);

  BDTCuts bdtCuts;
  float background_bdt_cut = bdtCuts.getCut(BDTType::BACKGROUND, -1);  

  // Retrieve P4 of selected objects
  ObjectP4 neutrino = {nullptr, "neutrino"};
  ObjectP4 lepton = {nullptr, "lepton"};
  ObjectP4 leptonic_B = {nullptr, "leptonic_B"};
  ObjectP4 hadronic_B = {nullptr, "hadronic_B"};
  ObjectP4 hadronic_first_jet = {nullptr, "hadronic_first_jet"};
  ObjectP4 hadronic_second_jet = {nullptr, "hadronic_second_jet"};

  ObjectP4 leptonic_W = {nullptr, "leptonic_W"};
  ObjectP4 hadronic_W = {nullptr, "hadronic_W"};

  ObjectP4 leptonic_T = {nullptr, "leptonic_T"};
  ObjectP4 hadronic_T = {nullptr, "hadronic_T"};

  ObjectP4 resonance = {nullptr, "resonance"};

  TClonesArray* chi2_lepton_p4 = nullptr;

  int numComb = 0;
  float discriminant = 0;
  if (kf) {
    setBranchAddress(mtt, "selectedNeutrinoP4_AfterKF", neutrino.p4);
    setBranchAddress(mtt, "selectedLeptonP4_AfterKF", lepton.p4);
    setBranchAddress(mtt, "selectedLeptonicBP4_AfterKF", leptonic_B.p4);
    setBranchAddress(mtt, "selectedHadronicBP4_AfterKF", hadronic_B.p4);
    setBranchAddress(mtt, "selectedFirstJetP4_AfterKF", hadronic_first_jet.p4);
    setBranchAddress(mtt, "selectedSecondJetP4_AfterKF", hadronic_second_jet.p4);
    setBranchAddress(mtt, "numComb_kf", numComb);
    setBranchAddress(mtt, "kf_chisquare", discriminant);
  } else if (chi2) {
    setBranchAddress(mtt, "selectedNeutrinoP4_AfterChi2", neutrino.p4);
    setBranchAddress(mtt, "selectedLeptonP4", chi2_lepton_p4);
    setBranchAddress(mtt, "selectedLeptonicBP4_AfterChi2", leptonic_B.p4);
    setBranchAddress(mtt, "selectedHadronicBP4_AfterChi2", hadronic_B.p4);
    setBranchAddress(mtt, "selectedFirstJetP4_AfterChi2", hadronic_first_jet.p4);
    setBranchAddress(mtt, "selectedSecondJetP4_AfterChi2", hadronic_second_jet.p4);
    setBranchAddress(mtt, "numComb_chi2", numComb);
    setBranchAddress(mtt, "bestSolChi2", discriminant);
  }

  float aplanarity;
  setBranchAddress(mtt, "aplanarity", aplanarity);

  float circularity;
  setBranchAddress(mtt, "circularity", circularity);

  float sphericity;
  setBranchAddress(mtt, "sphericity", sphericity);

  int lepton_index;
  setBranchAddress(mtt, "selectedLeptonIndex_in_array", lepton_index);

  int n_leptons;
  if (isSemiMu) {
    setBranchAddress(mtt, "nGoodMuons", n_leptons);
  } else {
    setBranchAddress(mtt, "nGoodElectrons", n_leptons);
  }

  float lepton_rel_iso_array[20];
  setBranchAddress(mtt, "muRelIso", lepton_rel_iso_array);

  float muonPt[20];
  float muonEta[20];
  float electronPt[20];
  float electronEta[20];
  if (isSemiMu) {
    setBranchAddress(mtt, "muonPt", muonPt);
    setBranchAddress(mtt, "muonEta", muonEta);
  } else {
    setBranchAddress(mtt, "electronPt", electronPt);
    setBranchAddress(mtt, "electronEta", electronEta);
  }

  float met;
  setBranchAddress(mtt, "MET", met);

  int32_t n_jets_mtt;
  setBranchAddress(mtt, "nJets", n_jets_mtt);

  float jetEta[100];
  setBranchAddress(mtt, "jetEta", jetEta);

  float pt_1stJet = 0;
  setBranchAddress(mtt, "1stjetpt", pt_1stJet);

  float pt_2ndJet = 0;
  setBranchAddress(mtt, "2ndjetpt", pt_2ndJet);

  float pt_3rdJet = 0;
  setBranchAddress(mtt, "3rdjetpt", pt_3rdJet);

  float pt_4thJet = 0;
  setBranchAddress(mtt, "4thjetpt", pt_4thJet);

  int n_btagged_jet = 0;
  if (btag >= 0)
    setBranchAddress(mtt, "nBtaggedJets_CSVM", n_btagged_jet);

  uint32_t n_jets;
  setBranchAddress(jets, "n_jets", n_jets);

  TClonesArray* jets_p4 = nullptr;
  setBranchAddress(jets, "jet_4vector", jets_p4);

  float jets_CSV_discriminant[100];
  setBranchAddress(jets, "jet_btag_CSV", jets_CSV_discriminant);

  std::vector<ObjectP4*> objects_p4 = {
    &neutrino,
    &lepton,
    &leptonic_B,
    &hadronic_B,
    &hadronic_first_jet,
    &hadronic_second_jet,

    &leptonic_W,
    &hadronic_W,

    &leptonic_T,
    &hadronic_T,

    &resonance
  };

  float lepton_weight = 0;
  setBranchAddress(mtt, "lepton_weight", lepton_weight);

  float btag_weight = 0;
  setBranchAddress(mtt, "btag_weight", btag_weight);

  float n_trueInteractions = 0;
  setBranchAddress(events, "nTrueInteractions", n_trueInteractions);
  
  uint32_t run = 0;
  setBranchAddress(events, "run", run);

  float generator_weight = 0;
  setBranchAddress(events, "generator_weight", generator_weight);

  uint32_t n_vertices = 0;
  setBranchAddress(vertices, "n_vertices", n_vertices);

  float output_weight = 0;

  PUReweighter* puReweigher = new PUReweighter(isSemiMu, puProfile, Systematic::NOMINAL);
  TopTriggerEfficiencyProvider::JES triggerJESSyst = TopTriggerEfficiencyProvider::NOMINAL;

  float lumi_run2012_A = 0;
  float lumi_run2012_B = 0;
  float lumi_run2012_C = 0;
  float lumi_run2012_D = 0;

  if (isSemiMu) {
    lumi_run2012_A = 0.876225;
    lumi_run2012_B = 4.412;
    lumi_run2012_C = 7.044;
    lumi_run2012_D = 7.368;
  } else {
    lumi_run2012_A = 0.876225;
    lumi_run2012_B = 4.399;
    lumi_run2012_C = 7.022;
    lumi_run2012_D = 7.369;
  }

  float lumi_run2012_AB = lumi_run2012_A + lumi_run2012_B;
  float lumi_run2012_CD = lumi_run2012_C + lumi_run2012_D;
  float lumi_total = lumi_run2012_AB + lumi_run2012_CD;

  float lumi_run2012_AB_over_total = lumi_run2012_AB / lumi_total;

  std::shared_ptr<TopTriggerEfficiencyProvider> m_trigger_efficiency_provider = std::make_shared<TopTriggerEfficiencyProvider>();
  TRandom2 random_generator;

  // Create output tree
  TFile* output = TFile::Open(outputFile.c_str(), "recreate");
  TTree* tree = new TTree(isSignal ? "signal" : "background", "");
  tree->SetAutoSave(0);

  createBranch(tree, "weight", weight);
  createBranch(tree, "aplanarity", aplanarity);
  createBranch(tree, "circularity", circularity);
  createBranch(tree, "sphericity", sphericity);
  createBranch(tree, "discriminant", discriminant);

  float mean_CSV = 0;
  createBranch(tree, "mean_csv", mean_CSV);
  
  createBranch(tree, "met", met);

  float HT = 0;
  createBranch(tree, "ht", HT);

  float ST = 0;
  createBranch(tree, "st", ST);

  float lepton_rel_iso = 0;
  createBranch(tree, "lepton_rel_iso", lepton_rel_iso);

  float theta_lepton = 0;
  createBranch(tree, "theta_lepton", theta_lepton);

  float cos_theta_lepton = 0;
  createBranch(tree, "cos_theta_lepton", cos_theta_lepton);

  float cos_theta_leading_top_resonance = 0;
  createBranch(tree, "cos_theta_leading_top_resonance", cos_theta_leading_top_resonance);

  const auto& createObjectBranches = [&](const std::string& prefix, ObjectP4* p4) -> Object* {
    Object* obj = new Object();
    obj->create(tree, prefix, p4);
    return obj;
  };

  std::vector<Object*> objects;
  for (auto& object_p4: objects_p4) {
    objects.push_back(createObjectBranches(object_p4->name, object_p4));
  }

  const auto& createRelationBranches = [&](const std::string& prefix1, const std::string& prefix2,
      ObjectP4* first_p4, ObjectP4* second_p4) -> Relation* {
    Relation* r = new Relation();
    r->create(tree, prefix1, prefix2, first_p4, second_p4);
    return r;
  };

  std::vector<Relation*> relations;
  for (auto i = objects_p4.begin(); i != objects_p4.end(); ++i) {
    for (auto j = i; ++j != objects_p4.end(); /* */) {
      relations.push_back(createRelationBranches((*i)->name, (*j)->name, *i, *j));
    }
  }
  
  uint64_t entries = mtt->GetEntries();
  if (maxEntries > 0)
    entries = maxEntries;

  ExtractorPostprocessing selection;

  for (uint64_t i = 0; i < entries; i++) {

    if ((i % 10000) == 0) {
      std::cout << "Processing entry " << i + 1 << " / " << entries << " (" << (float) i / entries * 100 << "%)" << std::endl;
    }

    mtt->GetEntry(i);
    jets->GetEntry(i);
    events->GetEntry(i);
    vertices->GetEntry(i);
    
    // Choose if we are run2012 A+B, or C+D

    bool isRun2012AB = false;
    double r = random_generator.Rndm();
    if (r < lumi_run2012_AB_over_total)
      isRun2012AB = true;

    // Selection

    if (numComb <= 0)
      continue;

    if (! selection.passJetsSel(pt_1stJet, pt_2ndJet, pt_3rdJet, pt_4thJet, isRun2012AB))
      continue;

    if (btag >= 0) {

      switch (btag) {
        case 0:
          if (n_btagged_jet != 0)
            continue;

          break;

        case 1:
          if (n_btagged_jet != 1)
            continue;

          break;

        case 2:
          if (n_btagged_jet < 2)
            continue;

          break;

        default:
          continue;
      }

    }

    float discriminant = bkgVsTTBDTReader.evaluate(i);
    if (discriminant <= background_bdt_cut)
      continue;

    double pt_lepton = 0;
    double eta_lepton = 0;
    if (isSemiMu) {
      pt_lepton = muonPt[lepton_index];
      eta_lepton = muonEta[lepton_index];
    } else {
      pt_lepton = electronPt[lepton_index];
      eta_lepton = electronEta[lepton_index];
    }

    output_weight = 1.;
    if (isRun2012AB) {
      m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunA, lumi_run2012_A);
      m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunB, lumi_run2012_B);
      m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunC, 0);
      m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunD, 0);
    } else {
      m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunA, 0);
      m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunB, 0);
      m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunC, lumi_run2012_C);
      m_trigger_efficiency_provider->setLumi(TopTriggerEfficiencyProvider::RunD, lumi_run2012_D);
    }

    // Compute trigger weight
    std::vector<double> triggerWeights = m_trigger_efficiency_provider->get_weight(pt_lepton, eta_lepton, pt_4thJet, jetEta[3], n_vertices, n_jets_mtt, isSemiMu, triggerJESSyst);
    double triggerWeight = triggerWeights[0];

    output_weight *= puReweigher->weight(n_trueInteractions) /* * generator_weight */ * weight * triggerWeight * lepton_weight * btag_weight;

    // Compute mean CSV value for all jets
    int n = 0;
    mean_CSV = 0;
    HT = 0;
    ST = 0;
    for (int i = 0; i < jets_p4->GetEntries(); i++) {
      TLorentzVector* p4 = (TLorentzVector*) (*jets_p4)[i];
      if (p4->Pt() > 30 && jets_CSV_discriminant[i] > 0.244) {
        n++;
        mean_CSV += jets_CSV_discriminant[i];
      }
      HT += p4->Pt();
    }

    if (n > 0)
      mean_CSV /= n;

    if (chi2_lepton_p4) {
      // Copy the quadrivector into lepton_p4
      TLorentzVector* p4 = (TLorentzVector*) (*chi2_lepton_p4)[0];
      lepton.p4 = new LorentzVector(p4->Pt(), p4->Eta(), p4->Phi(), p4->Energy());
    }

    ST = HT + lepton.p4->Pt();

    lepton_rel_iso = lepton_rel_iso_array[lepton_index];

    LorentzVector leptonic_W_p4 = *lepton.p4 + *neutrino.p4;
    leptonic_W.p4 = &leptonic_W_p4;

    LorentzVector hadronic_W_p4 = *hadronic_first_jet.p4 + *hadronic_second_jet.p4;
    hadronic_W.p4 = &hadronic_W_p4;

    LorentzVector leptonic_T_p4 = leptonic_W_p4 + *leptonic_B.p4;
    leptonic_T.p4 = &leptonic_T_p4;

    LorentzVector hadronic_T_p4 = hadronic_W_p4 + *hadronic_B.p4;
    hadronic_T.p4 = &hadronic_T_p4;

    LorentzVector resonance_p4 = hadronic_T_p4 + leptonic_T_p4;
    resonance.p4 = &resonance_p4;

    // Compute theta lepton
    // θ_l is the production angle of the positively charged lepton in the rest frame of its parent top
    // or anti-top, with respect to the direction of the parent top or anti-top in the tt rest frame
    //
    // Find the boost to go in tt rest frame, and apply it to the leptonic top to have the top direction
    ROOT::Math::Boost resonance_in_rest_frame_boost = ROOT::Math::Boost(resonance_p4.BoostToCM());
    LorentzVector leptonic_T_p4_in_resonance_rest_frame = resonance_in_rest_frame_boost(leptonic_T_p4);

    ROOT::Math::XYZVector leptonic_T_direction = leptonic_T_p4_in_resonance_rest_frame.Vect();
    
    // Find the boost to go in the top rest frame, and apply it to the lepton to have the lepton direction
    ROOT::Math::Boost leptonic_T_in_rest_frame_boost = ROOT::Math::Boost(leptonic_T_p4.BoostToCM());
    LorentzVector lepton_p4_in_resonance_rest_frame = leptonic_T_in_rest_frame_boost(*lepton.p4);

    ROOT::Math::XYZVector lepton_direction = lepton_p4_in_resonance_rest_frame.Vect();

    cos_theta_lepton = ROOT::Math::VectorUtil::CosTheta(leptonic_T_direction, lepton_direction);
    theta_lepton = std::acos(cos_theta_lepton);

    // Compute cos(highest_pt_top, resonance) in resonance rest frame
    LorentzVector& leading_top = (leptonic_T_p4.Pt() > hadronic_T_p4.Pt())
      ? leptonic_T_p4
      : hadronic_T_p4;

    cos_theta_leading_top_resonance = ROOT::Math::VectorUtil::CosTheta(leading_top.Vect(), resonance_p4.Vect());    

    for (auto& object: objects) {
      object->fill();
    }

    for (auto& relation: relations) {
      relation->fill();
    }

    if (chi2_lepton_p4) {
      delete lepton.p4;
    }

    tree->Fill();
  }

  tree->Write();
  output->Close();

  delete output;

  delete mtt;
  delete jets;
  delete events;
  delete vertices;
}

int main(int argc, char** argv) {

  try {
    TCLAP::CmdLine cmd("Create a tree from a skim suitable to train a TMVA MVA algorithm", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", true, "", "string", cmd);

    TCLAP::SwitchArg signalArg("", "signal", "Is this signal?", false);
    TCLAP::SwitchArg bkgArg("", "background", "Is this background?", false);

    cmd.xorAdd(signalArg, bkgArg);

    TCLAP::SwitchArg semimuArg("", "semimu", "Semi mu?", false);
    TCLAP::SwitchArg semieArg("", "semie", "Semi e?", false);

    cmd.xorAdd(semimuArg, semieArg);

    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets to select", true, -1, "int", cmd);

    TCLAP::ValueArg<std::string> bdtWeightsArg("", "bdt-weights", "XML file containing the BDT weighted created by TMVA", true, "", "string", cmd);

    //TCLAP::ValueArg<std::string> typeArg("", "type", "current inputfile type (semie or semimu)", true, "", "string", cmd);
    TCLAP::ValueArg<std::string> pileupArg("", "pileup", "PU profile used for MC production", false, "S10", "string", cmd);
    TCLAP::ValueArg<int> maxEntriesArg("n", "", "Maximal number of entries to process", false, -1, "int", cmd);
    TCLAP::ValueArg<double> weightArg("", "weight", "Sample weight", false, 1., "double", cmd);

    //TCLAP::ValueArg<std::string> pdfSystArg("", "pdf-syst", "PDF systematic to compute", false, "nominal", "string", cmd);
    //TCLAP::ValueArg<std::string> jecSystArg("", "jec-syst", "Computing trigger weight for this JEC up / down", false, "nominal", "string", cmd);
    //TCLAP::ValueArg<std::string> triggerSystArg("", "trigger-syst", "Computing trigger weight systematic", false, "nominal", "string", cmd);
    //TCLAP::ValueArg<std::string> pileupSystArg("", "pileup-syst", "PU profile to use for pileup reweigthing", false, "nominal", "string", cmd);
    //TCLAP::ValueArg<std::string> btagSystArg("", "btag-syst", "Compute btag weight systematic", false, "nominal", "string", cmd);
    //TCLAP::ValueArg<std::string> leptonSystArg("", "lepton-syst", "Compute lepton weight systematic", false, "nominal", "string", cmd);

    //TCLAP::SwitchArg skimArg("", "skim", "Run over a skimmed file", cmd, false);

    TCLAP::SwitchArg chi2Arg("", "chi2", "Use chi2 sorting algorithm", false);
    TCLAP::SwitchArg mvaArg("", "mva", "Use MVA instead of chi2", false);
    TCLAP::SwitchArg kfArg("", "kf", "Use KF instead of chi2", false);
    TCLAP::SwitchArg hybridArg("", "hybrid", "Use hybrid method for sorting algorithm", false);
    std::vector<TCLAP::Arg*>  xorlist;
    xorlist.push_back(&chi2Arg);
    xorlist.push_back(&mvaArg);
    xorlist.push_back(&kfArg);
    xorlist.push_back(&hybridArg);
    cmd.xorAdd( xorlist );

    //TCLAP::SwitchArg zprimeArg("", "zprime", "Is this Zprime analysis?", false);
    //TCLAP::SwitchArg higgsArg("", "higgs", "Is this Higgs analysis?", false);
    //cmd.xorAdd(zprimeArg, higgsArg);

    cmd.parse(argc, argv);

    std::string p = pileupArg.getValue();
    std::transform(p.begin(), p.end(), p.begin(), ::tolower);
    if (p == "s6")
      puProfile = PUProfile::S6;
    else if (p == "s7")
      puProfile = PUProfile::S7;
    else if (p == "s10")
      puProfile = PUProfile::S10;

    //std::string triggerSyst = triggerSystArg.getValue();
    //std::transform(triggerSyst.begin(), triggerSyst.end(), triggerSyst.begin(), ::tolower);
    //if (triggerSyst != "nominal" && triggerSyst != "up" && triggerSyst != "down") {
      //std::cerr << "--trigger-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      //exit(1);
    //}


    //std::string jecSyst = jecSystArg.getValue();
    //std::transform(jecSyst.begin(), jecSyst.end(), jecSyst.begin(), ::tolower);
    //if (jecSyst != "nominal" && jecSyst != "up" && jecSyst != "down") {
      //std::cerr << "--jec-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      //exit(1);
    //}

    //std::string puSyst = pileupSystArg.getValue();
    //std::transform(puSyst.begin(), puSyst.end(), puSyst.begin(), ::tolower);
    //if (puSyst != "nominal" && puSyst != "up" && puSyst != "down") {
      //std::cerr << "--pilup-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      //exit(1);
    //}

    //std::string pdfSyst = pdfSystArg.getValue();
    //std::transform(pdfSyst.begin(), pdfSyst.end(), pdfSyst.begin(), ::tolower);
    //if (pdfSyst != "nominal" && pdfSyst != "up" && pdfSyst != "down") {
      //std::cerr << "--pdf-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      //exit(1);
    //}

    //std::string leptonSyst = leptonSystArg.getValue();
    //std::transform(leptonSyst.begin(), leptonSyst.end(), leptonSyst.begin(), ::tolower);
    //if (leptonSyst != "nominal" && leptonSyst != "up" && leptonSyst != "down") {
      //std::cerr << "--lepton-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      //exit(1);
    //}

    //std::string btagSyst = btagSystArg.getValue();
    //std::transform(btagSyst.begin(), btagSyst.end(), btagSyst.begin(), ::tolower);
    //if (btagSyst != "nominal" && btagSyst != "up" && btagSyst != "down") {
      //std::cerr << "--btag-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      //exit(1);
    //}
    
    bool isSignal = signalArg.isSet();
    //bool isZprime = zprimeArg.isSet();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    createTree(inputFiles, outputFileArg.getValue(), maxEntriesArg.getValue(), weightArg.getValue(), isSignal, semimuArg.isSet(), btagArg.getValue(), mvaArg.getValue(), chi2Arg.getValue(), kfArg.getValue(), hybridArg.getValue(), bdtWeightsArg.getValue()); 

  } catch (TCLAP::ArgException& e) {
    std::cout << e.what() << std::endl;
  }

}
