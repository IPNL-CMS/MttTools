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

#include <tclap/CmdLine.h>

#include "../PUReweighting/PUReweighter.h"
#include "TopTriggerEfficiencyProvider.h"

#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
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

void createTree(const std::vector<std::string>& inputFiles, const std::string& outputFile, int64_t maxEntries, float weight, bool isSignal, bool mva, bool chi2, bool kf, bool hybrid) {
  TChain* mtt = loadChain(inputFiles, "Mtt");
  TChain* jets = loadChain(inputFiles, "jet_PF");

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
  if (kf) {
    setBranchAddress(mtt, "selectedNeutrinoP4_AfterKF", neutrino.p4);
    setBranchAddress(mtt, "selectedLeptonP4_AfterKF", lepton.p4);
    setBranchAddress(mtt, "selectedLeptonicBP4_AfterKF", leptonic_B.p4);
    setBranchAddress(mtt, "selectedHadronicBP4_AfterKF", hadronic_B.p4);
    setBranchAddress(mtt, "selectedFirstJetP4_AfterKF", hadronic_first_jet.p4);
    setBranchAddress(mtt, "selectedSecondJetP4_AfterKF", hadronic_second_jet.p4);
    setBranchAddress(mtt, "numComb_kf", numComb);
  } else if (chi2) {
    setBranchAddress(mtt, "selectedNeutrinoP4_AfterChi2", neutrino.p4);
    setBranchAddress(mtt, "selectedLeptonP4", chi2_lepton_p4);
    setBranchAddress(mtt, "selectedLeptonicBP4_AfterChi2", leptonic_B.p4);
    setBranchAddress(mtt, "selectedHadronicBP4_AfterChi2", hadronic_B.p4);
    setBranchAddress(mtt, "selectedFirstJetP4_AfterChi2", hadronic_first_jet.p4);
    setBranchAddress(mtt, "selectedSecondJetP4_AfterChi2", hadronic_second_jet.p4);
    setBranchAddress(mtt, "numComb_chi2", numComb);
  }

  float aplanarity;
  setBranchAddress(mtt, "aplanarity", aplanarity);

  float circularity;
  setBranchAddress(mtt, "circularity", circularity);

  float sphericity;
  setBranchAddress(mtt, "sphericity", sphericity);

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

  // Create output tree
  TFile* output = TFile::Open(outputFile.c_str(), "recreate");
  TTree* tree = new TTree(isSignal ? "signal" : "background", "");
  tree->SetAutoSave(0);

  createBranch(tree, "weight", weight);
  createBranch(tree, "aplanarity", aplanarity);
  createBranch(tree, "circularity", circularity);
  createBranch(tree, "sphericity", sphericity);

  float mean_CSV = 0;
  createBranch(tree, "mean_csv", mean_CSV);

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

  for (uint64_t i = 0; i < entries; i++) {

    if ((i % 10000) == 0) {
      std::cout << "Processing entry " << i + 1 << " / " << entries << " (" << (float) i / entries * 100 << "%)" << std::endl;
    }

    mtt->GetEntry(i);
    jets->GetEntry(i);

    if (numComb <= 0)
      continue;

    // Compute mean CSV value for all jets
    int n = 0;
    mean_CSV = 0;
    for (int i = 0; i < jets_p4->GetEntries(); i++) {
      TLorentzVector* p4 = (TLorentzVector*) (*jets_p4)[i];
      if (p4->Pt() > 30) {
        n++;
        mean_CSV += jets_CSV_discriminant[i];
      }
    }

    mean_CSV /= n;

    if (chi2_lepton_p4) {
      // Copy the quadrivector into lepton_p4
      TLorentzVector* p4 = (TLorentzVector*) (*chi2_lepton_p4)[0];
      lepton.p4 = new LorentzVector(p4->Pt(), p4->Eta(), p4->Phi(), p4->Energy());
    }

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

    createTree(inputFiles, outputFileArg.getValue(), maxEntriesArg.getValue(), weightArg.getValue(), isSignal, mvaArg.getValue(), chi2Arg.getValue(), kfArg.getValue(), hybridArg.getValue()); 

  } catch (TCLAP::ArgException& e) {
    std::cout << e.what() << std::endl;
  }

}
