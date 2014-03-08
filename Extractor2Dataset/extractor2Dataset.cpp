//#include <Riostream.h>
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
#include <TRandom2.h>

#include <tclap/CmdLine.h>

#include "ExtractorPostprocessing.h"

#include "../PUReweighting/PUReweighter.h"

#include "TopTriggerEfficiencyProvider.h"

bool OVERRIDE_TYPE;
std::string OVERRIDED_TYPE;
PUProfile puProfile;

ExtractorPostprocessing selection;

void reduce(TChain* mtt, TChain* event, TChain* vertices, const std::string& outputFile, bool isData, const std::string& type, int max, double lumi_weight, const std::string& puSyst, const std::string& pdfSyst, const std::string& jecSyst, const std::string& triggerSyst, const std::string& leptonSyst, const std::string& btagSyst, bool useMVA, bool runOnSkim) {

  TString outputFileFormated = OVERRIDE_TYPE ? outputFile : TString::Format(outputFile.c_str(), (isData) ? "" : type.c_str());

  std::map<int, TTree*> outputTrees {
    {0, new TTree("dataset_0btag", "dataset for 0 b-tagged jet") },
      {1, new TTree("dataset_1btag", "dataset for 1 b-tagged jet") },
      {2, new TTree("dataset_2btag", "dataset for at least 2 b-tagged jets") }
  };

  float mtt_afterReco, pt_1stJet, pt_2ndJet, pt_3rdJet, pt_4thJet, bestSolChi2;
  int isSel = 1, nBtaggedJets_CSVM;
  uint32_t run = 0;

  float lepton_weight = 1;
  float btag_weight = 1;

  float lepton_weight_error = 0;
  float btag_weight_error = 0;

  auto SetBranchAddress = [&](const std::string& param, void* address) {
    mtt->SetBranchStatus(param.c_str(), 1);
    mtt->SetBranchAddress(param.c_str(), address);
  };

  mtt->SetBranchStatus("*", 0);
  mtt->SetBranchStatus("1stjetpt", 1);
  mtt->SetBranchStatus("2ndjetpt", 1);
  mtt->SetBranchStatus("3rdjetpt", 1);
  mtt->SetBranchStatus("4thjetpt", 1);
  mtt->SetBranchStatus("isSel", 1);

  if (puProfile == PUProfile::S6)
    mtt->SetBranchStatus("nBtaggedJets_TCHET", 1);
  else
    mtt->SetBranchStatus("nBtaggedJets_CSVM", 1);
  if (! isData) {
    mtt->SetBranchStatus("lepton_weight", 1);
    mtt->SetBranchStatus("btag_weight", 1);

    if (leptonSyst == "up")
      SetBranchAddress("lepton_weight_error_high", &lepton_weight_error);

    if (leptonSyst == "down")
      SetBranchAddress("lepton_weight_error_low", &lepton_weight_error);

    if (btagSyst == "up")
      SetBranchAddress("btag_weight_error_high", &btag_weight_error);

    if (btagSyst == "down")
      SetBranchAddress("btag_weight_error_low", &btag_weight_error);
  }

  int numComb = 0;
  if (useMVA)
    SetBranchAddress("numComb_MVA", &numComb);
  else
    SetBranchAddress("numComb_chi2", &numComb);

  if (useMVA) {
    SetBranchAddress("mtt_AfterMVA", &mtt_afterReco);
  } else {
    SetBranchAddress("mtt_AfterChi2", &mtt_afterReco);
    SetBranchAddress("bestSolChi2", &bestSolChi2);
  }

  mtt->SetBranchAddress("1stjetpt", &pt_1stJet, NULL);
  mtt->SetBranchAddress("2ndjetpt", &pt_2ndJet, NULL);
  mtt->SetBranchAddress("3rdjetpt", &pt_3rdJet, NULL);
  mtt->SetBranchAddress("4thjetpt", &pt_4thJet, NULL);
  mtt->SetBranchAddress("isSel", &isSel, NULL);
  if (puProfile == PUProfile::S6)
    mtt->SetBranchAddress("nBtaggedJets_TCHET", &nBtaggedJets_CSVM, NULL);
  else
    mtt->SetBranchAddress("nBtaggedJets_CSVM", &nBtaggedJets_CSVM, NULL);
  if (!isData) {
    mtt->SetBranchAddress("lepton_weight", &lepton_weight);
    mtt->SetBranchAddress("btag_weight", &btag_weight);
  }

  bool isSemiMu = (type == "semimu");

  float muonPt[100];
  float muonEta[100];
  float electronPt[100];
  float electronEta[100];
  if (isSemiMu) {
    mtt->SetBranchStatus("muonPt", 1);
    mtt->SetBranchStatus("muonEta", 1);
    mtt->SetBranchAddress("muonPt", muonPt);
    mtt->SetBranchAddress("muonEta", muonEta);
  } else {
    mtt->SetBranchStatus("electronPt", 1);
    mtt->SetBranchStatus("electronEta", 1);
    mtt->SetBranchAddress("electronPt", electronPt);
    mtt->SetBranchAddress("electronEta", electronEta);
  }

  int32_t nJets;
  float jetEta[100];
  mtt->SetBranchStatus("nJets", 1);
  mtt->SetBranchStatus("jetEta", 1);
  mtt->SetBranchAddress("nJets", &nJets);
  mtt->SetBranchAddress("jetEta", jetEta);

  float n_trueInteractions;
  float generator_weight;

  event->SetBranchStatus("*", 0);
  event->SetBranchStatus("nTrueInteractions", 1);
  event->SetBranchStatus("run", 1);
  event->SetBranchStatus("generator_weight", 1);

  event->SetBranchAddress("nTrueInteractions", &n_trueInteractions, NULL);
  event->SetBranchAddress("run", &run, NULL);
  event->SetBranchAddress("generator_weight", &generator_weight, NULL);

  int32_t n_vertices = 0;
  vertices->SetBranchStatus("*", 0);
  vertices->SetBranchStatus("n_vertices", 1);
  vertices->SetBranchAddress("n_vertices", &n_vertices);

  float output_weight;
  int lepton = isSemiMu ? 13 : 11;
  for (auto& outputTree: outputTrees) {
    outputTree.second->Branch("mtt", &mtt_afterReco, "mtt/F");
    outputTree.second->Branch("weight", &output_weight, "weight/F");
    outputTree.second->Branch("lepton_type", &lepton, "lepton_type/I");
  }

  PUReweighter* puReweigher = NULL;
  if (! isData) {
    Systematic syst = Systematic::NOMINAL;
    if (puSyst == "up")
      syst = Systematic::UP;
    else if (puSyst == "down")
      syst = Systematic::DOWN;

    puReweigher = new PUReweighter(isSemiMu, puProfile, syst);
  }

  TopTriggerEfficiencyProvider::JES triggerJESSyst = TopTriggerEfficiencyProvider::NOMINAL;
  if (jecSyst == "up")
    triggerJESSyst = TopTriggerEfficiencyProvider::UP;
  else if (jecSyst == "down")
    triggerJESSyst = TopTriggerEfficiencyProvider::DOWN;

  std::vector<float>* pdfWeights = NULL;
  bool doPDFSyst = pdfSyst != "nominal";
  if (doPDFSyst) {
    mtt->SetBranchStatus("pdf_weights", 1);
    mtt->SetBranchAddress("pdf_weights", &pdfWeights);
  }

  int64_t entries = mtt->GetEntries();

  std::map<int, int64_t> selectedEntries;

   // 2012 luminosity
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

  for (int64_t i = 0; i < entries; i++) {
    mtt->GetEntry(i);
    event->GetEntry(i);
    vertices->GetEntry(i);

    if (i % 1000000 == 0) {
      std::cout << "Processing event #" << i + 1 << " over " << entries << " (" << (float) i / entries * 100 << " %)" << std::endl;
    }

    // Choose if we are run2012 A+B, or C+D
    bool isRun2012AB = false;
    if (! isData) {
      double r = random_generator.Rndm();
      if (r < lumi_run2012_AB_over_total)
        isRun2012AB = true;

      //if (mIsMC) {
      //hRunPeriod->Fill( isRun2012AB ? 0 : 1 );
    //}
    } else {
      isRun2012AB = (run <= 196531);
    }

    double ptLepton = 0;
    double etaLepton = 0;
    if (isSemiMu) {
      ptLepton = muonPt[0];
      etaLepton = muonEta[0];
    } else {
      ptLepton = electronPt[0];
      etaLepton = electronEta[0];
    }

    if (isSemiMu && !selection.passMuonSel(ptLepton, etaLepton))
      continue;

    if (!isSemiMu && !selection.passElectronSel(ptLepton, etaLepton))
      continue;

    if (! selection.passExtractorSel(isSel, numComb, mtt_afterReco))
      continue;

    if (! selection.passJetsSel(pt_1stJet, pt_2ndJet, pt_3rdJet, pt_4thJet, isRun2012AB))
      continue;

    if (useMVA) {

    } else {
      if (! selection.passChi2Sel(bestSolChi2))
        continue;
    }

    // Good event
    int index;
    if (nBtaggedJets_CSVM == 0)
      index = 0;
    else if (nBtaggedJets_CSVM == 1)
      index = 1;
    else
      index = 2;

    if (max > 0 && selectedEntries[index] >= max)
      continue;

    output_weight = 1.;
    if (! isData) {
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
      std::vector<double> triggerWeights = m_trigger_efficiency_provider->get_weight(ptLepton, etaLepton, pt_4thJet, jetEta[3], n_vertices, nJets, isSemiMu, triggerJESSyst);
      double triggerWeight = triggerWeights[0];
      if (triggerSyst == "up")
        triggerWeight = triggerWeight + triggerWeights[1];
      else if (triggerSyst == "down")
        triggerWeight = triggerWeight - triggerWeights[1];

      if (leptonSyst == "up")
        lepton_weight += lepton_weight_error;
      else if (leptonSyst == "down")
        lepton_weight -= lepton_weight_error;

      if (btagSyst == "up")
        btag_weight += btag_weight_error;
      else if (btagSyst == "down")
        btag_weight -= btag_weight_error;

      output_weight *= puReweigher->weight(n_trueInteractions) * generator_weight * lumi_weight * triggerWeight * lepton_weight * btag_weight;
    }

    if (doPDFSyst) {
      // 4? pdf systematics
      double sum = 0;
      for (unsigned int i = 0; i < (pdfWeights->size() / 2); i++) {
        int up_index = 2 * i;
        int down_index = up_index + 1;

        double up = (*pdfWeights)[up_index];
        double down = (*pdfWeights)[down_index];

        /*
           std::cout << "up weight: " << up << std::endl;
           std::cout << "down weight: " << down << std::endl;
           */

        if (pdfSyst == "up") {
          sum += pow(std::max(std::max(up - 1, down - 1), 0.), 2);
        } else {
          sum += pow(std::max(std::max(1 - up, 1 - down), 0.), 2);
        }
      }

      double pdf_weight = sqrt(sum);
      //std::cout << "event weight: " << pdf_weight << std::endl;
      if (pdfSyst == "down")
        pdf_weight *= -1;

      output_weight *= (pdf_weight / 1.645) + 1;
    }

    outputTrees[index]->Fill();
    selectedEntries[index]++;
  }

  TFile* output = TFile::Open(outputFileFormated, "recreate");

  for (auto& tree: outputTrees) {
    tree.second->SetDirectory(output);
    tree.second->Write();
  }

  output->Close();
  delete output;

  delete puReweigher;
}

void loadChain(const std::vector<std::string>& inputFiles, TChain*& mtt, TChain*& event, TChain*& vertices) {

  mtt = new TChain("Mtt");
  event = new TChain("event");
  vertices = new TChain("Vertices");

  for (const std::string& file: inputFiles) {
    mtt->Add(file.c_str());
    event->Add(file.c_str());
    vertices->Add(file.c_str());
  }

  event->SetCacheSize(30*1024*1024);
  mtt->SetCacheSize(30*1024*1024);
  vertices->SetCacheSize(30*1024*1024);
}

void reduce(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isData, const std::string& type, int max, double generator_weight, const std::string& puSyst, const std::string& pdfSyst, const std::string& jecSyst, const std::string& triggerSyst, const std::string& leptonSyst, const std::string& btagSyst, bool useMVA, bool runOnSkim) {

  TChain* mtt = NULL, *event = NULL, *vertices = NULL;

  loadChain(inputFiles, mtt, event, vertices);
  reduce(mtt, event, vertices, outputFile, isData, type, max, generator_weight, puSyst, pdfSyst, jecSyst, triggerSyst, leptonSyst, btagSyst, useMVA, runOnSkim);

  delete mtt;
  delete event;
}

void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

int main(int argc, char** argv)
{
  try {
    TCLAP::CmdLine cmd("reduce Zprime dataset", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", true, "", "string", cmd);

    TCLAP::SwitchArg dataArg("", "data", "Is this data?", false);
    TCLAP::SwitchArg mcArg("", "mc", "Is this mc?", false);

    cmd.xorAdd(dataArg, mcArg);

    TCLAP::ValueArg<std::string> typeArg("", "type", "current inputfile type (semie or semimu)", true, "", "string", cmd);
    TCLAP::ValueArg<std::string> pileupArg("", "pileup", "PU profile used for MC production", false, "S10", "string", cmd);
    TCLAP::ValueArg<int> maxEntriesArg("n", "", "Maximal number of entries to process", false, -1, "int", cmd);
    TCLAP::ValueArg<double> generatorWeightArg("", "weight", "MC generator weight", false, 1., "double", cmd);

    TCLAP::ValueArg<std::string> pdfSystArg("", "pdf-syst", "PDF systematic to compute", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> jecSystArg("", "jec-syst", "Computing trigger weight for this JEC up / down", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> triggerSystArg("", "trigger-syst", "Computing trigger weight systematic", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> pileupSystArg("", "pileup-syst", "PU profile to use for pileup reweigthing", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> btagSystArg("", "btag-syst", "Compute btag weight systematic", false, "nominal", "string", cmd);
    TCLAP::ValueArg<std::string> leptonSystArg("", "lepton-syst", "Compute lepton weight systematic", false, "nominal", "string", cmd);

    TCLAP::SwitchArg skimArg("", "skim", "Run over a skimmed file", cmd, false);
    TCLAP::SwitchArg mvaArg("", "mva", "Use MVA instead of chi2", cmd, false);

    cmd.parse(argc, argv);

    std::string p = pileupArg.getValue();
    std::transform(p.begin(), p.end(), p.begin(), ::tolower);
    if (p == "s6")
      puProfile = PUProfile::S6;
    else if (p == "s7")
      puProfile = PUProfile::S7;
    else if (p == "s10")
      puProfile = PUProfile::S10;

    std::string triggerSyst = triggerSystArg.getValue();
    std::transform(triggerSyst.begin(), triggerSyst.end(), triggerSyst.begin(), ::tolower);
    if (triggerSyst != "nominal" && triggerSyst != "up" && triggerSyst != "down") {
      std::cerr << "--trigger-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }


    std::string jecSyst = jecSystArg.getValue();
    std::transform(jecSyst.begin(), jecSyst.end(), jecSyst.begin(), ::tolower);
    if (jecSyst != "nominal" && jecSyst != "up" && jecSyst != "down") {
      std::cerr << "--jec-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    std::string puSyst = pileupSystArg.getValue();
    std::transform(puSyst.begin(), puSyst.end(), puSyst.begin(), ::tolower);
    if (puSyst != "nominal" && puSyst != "up" && puSyst != "down") {
      std::cerr << "--pilup-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    std::string pdfSyst = pdfSystArg.getValue();
    std::transform(pdfSyst.begin(), pdfSyst.end(), pdfSyst.begin(), ::tolower);
    if (pdfSyst != "nominal" && pdfSyst != "up" && pdfSyst != "down") {
      std::cerr << "--pdf-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    std::string leptonSyst = leptonSystArg.getValue();
    std::transform(leptonSyst.begin(), leptonSyst.end(), leptonSyst.begin(), ::tolower);
    if (leptonSyst != "nominal" && leptonSyst != "up" && leptonSyst != "down") {
      std::cerr << "--lepton-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }

    std::string btagSyst = btagSystArg.getValue();
    std::transform(btagSyst.begin(), btagSyst.end(), btagSyst.begin(), ::tolower);
    if (btagSyst != "nominal" && btagSyst != "up" && btagSyst != "down") {
      std::cerr << "--btag-syst can only be 'nominal', 'up' or 'down'" << std::endl;
      exit(1);
    }
    
    bool isData = dataArg.isSet();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    reduce(inputFiles, outputFileArg.getValue(), isData, typeArg.getValue(), maxEntriesArg.getValue(), generatorWeightArg.getValue(), puSyst, pdfSyst, jecSyst, triggerSyst, leptonSyst, btagSyst, mvaArg.getValue(), skimArg.getValue()); 

  } catch (TCLAP::ArgException& e) {
    std::cout << e.what() << std::endl;
  }

}
