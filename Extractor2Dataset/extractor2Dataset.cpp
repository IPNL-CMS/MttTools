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

#include "../PUReweighting/PUReweighter.h"

#include "TopTriggerEfficiencyProvider.h"

bool OVERRIDE_TYPE;
std::string OVERRIDED_TYPE;
PUProfile puProfile;

void reduce(TChain* mtt, TChain* event, TChain* vertices, const std::string& outputFile, bool isData, const std::string& type, int max, double lumi_weight, const std::string& puSyst, const std::string& pdfSyst) {

  TString outputFileFormated = OVERRIDE_TYPE ? outputFile : TString::Format(outputFile.c_str(), (isData) ? "" : type.c_str());

  std::map<int, TTree*> outputTrees {
    {0, new TTree("dataset_0btag", "dataset for 0 b-tagged jet") },
      {1, new TTree("dataset_1btag", "dataset for 1 b-tagged jet") },
      {2, new TTree("dataset_2btag", "dataset for at least 2 b-tagged jets") }
  };

  float mtt_afterChi2, pt_1stJet, pt_2ndJet, pt_3rdJet, pt_4thJet, bestSolChi2;
  int isSel, nBtaggedJets_CSVM;
  uint32_t run = 0;

  float lepton_weight = 1;
  float btag_weight = 1;

  mtt->SetBranchStatus("*", 0);
  mtt->SetBranchStatus("mtt_AfterChi2", 1);
  mtt->SetBranchStatus("1stjetpt", 1);
  mtt->SetBranchStatus("2ndjetpt", 1);
  mtt->SetBranchStatus("3rdjetpt", 1);
  mtt->SetBranchStatus("4thjetpt", 1);
  mtt->SetBranchStatus("bestSolChi2", 1);
  mtt->SetBranchStatus("isSel", 1);
  if (puProfile == PUProfile::S6)
    mtt->SetBranchStatus("nBtaggedJets_TCHET", 1);
  else
    mtt->SetBranchStatus("nBtaggedJets_CSVM", 1);
  if (! isData) {
    mtt->SetBranchStatus("lepton_weight", 1);
    mtt->SetBranchStatus("btag_weight", 1);
  }

  mtt->SetBranchAddress("mtt_AfterChi2", &mtt_afterChi2, NULL);
  mtt->SetBranchAddress("1stjetpt", &pt_1stJet, NULL);
  mtt->SetBranchAddress("2ndjetpt", &pt_2ndJet, NULL);
  mtt->SetBranchAddress("3rdjetpt", &pt_3rdJet, NULL);
  mtt->SetBranchAddress("4thjetpt", &pt_4thJet, NULL);
  mtt->SetBranchAddress("bestSolChi2", &bestSolChi2, NULL);
  mtt->SetBranchAddress("isSel", &isSel, NULL);
  if (puProfile == PUProfile::S6)
    mtt->SetBranchAddress("nBtaggedJets_TCHET", &nBtaggedJets_CSVM, NULL);
  else
    mtt->SetBranchAddress("nBtaggedJets_CSVM", &nBtaggedJets_CSVM, NULL);
  if (!isData) {
    mtt->SetBranchAddress("lepton_weight", &lepton_weight);
    mtt->SetBranchAddress("btag_weight", &btag_weight);
  }


  float muonPt[100];
  float muonEta[100];
  float electronPt[100];
  float electronEta[100];
  if (type == "semimu") {
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

  int32_t n_vertices;
  vertices->SetBranchStatus("*", 0);
  vertices->SetBranchStatus("n_vertices", 1);
  vertices->SetBranchAddress("n_vertices", &n_vertices);

  float output_weight;
  int lepton = (type == "semimu") ? 13 : 11;
  for (auto& outputTree: outputTrees) {
    outputTree.second->Branch("mtt", &mtt_afterChi2, "mtt/F");
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

    puReweigher = new PUReweighter(type == "semimu", puProfile, syst);
  }

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

  if (type == "semimu") {
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

    float firstJetCut = 0, secondJetCut = 0, thirdJetCut = 0;
    if (isRun2012AB) {
      firstJetCut = 45;
      secondJetCut = 45;
      thirdJetCut = 45;
    } else {
      firstJetCut = 55;
      secondJetCut = 45;
      thirdJetCut = 35;
    }

    // Lepton selection done on Extractor
    // Selection
    if (isSel == 1 && pt_1stJet > firstJetCut && pt_2ndJet > secondJetCut && pt_3rdJet > thirdJetCut && mtt_afterChi2 > 0. && bestSolChi2 < 500.) {
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

      double ptLepton = 0;
      double etaLepton = 0;
      if (type == "semimu") {
        ptLepton = muonPt[0];
        etaLepton = muonEta[0];
      } else {
        ptLepton = electronPt[0];
        etaLepton = electronEta[0];
        
        // The TOP reference selection exclude electron with
        // SuperCluster eta between 1.4442 and 1.5660
        // The TopTrigger efficiency does the same thing, but
        // using electron eta instead of SuperCluster eta.
        // Redo a cut here on electron eta
        // FIXME?
        if (fabs(etaLepton) >= 1.442 && fabs(etaLepton) < 1.5660)
          continue;
      }

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
        double triggerWeight = m_trigger_efficiency_provider->get_weight(ptLepton, etaLepton, pt_4thJet, jetEta[3], n_vertices, nJets, type == "semimu", TopTriggerEfficiencyProvider::NOMINAL)[0];
        output_weight *= puReweigher->weight(n_trueInteractions) * generator_weight * lumi_weight * triggerWeight * lepton_weight/* * btag_weight*/;
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
  vertices = new TChain("vertices");

  for (const std::string& file: inputFiles) {
    mtt->Add(file.c_str());
    event->Add(file.c_str());
    vertices->Add(file.c_str());
  }

  event->SetCacheSize(30*1024*1024);
  mtt->SetCacheSize(30*1024*1024);
  vertices->SetCacheSize(30*1024*1024);
}

void reduce(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isData, const std::string& type, int max, double generator_weight, const std::string& puSyst, const std::string& pdfSyst) {

  TChain* mtt = NULL, *event = NULL, *vertices = NULL;

  loadChain(inputFiles, mtt, event, vertices);
  reduce(mtt, event, vertices, outputFile, isData, type, max, generator_weight, puSyst, pdfSyst);

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
    TCLAP::ValueArg<std::string> pileupSystArg("", "pileup-syst", "PU profile to use for pileup reweigthing", false, "nominal", "string", cmd);
    TCLAP::ValueArg<int> maxEntriesArg("n", "", "Maximal number of entries to process", false, -1, "int", cmd);
    TCLAP::ValueArg<double> generatorWeightArg("", "weight", "MC generator weight", false, 1., "double", cmd);

    TCLAP::ValueArg<std::string> pdfSystArg("", "pdf-syst", "PDF systematic to compute", false, "nominal", "string", cmd);

    cmd.parse(argc, argv);

    std::string p = pileupArg.getValue();
    std::transform(p.begin(), p.end(), p.begin(), ::tolower);
    if (p == "s6")
      puProfile = PUProfile::S6;
    else if (p == "s7")
      puProfile = PUProfile::S7;
    else if (p == "s10")
      puProfile = PUProfile::S10;

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
    
    bool isData = dataArg.isSet();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    reduce(inputFiles, outputFileArg.getValue(), isData, typeArg.getValue(), maxEntriesArg.getValue(), generatorWeightArg.getValue(), puSyst, pdfSyst); 

  } catch (TCLAP::ArgException& e) {
    std::cout << e.what() << std::endl;
  }

}
