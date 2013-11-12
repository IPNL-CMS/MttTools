#include <TH2.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <TProfile.h>
#include <fstream>
#include <TLorentzVector.h>

#include "../PUReweighting/PUReweighter.h"
#include "tclap/CmdLine.h"

// Load libdpm at startup, on order to be sure that rfio files are working
#include <dlfcn.h>
struct Dummy
{
  Dummy()
  {
    dlopen("libdpm.so", RTLD_NOW|RTLD_GLOBAL);
  }
};
static Dummy foo;

void createBTaggingEffHistos(const std::vector<std::string>& files, const std::string& outputFile) {
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2(true);

  TH2* h_n_reco_bjets = nullptr;
  TH2* h_n_gen_bjets = nullptr;

  TH2* h_n_reco_bjets_in_cjets = nullptr;
  TH2* h_n_gen_cjets = nullptr;

  TH2* h_n_reco_bjets_in_lightjets = nullptr;
  TH2* h_n_gen_lightjets = nullptr;

  for (const std::string& file: files) {
    std::cout << "Opening '" << file << "'" << std::endl;
    TFile* f = TFile::Open(file.c_str());

    TH2* n_reco_bjets = static_cast<TH2*>(f->Get("number_of_reco_bjets"));
    n_reco_bjets->Sumw2();
    TH2* n_gen_bjets = static_cast<TH2*>(f->Get("number_of_gen_bjets"));
    n_gen_bjets->Sumw2();

    TH2* n_reco_bjets_in_cjets = static_cast<TH2*>(f->Get("number_of_reco_fake_bjets_among_cjets"));
    n_reco_bjets_in_cjets->Sumw2();
    TH2* n_gen_cjets = static_cast<TH2*>(f->Get("number_of_gen_cjets"));
    n_gen_cjets->Sumw2();

    TH2* n_reco_bjets_in_lightjets = static_cast<TH2*>(f->Get("number_of_reco_fake_bjets_among_lightjets"));
    n_reco_bjets_in_lightjets->Sumw2();
    TH2* n_gen_lightjets = static_cast<TH2*>(f->Get("number_of_gen_lightjets"));
    n_gen_lightjets->Sumw2();

    if (! h_n_reco_bjets) {
      h_n_reco_bjets = (TH2*) n_reco_bjets->Clone();
      h_n_reco_bjets_in_cjets = (TH2*) n_reco_bjets_in_cjets->Clone();
      h_n_reco_bjets_in_lightjets = (TH2*) n_reco_bjets_in_lightjets->Clone();

      h_n_gen_bjets = (TH2*) n_gen_bjets->Clone();
      h_n_gen_cjets = (TH2*) n_gen_cjets->Clone();
      h_n_gen_lightjets = (TH2*) n_gen_lightjets->Clone();

    } else {
      h_n_reco_bjets->Add(n_reco_bjets);
      h_n_reco_bjets_in_cjets->Add(n_reco_bjets_in_cjets);
      h_n_reco_bjets_in_lightjets->Add(n_reco_bjets_in_lightjets);

      h_n_gen_bjets->Add(n_gen_bjets);
      h_n_gen_cjets->Add(n_gen_cjets);
      h_n_gen_lightjets->Add(n_gen_lightjets);
    }

    f->Close();
    delete f;
  }

  TH2* btagging_efficiency = (TH2*) h_n_reco_bjets->Clone();
  btagging_efficiency->SetName("btagging_efficiency");

  TH2* c_fakerate = (TH2*) h_n_reco_bjets_in_cjets->Clone();
  c_fakerate->SetName("cjets_fakerate");

  TH2* light_fakerate = (TH2*) h_n_reco_bjets_in_lightjets->Clone();
  light_fakerate->SetName("lightjets_fakerate");

  btagging_efficiency->Divide(h_n_gen_bjets);
  c_fakerate->Divide(h_n_gen_cjets);
  light_fakerate->Divide(h_n_gen_lightjets);

  TFile* f = TFile::Open(outputFile.c_str(), "recreate");
  btagging_efficiency->Write();
  c_fakerate->Write();
  light_fakerate->Write();
  f->Close();
  delete f;
}

void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

int main(int argc, char** argv) {

  try {

    TCLAP::CmdLine cmd("Convert extractor tuples to histograms", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", true, "", "string", cmd);

    //TCLAP::SwitchArg dataArg("", "data", "Is this data?", false);
    //TCLAP::SwitchArg mcArg("", "mc", "Is this mc?", false);

    //cmd.xorAdd(dataArg, mcArg);

    //TCLAP::SwitchArg semimuArg("", "semimu", "Is this semi-mu channel?", false);
    //TCLAP::SwitchArg semieArg("", "semie", "Is this semi-e channel?", false);

    //cmd.xorAdd(semimuArg, semieArg);

    //TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jet to require", true, 2, "int", cmd);

    //TCLAP::ValueArg<std::string> pileupArg("", "pileup", "PU profile used for MC production", false, "S10", "string", cmd);

    //TCLAP::ValueArg<std::string> pileupSystArg("", "pileup-syst", "PU profile to use for pileup reweigthing", false, "nominal", "string", cmd);
    //TCLAP::ValueArg<std::string> pdfSystArg("", "pdf-syst", "PDF systematic to compute", false, "nominal", "string", cmd);


    //TCLAP::ValueArg<int> maxEntriesArg("n", "", "Maximal number of entries to process", false, -1, "int", cmd);
    //TCLAP::ValueArg<double> generatorWeightArg("", "weight", "MC generator weight", false, 1., "double", cmd);

    cmd.parse(argc, argv);

    //PUProfile puProfile;
    //std::string p = pileupArg.getValue();
    //std::transform(p.begin(), p.end(), p.begin(), ::tolower);
    //if (p == "s6")
      //puProfile = PUProfile::S6;
    //else if (p == "s7")
      //puProfile = PUProfile::S7;
    //else if (p == "s10")
      //puProfile = PUProfile::S10;

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
    
    //bool isData = dataArg.isSet();

    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }
    
    //Extractor2Histos convertor(inputFiles, outputFileArg.getValue(), semimuArg.isSet(), !isData, btagArg.getValue());
    //convertor.Loop();
    createBTaggingEffHistos(inputFiles, outputFileArg.getValue());

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
