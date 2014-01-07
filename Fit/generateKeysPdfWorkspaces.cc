#include <iostream>
#include <sys/wait.h>
#include <fstream>
#include <cmath>
#include <memory>
#include <chrono>
#include <unistd.h>
#include <fcntl.h>

#include <tclap/CmdLine.h>
#include <json/json.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <TCanvas.h>
#include <TString.h>
#include <TH1.h>
#include <TFile.h>
#include <TF1.h>
#include <TRandom.h>
#include <TRandom2.h>
#include <TParameter.h>
#include <TChain.h>

#include <RooWorkspace.h>
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooKeysPdf.h>
#include <RooMsgService.h>
#include <Roo1DTable.h>
#include <RooRandom.h>
#include <RooBinning.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>

#include "Utils.h"

#define NUM_ITER 150

std::string base_path = "";

void process(int mass, bool muonsOnly, int btag, const std::string& file, bool singleFile, const std::string& signalDatasetFile, __attribute__((unused)) bool fixBackground) {
  
  RooRandom::randomGenerator()->SetSeed(0);

  std::vector<pid_t> children; // For fork()

  std::string analysisName = getAnalysisName();
  base_path = "analysis/" + getAnalysisUUID();

  std::stringstream ss;
  ss << mass;
  std::string massStr = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  std::map<int, RooKeysPdf*> muon_pdf_rookeys;
  std::map<int, RooKeysPdf*> electron_pdf_rookeys;
  std::map<int, std::shared_ptr<TFile>> workspace_files;

  int btagMin = 0, btagMax = 0;
  if (btag <= 2) {
    btagMin = btagMax = btag;
  } else if (btag == 3) {
    btagMin = 1;
    btagMax = 2;
  } else {
    std::cerr << "ERROR: Unsupported number of b-tag" << std::endl;
    exit(1);
  }

  for (int i = btagMin; i <= btagMax; i++) {
    TString workspace_file = TString::Format("%s/frit/nominal-Zprime%d_%s_%d_btag_workspace.root", base_path.c_str(), mass, analysisName.c_str(), i);
    workspace_files[i].reset(TFile::Open(workspace_file, "read"));

    RooWorkspace* workspace = static_cast<RooWorkspace*>(workspace_files[i]->Get("w"));
    if (! workspace) {
      std::cerr << "ERROR: Workspace not found!" << std::endl;
      exit(1);
    }
    workspace->SetName(TString::Format("old_w_%d", i));

    muon_pdf_rookeys[i] = dynamic_cast<RooKeysPdf*>(workspace->pdf("signal_muon"));
    muon_pdf_rookeys[i]->SetName(TString::Format("old_signal_muon_%d", i));

    electron_pdf_rookeys[i] = dynamic_cast<RooKeysPdf*>(workspace->pdf("signal_electron"));
    electron_pdf_rookeys[i]->SetName(TString::Format("old_signal_electron_%d", i));
  }


  RooRealVar mtt("mtt", "mtt", 500, 2000, "GeV/c^2");
  RooRealVar weight("weight", "weight", 0, 100000);

  RooCategory lepton_type("lepton_type", "lepton_type");
  lepton_type.defineType("muon", 13);
  if (! muonsOnly) {
    lepton_type.defineType("electron", 11);
  }

  // Bin dataset
  int nBins = 1500 / 4;
  /*
  RooBinning binning(nBins, 500, 2000);
  mtt.setBinning(binning);
  */
  mtt.setRange(500, 2000);

  // Set binning to 4 GeV
  mtt.setBins(nBins);

  TParameter<int> n("iterations", NUM_ITER);

  // Load signal dataset
  std::map<int, std::shared_ptr<TChain>> chains;
  std::map<int, std::shared_ptr<RooDataSet>> signalDatasets;
  std::map<int, int> nMuons;
  std::map<int, int> nElectrons;
  std::map<int, std::shared_ptr<RooDataHist>> binnedSignalDataset_muons;
  std::map<int, std::shared_ptr<RooDataHist>> binnedSignalDataset_electrons;

  std::map<int, std::shared_ptr<RooHistPdf>> muon_pdf;
  std::map<int, std::shared_ptr<RooHistPdf>> electron_pdf;
  std::map<int, std::shared_ptr<RooAbsPdf::GenSpec>> genSpecs_muon;
  std::map<int, std::shared_ptr<RooAbsPdf::GenSpec>> genSpecs_electron;

  TString keysFile = TString::Format("%s/keyspdf_systematics_%d_%d_btag_pdf.root", base_path.c_str(), mass, btag);
  RooWorkspace keys_workspace("w", "Frit signal workspace");
  keys_workspace.import(mtt);
  keys_workspace.import(n);

  for (int i = btagMin; i <= btagMax; i++) {
    TString treeName = TString::Format("dataset_%dbtag", i);
    chains[i].reset(new TChain(treeName));
    chains[i]->Add(signalDatasetFile.c_str());

    signalDatasets[i].reset(new RooDataSet(TString::Format("dataset_%d", i), "dataset", RooArgSet(mtt, lepton_type, weight), RooFit::Import(*chains[i]), RooFit::WeightVar(weight)));
    // Keep 7000 events max
    //signalDatasets[i].reset(static_cast<RooDataSet*>(signalDatasets[i]->reduce(RooFit::EventRange(0, 7000), RooFit::SelectVars(RooArgSet(mtt, lepton_type)))));

    Roo1DTable* table = signalDatasets[i]->table(lepton_type);
    nMuons[i] = table->get("muon");
    nElectrons[i] = table->get("electron");
    delete table;

    std::cout << i << " b-tag: " << nMuons[i] << " muon events, " << nElectrons[i] << " electron events" << std::endl;

    RooDataHist* binnedSignalDataset = signalDatasets[i]->binnedClone();

    binnedSignalDataset_muons[i].reset(static_cast<RooDataHist*>(binnedSignalDataset->reduce(RooFit::SelectVars(RooArgSet(mtt)), RooFit::Cut("lepton_type == 13"), RooFit::Name(TString::Format("binned_dataset_muon_%d", i)))));
    binnedSignalDataset_electrons[i].reset(static_cast<RooDataHist*>(binnedSignalDataset->reduce(RooFit::SelectVars(RooArgSet(mtt)), RooFit::Cut("lepton_type == 11"), RooFit::Name(TString::Format("binned_dataset_electron_%d", i)))));

    delete binnedSignalDataset;

    // Create pdf from dataset
    muon_pdf[i].reset(new RooHistPdf(TString::Format("muon_hist_pdf_%d", i), "muon_hist_pdf", RooArgSet(mtt), *binnedSignalDataset_muons[i]));
    electron_pdf[i].reset(new RooHistPdf(TString::Format("electron_hist_pdf_%d", i), "electron_hist_pdf", RooArgSet(mtt), *binnedSignalDataset_electrons[i]));

    genSpecs_muon[i].reset(muon_pdf[i]->prepareMultiGen(RooArgSet(mtt), RooFit::NumEvents(nMuons[i]), RooFit::Extended(true), RooFit::AutoBinned(false)/*, RooFit::Verbose(true)*/));

    if (! muonsOnly) {
      genSpecs_electron[i].reset(electron_pdf[i]->prepareMultiGen(RooArgSet(mtt), RooFit::NumEvents(nElectrons[i]), RooFit::Extended(true), RooFit::AutoBinned(false)));
    }

    keys_workspace.import(*muon_pdf[i]);
    keys_workspace.import(*electron_pdf[i]);
    keys_workspace.import(*muon_pdf_rookeys[i]);
    keys_workspace.import(*electron_pdf_rookeys[i]);

  }

  for (int i = 0; i < NUM_ITER; i++) {

    std::cout << "Iteration #" << i + 1 << " over " << NUM_ITER << std::endl;


    std::map<int, std::shared_ptr<RooDataSet>> muon_toyData;
    std::map<int, std::shared_ptr<RooDataSet>> electron_toyData;

    std::map<int, std::shared_ptr<RooKeysPdf>> muon_toy_pdf;
    std::map<int, std::shared_ptr<RooKeysPdf>> electron_toy_pdf;

    // Create RooKeysPdf from toys
    for (int b = btagMin; b <= btagMax; b++) {

      std::cout << "Generating distribution for " << b << " b-tag..." << std::endl;
      muon_toyData[b].reset(muon_pdf[b]->generate(*genSpecs_muon[b]));
      if (!muonsOnly)
        electron_toyData[b].reset(electron_pdf[b]->generate(*genSpecs_electron[b]));
      
      std::cout << "Generating keys pdf for " << b << " b-tag ..." << std::endl;
      std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
      muon_toy_pdf[b].reset(new RooKeysPdf("signal_muon", "Keys pdf for signal", mtt, *muon_toyData[b], RooKeysPdf::MirrorBoth, 2));
      if (!muonsOnly) {
        electron_toy_pdf[b].reset(new RooKeysPdf("signal_electron", "Keys pdf for signal", mtt, *electron_toyData[b], RooKeysPdf::MirrorBoth, 2));
      }
      std::chrono::seconds t = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t0);
      std::cout << "Done in " << t.count() << "s." << std::endl;

      // Save the pdf in a workspace

      TString workspaceFile = TString::Format("%s/frit/temporary_Zprime%d_%s_%d_btag_workspace.root", base_path.c_str(), mass, analysisName.c_str(), b);
      RooWorkspace temp_workspace("w", "Frit signal workspace");

      temp_workspace.import(*muon_toy_pdf[b]);
      if (! muonsOnly)
        temp_workspace.import(*electron_toy_pdf[b]);

      temp_workspace.writeToFile(workspaceFile, true);

      TString datasetName = TString::Format("dataset_muon_%d_%d", b, i);
      muon_toyData[b]->SetName(datasetName);

      datasetName = TString::Format("dataset_electron_%d_%d", b, i);
      electron_toyData[b]->SetName(datasetName);

      keys_workspace.import(*muon_toyData[b]);
      keys_workspace.import(*electron_toyData[b]);
    }

    TString workspaceFile = TString::Format("%s/frit/temporary_Zprime%d_%s_%%d_btag_workspace.root", base_path.c_str(), mass, analysisName.c_str());
    TString outputCombineWorkspaceFile = TString::Format("%s/combine/zprime_%d_workspace_%d.root", base_path.c_str(), mass, i);
    std::cout << "Fitting..." << std::endl;
    pid_t child = fork();
    if (child == 0) {

      int fd = open("/dev/null", O_WRONLY);
      dup2(fd, STDOUT_FILENO); 
      dup2(fd, STDIN_FILENO);

      const std::string parameter = singleFile ? "-i" : "--input-list";
      if (execl("./createCombineWorkspace", "createCombineWorkspace", "-m", massStr.c_str(), parameter.c_str(), file.c_str(), "--b-tag", btagStr.c_str(), "-o", outputCombineWorkspaceFile.Data(), "--workspace", workspaceFile.Data(), NULL) < 0) {
        perror("Can't execute createCombineWorkspace");
      }
      exit(0);
    } else {
      waitpid(child, NULL, 0);
    }

    std::cout << "Done." << std::endl;

    for (int b = btagMin; b <= btagMax; b++) {
      TString pdfName = TString::Format("signal_muon_%d_%d", b, i);
      muon_toy_pdf[b]->SetName(pdfName);

      pdfName = TString::Format("signal_electron_%d_%d", b, i);
      electron_toy_pdf[b]->SetName(pdfName);

      keys_workspace.import(*muon_toy_pdf[b]);
      keys_workspace.import(*electron_toy_pdf[b]);
    }
  }

  keys_workspace.writeToFile(keysFile); 
}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Compute RooKeysPDF systematics", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> signalDatasetFileArg("", "signal", "The signal dataset input file", true, "", "string", cmd);
    TCLAP::SwitchArg muonsOnlyArg("", "muons-only", "Compute sigmaref using only semi-mu data", cmd);
    TCLAP::ValueArg<int> massArg("m", "mass", "Zprime mass", true, 0, "integer", cmd);
    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jets", false, 2, "integer", cmd);

    TCLAP::SwitchArg fixBackgroundArg("", "fix-background", "Fix background when fitting", cmd);

    cmd.parse(argc, argv);

    int mass = massArg.getValue();

    RooMsgService::instance().setStreamStatus(0, false);
    RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().addStream(RooFit::ERROR); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
    RooMsgService::instance().setSilentMode(true);

    process(mass, muonsOnlyArg.getValue(), btagArg.getValue(), inputFileArg.isSet() ? inputFileArg.getValue() : inputListArg.getValue(), inputFileArg.isSet(), signalDatasetFileArg.getValue(), fixBackgroundArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
