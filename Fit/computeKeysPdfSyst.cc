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

#define NUM_ITER 500

std::string base_path = "";

double getSigmaReference(int mass, int btag) {

  Json::Reader reader;
  Json::Value refRoot;
  std::ifstream file((base_path + "/sigma_reference.json").c_str());
  bool success = reader.parse(file, refRoot);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "sigma_reference.json" << "'. Exiting." << std::endl;
    std::cerr << "You may need to run fitMtt first." << std::endl;
    exit(1);
  }

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  refRoot = refRoot[getAnalysisUUID()];

  if (!refRoot.isMember(strMass)) {
    std::cerr << "ERROR: mass '" << mass << "' not found in JSON file. Exiting." << std::endl;
    exit(1);
  }

  double sigma_ref = refRoot[strMass][btagStr]["sigma"].asDouble();

  return sigma_ref;
}

double getNumberOfEventsReference(int mass, int btag) {

  Json::Reader reader;
  Json::Value refRoot;
  std::ifstream file((base_path + "/sigma_reference.json").c_str());
  bool success = reader.parse(file, refRoot);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse '" << "sigma_reference.json" << "'. Exiting." << std::endl;
    std::cerr << "You may need to run fitMtt first." << std::endl;
    exit(1);
  }

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  refRoot = refRoot[getAnalysisUUID()];

  if (!refRoot.isMember(strMass)) {
    std::cerr << "ERROR: mass '" << mass << "' not found in JSON file. Exiting." << std::endl;
    exit(1);
  }

  return refRoot[strMass][btagStr]["events"].asDouble();
}

void process(int mass, bool muonsOnly, int btag, const std::string& file, bool singleFile, const std::string& signalDatasetFile) {
  
  RooRandom::randomGenerator()->SetSeed(0);

  std::vector<pid_t> children; // For fork()
  int shmid;

  std::string analysisName = getAnalysisName();
  base_path = "analysis/" + getAnalysisUUID();

  TRandom* random = new TRandom2(0);

  // Create a shared memory area to store the cross sections
  key_t key = (key_t) random->Uniform(1024, 65535);
  if ((shmid = shmget(key, sizeof(SHMFitResults) * 8, IPC_CREAT | 0666)) < 0) {
    perror("Can't create shared memory area");
    exit(1);
  }

  std::stringstream ss;
  ss << mass;

  std::stringstream strKey;
  strKey << key;

  SHMFitResults* shm = NULL;
  if ((shm = static_cast<SHMFitResults*>(shmat(shmid, NULL, 0))) == (void *) -1) {
    perror("Can't map shared memory to local memory");
    exit(1);
  }

  double events_reference = getNumberOfEventsReference(mass, btag);

  /*
     pid_t child = fork();
     if (child == 0) {

     int fd = open("/dev/null", O_WRONLY);
     dup2(fd, STDOUT_FILENO); 
     dup2(fd, STDIN_FILENO);

     char** params = getSystCLParameters(ss.str(), file, singleFile, muonsOnly, btag, "--shared-memory", "--shm-key", strKey.str().c_str(), NULL);
     execv("./fitMtt", params);
     exit(0);
     } else {
     waitpid(child, NULL, 0);
     }

     std::cout << "Sigma ref: JSON = " << sigma_reference << "; C++ = " << shm->sigma << std::endl;
     */

  TString filename = TString::Format("keyspdf_systematics_%d_%d_btag.root", mass, btag);
  TH1* sigma = new TH1D("sigma", "sigma", 80, -3, 3);
  TH1* pull = new TH1D("pull", "pull", 80, -3, 3);
  TH1* events = new TH1D("events", "events", 80, events_reference - 100, events_reference + 100);
  TH1* residuals = new TH1D("residuals", "residuals", 80, -200, 200);

  TString workspace_file = TString::Format("%s/frit/nominal-Zprime%d_%s_%d_btag_workspace.root", base_path.c_str(), mass, analysisName.c_str(), btag);
  TFile *w = TFile::Open(workspace_file, "read");

  RooWorkspace* workspace = static_cast<RooWorkspace*>(w->Get("w"));
  if (! workspace) {
    std::cerr << "ERROR: Workspace not found!" << std::endl;
    exit(1);
  }
  workspace->SetName("old_w");

  RooKeysPdf* muon_pdf_rookeys = dynamic_cast<RooKeysPdf*>(workspace->pdf("signal_muon"));
  muon_pdf_rookeys->SetName("old_signal_muon");
  RooKeysPdf* electron_pdf_rookeys = dynamic_cast<RooKeysPdf*>(workspace->pdf("signal_electron"));
  electron_pdf_rookeys->SetName("old_signal_electron");

  RooRealVar mtt("mtt", "mtt", 500, 2000, "GeV/c^2");
  //RooRealVar weight("weight", "weight", 0, 100000);

  RooCategory lepton_type("lepton_type", "lepton_type");
  lepton_type.defineType("muon", 13);
  if (! muonsOnly) {
    lepton_type.defineType("electron", 11);
  }

  // Bin dataset
  //int nBins = 1500 / 5;
  int nBins = 50;
  /*
  RooBinning binning(nBins, 500, 2000);
  mtt.setBinning(binning);
  */
  mtt.setRange(500, 2000);

  // Set binning to 5 GeV
  mtt.setBins(nBins);


  // Load signal dataset
  TString treeName = TString::Format("dataset_%dbtag", btag);
  TChain* chain = new TChain(treeName);
  chain->Add(signalDatasetFile.c_str());

  RooDataSet signalDataset("dataset", "dataset", RooArgSet(mtt, lepton_type), RooFit::Import(*chain));

  Roo1DTable* table = signalDataset.table(lepton_type);
  int nMuons = table->get("muon");
  int nElectrons = table->get("electron");

  RooDataHist* binnedSignalDataset = signalDataset.binnedClone();

  RooDataHist* binnedSignalDataset_muons = static_cast<RooDataHist*>(binnedSignalDataset->reduce(RooFit::SelectVars(RooArgSet(mtt)), RooFit::Cut("lepton_type == 13"), RooFit::Name("binned_dataset_muon")));
  RooDataHist* binnedSignalDataset_electrons = static_cast<RooDataHist*>(binnedSignalDataset->reduce(RooFit::SelectVars(RooArgSet(mtt)), RooFit::Cut("lepton_type == 11"), RooFit::Name("binned_dataset_electron")));

  delete binnedSignalDataset;

  // Keys pdf from dataset
  //RooKeysPdf muon_pdf_binned_keys("muon_binned_keys_pdf", "muon_binned_keys_pdf", mtt, *binnedSignalDataset_muons, RooKeysPdf::MirrorBoth, 2);
  //RooKeysPdf electron_pdf_binned_keys("electron_binned_keys_pdf", "electron_binned_keys_pdf", mtt, *binnedSignalDataset_electrons, RooKeysPdf::MirrorBoth, 2);

  // Create pdf from dataset
  RooHistPdf muon_pdf("muon_hist_pdf", "muon_hist_pdf", RooArgSet(mtt), *binnedSignalDataset_muons);
  RooHistPdf electron_pdf("electron_hist_pdf", "electron_hist_pdf", RooArgSet(mtt), *binnedSignalDataset_electrons);

  std::shared_ptr<RooAbsPdf::GenSpec> genSpecs_muon = std::shared_ptr<RooAbsPdf::GenSpec>(
      muon_pdf.prepareMultiGen(RooArgSet(mtt), RooFit::NumEvents(nMuons), RooFit::Extended(true), RooFit::AutoBinned(false)/*, RooFit::Verbose(true)*/)
      );

  std::shared_ptr<RooAbsPdf::GenSpec> genSpecs_electron;
  if (! muonsOnly) {
    genSpecs_electron = std::shared_ptr<RooAbsPdf::GenSpec>(
        electron_pdf.prepareMultiGen(RooArgSet(mtt), RooFit::NumEvents(nElectrons), RooFit::Extended(true), RooFit::AutoBinned(false))
        );
  }

  TString keysFile = TString::Format("%s/keyspdf_systematics_%d_%d_btag_pdf.root", base_path.c_str(), mass, btag);
  RooWorkspace keys_workspace("w", "Frit signal workspace");

  keys_workspace.import(mtt);
  keys_workspace.import(muon_pdf);
  keys_workspace.import(electron_pdf);
  keys_workspace.import(*muon_pdf_rookeys);
  keys_workspace.import(*electron_pdf_rookeys);

  //keys_workspace.import(muon_pdf_binned_keys);
  //keys_workspace.import(electron_pdf_binned_keys);

  //keys_workspace.import(*binnedSignalDataset_muons);
  //keys_workspace.import(*binnedSignalDataset_electrons);

  TParameter<int> n("iterations", NUM_ITER);
  keys_workspace.import(n);

  for (int i = 0; i < NUM_ITER; i++) {

    std::cout << "Iteration #" << i + 1 << " over " << NUM_ITER << std::endl;

    std::cout << "Generating distribution..." << std::endl;
    RooDataSet* muon_toyData = nullptr;
    RooDataSet* electron_toyData = nullptr;

    muon_toyData = muon_pdf.generate(*genSpecs_muon);
    if (!muonsOnly)
      electron_toyData = electron_pdf.generate(*genSpecs_electron);

    // Create RooKeysPdf from toys
    std::cout << "Generating keys pdf..." << std::endl;
    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();

    RooKeysPdf* muon_toy_pdf = new RooKeysPdf("signal_muon", "Keys pdf for signal", mtt, *muon_toyData, RooKeysPdf::MirrorBoth, 2);
    RooKeysPdf* electron_toy_pdf = nullptr;
    if (!muonsOnly) {
      electron_toy_pdf = new RooKeysPdf("signal_electron", "Keys pdf for signal", mtt, *electron_toyData, RooKeysPdf::MirrorBoth, 2);
    }

    std::chrono::seconds t = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t0);
    std::cout << "Done in " << t.count() << "s." << std::endl;

    TString datasetName = TString::Format("dataset_muon_%d", i);
    muon_toyData->SetName(datasetName);

    datasetName = TString::Format("dataset_electron_%d", i);
    electron_toyData->SetName(datasetName);

    keys_workspace.import(*muon_toyData);
    keys_workspace.import(*electron_toyData);

    delete muon_toyData;
    delete electron_toyData;

    // Save the pdf in a workspace

    TString workspaceFile = TString::Format("%s/frit/temporary_Zprime%d_%s_%d_btag_workspace.root", base_path.c_str(), mass, analysisName.c_str(), btag);
    RooWorkspace temp_workspace("w", "Frit signal workspace");

    temp_workspace.import(*muon_toy_pdf);
    if (! muonsOnly)
      temp_workspace.import(*electron_toy_pdf);

    temp_workspace.writeToFile(workspaceFile, true);

    std::cout << "Fitting..." << std::endl;

    pid_t child = fork();
    if (child == 0) {

      int fd = open("/dev/null", O_WRONLY);
      dup2(fd, STDOUT_FILENO); 
      dup2(fd, STDIN_FILENO);

      char** params = getSystCLParameters(ss.str(), file, singleFile, muonsOnly, btag, "--shared-memory", "--shm-key", strKey.str().c_str(), "--workspace", workspaceFile.Data(), NULL);
      if (execv("./fitMtt", params) < 0) {
        perror("Can't execute fitMtt");
      }
      exit(0);
    } else {
      waitpid(child, NULL, 0);
    }

    std::cout << "Done. # events: " << shm->nSignalEvents << "; Cross-section: " << shm->sigma << " pb." << std::endl << std::endl;

    sigma->Fill(shm->sigma);
    events->Fill(shm->nSignalEvents);
    pull->Fill((shm->nSignalEvents - events_reference) / shm->nSignalEvents_error);
    residuals->Fill(shm->nSignalEvents - events_reference);

    TString pdfName = TString::Format("signal_muon_%d", i);
    muon_toy_pdf->SetName(pdfName);

    pdfName = TString::Format("signal_electron_%d", i);
    electron_toy_pdf->SetName(pdfName);

    keys_workspace.import(*muon_toy_pdf);
    keys_workspace.import(*electron_toy_pdf);

    delete muon_toy_pdf;
    delete electron_toy_pdf;
  }

  keys_workspace.writeToFile(keysFile); 

  TFile* f = TFile::Open(filename, "recreate");
  sigma->Write();
  pull->Write();
  events->Write();
  residuals->Write();
  f->Close();

  delete f;

  delete sigma;
  delete pull;

  delete w;

  delete binnedSignalDataset_muons;
  delete binnedSignalDataset_electrons;

  shmdt(shm);

  shmctl(shmid, IPC_RMID, NULL);
}

void saveSystematic(int mass, int btag, double syst) {

  FILE* lock = fopen((base_path + "/systematics.lock").c_str(), "w+");
  lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  Json::Reader reader;
  Json::Value root;
  std::ifstream file((base_path + "/systematics.json").c_str());
  reader.parse(file, root);
  file.close();

  std::stringstream ss;
  ss << mass;
  std::string strMass = ss.str();

  ss.clear(); ss.str(std::string());
  ss << btag;
  std::string btagStr = ss.str();

  root[getAnalysisUUID()][strMass][btagStr]["signal_pdf"] = syst;

  FILE* fd = fopen((base_path + "/systematics.json").c_str(), "w+");
  Json::StyledWriter writer;
  const std::string json = writer.write(root);
  fwrite(json.c_str(), json.length(), 1, fd);
  fclose(fd);

  fclose(lock);
}

void computeSyst(int mass, int btag) {

  double sigma_ref = getSigmaReference(mass, btag);

  std::cout << std::endl << std::endl;

  TString filename = TString::Format("%s/keyspdf_systematics_%d_%d_btag.root", base_path.c_str(), mass, btag);
  TFile* f = TFile::Open(filename, "read");

  TH1* sigma = static_cast<TH1*>(f->Get("sigma"));
  TF1* gaussian = new TF1("g", "gaus", -3, 3);
  sigma->Fit(gaussian, "Q");

  double syst = gaussian->GetParameter(2) / sigma_ref;
  saveSystematic(mass, btag, syst);

  TH1* pull = static_cast<TH1*>(f->Get("pull"));
  pull->Fit("gaus", "Q");

  TH1* events = static_cast<TH1*>(f->Get("events"));
  events->Fit("gaus", "Q");

  TH1* residuals = static_cast<TH1*>(f->Get("residuals"));
  residuals->Fit("gaus", "Q");

  filename = TString::Format("%s/keyspdf_systematics_%d_%d_btag_fits.root", base_path.c_str(), mass, btag);
  TFile* output = TFile::Open(filename, "recreate");
  sigma->Write();
  pull->Write();
  events->Write();
  residuals->Write();

  output->Close();
  delete output;

  f->Close();
  delete f;

  std::cout << "Signal PDF syst for MZ' = " << mass << " : " << syst << std::endl;
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
    TCLAP::SwitchArg extractArg("", "dont-extract", "Only compute systematic with previous results, don't run fitMtt", cmd);

    cmd.parse(argc, argv);

    int mass = massArg.getValue();

    RooMsgService::instance().setStreamStatus(0, false);
    RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().addStream(RooFit::ERROR); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
    RooMsgService::instance().setSilentMode(true);

    if (! extractArg.getValue())
      process(mass, muonsOnlyArg.getValue(), btagArg.getValue(), inputFileArg.isSet() ? inputFileArg.getValue() : inputListArg.getValue(), inputFileArg.isSet(), signalDatasetFileArg.getValue());

    computeSyst(mass, btagArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
