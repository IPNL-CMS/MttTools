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
#include <RooAbsPdf.h>
#include <RooFitResult.h>

#include "Functions.h"


void loadInputFiles(const std::string& filename, std::vector<std::string>& files) {

  ifstream ifs(filename.c_str());
  std::string line;

  while (getline(ifs, line))
    files.push_back(line);

  ifs.close();
}

TChain* loadChain(const std::vector<std::string>& files, const std::string& name) {
  TChain* c = new TChain(name.c_str());

  for (const std::string& file: files) {
    c->Add(file.c_str());
  }

  return c;
}

#include "Utils.h"

void process(const std::vector<std::string>& inputFiles, TFile* outputFile, bool usePDF, int btag) {

  TString treeName = TString::Format("dataset_%dbtag", btag);
  
  std::shared_ptr<TChain> chain(loadChain(inputFiles, treeName.Data()));

  double xMin = 550;
  double xMax = 2000;

  if (! usePDF) {
    xMin = 0;
    xMax = 5000;
  }

  RooRealVar mtt("mtt", "mtt", xMin, xMax, "GeV/c^2");
  RooRealVar weight("weight", "weight", 0, 10000);

  RooCategory lepton_type("lepton_type", "lepton_type");
  lepton_type.defineType("electron", 11);
  lepton_type.defineType("muon", 13);

  mtt.setBins((xMax -xMin) / 4.);

  std::shared_ptr<RooDataSet> dataset(new RooDataSet(treeName, "dataset", RooArgSet(mtt, lepton_type, weight), RooFit::Import(*chain), RooFit::WeightVar(weight)));

  // First, bin dataset
  Roo1DTable* table = dataset->table(lepton_type);
  int nMuons = table->get("muon");
  int nElectrons = table->get("electron");
  delete table;

  std::shared_ptr<RooDataHist> binnedSignalDataset(dataset->binnedClone());

  std::shared_ptr<RooDataHist> binnedSignalDataset_muons(static_cast<RooDataHist*>(binnedSignalDataset->reduce(RooFit::SelectVars(RooArgSet(mtt, weight)), RooFit::Cut("lepton_type == 13"), RooFit::Name(TString::Format("binned_dataset_muon_%d", btag)))));
  std::shared_ptr<RooDataHist> binnedSignalDataset_electrons(static_cast<RooDataHist*>(binnedSignalDataset->reduce(RooFit::SelectVars(RooArgSet(mtt, weight)), RooFit::Cut("lepton_type == 11"), RooFit::Name(TString::Format("binned_dataset_electron_%d", btag)))));
   
  // Create pdf from dataset
  std::shared_ptr<RooAbsPdf> muon_pdf;
  std::shared_ptr<RooAbsPdf> electron_pdf;

  // Used only if generating from PDF
  std::shared_ptr<BaseFunction> muon_pdf_base;
  std::shared_ptr<BaseFunction> electron_pdf_base;

  if (! usePDF) {
    // Generate using a histogram
    muon_pdf.reset(new RooHistPdf(TString::Format("muon_hist_pdf_%d", btag), "muon_hist_pdf", RooArgSet(mtt), *binnedSignalDataset_muons));
    electron_pdf.reset(new RooHistPdf(TString::Format("electron_hist_pdf_%d", btag), "electron_hist_pdf", RooArgSet(mtt), *binnedSignalDataset_electrons));
  } else {

    Json::Value pdf;
    pdf["name"] = "faltB";
    pdf["parameters"] = "faltB";

    muon_pdf_base = createPdf(pdf, "./fit_configuration", mtt, nullptr, 0, "background", "muon");
    electron_pdf_base = createPdf(pdf, "./fit_configuration", mtt, nullptr, 0, "background", "electron");
    
    muon_pdf = muon_pdf_base->getSharedPdf();
    electron_pdf = electron_pdf_base->getSharedPdf();

    RooFitResult* fitResult = muon_pdf->fitTo(*binnedSignalDataset_muons, RooFit::Save(), RooFit::Optimize(0));
    std::cout << "Muon fit" << std::endl;
    fitResult->Print("v");
    delete fitResult;

    fitResult = electron_pdf->fitTo(*binnedSignalDataset_electrons, RooFit::Save(), RooFit::Optimize(0));
    std::cout << "Electron fit" << std::endl;
    fitResult->Print("v");
    delete fitResult;
  }
  
  //genSpecs_muon[i].reset(muon_pdf[i]->prepareMultiGen(RooArgSet(mtt), RooFit::NumEvents(nMuons[i]), RooFit::Extended(true), RooFit::AutoBinned(false)[>, RooFit::Verbose(true)<]));

  //if (! muonsOnly) {
  //genSpecs_electron[i].reset(electron_pdf[i]->prepareMultiGen(RooArgSet(mtt), RooFit::NumEvents(nElectrons[i]), RooFit::Extended(true), RooFit::AutoBinned(false)));
  //}
  
  std::cout << "Generating #" << nMuons << " muons and #" << nElectrons << " electrons..." << std::endl;
  std::shared_ptr<RooDataSet> generated_dataset_muon(muon_pdf->generate(RooArgSet(mtt), RooFit::NumEvents(nMuons), RooFit::Extended(true), RooFit::AutoBinned(false)));
  std::shared_ptr<RooDataSet> generated_dataset_electron(electron_pdf->generate(RooArgSet(mtt), RooFit::NumEvents(nElectrons), RooFit::Extended(true), RooFit::AutoBinned(false)));
  std::cout << "done." << std::endl;

  std::shared_ptr<TTree> tree(new TTree(treeName, treeName));
  float _mtt = 0;
  float _weight = 1.;
  int _lepton_type = 11;

  tree->Branch("mtt", &_mtt, "mtt/F");
  tree->Branch("weight", &_weight, "weight/F");
  tree->Branch("lepton_type", &_lepton_type, "lepton_type/I");

  int electronEntries = generated_dataset_electron->sumEntries();
  for (int i = 0; i < electronEntries; i++) {

    const RooArgSet* vars = generated_dataset_electron->get(i);
    _mtt = static_cast<RooRealVar*>(vars->find("mtt"))->getVal();

    tree->Fill();
  }

  _lepton_type = 13;
  int muonEntries = generated_dataset_muon->sumEntries();
  for (int i = 0; i < muonEntries; i++) {

    const RooArgSet* vars = generated_dataset_muon->get(i);
    _mtt = static_cast<RooRealVar*>(vars->find("mtt"))->getVal();

    tree->Fill();
  }

  outputFile->cd();
  tree->Write();
}

void process(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool usePDF) {

  RooRandom::randomGenerator()->SetSeed(0);

  TFile* output = TFile::Open(outputFile.c_str(), "recreate");
  process(inputFiles, output, usePDF, 0);
  process(inputFiles, output, usePDF, 1);
  process(inputFiles, output, usePDF, 2);

  output->Close();
  delete output;
}


//void process(int mass, bool muonsOnly, int btag, const std::string& file, bool singleFile, const std::string& signalDatasetFile, bool fixBackground) {
  
  //RooRandom::randomGenerator()->SetSeed(0);

  //std::vector<pid_t> children; // For fork()
  //int shmid;

  //std::string analysisName = getAnalysisName();
  //base_path = "analysis/" + getAnalysisUUID();

  //TRandom* random = new TRandom2(0);

  //// Create a shared memory area to store the cross sections
  //key_t key = (key_t) random->Uniform(1024, 65535);
  //if ((shmid = shmget(key, sizeof(SHMFitResults) * 8, IPC_CREAT | 0666)) < 0) {
    //perror("Can't create shared memory area");
    //exit(1);
  //}
  
  //delete random;

  //std::stringstream ss;
  //ss << mass;

  //std::stringstream strKey;
  //strKey << key;

  //SHMFitResults* shm = NULL;
  //if ((shm = static_cast<SHMFitResults*>(shmat(shmid, NULL, 0))) == (void *) -1) {
    //perror("Can't map shared memory to local memory");
    //exit(1);
  //}

  //double sigma_reference = getSigmaReference(mass, btag);
  //double events_reference = getNumberOfEventsReference(mass, btag);

  //[>
     //pid_t child = fork();
     //if (child == 0) {

     //int fd = open("/dev/null", O_WRONLY);
     //dup2(fd, STDOUT_FILENO); 
     //dup2(fd, STDIN_FILENO);

     //char** params = getSystCLParameters(ss.str(), file, singleFile, muonsOnly, btag, "--shared-memory", "--shm-key", strKey.str().c_str(), NULL);
     //execv("./fitMtt", params);
     //exit(0);
     //} else {
     //waitpid(child, NULL, 0);
     //}

     //std::cout << "Sigma ref: JSON = " << sigma_reference << "; C++ = " << shm->sigma << std::endl;
     //*/

  //std::map<int, RooKeysPdf*> muon_pdf_rookeys;
  //std::map<int, RooKeysPdf*> electron_pdf_rookeys;
  //std::map<int, std::shared_ptr<TFile>> workspace_files;

  //TString filename = TString::Format("%s/keyspdf_systematics_%d_%d_btag.root", base_path.c_str(), mass, btag);
  //TH1* sigma = new TH1D("sigma", "sigma", 80, sigma_reference - 2, sigma_reference + 2);
  //TH1* pull = new TH1D("pull", "pull", 80, -3, 3);
  //TH1* events = new TH1D("events", "events", 80, events_reference - 100, events_reference + 100);
  //TH1* residuals = new TH1D("residuals", "residuals", 80, -200, 200);

  

  //int btagMin = 0, btagMax = 0;
  //if (btag <= 2) {
    //btagMin = btagMax = btag;
  //} else if (btag == 3) {
    //btagMin = 1;
    //btagMax = 2;
  //} else {
    //std::cerr << "ERROR: Unsupported number of b-tag" << std::endl;
    //exit(1);
  //}

  //for (int i = btagMin; i <= btagMax; i++) {
    //TString workspace_file = TString::Format("%s/frit/nominal-Zprime%d_%s_%d_btag_workspace.root", base_path.c_str(), mass, analysisName.c_str(), i);
    //workspace_files[i].reset(TFile::Open(workspace_file, "read"));

    //RooWorkspace* workspace = static_cast<RooWorkspace*>(workspace_files[i]->Get("w"));
    //if (! workspace) {
      //std::cerr << "ERROR: Workspace not found!" << std::endl;
      //exit(1);
    //}
    //workspace->SetName(TString::Format("old_w_%d", i));

    //muon_pdf_rookeys[i] = dynamic_cast<RooKeysPdf*>(workspace->pdf("signal_muon"));
    //muon_pdf_rookeys[i]->SetName(TString::Format("old_signal_muon_%d", i));

    //electron_pdf_rookeys[i] = dynamic_cast<RooKeysPdf*>(workspace->pdf("signal_electron"));
    //electron_pdf_rookeys[i]->SetName(TString::Format("old_signal_electron_%d", i));
  //}


  //RooRealVar mtt("mtt", "mtt", 500, 2000, "GeV/c^2");
  ////RooRealVar weight("weight", "weight", 0, 100000);

  //RooCategory lepton_type("lepton_type", "lepton_type");
  //lepton_type.defineType("muon", 13);
  //if (! muonsOnly) {
    //lepton_type.defineType("electron", 11);
  //}

  //// Bin dataset
  ////int nBins = 1500 / 5;
  //int nBins = 50;
  //[>
  //RooBinning binning(nBins, 500, 2000);
  //mtt.setBinning(binning);
  //*/
  //mtt.setRange(500, 2000);

  //// Set binning to 5 GeV
  //mtt.setBins(nBins);

  //TParameter<int> n("iterations", NUM_ITER);

  //// Load signal dataset
  //std::map<int, std::shared_ptr<TChain>> chains;
  //std::map<int, std::shared_ptr<RooDataSet>> signalDatasets;
  //std::map<int, int> nMuons;
  //std::map<int, int> nElectrons;
  //std::map<int, std::shared_ptr<RooDataHist>> binnedSignalDataset_muons;
  //std::map<int, std::shared_ptr<RooDataHist>> binnedSignalDataset_electrons;

  //std::map<int, std::shared_ptr<RooHistPdf>> muon_pdf;
  //std::map<int, std::shared_ptr<RooHistPdf>> electron_pdf;
  //std::map<int, std::shared_ptr<RooAbsPdf::GenSpec>> genSpecs_muon;
  //std::map<int, std::shared_ptr<RooAbsPdf::GenSpec>> genSpecs_electron;

  //TString keysFile = TString::Format("%s/keyspdf_systematics_%d_%d_btag_pdf.root", base_path.c_str(), mass, btag);
  //RooWorkspace keys_workspace("w", "Frit signal workspace");
  //keys_workspace.import(mtt);
  //keys_workspace.import(n);

  //for (int i = btagMin; i <= btagMax; i++) {
    //TString treeName = TString::Format("dataset_%dbtag", i);
    //chains[i].reset(new TChain(treeName));
    //chains[i]->Add(signalDatasetFile.c_str());

    //signalDatasets[i].reset(new RooDataSet(TString::Format("dataset_%d", i), "dataset", RooArgSet(mtt, lepton_type), RooFit::Import(*chains[i])));
    //// Keep 7000 events max
    ////signalDatasets[i].reset(static_cast<RooDataSet*>(signalDatasets[i]->reduce(RooFit::EventRange(0, 7000), RooFit::SelectVars(RooArgSet(mtt, lepton_type)))));

    //Roo1DTable* table = signalDatasets[i]->table(lepton_type);
    //nMuons[i] = table->get("muon");
    //nElectrons[i] = table->get("electron");
    //delete table;

    //std::cout << i << " b-tag: " << nMuons[i] << " muon events, " << nElectrons[i] << " electron events" << std::endl;

    //RooDataHist* binnedSignalDataset = signalDatasets[i]->binnedClone();

    //binnedSignalDataset_muons[i].reset(static_cast<RooDataHist*>(binnedSignalDataset->reduce(RooFit::SelectVars(RooArgSet(mtt)), RooFit::Cut("lepton_type == 13"), RooFit::Name(TString::Format("binned_dataset_muon_%d", i)))));
    //binnedSignalDataset_electrons[i].reset(static_cast<RooDataHist*>(binnedSignalDataset->reduce(RooFit::SelectVars(RooArgSet(mtt)), RooFit::Cut("lepton_type == 11"), RooFit::Name(TString::Format("binned_dataset_electron_%d", i)))));

    //delete binnedSignalDataset;

    //// Create pdf from dataset
    //muon_pdf[i].reset(new RooHistPdf(TString::Format("muon_hist_pdf_%d", i), "muon_hist_pdf", RooArgSet(mtt), *binnedSignalDataset_muons[i]));
    //electron_pdf[i].reset(new RooHistPdf(TString::Format("electron_hist_pdf_%d", i), "electron_hist_pdf", RooArgSet(mtt), *binnedSignalDataset_electrons[i]));

    //genSpecs_muon[i].reset(muon_pdf[i]->prepareMultiGen(RooArgSet(mtt), RooFit::NumEvents(nMuons[i]), RooFit::Extended(true), RooFit::AutoBinned(false)[>, RooFit::Verbose(true)<]));

    //if (! muonsOnly) {
      //genSpecs_electron[i].reset(electron_pdf[i]->prepareMultiGen(RooArgSet(mtt), RooFit::NumEvents(nElectrons[i]), RooFit::Extended(true), RooFit::AutoBinned(false)));
    //}

    //keys_workspace.import(*muon_pdf[i]);
    //keys_workspace.import(*electron_pdf[i]);
    //keys_workspace.import(*muon_pdf_rookeys[i]);
    //keys_workspace.import(*electron_pdf_rookeys[i]);

  //}

  //for (int i = 0; i < NUM_ITER; i++) {

    //std::cout << "Iteration #" << i + 1 << " over " << NUM_ITER << std::endl;


    //std::map<int, std::shared_ptr<RooDataSet>> muon_toyData;
    //std::map<int, std::shared_ptr<RooDataSet>> electron_toyData;

    //std::map<int, std::shared_ptr<RooKeysPdf>> muon_toy_pdf;
    //std::map<int, std::shared_ptr<RooKeysPdf>> electron_toy_pdf;

    //// Create RooKeysPdf from toys
    //for (int b = btagMin; b <= btagMax; b++) {

      //std::cout << "Generating distribution for " << b << " b-tag..." << std::endl;
      //muon_toyData[b].reset(muon_pdf[b]->generate(*genSpecs_muon[b]));
      //if (!muonsOnly)
        //electron_toyData[b].reset(electron_pdf[b]->generate(*genSpecs_electron[b]));
      
      //std::cout << "Generating keys pdf for " << b << " b-tag ..." << std::endl;
      //std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
      //muon_toy_pdf[b].reset(new RooKeysPdf("signal_muon", "Keys pdf for signal", mtt, *muon_toyData[b], RooKeysPdf::MirrorBoth, 2));
      //if (!muonsOnly) {
        //electron_toy_pdf[b].reset(new RooKeysPdf("signal_electron", "Keys pdf for signal", mtt, *electron_toyData[b], RooKeysPdf::MirrorBoth, 2));
      //}
      //std::chrono::seconds t = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - t0);
      //std::cout << "Done in " << t.count() << "s." << std::endl;

      //// Save the pdf in a workspace

      //TString workspaceFile = TString::Format("%s/frit/temporary_Zprime%d_%s_%d_btag_workspace.root", base_path.c_str(), mass, analysisName.c_str(), b);
      //RooWorkspace temp_workspace("w", "Frit signal workspace");

      //temp_workspace.import(*muon_toy_pdf[b]);
      //if (! muonsOnly)
        //temp_workspace.import(*electron_toy_pdf[b]);

      //temp_workspace.writeToFile(workspaceFile, true);

      //TString datasetName = TString::Format("dataset_muon_%d_%d", b, i);
      //muon_toyData[b]->SetName(datasetName);

      //datasetName = TString::Format("dataset_electron_%d_%d", b, i);
      //electron_toyData[b]->SetName(datasetName);

      //keys_workspace.import(*muon_toyData[b]);
      //keys_workspace.import(*electron_toyData[b]);
    //}

    //TString workspaceFile = TString::Format("%s/frit/temporary_Zprime%d_%s_%%d_btag_workspace.root", base_path.c_str(), mass, analysisName.c_str());
    //std::cout << "Fitting..." << std::endl;
    //pid_t child = fork();
    //if (child == 0) {

      //int fd = open("/dev/null", O_WRONLY);
      //dup2(fd, STDOUT_FILENO); 
      //dup2(fd, STDIN_FILENO);

      //char** params = getSystCLParameters(ss.str(), file, singleFile, fixBackground, muonsOnly, btag, "--shared-memory", "--shm-key", strKey.str().c_str(), "--workspace", workspaceFile.Data(), NULL);
      //if (execv("./fitMtt", params) < 0) {
        //perror("Can't execute fitMtt");
      //}
      //exit(0);
    //} else {
      //waitpid(child, NULL, 0);
    //}

    //std::cout << "Done. # events: " << shm->nSignalEvents << "; Cross-section: " << shm->sigma << " pb." << std::endl << std::endl;

    //sigma->Fill(shm->sigma);
    //events->Fill(shm->nSignalEvents);
    //pull->Fill((shm->nSignalEvents - events_reference) / shm->nSignalEvents_error);
    //residuals->Fill(shm->nSignalEvents - events_reference);

    //for (int b = btagMin; b <= btagMax; b++) {
      //TString pdfName = TString::Format("signal_muon_%d_%d", b, i);
      //muon_toy_pdf[b]->SetName(pdfName);

      //pdfName = TString::Format("signal_electron_%d_%d", b, i);
      //electron_toy_pdf[b]->SetName(pdfName);

      //keys_workspace.import(*muon_toy_pdf[b]);
      //keys_workspace.import(*electron_toy_pdf[b]);
    //}
  //}

  //keys_workspace.writeToFile(keysFile); 

  //TFile* f = TFile::Open(filename, "recreate");
  //sigma->Write();
  //pull->Write();
  //events->Write();
  //residuals->Write();
  //f->Close();

  //delete f;

  //delete sigma;
  //delete pull;

  //shmdt(shm);

  //shmctl(shmid, IPC_RMID, NULL);
//}

//void saveSystematic(int mass, int btag, double syst) {

  //FILE* lock = fopen((base_path + "/systematics.lock").c_str(), "w+");
  //lockf(fileno(lock), F_LOCK, 0); // This will block until we have the right to write in the file

  //Json::Reader reader;
  //Json::Value root;
  //std::ifstream file((base_path + "/systematics.json").c_str());
  //reader.parse(file, root);
  //file.close();

  //std::stringstream ss;
  //ss << mass;
  //std::string strMass = ss.str();

  //ss.clear(); ss.str(std::string());
  //ss << btag;
  //std::string btagStr = ss.str();

  //root[getAnalysisUUID()][strMass][btagStr]["signal_pdf"] = syst;

  //FILE* fd = fopen((base_path + "/systematics.json").c_str(), "w+");
  //Json::StyledWriter writer;
  //const std::string json = writer.write(root);
  //fwrite(json.c_str(), json.length(), 1, fd);
  //fclose(fd);

  //fclose(lock);
//}

int main(int argc, char** argv) {
  try {
    TCLAP::CmdLine cmd("Generate a toy dataset", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "A text file containing a list of input files", true, "", "string");
    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "The input file", true, "", "string");

    cmd.xorAdd(inputListArg, inputFileArg);

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "The output file", true, "", "string", cmd);

    TCLAP::SwitchArg usePDFArg("", "use-pdf", "Use the background PDF to generate toys", cmd);
    cmd.parse(argc, argv);

    RooMsgService::instance().setStreamStatus(0, false);
    RooMsgService::instance().setStreamStatus(1, false);
    RooMsgService::instance().addStream(RooFit::ERROR); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL
    RooMsgService::instance().setSilentMode(true);
    
    std::vector<std::string> inputFiles;
    if (inputFileArg.isSet()) {
      inputFiles.push_back(inputFileArg.getValue());
    } else {
      loadInputFiles(inputListArg.getValue(), inputFiles);
    }

    process(inputFiles, outputFileArg.getValue(), usePDFArg.getValue());

  } catch (TCLAP::ArgException& e) {
    std::cerr << e.error() << std::endl;
  }
}
