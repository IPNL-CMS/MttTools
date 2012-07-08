//#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <list>

#include <TChain.h>
#include <TFile.h>

#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooMsgService.h>
using namespace RooFit;

#include <json/json.h>
#include <tclap/CmdLine.h>

struct Dataset {
  std::string name;
  int mass;
};

void loadDatasets(std::vector<Dataset>& datasets) {

  Json::Reader reader;
  Json::Value root;
  std::ifstream file("parameters.json");
  bool success = reader.parse(file, root);
  file.close();
  if (! success) {
    std::cerr << "ERROR: Failed to parse 'parameters.json'" << std::endl;
    exit(1);
  }

  //FIXME: Need a check
  Json::Value params = root["parameters"]["datasets"]["mc"];
  Json::Value::Members members = params.getMemberNames();
  for (Json::Value::Members::iterator it = members.begin(); it != members.end(); ++it) {
    const std::string dataset = *it;
    const int mass = params[dataset].asInt();
    Dataset d = {dataset, mass};
    datasets.push_back(d);
  }
}

std::vector<RooDataSet*> reduceForSpecificMassAndChannel(TChain *chain) {
  RooRealVar mtt_AfterChi2("mtt_AfterChi2", "mtt_AfterChi2",0,3000, "GeV/c^2");
  RooRealVar mtt_1stjetpt("mtt_1stjetpt", "mtt_1stjetpt",0, 1000, "GeV/c");
  RooRealVar mtt_2ndjetpt("mtt_2ndjetpt", "mtt_2ndjetpt",0, 1000, "GeV/c");
  RooRealVar mtt_BestSolChi2("mtt_BestSolChi2", "mtt_BestSolChi2",-1, 500, "");
  RooRealVar mtt_isSel("mtt_isSel", "mtt_isSel",-1, 100, "");
  RooRealVar mtt_NBtaggedJets_TCHEL("mtt_NBtaggedJets_TCHEL", "mtt_NBtaggedJets_TCHEL",-1, 20, "");

  RooArgSet argSet(mtt_AfterChi2, mtt_isSel, mtt_BestSolChi2, mtt_1stjetpt, mtt_2ndjetpt, mtt_NBtaggedJets_TCHEL); 
  RooDataSet allData("allData", "chained events", chain, argSet);
  //allData.Print("v");

  // Reduce data set
  // signal region
  // TODO: Read selection from file?
  RooDataSet* reducedData = (RooDataSet*) allData.reduce(RooArgSet(mtt_AfterChi2, mtt_NBtaggedJets_TCHEL), "mtt_isSel == 1 && mtt_1stjetpt > 70. && mtt_2ndjetpt > 50. && mtt_AfterChi2 > 0. && mtt_BestSolChi2 < 500.");

  // mtt_NBtaggedJets_TCHEL > 1 &&
  RooDataSet* reducedData_0btag = (RooDataSet*) reducedData->reduce(RooArgSet(mtt_AfterChi2), "mtt_NBtaggedJets_TCHEL == 0");
  RooDataSet* reducedData_1btag = (RooDataSet*) reducedData->reduce(RooArgSet(mtt_AfterChi2), "mtt_NBtaggedJets_TCHEL == 1");
  RooDataSet* reducedData_2btag = (RooDataSet*) reducedData->reduce(RooArgSet(mtt_AfterChi2), "mtt_NBtaggedJets_TCHEL > 1");

  std::vector<RooDataSet*> reducedDatasets;
  reducedDatasets.push_back(reducedData_0btag);
  reducedDatasets.push_back(reducedData_1btag);
  reducedDatasets.push_back(reducedData_2btag);

  //RedData->Print("v");

  //TODO: Needed?!
  //TString outputfile = TString::Format("ds_%s-%s_%s.root", "Zprime_M1250GeV_W1250MeV", i"nominal", channel.c_str());
  //TFile* histoFile = new TFile(outputfile.Data(), "RECREATE");
  //RedData->Write();
  //histoFile->Close();

  return reducedDatasets;
}

bool reduceForSpecificDataset(const std::string& dataset, const std::string& jecType, const int mass, const std::string& rfioPath) {
  std::cout << "[" << getpid() << "] Reducing " << dataset << " for " << jecType << std::endl;
  
  RooRealVar lepton_type("lepton_type", "lepton_type", 11);

  std::cout << "[" << getpid() << "] Reducing for electrons ..." << std::endl;

  TString rfioDataset;

  if (jecType == "nominal") {
    rfioDataset = TString::Format("%s/mc/%%s/%s/*.root", rfioPath.c_str(), dataset.c_str());
  } else {
    rfioDataset = TString::Format("%s/mc/%%s-%s/%s-%s/*.root", rfioPath.c_str(), jecType.c_str(), dataset.c_str(), jecType.c_str());
  }

  TChain *chain = new TChain("Mtt");
  if (! chain->Add(TString::Format(rfioDataset, "semie")))
    return false;

  std::vector<RooDataSet*> reducedDataSetsElectron = reduceForSpecificMassAndChannel(chain); 

  delete chain;

  std::cout << "[" << getpid() << "] Reducing for muons ..." << std::endl;

  chain = new TChain("Mtt");
  if (! chain->Add(TString::Format(rfioDataset, "semimu")))
    return false;

  std::vector<RooDataSet*> reducedDataSetsMuon = reduceForSpecificMassAndChannel(chain);

  for (int i = 0; i < 3; i++) {
    lepton_type.setVal(11);
    reducedDataSetsElectron[i]->addColumn(lepton_type);
    lepton_type.setVal(13);
    reducedDataSetsMuon[i]->addColumn(lepton_type);

    // Make final dataset
    RooDataSet reducedDataSetFinal(*reducedDataSetsElectron[i], "reducedDataSet");
    reducedDataSetFinal.append(*reducedDataSetsMuon[i]);

    TString outputfile = TString::Format("ds_Zprime%d-%s_all_%d_btag.txt", mass, jecType.c_str(), i);
    reducedDataSetFinal.write(outputfile.Data());
  }


  std::cout << "[" << getpid() << "] Done." << std::endl;

  delete chain;

  return true;
}

int main(int argc, char** argv)
{
  try {
    TCLAP::CmdLine cmd("reduce Zprime dataset", ' ', "0.1");

    TCLAP::ValueArg<std::string> fileArg("f", "file", "The parameter file", false, "reduce_dataset.txt", "string");

    std::vector<std::string> allowed;
    allowed.push_back("mc");
    allowed.push_back("data");
    TCLAP::ValuesConstraint<std::string> allowedVals( allowed );

    TCLAP::ValueArg<std::string> rfiopathArg("p", "path", "Path where root file are located on rfio", false, "rfio:///dpm/in2p3.fr/home/cms/data/store/user/beaupere/Extractor_428_03jan12/", "string");

    cmd.add(fileArg);
    cmd.add(rfiopathArg);

    cmd.parse(argc, argv);

    std::string file = fileArg.getValue();
    std::string rfioPath = rfiopathArg.getValue();

    // Set RooFit verbosity
    RooMsgService::instance().setStreamStatus(0,false);
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().addStream(RooFit::ERROR); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL

    std::vector<std::string> JEC;
    JEC.push_back("nominal");
    JEC.push_back("JECup");
    JEC.push_back("JECdown");

    std::vector<Dataset> datasets;
    loadDatasets(datasets);
    for (std::vector<Dataset>::iterator it = datasets.begin(); it != datasets.end(); ++it) {

      std::vector<pid_t> children;
      std::vector<std::string>::const_iterator jecType = JEC.begin();

      std::cout << "Processing " << it->name << std::endl;

      for (; jecType != JEC.end(); ++jecType) {

        pid_t pid = fork();
        if (pid == 0) {
          reduceForSpecificDataset(it->name, *jecType, it->mass, rfioPath);
          exit(0);
        } else {
          children.push_back(pid);
        }
      }

      // Wait until all child are done before processing the next dataset
      for (std::vector<pid_t>::iterator child = children.begin(); child != children.end(); ++child) {
        waitpid(*child, NULL, 0);
      }

      std::cout << "All done for " << it->name << std::endl;
    }

    std::cout << "All done, now exiting ..." << std::endl;

  } catch (TCLAP::ArgException& e) {
  }

}
