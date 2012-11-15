//#include <Riostream.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <pthread.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <errno.h>

#include <list>
#include <map>

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include "RooMsgService.h"
using namespace RooFit;

#include <json/json.h>
#include <tclap/CmdLine.h>

std::vector<std::string> JEC;

struct thread_data {
  std::map<int, std::vector<std::string> > datasets;
  std::string rfioPath;
  std::string prefix;
  std::string jecType;
};

std::vector<RooDataSet*> reduce(TChain *chain) {
  RooRealVar mtt_AfterChi2("mtt_AfterChi2", "mtt_AfterChi2", 0, 3000, "GeV/c^2");
  RooRealVar mtt_1stjetpt("mtt_1stjetpt", "mtt_1stjetpt", 0, 1000, "GeV/c");
  RooRealVar mtt_2ndjetpt("mtt_2ndjetpt", "mtt_2ndjetpt", 0, 1000, "GeV/c");
  RooRealVar mtt_BestSolChi2("mtt_BestSolChi2", "mtt_BestSolChi2",-1, 500, "");
  RooRealVar mtt_isSel("mtt_isSel", "mtt_isSel",-1, 100, "");
  RooRealVar mtt_NBtaggedJets_TCHEL("mtt_NBtaggedJets_TCHEL", "mtt_NBtaggedJets_TCHEL",-1, 20, "");

  RooArgSet argSet(mtt_AfterChi2, mtt_isSel, mtt_BestSolChi2, mtt_1stjetpt, mtt_2ndjetpt, mtt_NBtaggedJets_TCHEL); 
  
  RooDataSet allData("allData", "chained events", chain, argSet);

  // Reduce data set
  // signal region
  // TODO: Read selection from file?
  RooDataSet* reducedData = (RooDataSet*) allData.reduce(RooArgSet(mtt_AfterChi2, mtt_NBtaggedJets_TCHEL), "mtt_isSel == 1 && mtt_1stjetpt > 70. && mtt_2ndjetpt > 50. &&  mtt_AfterChi2 > 0. && mtt_BestSolChi2 < 500.");

  RooDataSet* reducedData_0btag = (RooDataSet*) reducedData->reduce(RooArgSet(mtt_AfterChi2), "mtt_NBtaggedJets_TCHEL == 0");
  RooDataSet* reducedData_1btag = (RooDataSet*) reducedData->reduce(RooArgSet(mtt_AfterChi2), "mtt_NBtaggedJets_TCHEL == 1");
  RooDataSet* reducedData_2btag = (RooDataSet*) reducedData->reduce(RooArgSet(mtt_AfterChi2), "mtt_NBtaggedJets_TCHEL > 1");

  std::vector<RooDataSet*> reducedDatasets;
  reducedDatasets.push_back(reducedData_0btag);
  reducedDatasets.push_back(reducedData_1btag);
  reducedDatasets.push_back(reducedData_2btag);

  // mtt_NBtaggedJets_TCHEL > 1 &&

  return reducedDatasets;
}

std::vector<RooDataSet*> reduceForSpecificDataset(const std::string& dataset, const std::string& jecType, const std::string& rfioPath) {
  std::cout << "[" << getpid() << "] Reducing dataset for " << dataset << ":" << jecType << std::endl;

  TString rfioDataset;
  //TString datasetWithJEC = TString::Format(dataset.c_str(), jecType.c_str());

  if (jecType == "nominal") {
    rfioDataset = TString::Format("%s/data/%s/*.root/Mtt", rfioPath.c_str(), dataset.c_str());
  } else {
    rfioDataset = TString::Format("%s/data/%s-%s/*.root/Mtt", rfioPath.c_str(), dataset.c_str(), jecType.c_str());
  }

  // It's actually not needed, but still funny!
  TString randomName = TString::Format("%d%s%s", getpid(), dataset.c_str(), jecType.c_str());
  TChain *chain = new TChain(randomName);
  if (! chain->Add(rfioDataset)) {
    std::cerr << "[" << getpid() << "] ERROR: No file loaded" << std::endl;
    return std::vector<RooDataSet*>();
  }

  std::vector<RooDataSet*> reducedDataSets = reduce(chain); 


  std::cout << "[" << getpid() << "] Done" << std::endl;

  delete chain;

  return reducedDataSets;
}

void datasetToStream(RooDataSet* dataset, ostream& stream) {
  // From http://root.cern.ch/root/html/src/RooDataSet.cxx.html#PLxiPB
  for (int i = 0; i < dataset->numEntries(); i++) {
    RooArgList list(*(dataset->get(i)), "line");
    list.writeToStream(stream, kTRUE);
  }
}

typedef struct {
  int d[2]; // 0 is for reading, 1 for writing
  bool closed;
  int currentFileIndex;
  int remaingSizeBeforeEnd;
} descriptor_t;

void *reduceDatasets(void * ptr) {
  thread_data *data = static_cast<thread_data*>(ptr);

  RooRealVar leptonType("leptonType", "leptonType", 11);
  std::string jecType = data->jecType;
  // fork() to process muons and electrons side by side

  // Delete files
  for (int i = 0; i < 3; i++) {
    ofstream f(TString::Format("%s-%s_%d_btag.txt", data->prefix.c_str(), (jecType).c_str(), i), ofstream::out);
    f.close();
  }

  std::vector<descriptor_t> descriptors;
  std::map<int, std::vector<std::string> >::const_iterator it = data->datasets.begin();

  for (; it != data->datasets.end(); ++it) {

    descriptor_t descriptor;
    descriptor.closed = false;
    descriptor.currentFileIndex = -1;
    descriptor.remaingSizeBeforeEnd = -1;

    //FIXME: Check for error
    if (pipe(descriptor.d) < 0) {
      perror("pipe");
    }

    fcntl(descriptor.d[0], F_SETFL, O_NONBLOCK);

    pid_t pid = fork();

    if (pid == 0) {
      // Close unused read end
      close(descriptor.d[0]);

      leptonType.setVal(it->first);
      std::vector<std::string>::const_iterator it2 = it->second.begin();

      RooDataSet* leptonicDataset[3] = {NULL};

      for (; it2 != it->second.end(); ++it2) {

        std::vector<RooDataSet*> dataset;

        dataset = reduceForSpecificDataset(*it2, jecType, data->rfioPath);

        for (int i = 0; i < 3; i++) {
          dataset[i]->addColumn(leptonType);
          if (leptonicDataset[i] == NULL) {
            leptonicDataset[i] = new RooDataSet(*dataset[i], TString::Format("reducedDataset_%d_btag", i).Data());
          } else {
            leptonicDataset[i]->append(*dataset[i]);
          }

          delete dataset[i];
        }
      }

      if (leptonicDataset[0] == NULL) {
        close(descriptor.d[1]);
        exit(1);
      }


      for (int i = 0; i < 3; i++) {
        std::stringstream streamedDataset;
        datasetToStream(leptonicDataset[i], streamedDataset);

        //std::cout << streamedDataset.str() << std::endl;

        streamedDataset.seekg(0, std::ios::end);
        unsigned int streamSize = streamedDataset.tellg();
        streamedDataset.seekg(0, std::ios::beg);


        if (! write(descriptor.d[1], (void *) &streamSize, sizeof(unsigned int)))
        {
          perror("size write");
          close(descriptor.d[1]);
          exit(2);
        }

        char buffer[4096];
        int size = 0;
        int totalSize = 0;
        while (true) {
          streamedDataset.read(buffer, 4095);
          if ((size = streamedDataset.gcount()) <= 0)
            break;
          totalSize += size;

          int writtenTotal = 0;
          while (writtenTotal < size) {

            int written = write(descriptor.d[1], (void *) &buffer[writtenTotal], size - writtenTotal);
            if (written <= 0) {
              perror("write");
              close(descriptor.d[1]);
              exit(2);
            }
            writtenTotal += written;
          }
        }

        delete leptonicDataset[i];
      }

      close(descriptor.d[1]);

      std::cout << "[" << getpid() << "] All done" << std::endl;
      exit(0);

    } else {
      // Close unused write end
      close(descriptor.d[1]);
      descriptors.push_back(descriptor);
    }

    //break;
  }

  while (true) {

    fd_set set;
    FD_ZERO(&set);

    int size = 0;
    for (std::vector<descriptor_t>::iterator desc = descriptors.begin(); desc != descriptors.end(); ++desc) {
      if (! desc->closed) {
        FD_SET(desc->d[0], &set);
        size++;
      }
    }

    if (size == 0)
      break;

    int ready = select(FD_SETSIZE, &set, NULL, NULL, NULL);

    if (ready <= 0) {
      if (errno != EINTR) {
        perror("select");
        break;
      } else
        continue;
    }

    for (std::vector<descriptor_t>::iterator desc = descriptors.begin(); desc != descriptors.end(); ++desc) {
      if (FD_ISSET(desc->d[0], &set)) {
        char buffer[4096] = {0};
        int red = 0;

        if (desc->currentFileIndex >= 2) {
          close(desc->d[0]);
          desc->closed = true;
          break;
        }

        if (desc->remaingSizeBeforeEnd <= 0) {
          unsigned int redSize = 0;
          read(desc->d[0], (void*) &redSize, sizeof(unsigned int));
          desc->remaingSizeBeforeEnd = redSize;
          desc->currentFileIndex++;
        }

        ofstream f(TString::Format("%s-%s_%d_btag.txt", data->prefix.c_str(), (jecType).c_str(), desc->currentFileIndex), ofstream::app);

        do
        {
          int bufferSize = std::min(4095, desc->remaingSizeBeforeEnd);

          red = read(desc->d[0], (void *) buffer, bufferSize);
          if (red <= 0)
            break;

          desc->remaingSizeBeforeEnd -= red;

          f.write(buffer, red);
          
        } while (desc->remaingSizeBeforeEnd > 0);

        f.close();
      }
    }
  }

  delete data;

  return NULL;
}

void loadDatasets(std::map<int, std::vector<std::string> >& datasets) {

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
  Json::Value params = root["parameters"]["datasets"]["data"];

  Json::Value electrons = params["electrons"];
  std::vector<std::string> datasetsElectrons;
  for (unsigned int i = 0; i < electrons.size(); i++) {
    datasetsElectrons.push_back(electrons[i].asString());
  }
  datasets[11] = datasetsElectrons;

  Json::Value muons = params["muons"];
  std::vector<std::string> datasetsMuons;
  for (unsigned int i = 0; i < muons.size(); i++) {
    datasetsMuons.push_back(muons[i].asString());
  }
  datasets[13] = datasetsMuons;
}

int main(int argc, char** argv)
{
  try {
    TCLAP::CmdLine cmd("reduce Zprime dataset", ' ', "0.1");

    TCLAP::ValueArg<std::string> fileArg("f", "file", "The parameter file", false, "reduce_dataset_data.txt", "string");
    TCLAP::ValueArg<std::string> rfiopathArg("p", "path", "Path where root file are located on rfio", false, "rfio:///dpm/in2p3.fr/home/cms/data/store/user/beaupere/Extractor_428_03jan12/", "string");
    TCLAP::ValueArg<std::string> prefixArg("o", "prefix", "Output file prefix", false, "ds_data_2012", "string");

    cmd.add(fileArg);
    cmd.add(rfiopathArg);
    cmd.add(prefixArg);

    cmd.parse(argc, argv);

    std::string file = fileArg.getValue();
    std::string rfioPath = rfiopathArg.getValue();
    std::string prefix = prefixArg.getValue();
    std::map<int, std::vector<std::string> > datasets;
    loadDatasets(datasets);

    if (datasets.size() != 2) {
      std::cerr << "Warning: the number of differents types is invalid. 2 expected, " << datasets.size() << " given" << std::endl;
      //FIXME ?
      return 1;
    }

    JEC.push_back("nominal");
    JEC.push_back("JECup");
    JEC.push_back("JECdown");

    // Set RooFit verbosity
    RooMsgService::instance().setStreamStatus(0,false);
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().addStream(RooFit::ERROR); // DEBUG  INFO  PROGRESS  WARNING ERROR FATAL

    std::vector<pid_t> threads;

    std::vector<std::string>::const_iterator jecType = JEC.begin();
    for (; jecType != JEC.end(); ++jecType) {
      thread_data *data = new thread_data();
      data->datasets = datasets;
      data->rfioPath = rfioPath;
      data->prefix = prefix;
      data->jecType = *jecType;

      pid_t pid = fork();
      if (pid == 0) {
        reduceDatasets((void *) data);
        exit(0);
      } else {
        threads.push_back(pid);
      }

    }

    for (std::vector<pid_t>::iterator it = threads.begin(); it != threads.end(); ++it) {
      waitpid(*it, NULL, 0);
    }

  } catch (TCLAP::ArgException& e) {
    std::cout << e.error() << std::endl;
  }

}
