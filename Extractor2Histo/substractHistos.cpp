#include <TH1.h>
#include <TStyle.h>
#include <TList.h>
#include <TCollection.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TKey.h>
#include <vector>
#include <map>
#include <fstream>
#include <TLorentzVector.h>

#include "tclap/CmdLine.h"

#include "yaml-cpp/yaml.h"

#include <boost/filesystem.hpp>

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

bool parseYAML(const std::string& file, float& crossSection, uint64_t& nGeneratedEvents) {
  YAML::Node f = YAML::LoadFile(file);

  if (! f["cross-section"]) {
    std::cout << "You need to specify a 'cross-section' value in '" << file << "'" << std::endl;
    return false;
  }

  crossSection = f["cross-section"].as<float>();

  if (! f["generated-events"]) {
    std::cout << "You need to specify a 'generated-events' value in '" << file << "'" << std::endl;
    return false;
  }

  nGeneratedEvents = f["generated-events"].as<uint64_t>();

  return true;
}

TFile* loadHistograms(const std::string& file, std::map<std::string, TH1*>& histograms) {
  TFile* f = TFile::Open(file.c_str());

  TList *list = f->GetListOfKeys();
  TIter iter(list);

  TKey* key;
  TObject* obj;

  while ((key = static_cast<TKey*>(iter()))) {
    obj = key->ReadObj();
    if (! obj->InheritsFrom("TH1"))
      continue;

    TH1* hist = static_cast<TH1*>(obj);
    histograms[hist->GetName()] = hist;
  }

  return f;
}

int main(int argc, char** argv) {

  try {

    TCLAP::CmdLine cmd("Convert extractor tuples to histograms", ' ', "0.1");

    TCLAP::ValueArg<std::string> outputFileArg("o", "output-file", "output file", true, "", "string", cmd);

    TCLAP::SwitchArg dataArg("", "data", "Is this data?", false);
    TCLAP::SwitchArg mcArg("", "mc", "Is this mc?", false);

    cmd.xorAdd(dataArg, mcArg);

    TCLAP::SwitchArg semimuArg("", "semimu", "Is this semi-mu channel?", false);
    TCLAP::SwitchArg semieArg("", "semie", "Is this semi-e channel?", false);

    cmd.xorAdd(semimuArg, semieArg);

    TCLAP::ValueArg<int> btagArg("", "b-tag", "Number of b-tagged jet to require", true, 2, "int", cmd);

    TCLAP::ValueArg<float> lumiArg("", "lumi", "Luminosity", false, 19700, "int", cmd);

    TCLAP::ValueArg<std::string> configArg("", "config-path", "Path to the configuration files", false, "./config", "int", cmd);

    TCLAP::UnlabeledMultiArg<std::string> inputFilesArg("inputFiles", "input files", true, "string", cmd);

    cmd.parse(argc, argv);

    bool isData = dataArg.isSet();

    std::vector<std::string> inputFiles = inputFilesArg.getValue();
    if (inputFiles.size() != 2) {
      throw TCLAP::ArgException("Only two input files are supported");
    }

    boost::filesystem::path root(configArg.getValue());
    
    boost::filesystem::path f1 (inputFiles[0]);
    boost::filesystem::path f2 (inputFiles[1]);

    f1 = root / f1.filename();
    f2 = root / f2.filename();

    f1.replace_extension("yaml");
    f2.replace_extension("yaml");

    if (! boost::filesystem::exists(f1) || ! boost::filesystem::exists(f2)) {
      std::cout << "Error: Please make sure that " << f1 << " and " << f2 << " exist and contains cross-sections and number of generated events." << std::endl;
      return 1;
    }

    // Get scale factor for first file
    float first_crossSection;
    uint64_t first_nGeneratedEvents;
    if (! parseYAML(f1.string(), first_crossSection, first_nGeneratedEvents))
      return 1;

    // Get scale factor for second file
    float second_crossSection;
    uint64_t second_nGeneratedEvents;
    if (! parseYAML(f2.string(), second_crossSection, second_nGeneratedEvents))
      return 1;

    std::cout << "First file cross-section: " << first_crossSection << std::endl;
    std::cout << "Second file cross-section: " << second_crossSection << std::endl;

    float first_generatedLumi = first_nGeneratedEvents / first_crossSection;
    float second_generatedLumi = second_nGeneratedEvents / second_crossSection;

    float first_scale = lumiArg.getValue() / first_generatedLumi;
    float second_scale = lumiArg.getValue() / second_generatedLumi;

    std::cout << "First file scale: " << first_scale << std::endl;
    std::cout << "Second file scale: " << second_scale << std::endl;
    
    std::map<std::string, TH1*> first_histograms;
    TFile* first = loadHistograms(inputFiles[0], first_histograms);

    std::map<std::string, TH1*> second_histograms;
    TFile* second = loadHistograms(inputFiles[1], second_histograms);

    if (first_histograms.size() != second_histograms.size()) {
      std::cout << "Warning: the two files don't have the same number of histograms." << std::endl;
    }

    std::map<std::string, TH1*>& a = first_histograms;
    std::map<std::string, TH1*>& b = second_histograms;

    float a_scale = first_scale;
    float b_scale = second_scale;

    if (b.size() > a.size()) {
      a = second_histograms;
      b = first_histograms;

      std::swap(a_scale, b_scale);
    }
    
    TFile* output = TFile::Open(outputFileArg.getValue().c_str(), "recreate");

    // Do a - b for each histograms, and save them in output
    for (auto& hist: a) {

      std::cout << "Substracting '" << hist.first << "'" << std::endl;

      TH1* clone = static_cast<TH1*>(hist.second->Clone());
      clone->Scale(a_scale);
      clone->Add(b[hist.first], -1 * b_scale);

      clone->Write();
      delete clone;
    }

    output->Close();
    delete output;

    first->Close();
    delete first;

    second->Close();
    delete second;

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
