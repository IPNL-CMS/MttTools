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
#include <TGraphAsymmErrors.h>

#include "tclap/CmdLine.h"

#include "yaml-cpp/yaml.h"

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

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

struct Input {
  boost::filesystem::path filename;
  uint64_t generated_events;
  float cross_section;
  int32_t type; // 0: Nominal, 1: Up, -1: Down
  float top_pt_weight;

  std::shared_ptr<TFile> file;
};

std::shared_ptr<TH1> getHistogram(const std::string& name, const std::vector<Input>& inputs, int type) {

  TH1* h = nullptr;

  for (const auto& input: inputs) {
    if (input.type != type)
      continue;

    TH1* f = static_cast<TH1*>(input.file->Get(name.c_str()));
    f->Scale(input.cross_section / (input.generated_events * input.top_pt_weight));

    if (! h) {
      h = static_cast<TH1*>(f->Clone());
      h->SetDirectory(NULL);
    } else
      h->Add(f);
  }

  return std::shared_ptr<TH1>(h);
}

int main(int argc, char** argv) {

  try {
    TCLAP::CmdLine cmd("Create uncertainties histograms from systematic samples", ' ', "0.1");

    TCLAP::ValueArg<std::string> inputFileArg("i", "input-file", "Input file in YML format", true, "", "string", cmd);

    cmd.parse(argc, argv);

    if (inputFileArg.getValue().length() == 0)
      return 1;

    YAML::Node root = YAML::LoadFile(inputFileArg.getValue());

    // Load systematics, btag and type
    std::vector<std::string> systs = root["systs"].as<std::vector<std::string>>();
    std::vector<int> btags = root["btags"].as<std::vector<int>>();
    std::vector<std::string> types = root["types"].as<std::vector<std::string>>();

    for (const auto& it: root["data"]) {
      std::vector<Input> inputs;
      const YAML::Node& inputs_node = it["inputs"];
      for (const auto& input_node: inputs_node) {
        Input input;
        input.filename = input_node["file"].as<std::string>();
        std::string type = input_node["type"].as<std::string>();
        if (type == "nominal")
          input.type = 0;
        else if (type == "up")
          input.type = 1;
        else if (type == "down")
          input.type = -1;

        if (input_node["generated-events"])
          input.generated_events = input_node["generated-events"].as<uint64_t>();
        else
          input.generated_events = 1;

        if (input_node["cross-section"])
          input.cross_section = input_node["cross-section"].as<float>();
        else
          input.cross_section = 1.;

        input.top_pt_weight = 1.;
        inputs.push_back(input);
      }

      boost::filesystem::path output = it["output"].as<std::string>();

      bool manual = false;
      if (it["manual"])
        manual = it["manual"].as<bool>();

      // Loop over systematics, btag & type
      for (auto& btag: btags) {
        for (auto& type: types) {

          std::string btag_str = std::to_string(btag);
          if (btag < 0) {
            btag_str = "all";
          }

          // Create main output file, containing all systematic uncertainties merged
          boost::format mainFormatter(output.string());
          mainFormatter.exceptions(boost::io::all_error_bits ^ (boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );
          mainFormatter % btag_str % type % "total_syst_errors";

          std::shared_ptr<TFile> mainOutputFile;
          
          if (! manual) {
            mainOutputFile.reset(TFile::Open(mainFormatter.str().c_str(), "recreate"));
          }

          std::map<std::string, std::shared_ptr<TH1>> total_uncertainties_hists;

          std::vector<std::string> goodSysts = (manual) ? std::vector<std::string>({"dummy"}) : systs;
          for (auto& syst: goodSysts) {
            
            // Load input files
            for (auto& input: inputs) {
              boost::format formatter(input.filename.string());
              formatter.exceptions(boost::io::all_error_bits ^ (boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );

              std::string systSuffix;
              if (input.type == 1)
                systSuffix += "Up";
              else if (input.type == -1)
                systSuffix += "Down";

              // Hack for JEC & JER
              if (syst == "JEC" || syst == "JER")
                std::transform(systSuffix.begin(), systSuffix.end(), systSuffix.begin(), ::tolower);

              std::string fullSyst = syst + systSuffix;
              formatter % btag_str % type;

              if (! manual)
                formatter % fullSyst;

              input.file.reset(TFile::Open(formatter.str().c_str()));

              // Try to find if a '.info' file exists, containing the mean
              // value of the top pt weights
              boost::filesystem::path p(formatter.str());
              p.replace_extension("info");

              if (boost::filesystem::exists(p)) {
                std::ifstream info(p.string());
                info >> input.top_pt_weight;
              }
            }

            // List all histograms
            // For each, retrieve the nominal histogram, as well as the up & down variation

            // Retrieve just the names of the histograms
            std::vector<std::string> histograms;
            TList* list = inputs[0].file->GetListOfKeys();
            TIter iter(list);
            TKey* key;
            TObject* obj;
            while ((key = (TKey*) iter())) {
              obj = key->ReadObj();
              if (obj->InheritsFrom("TH1")) {
                histograms.push_back(obj->GetName());
              }
            }

            // Create output file
            boost::format formatter(output.string());
            formatter.exceptions(boost::io::all_error_bits ^ (boost::io::too_many_args_bit | boost::io::too_few_args_bit )  );
            formatter % btag_str % type % syst;

            std::shared_ptr<TFile> outputFile(TFile::Open(formatter.str().c_str(), "recreate"));
            std::vector<std::shared_ptr<TH1>> uncertainties_hists;

            for (const auto& histogram: histograms) {
              std::shared_ptr<TH1> nominal = getHistogram(histogram, inputs, 0);
              std::shared_ptr<TH1> up = getHistogram(histogram, inputs, 1);
              std::shared_ptr<TH1> down = getHistogram(histogram, inputs, -1);

              std::string name = nominal->GetName();
              nominal->SetName(TString::Format("%s_%s_%s", name.c_str(), syst.c_str(), "nominal"));
              up->SetName(TString::Format("%s_%s_%s", name.c_str(), syst.c_str(), "up"));
              down->SetName(TString::Format("%s_%s_%s", name.c_str(), syst.c_str(), "down"));
              
              outputFile->cd();
              nominal->Write();
              up->Write();
              down->Write();

              if (! manual) {
                mainOutputFile->cd();
                nominal->Write();
                up->Write();
                down->Write();
              }

              std::shared_ptr<TH1> uncertainties(static_cast<TH1*>(nominal->Clone()));
              uncertainties->SetDirectory(NULL);
              uncertainties->Reset(); // Keep binning but remove all events
              uncertainties->SetName(name.c_str());

              std::shared_ptr<TGraphAsymmErrors> uncertainties_asym = std::make_shared<TGraphAsymmErrors>();
              uncertainties_asym->SetName(TString::Format("%s_%s_asym", name.c_str(), syst.c_str()));

              std::shared_ptr<TH1>& total_uncertainties = total_uncertainties_hists[histogram];
              if (! manual) {
                if (! total_uncertainties.get()) {
                  total_uncertainties.reset(static_cast<TH1*>(uncertainties->Clone()));
                  total_uncertainties->SetDirectory(NULL);
                }
              }

              uint32_t index = 0;
              for (int i = 1; i <= uncertainties->GetNbinsX(); i++) {
                float nominal_x = nominal->GetBinCenter(i);
                float nominal_value = nominal->GetBinContent(i);

                float norm_plus = 0;
                float norm_minus = 0;
                if (up->GetBinError(i) / up->GetBinContent(i) < 0.4)
                  norm_plus = up->GetBinContent(i) - nominal_value;
                if (down->GetBinError(i) / down->GetBinContent(i) < 0.4)
                  norm_minus = nominal_value - down->GetBinContent(i);

                float delta = (fabs(norm_plus) + fabs(norm_minus)) / 2.;
                float delta_percent = (delta == 0 || nominal_value == 0) ? 0 : delta / nominal_value;

                uncertainties->SetBinContent(i, 0.);
                uncertainties->SetBinError(i, delta_percent);

                uncertainties_asym->SetPoint(index, nominal_x, 0.);
                uncertainties_asym->SetPointError(index, 0, 0, delta_percent, delta_percent);
                index++;

                if (! manual) {
                  total_uncertainties->SetBinContent(i, 0);
                  float total_error = total_uncertainties->GetBinError(i);
                  total_uncertainties->SetBinError(i, sqrt(delta_percent * delta_percent + total_error * total_error));
                }
              }

              uncertainties_hists.push_back(uncertainties);

              outputFile->cd();
              uncertainties_asym->Write();

              if (! manual) {
                mainOutputFile->cd();
                uncertainties_asym->Write();
              }
            }

            outputFile->cd();
            for (auto& h: uncertainties_hists) {
              h->Write();
            }
            outputFile->Close();
          }

          if (! manual) {
            mainOutputFile->cd();
            for (auto& it: total_uncertainties_hists) {
              it.second->Write();
            }
            mainOutputFile->Close();
          }
        }
      }
      

      
    }
    

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  return 0;
}
