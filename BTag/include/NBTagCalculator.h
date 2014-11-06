#pragma once

#include <vector>
#include <string>

#include <TChain.h>
#include <TBranch.h>
#include <TClonesArray.h>

#include <BTaggingEfficiencyProvider.h>
#include <BTagUtils.h>

class NBTagCalculator {

  public:
    NBTagCalculator(const std::vector<std::string>& inputFiles, SystVariation systVariation):
      m_systVariation(systVariation) {

        m_inputFiles = inputFiles;

        TChain* jets = loadChain("jet_PF");
        TChain* event = loadChain("event");

        setBranchAddress(jets, "n_jets", n_jets);
        setBranchAddress(jets, "jet_4vector", jet_p4);
        setBranchAddress(jets, "jet_isPFJetLoose", jet_isPFLoose);
        setBranchAddress(jets, "jet_algo_parton_flavor", jet_flavor);
        setBranchAddress(jets, "jet_btag_CSV", jet_CSV);
        setBranchAddress(jets, "jet_scaleFactor", jet_scaleFactor);

        setBranchAddress(event, "run", run);
        setBranchAddress(event, "evtID", eventNumber);

        // B-tagging efficiency provider for computing correct number of b-tagged jets
        m_btagging_efficiency_provider = std::make_shared<BTaggingEfficiencyProvider>("../BTag/TT_powheg_btagging_efficiency.root");
      }

    uint32_t getNumberOfBTaggedJets(uint32_t entry);

  private:
    TChain* loadChain(const std::string& treeName) {
      std::shared_ptr<TChain> output(new TChain(treeName.c_str()));

      for (const std::string& file: m_inputFiles) {
        output->Add(file.c_str());
      }

      output->SetCacheSize(500 * 1024 * 1024);
      output->SetBranchStatus("*", 0);

      m_chains.push_back(output);
      return output.get();
    }

    template<typename T>
      void setBranchAddress(TChain* chain, const std::string& name, T& address) {
        TBranch* branch = chain->FindBranch(name.c_str());
        if (!branch)
          return;

        activateBranch(branch);
        chain->AddBranchToCache(name.c_str(), true);
        chain->SetBranchAddress(name.c_str(), &address);
      }

    void activateBranch(TBranch* branch) {
      branch->SetStatus(1);
      TObjArray* objArray = branch->GetListOfBranches();
      for (int i = 0; i < objArray->GetEntries(); i++) {
        TBranch* b = static_cast<TBranch*> (objArray->At(i));
        if (b) {
          activateBranch(b);
        }
      }
    }

    std::vector<std::string> m_inputFiles;
    std::vector<std::shared_ptr<TChain>> m_chains;
    std::shared_ptr<BTaggingEfficiencyProvider> m_btagging_efficiency_provider;
    SystVariation m_systVariation;

    uint32_t n_jets = 0;
    int      jet_isPFLoose[100];
    int      jet_flavor[100];
    float    jet_CSV[100];
    std::vector<std::vector<double>>* jet_scaleFactor = nullptr;
    TClonesArray* jet_p4 = nullptr;

    uint32_t run = 0;
    uint32_t eventNumber = 0;
};
