#pragma once

#include <vector>
#include <string>
#include <memory>

#include <TBranch.h>
#include <TChain.h>
#include <TMVA/Reader.h>

class MVAReader {
  public:
    MVAReader(const std::vector<std::string>& inputFiles, bool isMC);

    virtual void initMVA(const std::string& weightFile);
    virtual float evaluate(uint64_t entry) final;

  protected:

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

    virtual void initTrees() = 0;
    virtual void setupVariables() = 0;

    virtual void computeVariables(uint64_t entry) {
      // Empty
    }

    std::shared_ptr<TMVA::Reader> m_reader;
    bool m_isMC;

  private:
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

    void createReader(const std::string& weightFile);

    std::vector<std::string> m_inputFiles;
    std::string m_weightFile;
    std::vector<std::shared_ptr<TChain>> m_chains;
};
