//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May  4 09:46:33 2012 by ROOT version 5.32/00
// from TTree Mtt/Analysis info
// found on file: /gridgroup/cms/brochet/MTT/CMSSW/CMSSW_5_2_3_patch4/src/Extractors/PatExtractor_2/analysis/tuples/MTT_TTJets_2012_v1_semimu.root
//////////////////////////////////////////////////////////

#pragma once

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include <iostream>

#include "../PUReweighting/PUReweighter.h"

class TBranch;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
//
#define MAX_ARRAY_SIZE 100

class createBDTTrees {
  public :
    TChain          *fMTT;   //!pointer to the analyzed TChain or TChain
    TChain          *fEvent;
    TChain          *fMuons;
    TChain          *fElectrons;
    TChain          *fMET;
    TChain          *fJet;
    Int_t           fCurrent; //!current Tree number in a TChain

    std::string     mDataset;
    bool            mIsMC;
    bool            mIsSemiMu;

    bool            eventIsAssociable;

    // Declaration of leaf types
    Int_t           MC_channel;
    Int_t           isSel;

    float           n_trueInteractions;

    int32_t         MC_leptonIndex;
    int32_t         MC_neutrinoIndex;

    int32_t         MC_leptonicBIndex;
    int32_t         MC_hadronicBIndex;
    int32_t         MC_hadronicFirstJetIndex;
    int32_t         MC_hadronicSecondJetIndex;

    std::string     mOutputFile;

    int             mBTag;
    PUProfile       mPUProfile;

    int             n_muons;
    TClonesArray*   muon_p4;

    int             n_electrons;
    TClonesArray*   electron_p4;

    TClonesArray*   jet_p4;
    uint32_t        n_jets;
    int32_t         jet_mc_index[100];

    TClonesArray*   met_p4;

    createBDTTrees(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isSemiMu, bool isMC, int btag, PUProfile puProfile);

    virtual ~createBDTTrees();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual void     Loop();
    virtual void     Init();
    void SetBranchAddress(TTree* t, const char* branchName, void* ptr, TBranch** branch);
    bool jetComesFromTTDecay(int mcIndex) const;
    bool jetPassSelection(uint32_t jetIndex);
    bool computeNeutrinoPz(const TLorentzVector& bJet, const TLorentzVector& lepton, TLorentzVector& neutrino) const;

    // BDT Trees
    float leptonic_B_Pt;
    float hadronic_B_Pt;
    float lightJet1p2_Pt;

    float leptonic_W_Pt;
    float hadronic_W_Pt;
    float leptonic_W_M;
    float hadronic_W_M;

    float leptonic_Top_Pt;
    float hadronic_Top_Pt;

    float leptonic_Top_M;
    float hadronic_Top_M;

    float delta_phi_tops;
    float delta_phi_lightjets;
    float delta_phi_W;

    float delta_R_tops;
    float delta_R_lightjets;
    float delta_R_W;

    float signal_weight;
    float background_weight;
};
