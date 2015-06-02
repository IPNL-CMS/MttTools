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

#include <set>
#include <iostream>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
//
#define MAX_ARRAY_SIZE 100

class PreSkim: public TObject {
  public :
    TChain          *fMTT;   //!pointer to the analyzed TChain or TChain
    TChain          *fVertices;
    TChain          *fEvent;
    TChain          *fJet;
    TChain          *fMet;

    TChain          *fMC;
    TChain          *fElectron;
    TChain          *fElectronLoose;
    TChain          *fMuon;
    TChain          *fMuonLoose;
    Int_t           fCurrent; //!current Tree number in a TChain

    std::string     mDataset;
    bool            mIsMC;
    bool            mDoPDFWeight;
    bool            mDoPDFWeightForSignal;
    bool            mIsSemiMu;

    Float_t         MC_boost_tt;

    float           MC_pt_tt;
    float           MC_eta_tt;
    //float           MC_pt_tt_com;
    //float           MC_eta_tt_com;

    //float           MC_top1Pt;
    //float           MC_top2Pt;
    //float           MC_top1Pt_com;
    //float           MC_top2Pt_com;

    //float           MC_top1Eta;
    //float           MC_top2Eta;
    //float           MC_top1Eta_com;
    //float           MC_top2Eta_com;

    // Declaration of leaf types
    Int_t           MC_channel;
    Float_t         MC_mtt;
    Int_t           nGoodMuons;
    Float_t         muonPt[MAX_ARRAY_SIZE];   //[nGoodMuons]
    Float_t         muonEta[MAX_ARRAY_SIZE];   //[nGoodMuons]
    Int_t           nGoodElectrons;
    Float_t         electronPt[MAX_ARRAY_SIZE];   //[nGoodElectrons]
    Float_t         electronEta[MAX_ARRAY_SIZE];   //[nGoodElectrons]
    Float_t         p_1stjetpt;
    Float_t         p_2ndjetpt;
    Float_t         p_3rdjetpt;
    Float_t         p_4thjetpt;
    Int_t           nJets;
    Float_t         jetEta[MAX_ARRAY_SIZE];   //[nJets]
    Float_t         jetPt[MAX_ARRAY_SIZE];   //[nJets]
    Int_t           nBtaggedJets_CSVM;
    Float_t         MET;
    Int_t           isSel;
    Float_t         bestSolChi2;
    Int_t           numComb_chi2;
    Int_t           numComb_MVA;
    Float_t         mLepTop_AfterChi2;
    Float_t         mHadTop_AfterChi2;
    Float_t         mtt_AfterChi2;
    Float_t         lepTopPt_AfterChi2;
    Float_t         lepTopEta_AfterChi2;
    Float_t         hadTopPt_AfterChi2;
    Float_t         hadTopEta_AfterChi2;
    Float_t         pt_tt_AfterChi2;
    Float_t         eta_tt_AfterChi2;
    Float_t         beta_tt_AfterChi2;
    int             selectedLeptonIndex_AfterChi2;

    Float_t         pdf_x1;
    Float_t         pdf_x2;
    Int_t           pdf_id1;
    Int_t           pdf_id2;
    Float_t         pdf_scale;
    int             m_Nmembers_pdf;
    int             m_N_error_pdf;
    float           m_pdf_weight_up[100];
    float           m_pdf_weight_down[100];
    float           m_alphas_weight_up;
    float           m_alphas_weight_down;

    bool            m_triggerPassed;

    TClonesArray* gen_top1_p4;
    TClonesArray* gen_top2_p4;

    TClonesArray* gen_lepton_p4;
    TClonesArray* gen_neutrino_p4;

    TClonesArray* gen_leptonic_B_p4;
    TClonesArray* gen_hadronic_B_p4;

    TClonesArray* gen_lightJet1_p4;
    TClonesArray* gen_lightJet2_p4;

    uint32_t             n_vertices;
    float           n_trueInteractions;
    float           m_lepton_weight;
    float           m_btag_weight;
    float           generator_weight;


    std::string     mOutputFile;

    TClonesArray*   jet_p4;
    TClonesArray*   met_p4;

    PreSkim(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isSemiMu, bool isMC, bool doPDFWeight, bool doPDFWeightForSignal);

    virtual ~PreSkim();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual void     Loop();
    virtual void     Init();

    std::set<TChain*> m_chains;

    template<typename T>
      void setBranchAddress(TChain* chain, const std::string& name, T& address) {

        if (! m_chains.count(chain))
          m_chains.insert(chain);

        TBranch* branch = chain->FindBranch(name.c_str());
        if (!branch)
          return;

        chain->SetBranchAddress(name.c_str(), &address);
      }

    void loadChain(const std::vector<std::string>& inputFiles, const std::string& treeName, TChain*& output);
};
