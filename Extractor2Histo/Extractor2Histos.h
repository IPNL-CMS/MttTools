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

#include <iostream>

class TBranch;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
//
#define MAX_ARRAY_SIZE 100

class Extractor2Histos {
  public :
    TChain          *fMTT;   //!pointer to the analyzed TChain or TChain
    TChain          *fVertices;
    TChain          *fEvent;
    TChain          *fLooseMuons;
    TChain          *fJet;
    Int_t           fCurrent; //!current Tree number in a TChain

    std::string     mDataset;
    bool            mIsMC;
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
    Int_t           MC_nPU;
    Int_t           nGoodMuons;
    Int_t           nLooseGoodMuons;
    Float_t         muonPt[MAX_ARRAY_SIZE];   //[nGoodMuons]
    Float_t         muonEta[MAX_ARRAY_SIZE];   //[nGoodMuons]
    Float_t         p_2DDrMin[MAX_ARRAY_SIZE];   //[nGoodMuons]
    Float_t         p_2DpTrel[MAX_ARRAY_SIZE];   //[nGoodMuons]
    Float_t         muRelIso[MAX_ARRAY_SIZE];   //[nGoodMuons]
    Int_t           nGoodElectrons;
    Float_t         electronPt[MAX_ARRAY_SIZE];   //[nGoodElectrons]
    Float_t         electronEta[MAX_ARRAY_SIZE];   //[nGoodElectrons]
    Float_t         elRelIso[MAX_ARRAY_SIZE];   //[nGoodElectrons]
    Int_t           hyperTight1MC[MAX_ARRAY_SIZE];   //[nGoodElectrons]
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
    Int_t           oneMatchedCombi;
    Float_t         bestSolChi2;
    Int_t           isBestSolMatched;
    Float_t         KFChi2;
    Int_t           numComb;
    Float_t         solChi2[MAX_ARRAY_SIZE];   //[numComb]
    Float_t         mLepTop_AfterChi2;
    Float_t         mHadTop_AfterChi2;
    Float_t         mtt_AfterChi2;
    Float_t         mLepTop_AfterChi2andKF;
    Float_t         mHadTop_AfterChi2andKF;
    Float_t         mtt_AfterChi2andKF;
    Float_t         lepTopPt_AfterChi2;
    Float_t         lepTopEta_AfterChi2;
    Float_t         hadTopPt_AfterChi2;
    Float_t         hadTopEta_AfterChi2;
    Float_t         pt_tt_AfterChi2;
    Float_t         eta_tt_AfterChi2;
    Float_t         beta_tt_AfterChi2;

    bool            m_triggerPassed;

    // List of branches
    TBranch        *b_MC_channel;   //!
    TBranch        *b_MC_mtt;   //!
    TBranch        *b_m_nPU;   //!
    TBranch        *b_nGoodMuons;   //!
    TBranch        *b_nLooseGoodMuons;   //!
    TBranch        *b_muonPt;   //!
    TBranch        *b_2DDrMin;   //!
    TBranch        *b_2DpTrel;   //!
    TBranch        *b_muRelIso;   //!
    TBranch        *b_nGoodElectrons;   //!
    TBranch        *b_electronPt;   //!
    TBranch        *b_elRelIso;   //!
    TBranch        *b_hyperTight1MC;   //!
    TBranch        *b_1stjetpt;   //!
    TBranch        *b_2ndjetpt;   //!
    TBranch        *b_3rdjetpt;   //!
    TBranch        *b_4thjetpt;   //!
    TBranch        *b_nJets;   //!
    TBranch        *b_jetEta;   //!
    TBranch        *b_jetPt;   //!
    TBranch        *b_nBtaggedJets_TCHEL;   //!
    TBranch        *b_nBtaggedJets_TCHEM;   //!
    TBranch        *b_nBtaggedJets_TCHET;   //!
    TBranch        *b_nBtaggedJets_TCHPL;   //!
    TBranch        *b_nBtaggedJets_TCHPM;   //!
    TBranch        *b_nBtaggedJets_TCHPT;   //!
    TBranch        *b_nBtaggedJets_SSVHEM;   //!
    TBranch        *b_nBtaggedJets_SSVHPT;   //!
    TBranch        *b_MET;   //!
    TBranch        *b_isSel;   //!
    TBranch        *b_oneMatchedCombi;   //!
    TBranch        *b_bestSolChi2;   //!
    TBranch        *b_isBestSolMatched;   //!
    TBranch        *b_KFChi2;   //!
    TBranch        *b_numComb;   //!
    TBranch        *b_solChi2;   //!
    TBranch        *b_mLepTop_AfterChi2;   //!
    TBranch        *b_mHadTop_AfterChi2;   //!
    TBranch        *b_mtt_AfterChi2;   //!
    TBranch        *b_mLepTop_AfterChi2andKF;   //!
    TBranch        *b_mHadTop_AfterChi2andKF;   //!
    TBranch        *b_mtt_AfterChi2andKF;   //!

    int             n_vertices;
    float           n_trueInteractions;
    float           m_weight;

    std::string     mOutputFile;

    int             mBTag;

    int             n_muons;
    float           muon_relIso[100];
    TClonesArray*   muon_p4;

    TClonesArray*   jet_p4;

    Extractor2Histos(const std::vector<std::string>& inputFiles, const std::string& outputFile, bool isSemiMu, bool isMC, int btag);

    virtual ~Extractor2Histos();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual void     Loop();
    virtual void     Init();
    void SetBranchAddress(TTree* t, const char* branchName, void* ptr, TBranch** branch);
};
