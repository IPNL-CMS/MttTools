#include "PUReweighter.h"

#include <TFile.h>
#include <iostream>

PUReweighter::PUReweighter(bool isSemiMu, const std::string& mcName):
  puHisto(NULL) {
    const std::string path = "../PUReweighting/";
    const std::string dataFileName = path + ((isSemiMu) ? "Mu_pileup_profile_2012.root" : "Ele_pileup_profile_2012.root");

    TString mcFileName = TString::Format((path + "summer12_computed_mc_%s_pu_truth_70bins.root").c_str(), mcName.c_str());

    TFile* dataFile = TFile::Open(dataFileName.c_str());
    TFile* mcFile = TFile::Open(mcFileName);

    if (! dataFile) {
      std::cerr << "Error: can't open " << dataFileName << ". No PU reweighting." << std::endl;
      return;
    }

    if (! mcFile) {
      std::cerr << "Error: can't open " << mcFileName << ". No PU reweighting." << std::endl;
      return;
    }

    TH1* dataHisto = static_cast<TH1*>(dataFile->Get("pileup"));
    TH1* mcHisto = static_cast<TH1*>(mcFile->Get("pileup"));

    //TODO: Check for NULL ptr

    // Normalize
    dataHisto->Scale(1.0 / dataHisto->Integral());
    mcHisto->Scale(1.0 / mcHisto->Integral());

    // MC * data / MC = data, so the weights are data/MC:
    puHisto = static_cast<TH1*>(dataHisto->Clone());
    puHisto->Divide(mcHisto);
    puHisto->SetDirectory(0); // "detach" the histo from the file

    std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;

    int NBins = puHisto->GetNbinsX();

    for (int ibin = 1; ibin < NBins + 1; ++ibin) {
      std::cout << "   " << ibin - 1 << " " << puHisto->GetBinContent(ibin) << std::endl;
    }

    dataFile->Close();
    mcFile->Close();

    delete dataFile;
    delete mcFile;
  }

double PUReweighter::weight(float interactions) const {
  if (!puHisto) {
    return 1.;
  }

  int bin = puHisto->GetXaxis()->FindBin(interactions);
  return puHisto->GetBinContent(bin);
} 
