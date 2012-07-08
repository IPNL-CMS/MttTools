void datasetToHisto(const std::string& file) {

  //fit region
  Float_t minmTT = 500;
  Float_t maxmTT = 2000;
  Int_t nBins = 100;

  RooRealVar Mtt_KF_reco("Mtt_KF_reco", "Mtt_KF_reco", minmTT, maxmTT, "GeV/c^2");
  RooCategory whichLepton("whichLepton", "whichLepton");
  whichLepton.defineType("electron", 11) ;
  whichLepton.defineType("muon", 13) ;

  RooDataSet *dataset = RooDataSet::read(file.c_str(), RooArgList(Mtt_KF_reco, whichLepton));

  TH1* histo_semimu = dataset->createHistogram("histo_semimu", Mtt_KF_reco, RooFit::Binning(nBins, minmTT, maxmTT), RooFit::Cut("whichLepton==13"));
  TH1* histo_semie = dataset->createHistogram("histo_semie", Mtt_KF_reco, RooFit::Binning(nBins, minmTT, maxmTT), RooFit::Cut("whichLepton==11"));

  TFile *f = TFile::Open(std::string(file + ".root").c_str(), "recreate");
  histo_semimu->Write();
  histo_semie->Write();
  f->Close();

  delete f;

  delete histo_semimu;
  delete histo_semie;
}
