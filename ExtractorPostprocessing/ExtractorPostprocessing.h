class ExtractorPostprocessing {
  public:
    bool passMuonSel(float pt, float eta);
    bool passElectronSel(float pt, float eta);

    bool passJetsSel(float firstJetPt, float secondJetPt, float thirdJetPt, float fourthJetPt, bool isRun2012AB);

    bool passExtractorSel(int isSel, int numComb, float mtt);

    bool passChi2Sel(float chi2);
    bool passMVASel(float mva);
};
