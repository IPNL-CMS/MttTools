#include <cmath>

#include <ExtractorPostprocessing.h>

bool ExtractorPostprocessing::passMuonSel(float pt, float eta) {
  return (pt > 26);
}
bool ExtractorPostprocessing::passElectronSel(float pt, float eta) {
  return ((pt > 30) && (!(fabs(eta) >= 1.442 && fabs(eta) < 1.5660)));
}

bool ExtractorPostprocessing::passJetsSel(float firstJetPt, float secondJetPt, float thirdJetPt, float fourthJetPt, bool isRun2012AB) {

  float firstJetCut = 0, secondJetCut = 0, thirdJetCut = 0, fourthJetCut = 20;
  if (isRun2012AB) {
    firstJetCut = 45;
    secondJetCut = 45;
    thirdJetCut = 45;
  } else {
    firstJetCut = 55;
    secondJetCut = 45;
    thirdJetCut = 35;
  }

  return (firstJetPt > firstJetCut) && (secondJetPt > secondJetCut) && (thirdJetPt > thirdJetCut) && (fourthJetPt > fourthJetCut);
}

bool ExtractorPostprocessing::passExtractorSel(int isSel, int numComb, float mtt) {
  return (isSel == 1) && (numComb > 0) && (mtt > 0);
}

bool ExtractorPostprocessing::passChi2Sel(float chi2) {
  return (chi2 < 500);
}

bool ExtractorPostprocessing::passMVASel(float mva) {
  return true;
}
