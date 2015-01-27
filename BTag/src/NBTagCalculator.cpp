#include "NBTagCalculator.h"

#include <TLorentzVector.h>

uint32_t NBTagCalculator::getNumberOfBTaggedJets(uint32_t entry) {

  for (auto& chain: m_chains)
    chain->GetEntry(entry);

  BTagUtils btagUtil(run, eventNumber, m_btagging_efficiency_provider);

  uint32_t numberOfBTaggedJets = 0;
  // Recompute number of b-tagged jets
  for (uint32_t j = 0; j < n_jets; j++) {

    TLorentzVector* p4 = static_cast<TLorentzVector*>((*jet_p4)[j]);
    Flavor flavor;
    switch (abs(jet_flavor[j])) {
      case 5:
        flavor = B;
        break;

      case 4:
        flavor = C;
        break;

      default:
        flavor = LIGHT;
        break;
    }

    float pt = p4->Pt();
    float eta = fabs(p4->Eta());

    if (pt < 30 || eta > 2.4)
      continue;

    bool isBTagged = jet_CSV[j] > 0.679;

    float scale_factor = (*jet_scaleFactor)[j][0];

    SystVariation btag_efficiency_syst_variation = NOMINAL;
    if (m_systVariation == UP) {
      scale_factor += (*jet_scaleFactor)[j][1];
      btag_efficiency_syst_variation = UP;
    } else if (m_systVariation == DOWN) {
      scale_factor -= (*jet_scaleFactor)[j][2];
      btag_efficiency_syst_variation = DOWN;
    }

    if (btagUtil.updateJetBTagStatus(isBTagged, pt, eta, flavor, scale_factor, btag_efficiency_syst_variation))
      numberOfBTaggedJets++;

  }

  return numberOfBTaggedJets;
}
