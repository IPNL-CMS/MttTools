#include <BTagUtils.h>

#include <tuple>
#include <iostream>

/**
 * Return true is this jet should be considered b-tagged, or false otherwise
 **/
bool BTagUtils::updateJetBTagStatus(bool isBTagged, float pt, float eta, Flavor flavor, float scale_factor, SystVariation syst) {

  std::tuple<double, double, double> eff = m_btagging_efficiency_provider->getEfficiency(flavor, pt, fabs(eta));

  float btag_efficiency = std::get<0>(eff);
  if (syst == UP)
    btag_efficiency += std::get<1>(eff);
  else if (syst == DOWN)
    btag_efficiency -= std::get<2>(eff);

  bool newBTag = isBTagged;

  float randomNumber = m_generator->Uniform(1.);

  if (scale_factor > 1) {

    if (! isBTagged) {
      float mistagPercent = (1. - scale_factor) / (1. - (1. / btag_efficiency) );

      if (randomNumber < mistagPercent) {
        newBTag = true;
      }
    }
  
  } else {

    if (isBTagged && randomNumber > scale_factor) {
      newBTag = false;
    }
  
  }

  return newBTag;
}
