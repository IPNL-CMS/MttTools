#pragma once

#include "BTaggingEfficiencyProvider.h"

#include <TRandom2.h>

enum SystVariation {
  NOMINAL = 0,
  UP = 1,
  DOWN = -1
};

class BTagUtils {
  public:
    BTagUtils(uint32_t runNumber, uint32_t eventNumber, std::shared_ptr<BTaggingEfficiencyProvider>& btag_efficiency_provider):
      m_btagging_efficiency_provider(btag_efficiency_provider) {
        // Create unique seed for each event, and reproductible on each reprocessing
        uint64_t seed = runNumber + eventNumber;

        // FIXME: Can result in precision lost. Is it a problem?
        m_generator = std::make_shared<TRandom2>((uint32_t) seed);
      }

    bool updateJetBTagStatus(bool isBTagged, float pt, float eta, Flavor flavor, float scale_factor, SystVariation syst);

  private:
    std::shared_ptr<TRandom2> m_generator;
    std::shared_ptr<BTaggingEfficiencyProvider> m_btagging_efficiency_provider;
};
