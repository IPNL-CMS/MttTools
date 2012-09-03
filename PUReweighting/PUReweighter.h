#include <string>
#include <TH1.h>

class PUReweighter {
  public:
    PUReweighter(bool isSemiMu, const std::string& mcName);

    ~PUReweighter() {
      delete puHisto;
    }

    double weight(float interactions) const;

  private:
    TH1* puHisto;
};

