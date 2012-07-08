#pragma once

#include <vector>

#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooBreitWigner.h>
#include <RooFFTConvPdf.h>

namespace Signal {

  class Function {
    public:
      Function(const std::string name):
        mMuPdf(NULL), mEPdf(NULL), mName(name) {
        };

      virtual ~Function() {
        delete mMuPdf;
        delete mEPdf;
      }

      std::string getName() const {
        return mName;
      };

      virtual RooAbsPdf& getMuPdf() {
        return *mMuPdf;
      }

      virtual RooAbsPdf& getEPdf() {
        return *mEPdf;
      }

      virtual void createPdfs(RooRealVar& observable, const int mass) = 0;


    protected:
      RooAbsPdf* mMuPdf;
      RooAbsPdf* mEPdf;

    private:
      std::string mName;
  };

  class CrystalBall : public Function {
    public:
      CrystalBall()
        : Function("CB"),
        mZ_mu(NULL), massShift(NULL), mZ_e(NULL), alphaCB_mu(NULL), alphaCB_e(NULL),
        sigmaMZ_mu(NULL), sigmaMZ_e(NULL), nCB_mu(NULL), nCB_e(NULL) {}

      virtual ~CrystalBall() {
        delete mZ_mu;
        delete massShift;
        delete mZ_e;
        delete alphaCB_mu;
        delete alphaCB_e;
        delete sigmaMZ_mu;
        delete sigmaMZ_e;
        delete nCB_mu;
        delete nCB_e;
      }

    protected:
      virtual void createPdfs(RooRealVar& observable, const int mass) {
        mZ_mu = new RooRealVar("mZ_mu", "Zprime mass", mass, observable.getBinning().lowBound(), observable.getBinning().highBound());
        massShift = new RooRealVar("massShiftEMu", "mass shift", 0.9, 0.6, 1.4);
        mZ_e = new RooFormulaVar("mZ_e", "Zprime mass", "mZ_mu*massShiftEMu", RooArgList(*mZ_mu, *massShift));
        alphaCB_mu = new RooRealVar("alphaCB_mu", "alpha Crystal Ball", 1.13, 0.1, 10.);
        alphaCB_e = new RooRealVar("alphaCB_e", "alpha Crystal Ball",  2.8, 0.1, 10.);
        sigmaMZ_mu = new RooRealVar("sigma_MZ_mu", "signal m_TT sigma1", 55., 50., 200.);
        sigmaMZ_e = new RooRealVar("sigma_MZ_e", "signal m_TT sigma1", 55., 50., 200.);
        nCB_mu = new RooRealVar("nCB_mu", "n Crystal Ball", 1., 0.5, 300.);
        nCB_e = new RooRealVar("nCB_e", "n Crystal Ball", 1., 0.5, 300.);

        mMuPdf = new RooCBShape("sigPdf_mu", "Zprime Crystal Ball", observable, *mZ_mu, *sigmaMZ_mu, *alphaCB_mu, *nCB_mu);
        mEPdf = new RooCBShape("sigPdf_e", "Zprime Crystal Ball", observable, *mZ_e, *sigmaMZ_e, *alphaCB_e, *nCB_e);
      };

    private:
      RooRealVar *mZ_mu;
      RooRealVar *massShift;
      RooFormulaVar* mZ_e;
      RooRealVar *alphaCB_mu;
      RooRealVar *alphaCB_e;
      RooRealVar *sigmaMZ_mu;
      RooRealVar *sigmaMZ_e;
      RooRealVar *nCB_mu;
      RooRealVar *nCB_e;
  };

  class BWCrystalBall : public Function {
    public:
      BWCrystalBall()
        : Function("CB"),
        CBmean_mu(NULL), CBmean_e(NULL), alphaCB_mu(NULL), alphaCB_e(NULL),
        sigmaMZ_mu(NULL), sigmaMZ_e(NULL), nCB_mu(NULL), nCB_e(NULL),
        sigmaBW_mu(NULL), sigmaBW_e(NULL), meanBW_mu(NULL), massShift(NULL),
        meanBW_e(NULL), BW_e(NULL), BW_mu(NULL) {}

      virtual ~BWCrystalBall() {
        delete CBmean_mu;
        delete CBmean_e;
        delete alphaCB_mu;
        delete alphaCB_e;
        delete sigmaMZ_mu;
        delete sigmaMZ_e;
        delete nCB_mu;
        delete nCB_e;
        delete CB_mu;
        delete CB_e;

        delete massShift;
        delete sigmaBW_mu;
        delete sigmaBW_e;
        delete meanBW_mu;
        delete meanBW_e;
        delete BW_e;
        delete BW_mu;
      }

    protected:
      virtual void createPdfs(RooRealVar& observable, const int mass) {
        CBmean_mu = new RooRealVar("CBmean_mu", "", 0.);
        CBmean_e = new RooRealVar("CBmean_e", "", 0.);

        alphaCB_mu = new RooRealVar("alphaCB_mu", "alpha Crystal Ball", 1, 0.4, 10.);
        alphaCB_e = new RooRealVar("alphaCB_e", "alpha Crystal Ball", 1, 0.4, 10.);
        sigmaMZ_mu = new RooRealVar("sigma_MZ_mu", "signal m_TT sigma1", 69., 1, 150.);
        sigmaMZ_e = new RooRealVar("sigma_MZ_e", "signal m_TT sigma1", 69., 1., 150.);
        nCB_mu = new RooRealVar("nCB_mu", "n Crystal Ball", 5., 0., 20.);
        nCB_e = new RooRealVar("nCB_e", "n Crystal Ball", 5., 0., 20.);

        CB_mu = new RooCBShape("CB_mu", "Zprime Crystal Ball", observable, *CBmean_mu, *sigmaMZ_mu, *alphaCB_mu, *nCB_mu);
        CB_e = new RooCBShape("CB_e", "Zprime Crystal Ball", observable, *CBmean_e, *sigmaMZ_e, *alphaCB_e, *nCB_e);

        sigmaBW_mu = new RooRealVar("sigmaBW_mu",	"BW sigma (GeV)",	40.,	1,	150);  
        meanBW_mu = new RooRealVar("meanBW_mu", "BW avg (GeV)", mass, observable.getBinning().lowBound(), observable.getBinning().highBound());

        BW_mu = new RooBreitWigner("BW_mu", "Breit Wigner theory", observable, *meanBW_mu, *sigmaBW_mu);

        sigmaBW_e = new RooRealVar("sigmaBW_e", "BW sigma (GeV)",	40.,	1,	150);  
        
        massShift = new RooRealVar("massShift", "mass shift", 0.9, 0.6, 1.4);
        meanBW_e = new RooFormulaVar("meanBW_e", "BW mean (GeV)", "meanBW_mu*massShift", RooArgList(*meanBW_mu, *massShift));

        BW_e = new RooBreitWigner("BW_e", "Breit Wigner theory", observable, *meanBW_e, *sigmaBW_e);

        // Convolution
        observable.setBins(5000, "fft");
        observable.setBins(5000, "cache");

        mMuPdf = new RooFFTConvPdf("sigPdf_mu", "FFT Conv CryBall and BW mu", observable, *BW_mu, *CB_mu);
        mEPdf = new RooFFTConvPdf("sigPdf_e", "FFT Conv CryBall and BW mu", observable, *BW_e, *CB_e);

        static_cast<RooFFTConvPdf*>(mMuPdf)->setBufferFraction(1);
        static_cast<RooFFTConvPdf*>(mEPdf)->setBufferFraction(1);
      };

    private:
      RooRealVar *CBmean_mu;
      RooRealVar *CBmean_e;
      RooRealVar *alphaCB_mu;
      RooRealVar *alphaCB_e;
      RooRealVar *sigmaMZ_mu;
      RooRealVar *sigmaMZ_e;
      RooRealVar *nCB_mu;
      RooRealVar *nCB_e;
      RooRealVar *sigmaBW_mu;
      RooRealVar *sigmaBW_e;
      RooRealVar *meanBW_mu;
      RooRealVar *massShift;
      RooFormulaVar *meanBW_e;

      RooAbsPdf* CB_mu;
      RooAbsPdf* CB_e;
      RooAbsPdf* BW_e;
      RooAbsPdf* BW_mu;
  };

  static Function* createFunction(const std::string& name, RooRealVar& observable, const int mass) {
    Function *fct = NULL;
    
    if (name == "crystalball") {
      fct = new CrystalBall();
    } else if (name == "bwcrystalball") {
      fct = new BWCrystalBall();
    }

    if (fct) {
      fct->createPdfs(observable, mass);
    }

    return fct;
  }
    
};
