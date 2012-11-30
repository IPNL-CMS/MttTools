#pragma once

#include "BaseFunction.h"

#include <RooCBShape.h>
#include <RooBreitWigner.h>
#include <RooFFTConvPdf.h>
#include <RooKeysPdf.h>

class KeysPdf : public BaseFunction {
  public:
    KeysPdf(const std::string& name)
      : BaseFunction(name) {
        mIsExtended = true;
      }

    virtual bool createPdf(RooRealVar& observable) {
      if (! mHasDataset)
        return false;

      mPdf.reset(new RooKeysPdf(mName.c_str(), "Keys pdf for signal", observable, *mDataset, RooKeysPdf::MirrorBoth, 2));
      return true;
    }

};

class CrystalBall : public BaseFunction {
  public:
    CrystalBall(const std::string& name)
      : BaseFunction(name) {}

    virtual bool createPdf(RooRealVar& observable) {
      if (mParameters.size() != 4)
        return false;

      mPdf.reset(new RooCBShape(mName.c_str(), "Zprime Crystal Ball", observable, *mParameters["mean"].get(), *mParameters["sigma"].get(), *mParameters["alpha"].get(), *mParameters["n"].get()));

      return true;
    };
};

class BWCrystalBall : public BaseFunction {
  public:
    BWCrystalBall(const std::string& name)
      : BaseFunction(name) {}

    virtual bool createPdf(RooRealVar& observable) {

      if (mParameters.size() != 6)
        return false;

      mCrystalBall.reset(new RooCBShape(std::string(mName + "_crystalball").c_str(), "Zprime Crystal Ball", observable, *mParameters["cb_mean"].get(), *mParameters["cb_sigma"].get(), *mParameters["cb_alpha"].get(), *mParameters["cb_n"].get()));


      mConvolution.reset(new RooBreitWigner(std::string(mName + "_breitwigner").c_str(), "Breit-Wigner", observable, *mParameters["bw_mean"].get(), *mParameters["bw_sigma"].get()));

      // Convolution
      observable.setBins(5000, "fft");
      observable.setBins(5000, "cache");

      mPdf.reset(new RooFFTConvPdf(std::string(mName + "_convolution").c_str(), "FFT Conv CryBall and BW", observable, *mCrystalBall.get(), *mConvolution.get()));

      static_cast<RooFFTConvPdf*>(mPdf.get())->setBufferFraction(1);

      return true;
    }

  private:
    std::shared_ptr<RooCBShape> mCrystalBall;
    std::shared_ptr<RooAbsPdf> mConvolution;
};
