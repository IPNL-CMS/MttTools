#pragma once

#include "BaseFunction.h"
#include "PDF.h"

#include <RooRealVar.h>
#include <RooExponential.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>

class FAltB : public BaseFunction {
  public:
    FAltB(const std::string& name)
      : BaseFunction(name) {}

  protected:
    virtual bool createPdf(RooRealVar& observable) {

      if (mParameters.size() != 3)
        return false;

      mPdf.reset(new FAltBPdf(mName.c_str(), mName.c_str(), observable, *mParameters["a"], *mParameters["b"], *mParameters["c"]));

      return true;
    }
};

class FAlt : public BaseFunction {
  public:
    FAlt(const std::string& name)
      : BaseFunction(name) {}

  protected:
    virtual bool createPdf(RooRealVar& observable) {

      if (mParameters.size() != 3)
        return false;

      mPdf.reset(new FAltPdf(mName.c_str(), mName.c_str(), observable, *mParameters["a"], *mParameters["b"], *mParameters["c"]));

      return true;
    }
};

class Exp : public BaseFunction {
  public:
    Exp(const std::string& name)
      : BaseFunction(name) {}

  protected:
    virtual bool createPdf(RooRealVar& observable) {

      if (mParameters.size() != 1)
        return false;

      mPdf.reset(new RooExponential(mName.c_str(), mName.c_str() , observable, *mParameters["c"]));

      return true;
    };
};
