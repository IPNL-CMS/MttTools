#pragma once

#include "BaseFunction.h"
#include "PDF.h"

#include <RooRealVar.h>
#include <RooExponential.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include <RooGenericPdf.h>

class FAltB : public BaseFunction {
  public:
    FAltB(const std::string& name)
      : BaseFunction(name) {}

  protected:
    virtual bool createPdf(RooRealVar& observable) {

      if (mParameters.size() != 3)
        return false;

      //mPdf.reset(new FAltBPdf(mName.c_str(), mName.c_str(), observable, *mParameters["a"], *mParameters["b"], *mParameters["c"]));
      mPdf.reset(new RooGenericPdf(mName.c_str(), mName.c_str(), "pow((1.-(@0/8000)+@3*pow((@0/8000),2)),@1)/pow(@0,@2)", RooArgList(observable, *mParameters["a"], *mParameters["b"], *mParameters["c"])));

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

      //mPdf.reset(new FAltPdf(mName.c_str(), mName.c_str(), observable, *mParameters["a"], *mParameters["b"], *mParameters["c"]));
      mPdf.reset(new RooGenericPdf(mName.c_str(), mName.c_str(), "pow((1. - (@0 / 8000.)), @1) / pow((@0 / 8000.), (@2 + @3 * log(@0 / 8000.)))", RooArgList(observable, *mParameters["a"], *mParameters["b"], *mParameters["c"])));

      return true;
    }
};

class FAltC : public BaseFunction {
  public:
    FAltC(const std::string& name)
      : BaseFunction(name) {}

  protected:
    virtual bool createPdf(RooRealVar& observable) {

      if (mParameters.size() != 2)
        return false;

      mPdf.reset(new FAltCPdf(mName.c_str(), mName.c_str(), observable, *mParameters["a"], *mParameters["b"]));

      return true;
    }
};

class UFO : public BaseFunction {
  public:
    UFO(const std::string& name)
      : BaseFunction(name) {}

  protected:
    virtual bool createPdf(RooRealVar& observable) {

      if (mParameters.size() != 2)
        return false;

      //mPdf.reset(new UFOPdf(mName.c_str(), mName.c_str(), observable, *mParameters["a"], *mParameters["b"]));
      mPdf.reset(new RooGenericPdf(mName.c_str(), mName.c_str(), "1 / (1. + exp((@0 / 8000. - @1) / @2))", RooArgList(observable, *mParameters["a"], *mParameters["b"])));

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

class GammaPlusLogNormal : public BaseFunction {
  public:
    GammaPlusLogNormal(const std::string& name)
      : BaseFunction(name) {}

  protected:
    virtual bool createPdf(RooRealVar& observable) {

      if (mParameters.size() != 5)
        return false;

      // alpha: @0
      // sigma: @1
      // theta: @2
      // shift: @3
      // sum: @4
      // mtt: @5

      TString constrained_mu = "log( (@0 - 1) * @2 ) + @1 * @1";
      TString lognormal = TString::Format("ROOT::Math::lognormal_pdf(@5, %s, @1, @3)", constrained_mu.Data());
      TString gamma = "ROOT::Math::gamma_pdf(@5, @0, @2, @3)";
      TString formula = TString::Format("@4 * %s + (1 - @4) * %s", lognormal.Data(), gamma.Data());

      // See GammaPlusLogNormalPdf code for a more clean implementation
      mPdf.reset(new RooGenericPdf(mName.c_str(), mName.c_str(),
            formula,
            RooArgList(*mParameters["alpha"], *mParameters["sigma"], *mParameters["theta"], *mParameters["shift"], *mParameters["sum"], observable)));

      return true;
    }
};

class Gamma : public BaseFunction {
  public:
    Gamma(const std::string& name)
      : BaseFunction(name) {}

  protected:
    virtual bool createPdf(RooRealVar& observable) {

      if (mParameters.size() != 3)
        return false;

      mPdf.reset(new GammaPdf(mName.c_str(), mName.c_str(), observable,
            // Gamma PDF
            *mParameters["alpha"], *mParameters["theta"],
            // Various
            *mParameters["shift"]));

      return true;
    }
};

