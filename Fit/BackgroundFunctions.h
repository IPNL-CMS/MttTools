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
