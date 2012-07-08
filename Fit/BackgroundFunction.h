#pragma once

#include <vector>
#include "PDF.h"

#include <RooRealVar.h>
#include <RooExponential.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>

namespace Background {

  class Function {
    public:
      Function(const std::string name):
        mMuPdf(NULL), mEPdf(NULL), mName(name) {
        };

      virtual ~Function() {
        delete mMuPdf;
        delete mEPdf;
        deleteVector(mMuParams);
        deleteVector(mEParams);
      }

      std::string getName() const {
        return mName;
      };

      virtual std::vector<RooRealVar*>& getMuParams() {
        return mMuParams;
      };

      virtual std::vector<RooRealVar*>& getEParams() {
        return mEParams;
      };

      virtual RooAbsPdf& getMuPdf() {
        return *mMuPdf;
      }

      virtual RooAbsPdf& getEPdf() {
        return *mEPdf;
      }

      virtual void createPdfs(RooRealVar& observable) = 0;


    protected:

      template<typename T>
        void deleteVector(std::vector<T*>& vector) {
          for (typename std::vector<T*>::iterator it = vector.begin(); it != vector.end(); ++it) {
            delete *it;
          }
        }

      std::vector<RooRealVar*> mMuParams;
      std::vector<RooRealVar*> mEParams;

      RooAbsPdf* mMuPdf;
      RooAbsPdf* mEPdf;

    private:
      std::string mName;
  };

  class Exp : public Function {
    public:
      Exp()
        : Function("Exp") {}

    protected:
      virtual void createPdfs(RooRealVar& observable) {
        mMuParams.push_back(new RooRealVar("coeff_mu", "coeff_mu", -0.0053, -0.1, 0.));
        mEParams.push_back(new RooRealVar("coeff_e", "coeff_e", -0.0044, -0.1, 0.));

        mMuPdf = new RooExponential("ExpBkg_mu","exponential tail for muon" , observable, *mMuParams[0]);
        mEPdf  = new RooExponential("ExpBkg_e" ,"exponential tail for electrons", observable, *mEParams[0]);
      };
  };

  class TwoExp : public Function {
    public:
      TwoExp()
        : Function("2Exp"), exp1_mu(NULL), exp2_mu(NULL), exp1_e(NULL), exp2_e(NULL) {}
      ~TwoExp() {
        delete exp1_mu;
        delete exp1_e;
        delete exp2_mu;
        delete exp2_e;
      }

    protected:
      virtual void createPdfs(RooRealVar& observable) {
        mMuParams.push_back(new RooRealVar("mixing_mu", "param1_mu", 0.6, 0.2, 1));
        mMuParams.push_back(new RooRealVar("coef1_mu", "param2_mu", -0.0053, -1., 0.));
        mMuParams.push_back(new RooRealVar("coef2_mu", "param3_mu", -0.053, -1., 0.));

        mEParams.push_back(new RooRealVar("mixing_e", "param1_e", 0.5, 0.2, 1));
        mEParams.push_back(new RooRealVar("coef1_e", "param2_e", -0.0053, -1, 0.));
        mEParams.push_back(new RooRealVar("coef2_e", "param3_e", -0.0053, -1, 0.));

        exp1_mu = new RooExponential("exp1_mu", "exp1_mu", observable, *mMuParams[1]);
        exp2_mu = new RooExponential("exp2_mu", "exp2_mu", observable, *mMuParams[2]);
        mMuPdf = new RooAddPdf("twoexp_mu", "2exp_mu", RooArgList(*exp1_mu, *exp2_mu), *mMuParams[0]);

        exp1_e = new RooExponential("exp1_e", "exp1_e", observable, *mEParams[1]);
        exp2_e = new RooExponential("exp2_e", "exp2_e", observable, *mEParams[2]);
        mEPdf = new RooAddPdf("twoexp_e", "2exp_e", RooArgList(*exp1_e, *exp2_e), *mEParams[0]);
      };

    private:
      RooAbsPdf* exp1_mu;
      RooAbsPdf* exp2_mu;
      RooAbsPdf* exp1_e;
      RooAbsPdf* exp2_e;
  };

  class FAlt : public Function {
    public:
      FAlt()
        : Function("FAlt") {}

    protected:
      virtual void createPdfs(RooRealVar& observable) {
        mMuParams.push_back(new RooRealVar("param1_mu", "param1_mu", 1., -50, 200));
        mMuParams.push_back(new RooRealVar("param2_mu", "param2_mu", 20., -50, 200));
        mMuParams.push_back(new RooRealVar("param3_mu", "param3_mu", 5., -50, 200));

        mEParams.push_back(new RooRealVar("param1_e", "param1_e", 1., -50, 200));
        mEParams.push_back(new RooRealVar("param2_e", "param2_e", 20., -50, 200));
        mEParams.push_back(new RooRealVar("param3_e", "param3_e", 5., -50, 200));

        mMuPdf = new FAltPdf("falt_mu", "falt_mu", observable, *mMuParams[0], *mMuParams[1], *mMuParams[2]);
        mEPdf = new FAltPdf("falt_e"," falt_e", observable, *mEParams[0], *mEParams[1], *mEParams[2]);
      };
  };

  class FAltB : public Function {
    public:
      FAltB()
        : Function("FAltB") {}

    protected:
      virtual void createPdfs(RooRealVar& observable) {
        mMuParams.push_back(new RooRealVar("param1_mu", "param1_mu", 80., -50, 200));
        mMuParams.push_back(new RooRealVar("param2_mu", "param2_mu", -4., -50, 200));
        mMuParams.push_back(new RooRealVar("param3_mu", "param3_mu", 1., 0, 100));

        mEParams.push_back(new RooRealVar("param1_e", "param1_e", 80., -50, 200));
        mEParams.push_back(new RooRealVar("param2_e", "param2_e", -4., -50., 200));
        mEParams.push_back(new RooRealVar("param3_e", "param3_e", 1., 0, 100));

        mMuPdf = new FAltBPdf("faltb_mu", "faltb_mu", observable, *mMuParams[0], *mMuParams[1], *mMuParams[2]);
        mEPdf = new FAltBPdf("faltb_e"," faltb_e", observable, *mEParams[0], *mEParams[1], *mEParams[2]);
      };
  };

  class FAltC : public Function {
    public:
      FAltC()
        : Function("FAltC") {}

    protected:
      virtual void createPdfs(RooRealVar& observable) {
        mMuParams.push_back(new RooRealVar("param1_mu", "param1_mu", 40., -50, 200));
        mMuParams.push_back(new RooRealVar("param2_mu", "param2_mu", -1., -50, 200));

        mEParams.push_back(new RooRealVar("param1_e", "param1_e", 40., 50, 200));
        mEParams.push_back(new RooRealVar("param2_e", "param2_e", -1., -50., 200));

        mMuPdf = new FAltCPdf("faltc_mu", "faltc_mu", observable, *mMuParams[0], *mMuParams[1]);
        mEPdf = new FAltCPdf("faltc_e"," faltc_e", observable, *mEParams[0], *mEParams[1]);
      };
  };

  class TwoPow : public Function {
    public:
      TwoPow()
        : Function("2Pow"), pow1_mu(NULL), pow2_mu(NULL), pow1_e(NULL), pow2_e(NULL) {}
      virtual ~TwoPow() {
        delete pow1_mu;
        delete pow1_e;
        delete pow2_mu;
        delete pow2_e;
      }

    protected:
      virtual void createPdfs(RooRealVar& observable) {
        mMuParams.push_back(new RooRealVar("param1_mu", "param1_mu", 0.5, 0, 10));
        mMuParams.push_back(new RooRealVar("param2_mu", "param2_mu", 1., -1., 100.));
        mMuParams.push_back(new RooRealVar("param3_mu", "param3_mu", 1., -1., 100.));

        mEParams.push_back(new RooRealVar("param1_e", "param1_e", 0.3, 0, 10));
        mEParams.push_back(new RooRealVar("param2_e", "param2_e", 1., -1, 100.));
        mEParams.push_back(new RooRealVar("param3_e", "param3_e", 1., -1, 100.));

        pow1_mu = new PowPdf("pow1_mu", "pow1_mu", observable, *mMuParams[1]);
        pow2_mu = new PowPdf("pow2_mu", "pow2_mu", observable, *mMuParams[2]);
        mMuPdf = new RooAddPdf("twopow_mu", "2pow_mu", RooArgList(*pow1_mu, *pow2_mu), *mMuParams[0]);

        pow1_e = new PowPdf("pow1_e", "pow1_e", observable, *mEParams[1]);
        pow2_e = new PowPdf("pow2_e", "pow2_e", observable, *mEParams[2]);
        mEPdf = new RooAddPdf("twopow_e", "2pow_e", RooArgList(*pow1_e, *pow2_e), *mEParams[0]);
      };

    private:
      RooAbsPdf* pow1_mu;
      RooAbsPdf* pow2_mu;
      RooAbsPdf* pow1_e;
      RooAbsPdf* pow2_e;
  };

  class FourPol : public Function {
    public:
      FourPol()
        : Function("4Pol") {}

    protected:
      virtual void createPdfs(RooRealVar& observable) {
        //mMuParams.push_back(new RooRealVar("param1_mu", "param1_mu", 285, 0., 500.));
        //mMuParams.push_back(new RooRealVar("param2_mu", "param2_mu", -0.75,  -2., 2.));
        //mMuParams.push_back(new RooRealVar("param3_mu", "param3_mu", 7.5e-4, 0., 1.));
        //mMuParams.push_back(new RooRealVar("param4_mu", "param4_mu", -3.33e-7, -1., 0.));
        //mMuParams.push_back(new RooRealVar("param5_mu", "param5_mu", 5.46e-11, 0., 1.));
        mMuParams.push_back(new RooRealVar("param1_mu", "param1_mu", 285, 200, 500));
        mMuParams.push_back(new RooRealVar("param2_mu", "param2_mu", -0.75, -2., 0.));
        mMuParams.push_back(new RooRealVar("param3_mu", "param3_mu", 7.5e-4, 0., 1.));
        mMuParams.push_back(new RooRealVar("param4_mu", "param4_mu", -3.33e-7, -1e-4, 0.));
        mMuParams.push_back(new RooRealVar("param5_mu", "param5_mu", 5.65e-11, 0., 1e-9));
        mMuParams.push_back(new RooRealVar("param6_mu", "param6_mu", 5.65e-13, 0., 1e-9));

        mEParams.push_back(new RooRealVar("param1_e", "param1_e", 285, 0., 500.));
        mEParams.push_back(new RooRealVar("param2_e", "param2_e", -0.75,  -2., 2.));
        mEParams.push_back(new RooRealVar("param3_e", "param3_e", 7.5e-4, 0., 1.));
        mEParams.push_back(new RooRealVar("param4_e", "param4_e", -3.33e-7, -1., 0.));
        mEParams.push_back(new RooRealVar("param5_e", "param5_e", 5.46e-11, 0., 1.));

        mMuPdf = new RooPolynomial("fourpol_mu", "4pol_mu", observable, RooArgList(*mMuParams[0], *mMuParams[1], *mMuParams[2], *mMuParams[3], *mMuParams[4], *mMuParams[5]), 0);
        mEPdf = new RooPolynomial("fourpol_e"," 4pol_e", observable, RooArgList(*mEParams[0], *mEParams[1], *mEParams[2], *mEParams[3], *mEParams[4]), 0);
      };
  };



  class FourLaurent : public Function {
    public:
      FourLaurent()
        : Function("FourLaurent") {}

    protected:
      virtual void createPdfs(RooRealVar& observable) {
        mMuParams.push_back(new RooRealVar("param1_mu", "param1_mu", -1.21, 0., 100.));
        mMuParams.push_back(new RooRealVar("param2_mu", "param2_mu", 2.4395e3, 0., 10000.));
        mMuParams.push_back(new RooRealVar("param3_mu", "param3_mu", -7.43e5, -1e7, 0.));

        mEParams.push_back(new RooRealVar("param1_e", "param1_e", -1.21, -100., 100.));
        mEParams.push_back(new RooRealVar("param2_e", "param2_e", 2.4395e3, 0., 10000.));
        mEParams.push_back(new RooRealVar("param3_e", "param3_e", -7.43e5, -1e7, 0.));

        mMuPdf = new LaurentPdf("fourlaurent_mu", "fourlaurent_mu", -3, -6, observable, RooArgList(*mMuParams[0], *mMuParams[1], *mMuParams[2]));
        mEPdf = new LaurentPdf("fourlaurent_e", "fourlaurent_e", -3, -6, observable, RooArgList(*mEParams[0], *mEParams[1], *mEParams[2]));
      };
  };

  static Function* createFunction(const std::string& name, RooRealVar& observable) {

    Function * fct = NULL;

    if (name == "exp") {
      fct = new Exp();
    } else if (name == "falt") {
      fct = new FAlt();
    } else if (name == "faltB") {
      fct = new FAltB();
    } else if (name == "faltC") {
      fct = new FAltC();
    } else if (name == "2Pow") {
      fct = new TwoPow();
    } else if (name == "4Pol") {
      fct = new FourPol();
    } else if (name == "4Lau") {
      fct = new FourLaurent();
    } else if (name == "2Exp") {
      fct = new TwoExp();
    }

    if (fct != NULL) {
      fct->createPdfs(observable);
    }

    return fct;
  }

};
