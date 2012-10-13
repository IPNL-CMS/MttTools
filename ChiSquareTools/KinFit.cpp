/* **********************************************************************
 *                                                                      *
 * KinFit : kinematical fit for semilept tt bar events             *
 *               based on                                               *
 *               - an external  full jet parametrisation                *
 *               - chi2 minisation                                      *
 *               - computed neutrino Pz (for initialisation)            *
 *               - 2 mtop + mw constraint                               *
 * Djamel BOUMEDIENE   Djamel.Boumedien@cern.ch                         *
 *                                                                      *
 * courtesy for help and code structure :                               *
 *                           Diane Cinca diane.cinca@cern.ch            *
 * TODO : parametrization of Neutrino + Muon + electron from .dat file  *
 * internal tracking :                                                  *
 * $Id: KinFit.cc,v 1.1 2011/11/10 15:09:01 sviret Exp $*
 ************************************************************************/

#include "KinFit.h"

double KinFit::PzNeutrino(const TLorentzVector& lept, TLorentzVector& neut, const TLorentzVector& bJet)
{
  double m_top = 173.5;
  double m_w = 80.385;
  if(!lept.E()) return 0;

  double x = (m_w * m_w - lept.M() * lept.M() + 2. * (neut.Px() * lept.Px() + neut.Py() * lept.Py())) / (2 * lept.E());
  double a = 1 - (lept.Pz() * lept.Pz()) / (lept.E() * lept.E());
  double b = -2. * (lept.Pz() / lept.E()) * x;
  double c = neut.Pt() * neut.Pt() - x * x;

  if (!a && !b) return 0;

  if(!a)
  {     
    neut.SetPz(-1 * c / b);
    neut.SetE(sqrt(neut.Px() * neut.Px() + neut.Py() * neut.Py() + neut.Pz() * neut.Pz()));
    return 1;
  }


  double delta = b * b - 4 * a *c;


  if (delta < 0)     // No solution, try to correct MET
  {
    double rat = neut.Py() / neut.Px();

    double u = 4./(lept.E()*lept.E())*((lept.Px()+rat*lept.Py())*(lept.Px()+rat*lept.Py())/(1+rat*rat)
        -(lept.E()*lept.E())+(lept.Pz()*lept.Pz()));

    double v = 4./(lept.E()*lept.E())*(m_w*m_w-lept.M()*lept.M())
      *(lept.Px()+rat*lept.Py())/sqrt(1+rat*rat);

    double w = (m_w*m_w-lept.M()*lept.M())*(m_w*m_w-lept.M()*lept.M())/(lept.E()*lept.E());

    double deltan = v * v - 4 * u * w;

    if (deltan < 0)
      return 0; // Hopeless, MET can't be corrected

    double pt      = 0.;
    double corfact = 0.;

    if(u==0)
    {
      pt = -w/v;
      if (pt <= 0) return 0; // There is no way out...

      corfact = pt / neut.Pt();
    }
    else // Deltan>=0 and u!=0
    {
      double pt2 = (v-(sqrt(deltan)))/(2*u);
      double pt1 = (v+(sqrt(deltan)))/(2*u);

      // Pas de correction car negative
      if(pt1<=0 && pt2<=0) return 0;

      if(pt1>0 && pt2<0) pt=pt1;
      if(pt2>0 && pt1<0) pt=pt2;
      if(pt1>0 && pt2>0)
      {
        (fabs(pt1-neut.Pt())<=fabs(pt2-neut.Pt()))
          ? pt=pt1
          : pt=pt2;     
      }

      corfact = pt/neut.Pt();
    }

    // Now we have the correction factor

    neut.SetPx(corfact*neut.Px());
    neut.SetPy(corfact*neut.Py());

    // Recompute the new parameters

    x = (m_w*m_w-lept.M()*lept.M()+2.*(neut.Px()*lept.Px()+neut.Py()*lept.Py()))/(2*lept.E());
    a = 1-(lept.Pz()*lept.Pz())/(lept.E()*lept.E());
    b = -2.*(lept.Pz()/lept.E())*x;
    c = neut.Px()*neut.Px()+neut.Py()*neut.Py()-x*x;

//         std::cout << "We have rescaled the MET " << lept->E() << " / " << corfact <<  " , now delta should be null:" << std::endl;
//         std::cout << "Previous delta: " << delta<< std::endl;


    delta= b*b-4*a*c;

    if (fabs(delta)<0.000001) delta=0.;

//    std::cout << "New delta     : " << delta << std::endl;

    if (delta != 0) return 0; // This should not happen, but who knows...
  }


  // We can go back to the normal path: 

  TLorentzVector TopCand1 = lept + bJet;
  TLorentzVector TopCand2 = lept + bJet;

  neut.SetPz((-b-(sqrt(delta)))/(2*a));
  neut.SetE(sqrt(neut.Px()*neut.Px()+neut.Py()*neut.Py()+neut.Pz()*neut.Pz()));
  TopCand1 += neut;

  neut.SetPz((-b+(sqrt(delta)))/(2*a));
  neut.SetE(sqrt(neut.Px()*neut.Px()+neut.Py()*neut.Py()+neut.Pz()*neut.Pz()));
  TopCand2 += neut;

  double mtt_1 = sqrt(std::max(0.,TopCand1.M2()));
  double mtt_2 = sqrt(std::max(0.,TopCand2.M2()));

  if(fabs(mtt_1-m_top) <= fabs(mtt_2-m_top)) // Otherwise it's OK
  {
    neut.SetPz((-b-(sqrt(delta)))/(2*a));
    neut.SetE(sqrt(neut.Px()*neut.Px()+neut.Py()*neut.Py()+neut.Pz()*neut.Pz()));
  }

  return 1;
} 

