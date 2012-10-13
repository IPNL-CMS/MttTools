#ifndef KINFIT_H
#define KINFIT_H

/// ROOT includes
#include "TLorentzVector.h"
#include "TMinuit.h"
#include <TROOT.h>
#include <TMath.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>

#include <memory> // for std::shared_ptr

/// modif BDE, reduction nombre de parametre car mtop fixee
static const int ParamNber = 16;
static const int MaxEtaBins=4;



class KinFit
{ 
public:
  static double PzNeutrino(const TLorentzVector& lept, TLorentzVector& neut, const TLorentzVector& bJet);
  
};


#endif

