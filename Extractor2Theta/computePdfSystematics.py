#! /usr/bin/env python

from ROOT import TFile, gDirectory

def GetKeyNames( self, dir = "" ):
    self.cd(dir)
    return [key.GetName() for key in gDirectory.GetListOfKeys()]
TFile.GetKeyNames = GetKeyNames

import os
from optparse import OptionParser
from math import sqrt

parser = OptionParser()
parser.add_option("-p", "--path", dest="root",
                  help="root directory where ROOT file are located. A directory structure with 'semimu', 'semie', 'data' is expected")

(options, args) = parser.parse_args()

if not options.root:
    parser.error("--path argument is required")

systs = [

        {
            "nominal": "MC_TT_powheg_theta_nominal.root",
            "variation_up": "Systematics/MC_TT_powheg_theta_pdfOnlyUp.root",
            "variation_down": "Systematics/MC_TT_powheg_theta_pdfOnlyDown.root",
            "output_up": "Systematics/MC_TT_powheg_theta_pdfUp.root",
            "output_down": "Systematics/MC_TT_powheg_theta_pdfDown.root",
            "sets": { # associate to each set the factor Cx by which we have to divide 
                # an uncertainty interval with a confidence level equal to X%,
                # to obtain the corresponding 68%C.L. value.
                # see http://www.hep.ucl.ac.uk/pdf4lhc/PDF4LHC_practical_guide.pdf
                "CT10nlo": 1.64485, 
                "MSTW2008nlo68cl": 1., 
                "MSTW2008nnlo68cl": 1.
                }

            },
#        {
            #"nominal": "MC_TT_madgraph_theta_nominal.root",
            #"variation_up": "Systematics/MC_TT_madgraph_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/MC_TT_madgraph_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/MC_TT_madgraph_theta_pdfUp.root",
            #"output_down": "Systematics/MC_TT_madgraph_theta_pdfDown.root"
        #},
#        {
            #"nominal": "Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_pdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_pdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_pdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_pdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_pdfDown.root"
        #},

        #{
            #"nominal": "Signal_S0_S_i_M400_cpl1_scalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_pdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M500_cpl1_scalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_pdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M600_cpl1_scalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_pdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M700_cpl1_scalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_pdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M800_cpl1_scalar_theta_nominal.root",
            #"variation_up": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_pdfOnlyUp.root",
            #"variation_down": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_pdfOnlyDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_pdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_pdfDown.root"
        #},

]

for type in ["semie", "semimu"]:
#for type in ["semimu"]:
    for syst in systs:

        # Open nominal file, as well as up and down variations
        nominal = TFile.Open(os.path.join(options.root, type, syst["nominal"]))
        variation_up = TFile.Open(os.path.join(options.root, type, syst["variation_up"]))
        variation_down = TFile.Open(os.path.join(options.root, type, syst["variation_down"]))

        print("Creating %s" % os.path.join(options.root, type, syst["output_up"]))
        output_up = TFile.Open(os.path.join(options.root, type, syst["output_up"]), "recreate")

        print("Creating %s" % os.path.join(options.root, type, syst["output_down"]))
        output_down = TFile.Open(os.path.join(options.root, type, syst["output_down"]), "recreate")

        for hist in nominal.GetKeyNames():
            h_nominal = nominal.Get(hist)

            h_varied_nominal_up = h_nominal.Clone()
            h_varied_nominal_up.SetDirectory(0)

            h_varied_nominal_down = h_nominal.Clone()
            h_varied_nominal_down.SetDirectory(0)

            for set, C_ in syst["sets"].items():
                list_h_variation_up = []
                list_h_variation_down = []
                h_varied_nominal_down.Reset()
                h_varied_nominal_up.Reset()
                histoName = h_nominal.GetName() + "_" + set
                h_varied_nominal_down.SetNameTitle(histoName, histoName)
                h_varied_nominal_up.SetNameTitle(histoName, histoName)

                for var in variation_up.GetKeyNames():
                    if hist in var and set in var:
                        h_variation = variation_up.Get(var)
                        assert(h_nominal.GetNbinsX() == h_variation.GetNbinsX())
                        list_h_variation_up.append(h_variation)

                for var in variation_down.GetKeyNames():
                    if hist in var and set in var:
                        h_variation = variation_down.Get(var)
                        assert(h_nominal.GetNbinsX() == h_variation.GetNbinsX())
                        list_h_variation_down.append(h_variation)

                for i in range(1, h_nominal.GetNbinsX() + 1):
                    nominal_value = h_nominal.GetBinContent(i)
                    assert(len(list_h_variation_up) == len(list_h_variation_down))
                    sum_up = 0
                    sum_down = 0
                    for j,h_var in enumerate(list_h_variation_up):
                        variation_up_value = list_h_variation_up[j].GetBinContent(i);
                        variation_down_value = list_h_variation_down[j].GetBinContent(i);
                        sum_up += pow(max(variation_up_value - nominal_value, variation_down_value - nominal_value, 0), 2)
                        sum_down += pow(max(nominal_value - variation_up_value, nominal_value - variation_down_value, 0), 2)
                    h_varied_nominal_up.SetBinContent(i, h_nominal.GetBinContent(i) + 1/C_ * sqrt(sum_up))
                    h_varied_nominal_down.SetBinContent(i, h_nominal.GetBinContent(i) - 1/C_ * sqrt(sum_down))

                output_up.cd()
                h_varied_nominal_up.Write()
                output_down.cd()
                h_varied_nominal_down.Write()
        output_up.Close()
        output_down.Close()

