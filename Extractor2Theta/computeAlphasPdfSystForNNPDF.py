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
            "output_up": "Systematics/MC_TT_powheg_theta_combinedAlphasPdfUp.root",
            "output_down": "Systematics/MC_TT_powheg_theta_combinedAlphasPdfDown.root",
        },

]

set = "NNPDF23_nlo"

for type in ["semie", "semimu"]:
#for type in ["semimu"]:
    for syst in systs:

        # Open nominal file, as well as up and down variations
        nominal = TFile.Open(os.path.join(options.root, type, syst["nominal"]))
        variation_up = TFile.Open(os.path.join(options.root, type, syst["variation_up"]))
        variation_down = TFile.Open(os.path.join(options.root, type, syst["variation_down"]))

        print("Adding histogramms to %s" % os.path.join(options.root, type, syst["output_up"]))
        output_up = TFile.Open(os.path.join(options.root, type, syst["output_up"]), "UPDATE")

        print("Adding histogramms to %s" % os.path.join(options.root, type, syst["output_down"]))
        output_down = TFile.Open(os.path.join(options.root, type, syst["output_down"]), "UPDATE")

        for hist in nominal.GetKeyNames():
            h_nominal = nominal.Get(hist)

            h_varied_nominal_up = h_nominal.Clone()
            h_varied_nominal_up.SetDirectory(0)

            h_varied_nominal_down = h_nominal.Clone()
            h_varied_nominal_down.SetDirectory(0)

            histoName = h_nominal.GetName() + "_" + set
            h_varied_nominal_down.SetNameTitle(histoName, histoName)
            h_varied_nominal_up.SetNameTitle(histoName, histoName)

            h_nominal_NNPDF = h_nominal.Clone()
            h_nominal_NNPDF.SetDirectory(0)
            histoName = h_nominal.GetName() + "_nominal_" + set
            h_nominal_NNPDF.SetNameTitle(histoName, histoName)

            list_h_variation_up = []
            list_h_variation_down = []

            for var in variation_up.GetKeyNames():
                if hist in var:
                    if set in var:
                        h_variation = variation_up.Get(var)
                        assert(h_nominal.GetNbinsX() == h_variation.GetNbinsX())
                        list_h_variation_up.append(h_variation)

            for var in variation_down.GetKeyNames():
                if hist in var:
                    if set in var:
                        h_variation = variation_down.Get(var)
                        assert(h_nominal.GetNbinsX() == h_variation.GetNbinsX())
                        list_h_variation_down.append(h_variation)

            for i in range(1, h_nominal.GetNbinsX() + 1):
                assert(len(list_h_variation_up) == len(list_h_variation_down))
                mean = 0
                sigma = 0
                # first compute mean value of mtt per bin
                for j,h_var in enumerate(list_h_variation_up):
                    variation_up_value = list_h_variation_up[j].GetBinContent(i);
                    variation_down_value = list_h_variation_down[j].GetBinContent(i);
                    mean += variation_up_value + variation_down_value
                mean /= (len(list_h_variation_up) * 2)
                h_nominal_NNPDF.SetBinContent(i, mean)
                # then compute sigma i.e. the variation wrt the mean value
                for j,h_var in enumerate(list_h_variation_up):
                    variation_up_value = list_h_variation_up[j].GetBinContent(i);
                    variation_down_value = list_h_variation_down[j].GetBinContent(i);
                    sigma += pow(variation_up_value - mean, 2)
                    sigma += pow(variation_down_value - mean, 2)
                sigma /= ((len(list_h_variation_up) * 2) - 1)
                sigma = sqrt(sigma)

                h_varied_nominal_up.SetBinContent(i, h_nominal_NNPDF.GetBinContent(i) + sigma)
                h_varied_nominal_down.SetBinContent(i, h_nominal_NNPDF.GetBinContent(i) - sigma)

            output_up.cd()
            h_varied_nominal_up.Write()
            h_nominal_NNPDF.Write()
            output_down.cd()
            h_varied_nominal_down.Write()
            h_nominal_NNPDF.Write()
        output_up.Close()
        output_down.Close()



