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
            "variation_alphas_up": "Systematics/MC_TT_powheg_theta_alphasUp.root",
            "variation_alphas_down": "Systematics/MC_TT_powheg_theta_alphasDown.root",
            "variation_pdf_up": "Systematics/MC_TT_powheg_theta_pdfUp.root",
            "variation_pdf_down": "Systematics/MC_TT_powheg_theta_pdfDown.root",
            "output_up": "Systematics/MC_TT_powheg_theta_alphasPdfUp.root",
            "output_down": "Systematics/MC_TT_powheg_theta_alphasPdfDown.root"
        },

        {
            "nominal": "Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_nominal.root",
            "variation_alphas_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasUp.root",
            "variation_alphas_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasDown.root",
            "variation_pdf_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_pdfUp.root",
            "variation_pdf_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_pdfDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        },
        {
            "nominal": "Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_nominal.root",
            "variation_alphas_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasUp.root",
            "variation_alphas_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasDown.root",
            "variation_pdf_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_pdfUp.root",
            "variation_pdf_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_pdfDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        },
        {
            "nominal": "Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_nominal.root",
            "variation_alphas_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasUp.root",
            "variation_alphas_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasDown.root",
            "variation_pdf_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_pdfUp.root",
            "variation_pdf_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_pdfDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        },
        {
            "nominal": "Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_nominal.root",
            "variation_alphas_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasUp.root",
            "variation_alphas_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasDown.root",
            "variation_pdf_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_pdfUp.root",
            "variation_pdf_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_pdfDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        },
        {
            "nominal": "Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_nominal.root",
            "variation_alphas_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasUp.root",
            "variation_alphas_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasDown.root",
            "variation_pdf_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_pdfUp.root",
            "variation_pdf_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_pdfDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        },

]

type = ["semie", "semimu"]

for type in ["semie", "semimu"]:
    for syst in systs:

        # Open nominal file, as well as up and down variations
        nominal = TFile.Open(os.path.join(options.root, type, syst["nominal"]))
        variation_alphas_up = TFile.Open(os.path.join(options.root, type, syst["variation_alphas_up"]))
        variation_alphas_down = TFile.Open(os.path.join(options.root, type, syst["variation_alphas_down"]))
        variation_pdf_up = TFile.Open(os.path.join(options.root, type, syst["variation_pdf_up"]))
        variation_pdf_down = TFile.Open(os.path.join(options.root, type, syst["variation_pdf_down"]))

        print("Creating %s" % os.path.join(options.root, type, syst["output_up"]))
        output_up = TFile.Open(os.path.join(options.root, type, syst["output_up"]), "recreate")

        print("Creating %s" % os.path.join(options.root, type, syst["output_down"]))
        output_down = TFile.Open(os.path.join(options.root, type, syst["output_down"]), "recreate")

        for hist in nominal.GetKeyNames():
            h_nominal = nominal.Get(hist)
            h_variation_alphas_up = variation_alphas_up.Get(hist)
            h_variation_alphas_down = variation_alphas_down.Get(hist)
            h_variation_pdf_up = variation_pdf_up.Get(hist)
            h_variation_pdf_down = variation_pdf_down.Get(hist)

            h_varied_nominal_up = h_nominal.Clone()
            h_varied_nominal_up.SetDirectory(0)

            h_varied_nominal_down = h_nominal.Clone()
            h_varied_nominal_down.SetDirectory(0)

            assert(h_nominal.GetNbinsX() == h_variation_alphas_up.GetNbinsX() == h_variation_alphas_down.GetNbinsX() == h_variation_pdf_up.GetNbinsX() == h_variation_pdf_down.GetNbinsX())

            for i in range(1, h_nominal.GetNbinsX() + 1):
                nominal_value = h_nominal.GetBinContent(i)
                sigma_alphas_value_up = h_variation_alphas_up.GetBinContent(i) - nominal_value;
                sigma_alphas_value_down = nominal_value - h_variation_alphas_down.GetBinContent(i);
                sigma_pdf_value_up = h_variation_pdf_up.GetBinContent(i) - nominal_value;
                sigma_pdf_value_down = nominal_value - h_variation_pdf_down.GetBinContent(i);

                sigma_up = sqrt( pow(sigma_alphas_value_up, 2) + pow(sigma_pdf_value_up, 2) )
                sigma_down = sqrt( pow(sigma_alphas_value_down, 2) + pow(sigma_pdf_value_down, 2) )

                h_varied_nominal_up.SetBinContent(i, h_varied_nominal_up.GetBinContent(i) + sigma_up)
                h_varied_nominal_down.SetBinContent(i, h_varied_nominal_down.GetBinContent(i) - sigma_down)

            output_up.cd()
            h_varied_nominal_up.Write()
            output_down.cd()
            h_varied_nominal_down.Write()

        output_up.Close()
        output_down.Close()


