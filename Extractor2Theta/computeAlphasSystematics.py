#! /usr/bin/env python

from ROOT import TFile, gDirectory

def GetKeyNames( self, dir = "" ):
    self.cd(dir)
    return [key.GetName() for key in gDirectory.GetListOfKeys()]
TFile.GetKeyNames = GetKeyNames

import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--path", dest="root",
                  help="root directory where ROOT file are located. A directory structure with 'semimu', 'semie', 'data' is expected")

(options, args) = parser.parse_args()

if not options.root:
    parser.error("--path argument is required")

systs = [
        
        {
            "nominal": "MC_TT_powheg_theta_nominal.root",
            "variation_up": "Systematics/MC_TT_powheg_theta_alphasOnlyUp.root",
            "variation_down": "Systematics/MC_TT_powheg_theta_alphasOnlyDown.root",
            "output_up": "Systematics/MC_TT_powheg_theta_alphasUp.root",
            "output_down": "Systematics/MC_TT_powheg_theta_alphasDown.root"
        },
        {
            "nominal": "Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_nominal.root",
            "variation_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasOnlyUp.root",
            "variation_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasOnlyDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasDown.root"
        },
        {
            "nominal": "Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_nominal.root",
            "variation_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasOnlyUp.root",
            "variation_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasOnlyDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasDown.root"
        },
        {
            "nominal": "Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_nominal.root",
            "variation_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasOnlyUp.root",
            "variation_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasOnlyDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasDown.root"
        },
        {
            "nominal": "Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_nominal.root",
            "variation_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasOnlyUp.root",
            "variation_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasOnlyDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasDown.root"
        },
        {
            "nominal": "Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_nominal.root",
            "variation_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasOnlyUp.root",
            "variation_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasOnlyDown.root",
            "output_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasUp.root",
            "output_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasDown.root"
        },

]

type = ["semie", "semimu"]

for type in ["semie", "semimu"]:
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
            h_variation_up = variation_up.Get(hist)
            h_variation_down = variation_down.Get(hist)

            h_varied_nominal_up = h_nominal.Clone()
            h_varied_nominal_up.SetDirectory(0)

            h_varied_nominal_down = h_nominal.Clone()
            h_varied_nominal_down.SetDirectory(0)

            assert(h_nominal.GetNbinsX() == h_variation_up.GetNbinsX() == h_variation_down.GetNbinsX())

            for i in range(1, h_nominal.GetNbinsX() + 1):
                nominal_value = h_nominal.GetBinContent(i)
                variation_value_up = h_variation_up.GetBinContent(i);
                variation_value_down = h_variation_down.GetBinContent(i);

                sigma_up = (variation_value_up - nominal_value) * 5 / 6
                sigma_down = (variation_value_down - nominal_value) * 5 / 6
                h_varied_nominal_up.SetBinContent(i, h_varied_nominal_up.GetBinContent(i) + sigma_up)
                h_varied_nominal_down.SetBinContent(i, h_varied_nominal_down.GetBinContent(i) + sigma_down)

            output_up.cd()
            h_varied_nominal_up.Write()
            output_down.cd()
            h_varied_nominal_down.Write()

        output_up.Close()
        output_down.Close()

