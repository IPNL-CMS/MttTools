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
            "variation_up": "Systematics/MC_TT_powheg_theta_combinedAlphasPdfUp.root",
            "variation_down": "Systematics/MC_TT_powheg_theta_combinedAlphasPdfDown.root",
            "output_up": "Systematics/MC_TT_powheg_theta_alphasPdfUp.root",
            "output_down": "Systematics/MC_TT_powheg_theta_alphasPdfDown.root",
            "sets": ["CT10nlo", "MSTW2008nlo68cl", "MSTW2008nnlo68cl", "NNPDF23_nlo"]
        },
#        {
            #"nominal": "MC_TT_madgraph_theta_nominal.root",
            #"variation_alphas_up": "Systematics/MC_TT_madgraph_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/MC_TT_madgraph_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/MC_TT_madgraph_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/MC_TT_madgraph_theta_pdfDown.root",
            #"output_up": "Systematics/MC_TT_madgraph_theta_alphasPdfUp.root",
            #"output_down": "Systematics/MC_TT_madgraph_theta_alphasPdfDown.root"
        #},
#        {
            #"nominal": "Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_alphasPdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M400_cpl1_scalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M400_cpl1_scalar_theta_alphasPdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M500_cpl1_scalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M500_cpl1_scalar_theta_alphasPdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M600_cpl1_scalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M600_cpl1_scalar_theta_alphasPdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M700_cpl1_scalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M700_cpl1_scalar_theta_alphasPdfDown.root"
        #},
        #{
            #"nominal": "Signal_S0_S_i_M800_cpl1_scalar_theta_nominal.root",
            #"variation_alphas_up": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_alphasUp.root",
            #"variation_alphas_down": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_alphasDown.root",
            #"variation_pdf_up": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_pdfUp.root",
            #"variation_pdf_down": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_pdfDown.root",
            #"output_up": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_alphasPdfUp.root",
            #"output_down": "Systematics/Signal_S0_S_i_M800_cpl1_scalar_theta_alphasPdfDown.root"
        #},

]

def retrieve_var_histo(hist, set, variation, h_nominal):
    #h_var = 0
    for var in variation.GetKeyNames():
        if hist in var and set in var:
            h_var = variation.Get(var)
            assert(h_nominal.GetNbinsX() == h_var.GetNbinsX())
    return h_var

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

            list_h_variation_up = []
            list_h_variation_down = []
            h_varied_nominal_down.Reset()
            h_varied_nominal_up.Reset()

            for var in variation_up.GetKeyNames():
                if hist in var:
                    if not "nominal" in var:
                        h_variation = variation_up.Get(var)
                        assert(h_nominal.GetNbinsX() == h_variation.GetNbinsX())
                        list_h_variation_up.append(h_variation)
                    else:
                        h_nominal_NNPDF = variation_up.Get(var)
                        assert(h_nominal.GetNbinsX() == h_nominal_NNPDF.GetNbinsX())

            for var in variation_down.GetKeyNames():
                if hist in var:
                    if not "nominal" in var:
                        h_variation = variation_down.Get(var)
                        assert(h_nominal.GetNbinsX() == h_variation.GetNbinsX())
                        list_h_variation_down.append(h_variation)

            for i in range(1, h_nominal.GetNbinsX() + 1):
                nominal_value = h_nominal.GetBinContent(i)
                nominal_NNPDF_value = h_nominal_NNPDF.GetBinContent(i)
                assert(len(list_h_variation_up) == len(list_h_variation_down))
                list_up_value = []
                list_down_value = []
                sigma_env_MSTW_nlo_up = 0
                sigma_env_MSTW_nlo_down = 0
                sigma_MSTW_nnlo_up = 0
                sigma_MSTW_nnlo_down = 0
                if nominal_value:
                    for j,h_var in enumerate(list_h_variation_up):
                        variation_up_value = list_h_variation_up[j].GetBinContent(i);
                        variation_down_value = list_h_variation_down[j].GetBinContent(i);

                        if "MSTW2008nlo" in h_var.GetName():
                            sigma_env_MSTW_nlo_up = (variation_up_value-nominal_value)/nominal_value;
                            sigma_env_MSTW_nlo_down = (nominal_value - variation_down_value)/nominal_value;

                        if "MSTW2008nnlo" in h_var.GetName():
                            sigma_MSTW_nnlo_up = variation_up_value - nominal_value;
                            sigma_MSTW_nnlo_down = nominal_value - variation_down_value;

                        if not "nnlo" in h_var.GetName():
                            list_up_value.append(variation_up_value)
                            list_down_value.append(variation_down_value)
   
                    # first compute the NLO-QCD envelope using CT10, MSTW and NNPDF nlo
                    # U: upper edge of the envelope
                    # L: lower edge of the envelope
                    # M: mid-point of the envelope
                    U = max(list_up_value)
                    L = min(list_down_value)
                    M = (U+L)/2.

                    sigma_env = 0
                    if M:
                        sigma_env = (U-M)/M

                    # Compute the Rplus/minus factor
                    Rplus = sigma_env / sigma_env_MSTW_nlo_up
                    Rminus = sigma_env / sigma_env_MSTW_nlo_down

                    # Compute the total envelope at NNLO-QCD which is given by the MSTW-NNLO alphas+PDF band, 
                    # multiplied by the rescaling factor Rplus/minus
                    sigma_up = sigma_MSTW_nnlo_up * Rplus
                    sigma_down = sigma_MSTW_nnlo_down * Rminus

                    h_varied_nominal_up.SetBinContent(i, h_nominal.GetBinContent(i) + sigma_up)
                    h_varied_nominal_down.SetBinContent(i, h_nominal.GetBinContent(i) - sigma_down)
                else:
                    h_varied_nominal_up.SetBinContent(i, h_nominal.GetBinContent(i))
                    h_varied_nominal_down.SetBinContent(i, h_nominal.GetBinContent(i))

            output_up.cd()
            h_varied_nominal_up.Write()
            output_down.cd()
            h_varied_nominal_down.Write()
        output_up.Close()
        output_down.Close()
