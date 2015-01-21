#! /usr/bin/env python

from ROOT import TFile, gDirectory

def GetKeyNames( self, dir = "" ):
    self.cd(dir)
    return [key.GetName() for key in gDirectory.GetListOfKeys()]
TFile.GetKeyNames = GetKeyNames

import os, datetime
from optparse import OptionParser
from math import sqrt

d = datetime.datetime.now().strftime("%d%b%y")

files = [

        ["MC_TT_powheg_histos_nominal.root", "Systematics/MC_TT_powheg_histos_%s.root"],
        ["Signal_S0_S_i_M400_cpl1_pseudoscalar_histos_nominal.root", "Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_histos_%s.root"],
        ["Signal_S0_S_i_M500_cpl1_pseudoscalar_histos_nominal.root", "Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_histos_%s.root"],
        ["Signal_S0_S_i_M600_cpl1_pseudoscalar_histos_nominal.root", "Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_histos_%s.root"],
        ["Signal_S0_S_i_M700_cpl1_pseudoscalar_histos_nominal.root", "Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_histos_%s.root"],
        ["Signal_S0_S_i_M800_cpl1_pseudoscalar_histos_nominal.root", "Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_histos_%s.root"],

]

#btags = [-1, 0, 1, 2]
btags = [1, 2]

def getBFolderName(btag):
    b = "%d-btag" % btag
    if btag == -1:
        b = "all-btag"
    return b

type = ["semie", "semimu"]

for file in files:
    for btag in btags:
    #for btag in [2]:
        for type in ["semie", "semimu"]:
            path = "plots/%s/%s/%s" % (d, getBFolderName(btag), type)

            # Open nominal file, as well as up and down variations
            nominal = TFile.Open(os.path.join(path, file[0]))
            variation_alphas_up = TFile.Open(os.path.join(path, file[1] % "alphasUp"))
            variation_alphas_down = TFile.Open(os.path.join(path, file[1] % "alphasDown"))
            variation_pdf_up = TFile.Open(os.path.join(path, file[1] % "pdfUp"))
            variation_pdf_down = TFile.Open(os.path.join(path, file[1] % "pdfDown"))

            print("Creating %s" % os.path.join(path, file[1] % "alphasPdfUp"))
            output_up = TFile.Open(os.path.join(path, file[1] % "alphasPdfUp"), "recreate")

            print("Creating %s" % os.path.join(path, file[1] % "alphasPdfDown"))
            output_down = TFile.Open(os.path.join(path, file[1] % "alphasPdfDown"), "recreate")

            histoList = ["mttSelected_btag_sel_reco_fullsel_binning15GeV"]
            #for hist in nominal.GetKeyNames():
            for hist in histoList:
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


