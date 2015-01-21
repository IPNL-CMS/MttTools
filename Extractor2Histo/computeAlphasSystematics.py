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
            variation_up = TFile.Open(os.path.join(path, file[1] % "alphasOnlyUp"))
            variation_down = TFile.Open(os.path.join(path, file[1] % "alphasOnlyDown"))

            print("Creating %s" % os.path.join(path, file[1] % "alphasUp"))
            output_up = TFile.Open(os.path.join(path, file[1] % "alphasUp"), "recreate")

            print("Creating %s" % os.path.join(path, file[1] % "alphasDown"))
            output_down = TFile.Open(os.path.join(path, file[1] % "alphasDown"), "recreate")


            histoList = ["mttSelected_btag_sel_reco_fullsel_binning15GeV"]
            #for hist in nominal.GetKeyNames():
            for hist in histoList:
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

