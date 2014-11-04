#! /bin/env python

import os, subprocess, datetime, tempfile, sys, math
from optparse import OptionParser

from ROOT import TH1, TFile

d = datetime.datetime.now().strftime("%d%b%y")

parser = OptionParser()
parser.add_option("-i", "--input", action="store", dest="input", help="Input folder where root files are located")
(option, args) = parser.parse_args()

if not option.input:
    parser.error("--input argument is required")

filename = "MC_TT_powheg_histos_nominal.root"
filename_syst = "MC_TT_powheg_histos_%s.root"

categories = ["semimu", "semie"]
btags = [1, 2]
def getBFolderName(btag):
    b = "%d-btag" % btag
    if btag == -1:
        b = "all-btag"

    return b

scale_factors = {
        "lepton": ["lepton_weight_fullsel", "lept"],
        "trigger": ["trigger_weight_fullsel", "trig"],
        "btag": ["btag_weight_fullsel", "btag"],
        }

for category in categories:
    for btag in btags:
        nominal_file = os.path.join(option.input, getBFolderName(btag), category, filename)
        nominal_file_ = TFile.Open(nominal_file)

        for scale_factor_name, scale_factor_data in scale_factors.items():

            hist = nominal_file_.Get(scale_factor_data[0])
            scale_factor_nominal = hist.GetMean()

            # Open scale factor up and down variation
            up_file_ = TFile.Open(os.path.join(option.input, getBFolderName(btag), category, "Systematics", filename_syst % (scale_factor_data[1] + "Up")))
            down_file_ = TFile.Open(os.path.join(option.input, getBFolderName(btag), category, "Systematics", filename_syst % (scale_factor_data[1] + "Down")))

            up_hist = up_file_.Get(scale_factor_data[0])
            down_hist = down_file_.Get(scale_factor_data[0])

            error = (math.fabs(up_hist.GetMean() - scale_factor_nominal) + math.fabs(scale_factor_nominal - down_hist.GetMean())) / 2

            #print("Nominal: %.3f, Up: %.3f, Down: %.3f, Error: %.3f" % (scale_factor_nominal, up_hist.GetMean(), down_hist.GetMean(), error))

            print("%s - %s - %s: %.3f \pm %.3f" % (category, getBFolderName(btag), scale_factor_name, scale_factor_nominal, error))


        nominal_file_.Close()
