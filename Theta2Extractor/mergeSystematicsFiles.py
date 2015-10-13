#! /bin/env python

import os, datetime, subprocess
from ROOT import TH1, TFile, TObject
from math import fabs, sqrt
from optparse import OptionParser

d = datetime.datetime.now().strftime("%d%b%y")

parser = OptionParser()
parser.add_option("-p", "--path", dest="root",
                  help="root directory where ROOT file are located.")
parser.add_option("-s", "--signal",
                  dest="signal", help="MC histograms after postfit with which signal? Precise!")

(options, args) = parser.parse_args()

if not options.root:
    parser.error("--path argument is required")

if not options.signal:
    parser.error("--signal argument is required")

def getBFolderName(btag):
    b = "%d-btag" % btag
    if btag == -1:
        b = "all-btag"

    return b

#btags = [0, 1, 2]
btags = [1, 2]

MC = {
        # Background
        "TT_powheg": "TT",
        #"TT_madgraph": "TT",
        "T": "single_top",
        "Tbar": "single_antitop",
        "DYJetsToLL_M-50": "zjets",
        "WJetsToLNu": "wjets",
        "dibosons": "dibosons",
        #"QCD_EMEnriched": "QCD",
        #"QCD_MuEnriched": "QCD",
     }

theta_syst = ["alphaspdf", "btag", "diboson_rate", "jec", "jer", "lept_e", "lept_mu", "lumi", "pu", "sat_rate", "scale", "st_rate", "trigger", "ttbar_rate", "w_rate", "z_rate"]

root = options.root

# Build output tree structure
for btag in btags:
    for type in ["semie", "semimu"]:
        path = "plots/%s/%s/%s" % (d, getBFolderName(btag), type)
        try:
            os.makedirs(path)
        except:
            pass

for btag in btags:
    for type in ["semie", "semimu"]:
        path = "plots/%s/%s/%s/Systematics" % (d, getBFolderName(btag), type)
        try:
            os.makedirs(path)
        except:
            pass


# systematics post fit
for syst in theta_syst:
    f_nominal = os.path.join(root, "%s" % (options.signal), "histos-mle-%s-fitMethod1-%s-nominal.root" % (options.signal, syst))
    if not os.path.exists(f_nominal):
        print("Warning input file '%s' does not exist. Skipping job." % f_nominal)
        continue
    f_nominal_ = TFile.Open(f_nominal)

    f_up = os.path.join(root, "%s" % (options.signal), "histos-mle-%s-fitMethod1-%s-up.root" % (options.signal, syst))
    if not os.path.exists(f_up):
        print("Warning input file '%s' does not exist. Skipping job." % f_up)
        continue
    f_up_ = TFile.Open(f_up)

    f_down = os.path.join(root, "%s" % (options.signal), "histos-mle-%s-fitMethod1-%s-down.root" % (options.signal, syst))
    if not os.path.exists(f_down):
        print("Warning input file '%s' does not exist. Skipping job." % f_down)
        continue
    f_down_ = TFile.Open(f_down)

    for type in ["semie", "semimu"]:
        tag = "e" if "semie" in type else "mu"
        # Special treatment for lept syst 
        if syst == "lept_mu" and type == "semie":
            continue
        if syst == "lept_e" and type == "semimu":
            continue

        for btag in btags:
            for filename, title in MC.items():
                hist_name = "mtt_%s_%dbtag__%s" % (tag, btag, title)
                hist_nominal = f_nominal_.Get(hist_name).Clone()
                hist_nominal.SetDirectory(0)
                hist_nominal.SetName("mttSelected_btag_sel_reco_fullsel_binning15GeV")
                hist_up = f_up_.Get(hist_name).Clone()
                hist_up.SetDirectory(0)
                hist_up.SetName("mttSelected_btag_sel_reco_fullsel_binning15GeV")
                hist_down = f_down_.Get(hist_name).Clone()
                hist_down.SetDirectory(0)
                hist_down.SetName("mttSelected_btag_sel_reco_fullsel_binning15GeV")

                path = "plots/%s/%s/%s/Systematics" % (d, getBFolderName(btag), type)
                out_name_nominal = os.path.join(path, "MC_%s_histos_%s_nominal_postfit.root" % (filename, syst))
                out_name_up = os.path.join(path, "MC_%s_histos_%s_up_postfit.root" % (filename, syst))
                out_name_down = os.path.join(path, "MC_%s_histos_%s_down_postfit.root" % (filename, syst))

                print("Creating %s" % out_name_nominal)
                output_nominal = TFile.Open(out_name_nominal, "recreate")
                output_nominal.cd()
                hist_nominal.Write()
                output_nominal.Close()

                print("Creating %s" % out_name_up)
                output_up = TFile.Open(out_name_up, "recreate")
                output_up.cd()
                hist_up.Write()
                output_up.Close()

                print("Creating %s" % out_name_down)
                output_down = TFile.Open(out_name_down, "recreate")
                output_down.cd()
                hist_down.Write()
                output_down.Close()

print("Merging MC...")
for type in ["semie", "semimu"]:
    tag = "e" if "semie" in type else "mu"
    for btag in btags:
        print("Category %s %d-btag" % (type, btag))
        path = "plots/%s/%s/%s/Systematics" % (d, getBFolderName(btag), type)
        for syst in theta_syst:
            if syst == "lept_mu" and type == "semie":
                continue
            if syst == "lept_e" and type == "semimu":
                continue
            args_nominal = ["hadd","-f", os.path.join(path, "MC_total_histos_%s_nominal_postfit.root" % (syst))]
            for filename, title in MC.items():
                args_nominal.append(os.path.join(path, "MC_%s_histos_%s_nominal_postfit.root" % (filename, syst)))
            subprocess.call(args_nominal)
            #print args_nominal
            args_up = ["hadd","-f", os.path.join(path, "MC_total_histos_%s_up_postfit.root" % (syst))]
            for filename, title in MC.items():
                args_up.append(os.path.join(path, "MC_%s_histos_%s_up_postfit.root" % (filename, syst)))
            subprocess.call(args_up)
            #print args_up
            args_down = ["hadd","-f", os.path.join(path, "MC_total_histos_%s_down_postfit.root" % (syst))]
            for filename, title in MC.items():
                args_down.append(os.path.join(path, "MC_%s_histos_%s_down_postfit.root" % (filename, syst)))
            subprocess.call(args_down)
            #print args_down


