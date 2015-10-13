#! /bin/env python

import os, datetime
from ROOT import TH1, TFile, TObject
from math import fabs, sqrt
from optparse import OptionParser

d = datetime.datetime.now().strftime("%d%b%y")

parser = OptionParser()
parser.add_option("-p", "--path", dest="root",
                  help="root directory where ROOT file are located.")
#parser.add_option("-o", "--output",
                  #dest="output", help="Name of output file")

(options, args) = parser.parse_args()

if not options.root:
    parser.error("--path argument is required")

#if not options.output:
    #parser.error("--output argument is required")

def getBFolderName(btag):
    b = "%d-btag" % btag
    if btag == -1:
        b = "all-btag"

    return b

#btags = [0, 1, 2]
btags = [1, 2]

MC = {
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

print("Processing...")

theta_syst = ["alphaspdf", "btag", "diboson_rate", "jec", "jer", "lept_e", "lept_mu", "lumi", "pu", "sat_rate", "scale", "st_rate", "trigger", "ttbar_rate", "w_rate", "z_rate"]

# nominal histograms post fit
f = os.path.join(root, "bkgOnly", "histos-mle-bkgOnly-fitMethod1.root")
if not os.path.exists(f):
    print("Warning input file '%s' does not exist. Skipping job." % f)
f_ = TFile.Open(f)

for type in ["semie", "semimu"]:
    tag = "e" if "semie" in type else "mu"
    for btag in btags:
        for filename, title in MC.items():
            hist_name = "mtt_%s_%dbtag__%s" % (tag, btag, title)
            hist = f_.Get(hist_name).Clone()
            hist.SetDirectory(0)
            hist.SetName("mttSelected_btag_sel_reco_fullsel_binning15GeV")

            path = "plots/%s/%s/%s" % (d, getBFolderName(btag), type)
            out_name = os.path.join(path, "MC_%s_histos_nominal_postfit.root" % filename)

            print("Creating %s" % out_name)
            output = TFile.Open(out_name, "recreate")
            output.cd()
            hist.Write()
            output.Close()

# systematics post fit
for syst in theta_syst:
    f_nominal = os.path.join(root, "bkgOnly", "histos-mle-bkgOnly-fitMethod1-%s-nominal.root" % (syst))
    if not os.path.exists(f_nominal):
        print("Warning input file '%s' does not exist. Skipping job." % f_nominal)
        continue
    f_nominal_ = TFile.Open(f_nominal)

    f_up = os.path.join(root, "bkgOnly", "histos-mle-bkgOnly-fitMethod1-%s-up.root" % (syst))
    if not os.path.exists(f_up):
        print("Warning input file '%s' does not exist. Skipping job." % f_up)
        continue
    f_up_ = TFile.Open(f_up)

    f_down = os.path.join(root, "bkgOnly", "histos-mle-bkgOnly-fitMethod1-%s-down.root" % (syst))
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
                hist_up = f_up_.Get(hist_name).Clone()
                hist_up.SetDirectory(0)
                hist_down = f_down_.Get(hist_name).Clone()
                hist_down.SetDirectory(0)

                path = "plots/%s/%s/%s/Systematics" % (d, getBFolderName(btag), type)
                out_name = os.path.join(path, "MC_%s_histos_total_syst_errors_postfit.root" % filename)
                out_name_per_syst = os.path.join(path, "MC_%s_histos_%s_syst_errors_postfit.root" % (filename, syst))

                print("Creating %s" % out_name_per_syst)
                output_per_syst = TFile.Open(out_name_per_syst, "recreate")
                uncertainties = hist_nominal.Clone()
                uncertainties.SetDirectory(None)
                uncertainties.Reset() # Keep binning but remove all events
                uncertainties.SetName("mttSelected_btag_sel_reco_fullsel_binning15GeV")
                if not os.path.exists(out_name):
                    print("Creating %s" % out_name)
                    output = TFile.Open(out_name, "recreate")
                    total_uncertainties = hist_nominal.Clone()
                    total_uncertainties.SetDirectory(None)
                    total_uncertainties.Reset() # Keep binning but remove all events
                    total_uncertainties.SetName("mttSelected_btag_sel_reco_fullsel_binning15GeV")
                    total_uncertainties.Write()
                else:
                    output = TFile.Open(out_name, "UPDATE")

                total_uncertainties = output.Get("mttSelected_btag_sel_reco_fullsel_binning15GeV")

                for i in range(1, hist_nominal.GetNbinsX()+1):

                    nominal_x = hist_nominal.GetBinCenter(i);
                    nominal_value = hist_nominal.GetBinContent(i);

                    norm_plus = 0;
                    norm_minus = 0;
                    if (hist_up.GetBinError(i) / hist_up.GetBinContent(i) < 0.4):
                        norm_plus = hist_up.GetBinContent(i) - nominal_value
                    if (hist_down.GetBinError(i) / hist_down.GetBinContent(i) < 0.4):
                        norm_minus = nominal_value - hist_down.GetBinContent(i)

                    delta = (fabs(norm_plus) + fabs(norm_minus)) / 2.
                    delta_percent = 0 if (delta == 0 or nominal_value == 0) else delta / nominal_value

                    uncertainties.SetBinContent(i, 0.)
                    uncertainties.SetBinError(i, delta_percent)

                    total_uncertainties.SetBinContent(i, 0)
                    total_error = total_uncertainties.GetBinError(i)
                    total_uncertainties.SetBinError(i, sqrt(delta_percent * delta_percent + total_error * total_error))

                output.cd()
                total_uncertainties.Write("",TObject.kOverwrite)
                output.Close()

                output_per_syst.cd()
                uncertainties.Write()
                output_per_syst.Close()



