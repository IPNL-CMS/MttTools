#! /bin/env python

import os, datetime
import json, re
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

        # Signal
        #"S0_S_i_M400_cpl1_scalar": "H400_scalar",
        #"S0_S_i_M500_cpl1_scalar": "H500_scalar",
        #"S0_S_i_M600_cpl1_scalar": "H600_scalar",
        #"S0_S_i_M700_cpl1_scalar": "H700_scalar",
        #"S0_S_i_M800_cpl1_scalar": "H800_scalar",
        "S0_S_i_M400_cpl1_pseudoscalar": "H400_pseudoscalar",
        "S0_S_i_M500_cpl1_pseudoscalar": "H500_pseudoscalar",
        "S0_S_i_M600_cpl1_pseudoscalar": "H600_pseudoscalar",
        "S0_S_i_M700_cpl1_pseudoscalar": "H700_pseudoscalar",
        "S0_S_i_M800_cpl1_pseudoscalar": "H800_pseudoscalar",
     }

root = options.root

input_cov = os.path.join(root, "cov_matrix_allParameters_withSignal_%s.json" % options.signal)
with open(input_cov) as json_data:
    cov = json.load(json_data)

corr = cov["cov"]
covv = cov["cov"]
#print corr
#print covv
for i in range(len(corr)):
    for j in range(len(corr[i])):
        if not (sqrt(covv[i][i] * covv[j][j]) == 0):
            #corr[i][j] = covv[i][j] / (sqrt(covv[i][i]) * sqrt(covv[j][j]))
            corr[i][j] = covv[i][j] / (sqrt(covv[i][i] * covv[j][j]))
        else:
            print corr[i][j]

#print corr

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

#theta_syst = ["alphaspdf", "btag", "diboson_rate", "jec", "jer", "lept_e", "lept_mu", "lumi", "pu", "sat_rate", "scale", "st_rate", "trigger", "ttbar_rate", "w_rate", "z_rate"]
theta_syst = cov["parameters"]
print("Theta systematics: ")
print theta_syst

# nominal histograms post fit
f = os.path.join(root, "%s" % (options.signal), "histos-mle-%s-fitMethod1.root" % options.signal)
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
            if "scalar" in title:
                out_name = os.path.join(path, "Signal_%s_histos_nominal_postfit.root" % filename)
            else:
                out_name = os.path.join(path, "MC_%s_histos_nominal_postfit.root" % filename)

            print("Creating %s" % out_name)
            output = TFile.Open(out_name, "recreate")
            output.cd()
            hist.Write()
            output.Close()


# systematics post fit
for type in ["semie", "semimu"]:
    tag = "e" if "semie" in type else "mu"
    for btag in btags:
        for filename, title in MC.items():
            hist_name = "mtt_%s_%dbtag__%s" % (tag, btag, title)
            sigma_x = {} # store sigma per mtt bin per syst 
            for syst in theta_syst:
                if syst == "beta_signal":
                    continue
                # Special treatment for lept syst 
                if syst == "lept_mu" and type == "semie":
                    continue
                if syst == "lept_e" and type == "semimu":
                    continue

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
                
                hist_nominal = f_nominal_.Get(hist_name).Clone()
                hist_nominal.SetDirectory(0)
                hist_up = f_up_.Get(hist_name).Clone()
                hist_up.SetDirectory(0)
                hist_down = f_down_.Get(hist_name).Clone()
                hist_down.SetDirectory(0)

                path = "plots/%s/%s/%s/Systematics" % (d, getBFolderName(btag), type)
                if "scalar" in title:
                    out_name_per_syst = os.path.join(path, "Signal_%s_histos_%s_syst_errors_postfit.root" % (filename, syst))
                else:
                    out_name_per_syst = os.path.join(path, "MC_%s_histos_%s_syst_errors_postfit.root" % (filename, syst))

                print("Creating %s" % out_name_per_syst)
                output_per_syst = TFile.Open(out_name_per_syst, "recreate")
                uncertainties_percent = hist_nominal.Clone()
                uncertainties_percent.SetDirectory(None)
                uncertainties_percent.Reset() # Keep binning but remove all events
                uncertainties_percent.SetName("mttSelected_btag_sel_reco_fullsel_binning15GeV")

                uncertainties = hist_nominal.Clone()
                uncertainties.SetDirectory(None)
                uncertainties.Reset() # Keep binning but remove all events
                uncertainties.SetName("mttSelected_btag_sel_reco_fullsel_binning15GeV_absolute_sigma")


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

                    uncertainties_percent.SetBinContent(i, 0.)
                    uncertainties_percent.SetBinError(i, delta_percent)

                    uncertainties.SetBinContent(i, 0.)
                    uncertainties.SetBinError(i, delta)
                output_per_syst.cd()
                uncertainties_percent.Write()
                uncertainties.Write()
                sigma_x[syst] = uncertainties
                output_per_syst.Close()
                f_nominal_.Close()
                f_up_.Close()
                f_down_.Close()


            path = "plots/%s/%s/%s/Systematics" % (d, getBFolderName(btag), type)
            if "scalar" in title:
                out_name = os.path.join(path, "Signal_%s_histos_total_syst_errors_postfit.root" % filename)
            else:
                out_name = os.path.join(path, "MC_%s_histos_total_syst_errors_postfit.root" % filename)
            print("Creating %s" % out_name)
            output = TFile.Open(out_name, "recreate")
            total_uncertainties = f_.Get(hist_name).Clone()
            total_uncertainties.SetDirectory(0)
            #total_uncertainties.Reset() # Keep binning but remove all events
            total_uncertainties.SetName("mttSelected_btag_sel_reco_fullsel_binning15GeV")
            for k in range(1, hist_nominal.GetNbinsX()+1):
                #total_uncertainties.SetBinContent(k, 0)
                delta_2 = 0.
                for i, syst_i in enumerate(theta_syst):
                    if syst_i == "beta_signal":
                        continue
                    if syst_i == "lept_mu" and type == "semie":
                        continue
                    if syst_i == "lept_e" and type == "semimu":
                        continue
                    for j, syst_j in enumerate(theta_syst):
                        if syst_j == "beta_signal":
                            continue
                        if syst_j == "lept_mu" and type == "semie":
                            continue
                        if syst_j == "lept_e" and type == "semimu":
                           continue
                        sigma_i = sigma_x[syst_i].GetBinError(k)
                        sigma_j = sigma_x[syst_j].GetBinError(k)
                        delta_2 = delta_2 + corr[i][j] * sigma_i * sigma_j # if i = j, corr[i][j] = 1
                total_uncertainties.SetBinError(k, sqrt(delta_2))
            output.cd()
            total_uncertainties.Write()
            output.Close()



