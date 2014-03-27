#! /bin/env python

import os
from ROOT import TH1, TFile
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--path", dest="root",
                  help="root directory where ROOT file are located. A directory structure with 'semimu', 'semie', 'data' is expected")
parser.add_option("-o", "--output",
                  dest="output", help="Name of output file")

(options, args) = parser.parse_args()

if not options.root:
    parser.error("--path argument is required")

if not options.output:
    parser.error("--output argument is required")

MC = {"TT_powheg": "TT", "T_tW-channel": "T_tW", "T_s-channel": "T_s", "T_t-channel": "T_t", "Tbar_tW-channel": "Tbar_tW", "Tbar_s-channel": "Tbar_s",
      "Tbar_t-channel": "Tbar_t", "DY1JetsToLL_M-50": "Z1Jets", "DY2JetsToLL_M-50": "Z2Jets", "DY3JetsToLL_M-50": "Z3Jets", "DY4JetsToLL_M-50": "Z4Jets",
      "W1JetsToLNu": "W1Jets", "W2JetsToLNu": "W2Jets", "W3JetsToLNu": "W3Jets", "W4JetsToLNu": "W4Jets"  
     }

Signal = {
            "S0_S_i_M400_cpl1_scalar": "H400_scalar",
            "S0_S_i_M500_cpl1_scalar": "H500_scalar",
            "S0_S_i_M600_cpl1_scalar": "H600_scalar",
            "S0_S_i_M700_cpl1_scalar": "H700_scalar",
            "S0_S_i_M800_cpl1_scalar": "H800_scalar",
            "S0_S_i_M400_cpl1_pseudoscalar": "H400_pseudoscalar",
            "S0_S_i_M500_cpl1_pseudoscalar": "H500_pseudoscalar",
            "S0_S_i_M600_cpl1_pseudoscalar": "H600_pseudoscalar",
            "S0_S_i_M700_cpl1_pseudoscalar": "H700_pseudoscalar",
            "S0_S_i_M800_cpl1_pseudoscalar": "H800_pseudoscalar",
        }

systs = []

if True:
    systs = ["JECup", "JECdown", "JERup", "JERdown", "puUp", "puDown"]

if True:
    systs += ["trigUp", "trigDown"]

if True:
    systs += ["leptUp", "leptDown"]

if True:
    systs += ["btagUp", "btagDown"]

if False:
    systs += ["pdfUp", "pdfDown"]

mergeHistos = False

toMerge = {
        "singletop": ["T_t", "T_tW", "T_s", "Tbar_t", "Tbar_tW", "Tbar_s"],
        "wjets": ["W1Jets", "W2Jets", "W3Jets", "W4Jets"],
        "zjets": ["Z1Jets", "Z2Jets", "Z3Jets", "Z4Jets"]
        }

mergedHistos = {}

def shouldBeMerged(name):
    if not mergeHistos:
        return False, None

    for title, files in toMerge.items():
        if name in files:
            return True, title

    return False, None

# Load MC files

histos = []

root = options.root

print("Processing...")

for type in ["semie", "semimu"]:
    tag = "e" if "semie" in type else "mu"
    for filename, title in MC.items():
        f = os.path.join(root, type, "MC_%s_theta_nominal.root" % filename)
        f_ = TFile.Open(f)
        for btag in [1, 2]:
            hist_name = "mtt_%dbtag" % btag
            hist = f_.Get(hist_name).Clone()
            hist.SetDirectory(0)
            
            ret, mergedTitle = shouldBeMerged(title)
            if ret:
                if mergedTitle + str(btag) in mergedHistos:
                    mergedHistos[mergedTitle + str(btag)].Add(hist)
                else:
                    hist.SetName("mtt_%s_%dbtag__%s" % (tag, btag, mergedTitle))
                    mergedHistos[mergedTitle + str(btag)] = hist.Clone()
                    mergedHistos[mergedTitle + str(btag)].SetDirectory(0)
            else:
                hist.SetName("mtt_%s_%dbtag__%s" % (tag, btag, title))
                histos.append(hist)
        f_.Close()

        for syst in systs:
            sign = "up" if "up" in syst.lower() else "down"
            theta_syst = syst.lower().replace("up", "").replace("down", "")

            f = os.path.join(root, type, "Systematics", "MC_%s_theta_%s.root" % (filename, syst))
            f_ = TFile.Open(f)
            for btag in [1, 2]:
                hist_name = "mtt_%dbtag" % btag
                hist = f_.Get(hist_name).Clone()
                hist.SetDirectory(0)

                ret, mergedTitle = shouldBeMerged(title)
                if ret:
                    hist.SetName("mtt_%s_%dbtag__%s__%s__%s" % (tag, btag, mergedTitle, theta_syst, sign))
                    mergedTitle = "%s__%s__%s" % (mergedTitle, theta_syst, sign)
                    if mergedTitle + str(btag) in mergedHistos:
                        mergedHistos[mergedTitle + str(btag)].Add(hist)
                    else:
                        hist.SetName("mtt_%s_%dbtag__%s" % (tag, btag, mergedTitle))
                        mergedHistos[mergedTitle + str(btag)] = hist.Clone()
                        mergedHistos[mergedTitle + str(btag)].SetDirectory(0)
                else:
                    hist.SetName("mtt_%s_%dbtag__%s__%s__%s" % (tag, btag, title, theta_syst, sign))
                    histos.append(hist)
            f_.Close()

    for title, histo in mergedHistos.items():
        histos.append(histo)

    mergedHistos.clear()

# Process data
for type in ["semie", "semimu"]:
    tag = "e" if "semie" in type else "mu"
    f = os.path.join(root, "data", "MTT_Data_%s.root" % type)
    f_ = TFile.Open(f)
    for btag in [1, 2]:
        hist_name = "mtt_%dbtag" % btag
        hist = f_.Get(hist_name).Clone()
        hist.SetDirectory(0)
        hist.SetName("mtt_%s_%dbtag__DATA" % (tag, btag))
        histos.append(hist)
    f_.Close()

# Process signal

for type in ["semie", "semimu"]:
    tag = "e" if "semie" in type else "mu"
    for filename, title in Signal.items():
        f = os.path.join(root, type, "Signal_%s_theta_nominal.root" % filename)
        f_ = TFile.Open(f)
        for btag in [1, 2]:
            hist_name = "mtt_%dbtag" % btag
            hist = f_.Get(hist_name).Clone()
            hist.SetDirectory(0)
            hist.SetName("mtt_%s_%dbtag__%s" % (tag, btag, title))
            histos.append(hist)
        f_.Close()

        for syst in systs:
            sign = "up" if "up" in syst.lower() else "down"
            theta_syst = syst.lower().replace("up", "").replace("down", "")

            f = os.path.join(root, type, "Systematics", "Signal_%s_theta_%s.root" % (filename, syst))
            f_ = TFile.Open(f)
            for btag in [1, 2]:
                hist_name = "mtt_%dbtag" % btag
                hist = f_.Get(hist_name).Clone()
                hist.SetDirectory(0)
                hist.SetName("mtt_%s_%dbtag__%s__%s__%s" % (tag, btag, title, theta_syst, sign))
                histos.append(hist)
            f_.Close()

output = TFile.Open(options.output, "recreate")
for hist in histos:
    hist.Write()

output.Close()
