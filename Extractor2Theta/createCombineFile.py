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

#MC = {"TT_powheg": "TT", "T_tW-channel": "T_tW", "T_s-channel": "T_s", "T_t-channel": "T_t", "Tbar_tW-channel": "Tbar_tW", "Tbar_s-channel": "Tbar_s",
      #"Tbar_t-channel": "Tbar_t", "DY1JetsToLL_M-50": "Z1Jets", "DY2JetsToLL_M-50": "Z2Jets", "DY3JetsToLL_M-50": "Z3Jets", "DY4JetsToLL_M-50": "Z4Jets",
      #"W1JetsToLNu": "W1Jets", "W2JetsToLNu": "W2Jets", "W3JetsToLNu": "W3Jets", "W4JetsToLNu": "W4Jets",
      #"ZZ": "ZZ", "WW": "WW", "WZ": "WZ"
     #}

MC = {
        "TT_powheg": "TT",
        "T": "single_top",
        "Tbar": "single_antitop",
        "DYJetsToLL_M-50": "zjets",
        "WJetsToLNu": "wjets",
        "dibosons": "dibosons"
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

            "ZPrimeToTTJets_M500GeV_W5GeV": "zp500_narrow",
            "ZPrimeToTTJets_M750GeV_W7p5GeV": "zp750_narrow",
            "ZPrimeToTTJets_M1000GeV_W10GeV": "zp1000_narrow",
            "ZPrimeToTTJets_M1250GeV_W12p5GeV": "zp1250_narrow",
            "ZPrimeToTTJets_M1500GeV_W15GeV": "zp1500_narrow",
            "ZPrimeToTTJets_M2000GeV_W20GeV": "zp2000_narrow",

            "ZPrimeToTTJets_M500GeV_W50GeV": "zp500_large",
            "ZPrimeToTTJets_M750GeV_W75GeV": "zp750_large",
            "ZPrimeToTTJets_M1000GeV_W100GeV": "zp1000_large",
            "ZPrimeToTTJets_M1250GeV_W125GeV": "zp1250_large",
            "ZPrimeToTTJets_M1500GeV_W150GeV": "zp1500_large",
            "ZPrimeToTTJets_M2000GeV_W200GeV": "zp2000_large",
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

if True:
    systs += ["matchingup", "matchingdown", "scaleup", "scaledown"]
    

if False:
    systs += ["pdfUp", "pdfDown"]

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
            
            hist.SetName("mtt_%s_%dbtag__%s" % (tag, btag, title))
            histos.append(hist)
        f_.Close()

        for syst in systs:
            sign = "Up" if "up" in syst.lower() else "Down"
            theta_syst = syst.lower().replace("up", "").replace("down", "")

            f = os.path.join(root, type, "Systematics", "MC_%s_theta_%s.root" % (filename, syst))
            if not os.path.exists(f):
                continue
            f_ = TFile.Open(f)
            for btag in [1, 2]:
                hist_name = "mtt_%dbtag" % btag
                hist = f_.Get(hist_name).Clone()
                hist.SetDirectory(0)

                hist.SetName("mtt_%s_%dbtag__%s__%s%s" % (tag, btag, title, theta_syst, sign))
                histos.append(hist)
            f_.Close()

# Process data
for type in ["semie", "semimu"]:
    tag = "e" if "semie" in type else "mu"
    f = os.path.join(root, "data", "MTT_Data_%s.root" % type)
    f_ = TFile.Open(f)
    for btag in [1, 2]:
        hist_name = "mtt_%dbtag" % btag
        hist = f_.Get(hist_name).Clone()
        hist.SetDirectory(0)
        hist.SetName("mtt_%s_%dbtag__data_obs" % (tag, btag))
        histos.append(hist)
    f_.Close()

# Process signal


for filename, title in Signal.items():
    signal_histos = []
    for type in ["semie", "semimu"]:
        tag = "e" if "semie" in type else "mu"
        f = os.path.join(root, type, "Signal_%s_theta_nominal.root" % filename)
        f_ = TFile.Open(f)
        for btag in [1, 2]:
            hist_name = "mtt_%dbtag" % btag
            hist = f_.Get(hist_name).Clone()
            hist.SetDirectory(0)
            hist.SetName("mtt_%s_%dbtag__%s" % (tag, btag, "signal"))
            signal_histos.append(hist)
        f_.Close()

        for syst in systs:
            sign = "Up" if "up" in syst.lower() else "Down"
            theta_syst = syst.lower().replace("up", "").replace("down", "")

            f = os.path.join(root, type, "Systematics", "Signal_%s_theta_%s.root" % (filename, syst))
            if not os.path.exists(f):
                continue
            f_ = TFile.Open(f)
            for btag in [1, 2]:
                hist_name = "mtt_%dbtag" % btag
                hist = f_.Get(hist_name).Clone()
                hist.SetDirectory(0)
                hist.SetName("mtt_%s_%dbtag__%s__%s%s" % (tag, btag, "signal", theta_syst, sign))
                signal_histos.append(hist)
            f_.Close()


    output_filename = options.output % title
    output = TFile.Open(output_filename, "recreate")
    for hist in histos:
        hist.Write()
    for sig in signal_histos:
        sig.Write()

    output.Close()
