#! /usr/bin/env python

import os
import subprocess
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--path", dest="root",
                  help="root directory where ROOT file are located. A directory structure with 'semimu', 'semie', 'data' is expected")

(options, args) = parser.parse_args()

if not options.root:
    parser.error("--path argument is required")

def formatFileName(token):
    return os.path.join(root, "MC_%s_theta_nominal.root" % token)

toMerge = [
        {
            "inputs": ["W1JetsToLNu", "W2JetsToLNu", "W3JetsToLNu", "W4JetsToLNu"],
            "output": "WJetsToLNu"
        },

        {
            "inputs": ["DY1JetsToLL_M-50", "DY2JetsToLL_M-50", "DY3JetsToLL_M-50", "DY4JetsToLL_M-50"],
            "output": "DYJetsToLL_M-50"
        },

        {
            "inputs": ["T_s-channel", "T_t-channel", "T_tW-channel"],
            "output": "T"
        },

        {
            "inputs": ["Tbar_s-channel", "Tbar_t-channel", "Tbar_tW-channel"],
            "output": "Tbar"
        },

        {
            "inputs": ["WW", "WZ", "ZZ"],
            "output": "dibosons"
        },

        {
            "inputs": ["QCD_pt15to30_bEnriched_MuEnrichedPt14", "QCD_pt30to50_bEnriched_MuEnrichedPt14", "QCD_pt50to150_bEnriched_MuEnrichedPt14", "QCD_pt150_bEnriched_MuEnrichedPt14"],
            "output": "QCD_MuEnriched"
        },

        {
            "inputs": ["QCD_Pt_20_30_EMEnriched", "QCD_Pt_30_80_EMEnriched", "QCD_Pt_80_170_EMEnriched", "QCD_Pt_170_250_EMEnriched", "QCD_Pt_250_350_EMEnriched", "QCD_Pt_350_EMEnriched"],
            "output": "QCD_EMEnriched"
        }
        ]

for type in ["semie", "semimu"]:
    root = os.path.join(options.root, type)

    for data in toMerge:
        filenames = map(formatFileName, data["inputs"])
        skip_job = False
        for file in filenames:
            if not os.path.exists(file):
                print("Warning input file '%s' does not exist. Skipping job." % file)
                skip_job = True
        if skip_job :
            continue
        args = ["hadd", "-f", formatFileName(data["output"])]
        args += filenames
        #print args
        subprocess.call(args)

systs = ["JECup", "JECdown", "JERup", "JERdown", "puUp", "puDown", "trigUp", "trigDown", "leptUp", "leptDown", "btagUp", "btagDown"]
for type in ["semie", "semimu"]:
    root = os.path.join(options.root, type, "Systematics")

    for syst in systs:
        for data in toMerge:
            filenames = map(formatFileName, data["inputs"])
            skip_job = False
            for file in filenames:
                if not os.path.exists(file.replace("nominal", syst)):
                    print("Warning input file '%s' does not exist. Skipping job." % file)
                    skip_job = True
            if skip_job :
                continue
            args = ["hadd", "-f", formatFileName(data["output"]).replace("nominal", syst)]
            args += [f.replace("nominal", syst) for f in filenames]

            subprocess.call(args)
            #print args
