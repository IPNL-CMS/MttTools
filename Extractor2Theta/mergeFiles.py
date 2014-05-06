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
        }
        ]

for type in ["semie", "semimu"]:
    root = os.path.join(options.root, type)

    for data in toMerge:
        filenames = map(formatFileName, data["inputs"])
        args = ["hadd", "-f", formatFileName(data["output"])]
        args += filenames

        subprocess.call(args)

systs = ["JECup", "JECdown", "JERup", "JERdown", "puUp", "puDown", "trigUp", "trigDown", "leptUp", "leptDown", "btagUp", "btagDown"]
for type in ["semie", "semimu"]:
    root = os.path.join(options.root, type, "Systematics")

    for syst in systs:
        for data in toMerge:
            filenames = map(formatFileName, data["inputs"])
            args = ["hadd", "-f", formatFileName(data["output"]).replace("nominal", syst)]
            args += [f.replace("nominal", syst) for f in filenames]

            subprocess.call(args)
            #print args
