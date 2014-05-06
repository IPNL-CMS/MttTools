#! /usr/bin/env python

from ROOT import TFile, gDirectory

def GetKeyNames( self, dir = "" ):
    self.cd(dir)
    return [key.GetName() for key in gDirectory.GetListOfKeys()]
TFile.GetKeyNames = GetKeyNames

import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--path", dest="root",
                  help="root directory where ROOT file are located. A directory structure with 'semimu', 'semie', 'data' is expected")

(options, args) = parser.parse_args()

if not options.root:
    parser.error("--path argument is required")

systs = [
        
        {
            "nominal": "MC_TT_madgraph_theta_nominal.root",
            "variation": "Systematics/MC_TT_madgraph_theta_matchingup.root",
            "apply-to": "MC_TT_powheg_theta_nominal.root",
            "output": "Systematics/MC_TT_powheg_theta_matchingup.root"
        },

        {
            "nominal": "MC_TT_madgraph_theta_nominal.root",
            "variation": "Systematics/MC_TT_madgraph_theta_matchingdown.root",
            "apply-to": "MC_TT_powheg_theta_nominal.root",
            "output": "Systematics/MC_TT_powheg_theta_matchingdown.root"
        },
        
        {
            "nominal": "MC_TT_madgraph_theta_nominal.root",
            "variation": "Systematics/MC_TT_madgraph_theta_scaleup.root",
            "apply-to": "MC_TT_powheg_theta_nominal.root",
            "output": "Systematics/MC_TT_powheg_theta_scaleup.root"
        },

        {
            "nominal": "MC_TT_madgraph_theta_nominal.root",
            "variation": "Systematics/MC_TT_madgraph_theta_scaledown.root",
            "apply-to": "MC_TT_powheg_theta_nominal.root",
            "output": "Systematics/MC_TT_powheg_theta_scaledown.root"
        },
]

type = ["semie", "semimu"]

for type in ["semie", "semimu"]:
    for syst in systs:

        # Open nominal file, as well as up and down variations
        nominal = TFile.Open(os.path.join(options.root, type, syst["nominal"]))
        variation = TFile.Open(os.path.join(options.root, type, syst["variation"]))

        reference = TFile.Open(os.path.join(options.root, type, syst["apply-to"]))

        print("Creating %s" % os.path.join(options.root, type, syst["output"]))
        output = TFile.Open(os.path.join(options.root, type, syst["output"]), "recreate")

        for hist in nominal.GetKeyNames():
            h_nominal = nominal.Get(hist)
            h_variation = variation.Get(hist)

            h_reference = reference.Get(hist)
            h_varied_reference = h_reference.Clone()
            h_varied_reference.SetDirectory(0)

            assert(h_nominal.GetNbinsX() == h_variation.GetNbinsX() == h_reference.GetNbinsX())

            for i in range(1, h_nominal.GetNbinsX() + 1):
                nominal_value = h_nominal.GetBinContent(i)
                variation_value = h_variation.GetBinContent(i);

                uncertainty = variation_value / nominal_value if nominal_value != 0 else 0
                h_varied_reference.SetBinContent(i, h_reference.GetBinContent(i) * uncertainty)

            output.cd()
            h_varied_reference.Write()
        output.Close()
