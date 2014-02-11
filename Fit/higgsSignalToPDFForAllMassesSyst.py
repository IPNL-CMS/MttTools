#!/bin/env python

import subprocess

masses = range(400, 801, 100)
btag = [1, 2]

systs = {
        "JECup": "--jec JECup",
        "JECdown": "--jec JECdown",
        "JERup": "--jer up",
        "JERdown": "--jer down",
        "puUp": "--pileup up",
        "puDown": "--pileup down"
        }

for mass in masses:
    for b in btag:
        for syst, param in systs.items():
            input = "../Extractor2Dataset/datasets/Latest/Systematics/Signal_S0_S_i_M%d_cpl1_pseudoscalar_histos_%s_merged.root" % (mass, syst)
            args = ["./higgsSignalToPDF", "-m", str(mass), "--b-tag", str(b), "-i", input, param]
            subprocess.call(args)
