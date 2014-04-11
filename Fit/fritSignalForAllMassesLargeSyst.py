#!/bin/env python

import subprocess
import json, uuid, sys, os, datetime

masses = {
        500: "50",
        750: "75",
        1000: "100",
        1250: "125",
        1500: "150",
        2000: "200"
        }

systs = {
        "JECup": "--jec JECup",
        "JECdown": "--jec JECdown",
        "JERup": "--jer up",
        "JERdown": "--jer down",
        "puUp": "--pileup up",
        "puDown": "--pileup down"
        }

btag = [1, 2]

for mass, width in masses.items():
        for b in btag:
            for syst, param in systs.items():
                inputFile = "../Extractor2Dataset/datasets/Latest/Systematics/Signal_ZPrimeToTTJets_M%dGeV_W%sGeV_dataset_%s_merged.root" % (mass, width, syst)
                args = ["./fritSignal", "-m", str(mass), "--b-tag", str(b), "-i", inputFile, param]
                subprocess.call(args)



