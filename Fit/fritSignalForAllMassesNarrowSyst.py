#!/bin/env python

import subprocess
import json, uuid, sys, os, datetime

masses = {
        500: "5",
        750: "7p5",
        1000: "10",
        1250: "12p5",
        1500: "15",
        2000: "20"
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


