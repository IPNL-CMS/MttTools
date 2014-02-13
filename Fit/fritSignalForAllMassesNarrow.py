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

btag = [1, 2]

for mass, width in masses.items():
    inputFile = "../Extractor2Dataset/datasets/Latest/Signal_ZPrimeToTTJets_M%dGeV_W%sGeV_dataset_nominal_merged.root" % (mass, width)
    for b in btag:
        args = ["./fritSignal", "-m", str(mass), "--b-tag", str(b), "-i", inputFile]
        subprocess.call(args)

