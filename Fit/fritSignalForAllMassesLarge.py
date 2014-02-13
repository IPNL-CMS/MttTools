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

btag = [1, 2]

for mass, width in masses.items():
    inputFile = "../Extractor2Dataset/datasets/Latest/Signal_ZPrimeToTTJets_M%dGeV_W%sGeV_dataset_nominal_merged.root" % (mass, width)
    for b in btag:
        args = ["./fritSignal", "-m", str(mass), "--b-tag", str(b), "-i", inputFile]
        subprocess.call(args)

