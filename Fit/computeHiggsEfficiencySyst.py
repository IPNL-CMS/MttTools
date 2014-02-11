#!/bin/env python

import subprocess

masses = range(400, 801, 100)
btag = [1, 2]

systs = ["jec", "jer", "pu"]

for mass in masses:
    for b in btag:
        for syst in systs:
            for sign in ["up", "down"]:
                args = ["./computeEfficiencies", "-m", str(mass), "--b-tag", str(b), "--gen-file", "data/s0_pseudoscalar_gen.json", "--%s" % syst, sign]
                subprocess.call(args)
