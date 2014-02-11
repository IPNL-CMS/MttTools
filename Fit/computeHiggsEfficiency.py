#!/bin/env python

import subprocess

masses = range(400, 801, 100)
btag = [1, 2]

for mass in masses:
    for b in btag:
        args = ["./computeEfficiencies", "-m", str(mass), "--b-tag", str(b), "--gen-file", "data/s0_pseudoscalar_gen.json"]
        subprocess.call(args)
