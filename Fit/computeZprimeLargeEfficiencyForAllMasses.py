#!/bin/env python

import subprocess

masses = [500, 750, 1000, 1250, 1500, 2000]
btag = [1, 2]

for mass in masses:
    for b in btag:
        args = ["./computeEfficiencies", "-m", str(mass), "--b-tag", str(b), "--gen-file", "data/zprime_large_gen.json"]
        subprocess.call(args)




