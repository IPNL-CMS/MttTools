#!/bin/env python

import subprocess

masses = [500, 750, 1000, 1250, 1500, 2000]

for mass in masses:
    args = ["./fitMtt", "-m", str(mass), "-i", "../Extractor2Dataset/datasets/Latest/data/MTT_Data_merged.root", "--b-tag", "3", "--combine"]
    subprocess.call(args)

