#!/usr/bin/env python

import os, subprocess, datetime

d = datetime.datetime.now().strftime("%d%b%y")

inputs = [
        ["MTT_Muons_2012_%s.root" % d, "lists/MTT_SingleMu_Run2012.list", "semimu"],
        ["MTT_Electrons_2012_%s.root" % d, "lists/MTT_SingleElectron_Run2012.list", "semie"]
        ]

jobs = []

def launch(input, output, type, btag):
    args = ["./extractorToHisto", "--input-list", input, "-o", output, "--data", "--%s" % type, "--b-tag", str(btag)]
    p = subprocess.Popen(args)
    jobs.append(p)

# Build output tree structure
for btag in [1, 2]:
    path = "plots/%s/%d-btag/data" % (d, btag)
    try:
        os.makedirs(path)
    except:
        pass

print("Extracting dataset ...")
for input in inputs:
    for btag in [1, 2]:
        path = "plots/%s/%d-btag/data" % (d, btag)
        launch(input[1], os.path.join(path, input[0]), input[2], btag)

for job in jobs:
    job.wait()

## All is done, merge

#print("Merging ...")
#for btag in [1, 2]:
    #args = ["hadd", "MTT_Data_merged_2012_%d_btag_%s.root" % (btag, d)]
    #for output in inputs:
        #args.append(output[0])

#subprocess.call(args)

#for output in inputs:
#  os.remove(output[0])
