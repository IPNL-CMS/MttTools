#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
        # Background
        ["HTT_MC_TT_histos_%s_nominal.root", "MC/HTT_TT_%s.list"],

        # Background + Signal
        ["HTT_BS_TT_plus_S0_M400_coupling_5_histos_%s_nominal.root", "MC/HTT_TT_plus_S0_M400_coupling_5_%s.list"]
        ]

def launch(input, output, btag):
    args = ["./extractorToHisto", "--input-list", input, "-o", output, "--mc", "--b-tag", str(btag)]
    if "semie" in input:
        args.append("--semie")
    elif "semimu" in input:
        args.append("--semimu")

    #args.append("--weight %.15f" % weight)

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

# Build output tree structure
for btag in [1, 2]:
    for type in ["semie", "semimu"]:
        path = "plots/%s/%d-btag/%s" % (d, btag, type)
        try:
            os.makedirs(path)
        except:
            pass

print("Extracting datasets...")

for file in files:
    for btag in [1, 2]:
        for type in ["semie", "semimu"]:
            path = "plots/%s/%d-btag/%s" % (d, btag, type)
            tmpfile.write(launch(file[1] % type, "%s/%s" % (path, file[0] % type), btag) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 
subprocess.call(args)
