#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

from optparse import OptionParser

parser = OptionParser()
parser.add_option("", "--from", dest="from_", type="string", help="The directory where are stored the ROOT files")

(options, args) = parser.parse_args()

if options.from_ is None or not os.path.isdir(options.from_):
  parser.error("You must specify a valid path using --from")

d = datetime.datetime.now().strftime("%d%b%y")

background = "HTT_MC_TT_histos_%s_nominal.root"

files = [
        ["HTT_TT_plus_S0_minus_TT_M400_coupling_5_histos_%s_nominal.root", "HTT_BS_TT_plus_S0_M400_coupling_5_histos_%s_nominal.root"]
        ]

def launch(file1, file2, output, btag):
    args = ["./substractHistos", "-o", output, "--mc", "--b-tag", str(btag)]
    if "semie" in input:
        args.append("--semie")
    elif "semimu" in input:
        args.append("--semimu")

    args.append(file1)
    args.append(file2)

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

print("Substracting histograms...")

for file in files:
    for btag in [1, 2]:
        for type in ["semie", "semimu"]:
            input = os.path.join(options.from_, "%d-btag/%s" % (btag, type))
            tmpfile.write(launch(os.path.join(input, file[1] % type), os.path.join(input, background % type), os.path.join(input, file[0] % type), btag) + "\n");

tmpfile.flush()

#print tmpfile.name
args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 
subprocess.call(args)
