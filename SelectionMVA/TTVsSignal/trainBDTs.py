#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime
from optparse import OptionParser

d = datetime.datetime.now().strftime("%d%b%y")

parser = OptionParser()
parser.add_option("-i", "--input-file", action="store", dest="input", default=False, help="Input file")

parser.add_option("", "--chi2", action="store_true", dest="chi2", default=False, help="Use Chi2 sorting algorithm")
parser.add_option("", "--kf", action="store_true", dest="kf", default=False, help="Use KF sorting algorithm")
(option, args) = parser.parse_args()

outputPath = os.path.join("bdt_trained", d)

sortingAlgoArg = ""
if option.kf:
    sortingAlgoArg = "--kf"
elif option.chi2:
    sortingAlgoArg = "--chi2"

if len(sortingAlgoArg) == 0:
    raise ValueError("Sorting algorithm is not specified")

def launch(input, output_path, output, name):
    if not os.path.exists(input):
        print("Warning input file '%s' does not exist. Skipping job." % input)
        return ""

    args = ["./trainMVA", "-i", input, "--output-path", output_path, "-o", output, "--name", name, "--bdt", sortingAlgoArg, "&>", os.path.join(output_path, "output.log")]

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/tmp/', delete = False)

# Build output tree structure
for btag in [-1, 0, 1, 2]:
    b = "%d-btag" % btag
    if btag == -1:
        b = "all-btag"
    path = os.path.join(outputPath, b)
    try:
        os.makedirs(path)
    except:
        pass

for btag in [-1, 0, 1, 2]:
    b = "%d-btag" % btag
    if btag == -1:
        b = "all-btag"
    path = os.path.join(outputPath, b)
    inputPath = os.path.join(option.input, b)
    outputFile = os.path.join(path, "BDT_trained.root")
    tmpfile.write(launch(os.path.join(inputPath, "mva_trees.root"), path, outputFile, "BDT_%s" % b) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "15"]
subprocess.call(args)
