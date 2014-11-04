#! /bin/env python

import os, subprocess, datetime, tempfile, sys, math, json
from optparse import OptionParser
from ROOT import TFile

def check_output(*popenargs, **kwargs):
    import subprocess
    r"""Run command with arguments and return its output as a byte string.
 
    Backported from Python 2.7 as it's implemented as pure python on stdlib.
 
    >>> check_output(['/usr/bin/python', '--version'])
    Python 2.6.2
    """

    env = os.environ.copy()
    if 'TERM' in env:
        del env['TERM']

    process = subprocess.Popen(stdout=subprocess.PIPE, env=env, *popenargs, **kwargs)
    output, unused_err = process.communicate()
    retcode = process.poll()
    if retcode:
        cmd = kwargs.get("args")
        if cmd is None:
            cmd = popenargs[0]
        error = subprocess.CalledProcessError(retcode, cmd)
        error.output = output
        raise error
    return output

parser = OptionParser()
parser.add_option("-i", "--input", action="store", dest="input", help="Input root file")
(option, args) = parser.parse_args()

if not option.input:
    parser.error("--input argument is required")

backgrounds = ["DATA", "TT", "single_antitop", "single_top", "zjets", "wjets", "dibosons", "H400_scalar", "H500_scalar", "H600_scalar", "H700_scalar", "H800_scalar", "H400_pseudoscalar", "H500_pseudoscalar", "H600_pseudoscalar", "H700_pseudoscalar", "H800_pseudoscalar"]

toMerge = [
        {
            "inputs": ["single_top", "single_antitop"],
            "output": "st"
        }
        ]

final_backgrounds = ["H400_scalar", "H500_scalar", "H600_scalar", "H700_scalar", "H800_scalar", "H400_pseudoscalar", "H500_pseudoscalar", "H600_pseudoscalar", "H700_pseudoscalar", "H800_pseudoscalar", "LINE", "dibosons", "st", "wjets", "zjets", "TT"]

categories = ["e", "mu"]
btags = [1, 2]

def getBTagName(btag):
    b = "%dbtag" % btag
    if btag == -1:
        b = "allbtag"

    return b

def getChannel(category, btag, background):
    return "mtt_%s_%s__%s" % (category, getBTagName(btag), background)

def getEvents(h):
    entries_pos = 0
    entries_neg = 0
    for i in range(1, h.GetNbinsX() + 1):
        c = h.GetBinContent(i)
        if c > 0:
            entries_pos += c
        else:
            entries_neg += c

    return abs(entries_pos) + abs(entries_neg)

f = TFile.Open(option.input)

events = {}
for category in categories:
    for btag in btags:
        for background in backgrounds:
            h = f.Get(getChannel(category, btag, background))

            # Get number of events in mtt histo
            event = getEvents(h)
            events[(category, btag, background)] = event


for mergeData in toMerge:
    for category in categories:
        for btag in btags:
            event = 0
            for input in mergeData["inputs"]:
                event += events[(category, btag, input)]
                del events[(category, btag, input)]

            events[(category, btag, mergeData["output"])] = event

def getPrettyName(name):

    if name[0] is 'H':
        mass = name[1:4]
        type = name[5:]
        return "\\sz (m = %s\\,\\GeV, %s)" % (
          mass,
          "scalar" if type == "scalar" else "pseudo-scalar"
        )

    if name is "TT":
        return "\\ttbar"
    elif name is "st":
        return "single top"
    elif name is "wjets":
        return "W + jets"
    elif name is "zjets":
        return "Z + jets"
    elif name is "dibosons":
        return "di-bosons"

    return name

# Print everything
events[("e", 1)] = 0
events[("e", 2)] = 0
events[("mu", 1)] = 0
events[("mu", 2)] = 0
for bkg in final_backgrounds:
    if bkg is 'LINE':
        print("\\midrule")
        continue

    print("%s & %.2f & %.2f & %.2f & %.2f \\\\" % (
        getPrettyName(bkg),
        events[("e", 1, bkg)],
        events[("e", 2, bkg)],
        events[("mu", 1, bkg)],
        events[("mu", 2, bkg)]
        ))

    events[("e", 1)] += events[("e", 1, bkg)]
    events[("e", 2)] += events[("e", 2, bkg)]
    events[("mu", 1)] += events[("mu", 1, bkg)]
    events[("mu", 2)] += events[("mu", 2, bkg)]
    
print("\\midrule")
print("Total & %.2f & %.2f & %.2f & %.2f \\\\" % (
    events[("e", 1)],
    events[("e", 2)],
    events[("mu", 1)],
    events[("mu", 2)]
    ))

print("\\midrule")
print("Data & %.2f & %.2f & %.2f & %.2f \\\\" % (
    events[("e", 1, "DATA")],
    events[("e", 2, "DATA")],
    events[("mu", 1, "DATA")],
    events[("mu", 2, "DATA")]
    ))
print("\\bottomrule")
