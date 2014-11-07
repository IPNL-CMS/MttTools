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

backgrounds = ["DATA", "TT", "single_antitop", "single_top", "zjets", "wjets", "dibosons"]

toMerge = [
        {
            "inputs": ["single_top", "single_antitop"],
            "output": "st"
        }
        ]

final_backgrounds = ["dibosons", "st", "wjets", "zjets", "TT"]

scales = {
        'TT': 1,
        'single_top': 1,
        'single_antitop': 1,
        'wjets': 1,
        'zjets': 1,
        'dibosons': 1
        }

categories = ["e", "mu"]
btags = [1, 2]

systs = ["jec", "jer", "pu", "btag", "lept", "trig", "matching", "scale"]

def getBTagName(btag):
    b = "%dbtag" % btag
    if btag == -1:
        b = "allbtag"

    return b

def getChannel(category, btag, background):
    return "mtt_%s_%s__%s" % (category, getBTagName(btag), background)

def getEvents(background, h):
    entries_pos = 0
    entries_neg = 0
    for i in range(1, h.GetNbinsX() + 1):
        c = h.GetBinContent(i)
        if c > 0:
            entries_pos += c
        else:
            entries_neg += c

    n = abs(entries_pos) + abs(entries_neg)

    if background in scales:
        n *= scales[background]

    return n

def getEventsSyst(f, category, btag, background, syst, syst_variation):
    channel = getChannel(category, btag, background)

    # FIXME: Is this really what we want?
    if (background == "wjets" or background == "zjets") and (syst == "matching" or syst == "scale"):
        return 0

    syst_channel = "%s__%s__%s" % (channel, syst, syst_variation)
    
    h = f.Get(syst_channel)
    if not h:
        return 0

    return getEvents(background, h)

f = TFile.Open(option.input)

events = {}
sys_errors = {}
for category in categories:
    for btag in btags:
        for background in backgrounds:
            h = f.Get(getChannel(category, btag, background))

            # Get number of events in mtt histo
            event = getEvents(background, h)
            events[(category, btag, background)] = event

            sys_error = 0
            for syst in systs:
                up = getEventsSyst(f, category, btag, background, syst, "up")
                down = getEventsSyst(f, category, btag, background, syst, "down")

                if up == 0 or down == 0:
                    continue

                error = (math.fabs(up - event) + math.fabs(event - down)) / 2
                sys_errors[(category, btag, background, syst)] = error

                sys_error = math.sqrt( sys_error**2 + error**2 )

            sys_errors[(category, btag, background)] = sys_error


for mergeData in toMerge:
    for category in categories:
        for btag in btags:
            event = 0
            sys_error = 0
            for input in mergeData["inputs"]:
                event += events[(category, btag, input)]
                del events[(category, btag, input)]

                # Merge systematic errors
                sys_error = math.sqrt(sys_error**2 + sys_errors[(category, btag, input)]**2)
                del sys_errors[(category, btag, input)]

            events[(category, btag, mergeData["output"])] = event
            sys_errors[(category, btag, mergeData["output"])] = sys_error

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

def getPrettySysName(name):
    if name == "jec":
        return "Jet Energy Scale"
    elif name == "jer":
        return "Jet Energy Resolution"
    elif name == "pu":
        return "Pileup"
    elif name == "btag":
        return "b-tagging"
    elif name == "lept":
        return "Lepton id and iso"
    elif name == "scale":
        return "Scale"
    elif name == "matching":
        return "Matching"
    elif name == "trig":
        return "Trigger"

    return name

# Print everything
events[("e", 1)] = 0
sys_errors[("e", 1)] = 0
events[("e", 2)] = 0
sys_errors[("e", 2)] = 0
events[("mu", 1)] = 0
sys_errors[("mu", 1)] = 0
events[("mu", 2)] = 0
sys_errors[("mu", 2)] = 0

for syst in systs:
    sys.stdout.write(' %s & ' % getPrettySysName(syst))
    for category in categories:
        for btag in btags:
            syst_percent = sys_errors[(category, btag, 'TT', syst)] / events[(category, btag, 'TT')] * 100
            sys.stdout.write("$%.1f\\,\\%%$" % (syst_percent))

            if category != "mu" or btag != 2:
                sys.stdout.write(" & ")


    print(" \\\\ ")

print("\midrule")
sys.stdout.write(" Total: & ")
for category in categories:
    for btag in btags:
        syst_percent = sys_errors[(category, btag, 'TT')] / events[(category, btag, 'TT')] * 100
        sys.stdout.write('$%.1f\\,\\%%$' % syst_percent)

        if category != "mu" or btag != 2:
            sys.stdout.write(" & ")

print(" \\\\ ")
print("\\bottomrule")
