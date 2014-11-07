#! /bin/env python

import os, subprocess, datetime, tempfile, sys, math, json
from optparse import OptionParser

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

import sys
sys.argv.append( '-b' )

from ROOT import TFile, TCanvas, TColor, TLegend

backgrounds = ["DATA", "TT", "single_antitop", "single_top", "zjets", "wjets", "dibosons", "H400_scalar", "H500_scalar", "H600_scalar", "H700_scalar", "H800_scalar", "H400_pseudoscalar", "H500_pseudoscalar", "H600_pseudoscalar", "H700_pseudoscalar", "H800_pseudoscalar"]

toMerge = [
        {
            "inputs": ["single_top", "single_antitop"],
            "output": "st"
        }
        ]

final_backgrounds = ["H400_scalar", "H500_scalar", "H600_scalar", "H700_scalar", "H800_scalar", "H400_pseudoscalar", "H500_pseudoscalar", "H600_pseudoscalar", "H700_pseudoscalar", "H800_pseudoscalar", "LINE", "dibosons", "st", "wjets", "zjets", "TT"]

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

def getChannelSyst(category, btag, background, syst, syst_variation):
    channel = getChannel(category, btag, background)

    # FIXME: Is this really what we want?
    if (background == "wjets" or background == "zjets") and (syst == "matching" or syst == "scale"):
        return ""

    syst_channel = "%s__%s__%s" % (channel, syst, syst_variation)
    return syst_channel

def drawRatio(nominal, up, down, syst, filename):
    c = TCanvas("c", "c", 800, 800)
    c.SetLeftMargin(0.15)

    ratio_up = up.Clone()
    ratio_up.Divide(nominal)
    ratio_up.SetStats(False)
    ratio_up.SetTitle("")

    ratio_down = down.Clone()
    ratio_down.Divide(nominal)

    ratio_up.SetLineColor(TColor.GetColor('#C02942'))
    ratio_down.SetLineColor(TColor.GetColor('#53777A'))

    ratio_up.Draw("hist L")
    ratio_down.Draw("hist L same")

    max = ratio_up.GetMaximum() * 1.02
    min = ratio_down.GetMinimum() * 0.98

    ratio_up.GetYaxis().SetTitle("ratio")
    ratio_up.GetYaxis().SetTitleOffset(2)
    ratio_up.GetXaxis().SetNdivisions(505)
    ratio_up.GetXaxis().SetTitle("m_{t#bar{t}}")
    ratio_up.SetMaximum(max)
    ratio_up.SetMinimum(min)

    l = TLegend(0.67, 0.79, 0.87, 0.89)
    l.SetBorderSize(0)
    l.SetTextFont(42)
    l.SetFillColor(0)

    l.AddEntry(ratio_up, "%s up" % syst, "L")
    l.AddEntry(ratio_down, "%s down" % syst, "L")

    l.Draw()

    c.Print(filename)
    pass

f = TFile.Open(option.input)

histos = {}
for category in categories:
    for btag in btags:
        for background in backgrounds:
            h = f.Get(getChannel(category, btag, background))

            histos[(category, btag, background)] = h

            for syst in systs:
                up = f.Get(getChannelSyst(category, btag, background, syst, "up"))
                down = f.Get(getChannelSyst(category, btag, background, syst, "down"))

                if up == None or down == None:
                    continue

                histos[(category, btag, background, syst, "up")] = up
                histos[(category, btag, background, syst, "down")] = down

for mergeData in toMerge:
    for category in categories:
        for btag in btags:
            for input in mergeData["inputs"]:

                if (category, btag, mergeData["output"]) in histos:
                    histos[(category, btag, mergeData["output"])].Add(histos[(category, btag, background)])
                else:
                    histos[(category, btag, mergeData["output"])] = histos[(category, btag, background)].Clone()

                for syst in systs:

                    if not (category, btag, background, syst, "up") in histos:
                        continue

                    if (category, btag, mergeData["output"], syst, "up") in histos:
                        histos[(category, btag, mergeData["output"], syst, "up")].Add(histos[(category, btag, background, syst, "up")])
                    else:
                        histos[(category, btag, mergeData["output"], syst, "up")] = histos[(category, btag, background, syst, "up")].Clone()

                    if (category, btag, mergeData["output"], syst, "down") in histos:
                        histos[(category, btag, mergeData["output"], syst, "down")].Add(histos[(category, btag, background, syst, "down")])
                    else:
                        histos[(category, btag, mergeData["output"], syst, "down")] = histos[(category, btag, background, syst, "down")].Clone()


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


for category in categories:
    for btag in btags:
        #for background in final_backgrounds:
        for background in ['TT']:
            for syst in systs:

                if not (category, btag, background, syst, "up") in histos:
                    continue

                filename = "%s_%s_%db_%s.pdf" % (background, category, btag, syst)
                drawRatio(histos[(category, btag, background)], histos[(category, btag, background, syst, "up")], histos[(category, btag, background, syst, "down")], syst, filename)
