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

from ROOT import TFile, TCanvas, TColor, TLegend, TLine, gPad, gROOT, TPaveText
gROOT.ProcessLine(".L tdrstyle.cc")
gROOT.ProcessLine("setTDRStyle();")

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

def getMax(a, b):
    return max(a.GetMaximum(), b.GetMaximum())

def getMin(a, b):
    return min(a.GetMinimum(), b.GetMinimum())

def drawRatio(nominal, up, down, syst, filename):

    # Clear first 3 bins
    for bin in range(1, 4):
        nominal.SetBinContent(bin, 1)
        up.SetBinContent(bin, 1)
        down.SetBinContent(bin, 1)

    c = TCanvas("c", "c", 800, 800)
    #c.SetLeftMargin(0.15)

    ratio_up = up.Clone()
    ratio_up.Divide(nominal)
    ratio_up.SetStats(False)
    ratio_up.SetTitle("")

    ratio_down = down.Clone()
    ratio_down.Divide(nominal)

    ratio_up.SetLineColor(TColor.GetColor('#C02942'))
    ratio_up.SetLineWidth(1)
    ratio_down.SetLineColor(TColor.GetColor('#031634'))
    ratio_down.SetLineStyle(2);
    ratio_down.SetLineWidth(1)

    ratio_up.Draw("hist")

    line = TLine(ratio_up.GetBinLowEdge(1), 1, ratio_up.GetXaxis().GetBinUpEdge(ratio_up.GetNbinsX()), 1)
    line.SetLineColor(TColor.GetColor("#bbbbbb"))
    line.SetLineWidth(1)
    line.SetLineStyle(1)
    line.Draw()

    ratio_up.Draw("hist same")
    ratio_down.Draw("hist same")

    max_ = getMax(ratio_up, ratio_down) * 1.02
    min_ = getMin(ratio_up, ratio_down) * 0.98

    ratio_up.GetYaxis().SetTitle("ratio (#pm 1 #sigma / nominal)")
    ratio_up.GetYaxis().SetTitleOffset(2)
    ratio_up.GetYaxis().SetLabelOffset(0.02)
    ratio_up.GetXaxis().SetTitleOffset(1.5)
    ratio_up.GetXaxis().SetNdivisions(505)
    ratio_up.GetXaxis().SetTitle("m_{t#bar{t}}")
    ratio_up.SetMaximum(max_)
    ratio_up.SetMinimum(min_)
    
    gPad.RedrawAxis()

    l = TLegend(0.60, 0.79, 0.87, 0.89)
    l.SetBorderSize(0)
    l.SetTextFont(42)
    l.SetFillColor(0)
    l.SetTextSize(0.03)

    l.AddEntry(ratio_up, "%s up" % getPrettySysName(syst), "L")
    l.AddEntry(ratio_down, "%s down" % getPrettySysName(syst), "L")

    l.Draw()

    lumi = TPaveText(c.GetLeftMargin(), 1 - 0.5 * c.GetTopMargin(), 1 - c.GetRightMargin(), 1, "brNDC")
    lumi.SetFillStyle(0)
    lumi.SetBorderSize(0)
    lumi.SetMargin(0)
    lumi.SetTextSize(0.6 * c.GetTopMargin())
    lumi.SetTextFont(42)
    lumi.SetTextAlign(33)
    lumi.AddText("19.67 fb^{-1} (8 TeV)")
    lumi.Draw()

    title = TPaveText(c.GetLeftMargin(), 1 - 0.5 * c.GetTopMargin(), 1 - c.GetRightMargin(), 1, "brNDC")
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    title.SetMargin(0)
    title.SetTextSize(0.75 * c.GetTopMargin())
    title.SetTextFont(62)
    title.SetTextAlign(13)
    title.AddText("CMS #font[52]{#scale[0.76]{Preliminary}}")
    title.Draw()

    c.Print(filename)
    pass

f = TFile.Open(option.input)

try:
    os.makedirs("systematics_ratios")
except:
    pass
    
histos = {}
for category in categories:
    for btag in btags:
        for background in backgrounds:
            h = f.Get(getChannel(category, btag, background))

            histos[(category, btag, background)] = h

            for syst in systs:

                if syst == "lept":
                    syst = "lept_" + category

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

                    if syst == "lept":
                        syst = "lept_" + category

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

def getPrettySysName(name):
    if name == "jec":
        return "JES"
    elif name == "jer":
        return "JER"
    elif name == "pu":
        return "Pile-up"
    elif name == "btag":
        return "b-tagging"
    elif "lept" in name:
        return "Lepton id"
    elif name == "scale":
        return "Scale"
    elif name == "matching":
        return "Matching"
    elif name == "trig":
        return "Trigger"

    return name


for category in categories:
    for btag in btags:
        #for background in final_backgrounds:
        for background in ['TT']:
            for syst in systs:

                if syst == "lept":
                    syst = "lept_" + category

                if not (category, btag, background, syst, "up") in histos:
                    print("Systematic '%s' not found for category %s:%d" % (syst, category, btag))
                    continue

                filename = "%s_%s_%db_%s.pdf" % (background, category, btag, syst)
                drawRatio(histos[(category, btag, background)], histos[(category, btag, background, syst, "up")], histos[(category, btag, background, syst, "down")], syst, os.path.join("systematics_ratios", filename))
