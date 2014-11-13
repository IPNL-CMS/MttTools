#! /usr/bin/env python

from optparse import OptionParser
parser = OptionParser()

parser.add_option("-i", "--input", dest="input")

(options, args) = parser.parse_args()

import scipy, scipy.stats
def Z_to_p(z_value):
    """
    Convert a Z-value to a p-value; inverse of :meth:`p_to_Z`
    """
    return scipy.stats.norm.sf(z_value)

# Standard Imports and calculators
import json, re
import ROOT
import array, sys, numpy
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
ROOT.gROOT.ProcessLine("setTDRStyle();")

with open(options.input) as json_data:
    dnll_ratios = json.load(json_data)

expected = dnll_ratios["expected_deltall_ratios"]
toy_ratios = dnll_ratios["expected_toy_1"]

expected_sorted = expected
expected_sorted.sort()

toy_ratios_sorted = toy_ratios
toy_ratios_sorted.sort()

min_ = min(expected)
max = max(expected)

max_95 = expected_sorted[int(0.998 * len(expected))]

h = ROOT.TH1F("expected", "expected", 50, min_, max_95)
h.SetLineColor(ROOT.TColor.GetColor("#542437"))
h.GetXaxis().SetTitle("q")
h.GetYaxis().SetTitleOffset(1.5)
h.GetYaxis().SetTitle("Number of pseudo-experiments")

q_data = dnll_ratios["observed_deltall_ratios"][0]

min_toy = min(toy_ratios)
max_95_toy = toy_ratios_sorted[int(0.998 * len(toy_ratios))]
toy = ROOT.TH1F("toy", "toy", 50, min_toy, max_95_toy)
toy.SetLineColor(ROOT.TColor.GetColor("#542437"))
toy.GetXaxis().SetTitle("q")
toy.GetYaxis().SetTitleOffset(1.5)
toy.GetYaxis().SetTitle("Number of pseudo-experiments")

for dnll in expected:
    h.Fill(dnll)

for dnll in toy_ratios:
    toy.Fill(dnll)

areaq = array.array('d', [0.5, 0.16, 0.84, 0.025, 0.975])
medianq = array.array('d', [0] * len(areaq))
h.GetQuantiles(len(areaq), medianq, areaq);

x = h.GetXaxis().GetBinUpEdge(h.FindBin(medianq[0]))
line = ROOT.TLine(x, 0, x, h.GetBinContent(h.FindBin(medianq[0])))
line.SetLineColor(ROOT.TColor.GetColor("#3B8686"))
line.SetLineWidth(3)

h_68 = h.Clone()
h_95 = h.Clone()
for i in range(1, h_68.GetNbinsX() + 1):
    low = h_68.GetXaxis().GetBinLowEdge(i)
    up = h_68.GetXaxis().GetBinUpEdge(i)

    if (low < medianq[1]) or (up > medianq[2]):
        h_68.SetBinContent(i, 0)
        h_68.SetBinError(i, 0)

    if (low < medianq[3]) or (up > medianq[4]):
        h_95.SetBinContent(i, 0)
        h_95.SetBinError(i, 0)

h_68.SetFillStyle(1001)
h_68.SetFillColor(ROOT.TColor.GetColor("#79BD9A"))
h_68.SetLineColor(h_68.GetFillColor())

h_95.SetFillStyle(1001)
h_95.SetFillColor(ROOT.TColor.GetColor("#A8DBA8"))
h_95.SetLineColor(h_68.GetFillColor())

c = ROOT.TCanvas("c", "c", 800, 800)

lumi = ROOT.TPaveText(c.GetLeftMargin(), 1 - 0.5 * c.GetTopMargin(), 1 - c.GetRightMargin(), 1, "brNDC")
lumi.SetFillStyle(0)
lumi.SetBorderSize(0)
lumi.SetMargin(0)
lumi.SetTextSize(0.6 * c.GetTopMargin())
lumi.SetTextFont(42)
lumi.SetTextAlign(33)
lumi.AddText("19.67 fb^{-1} (8 TeV)")

title = ROOT.TPaveText(c.GetLeftMargin(), 1 - 0.5 * c.GetTopMargin(), 1 - c.GetRightMargin(), 1, "brNDC")
title.SetFillStyle(0)
title.SetBorderSize(0)
title.SetMargin(0)
title.SetTextSize(0.75 * c.GetTopMargin())
title.SetTextFont(62)
title.SetTextAlign(13)
title.AddText("CMS #font[52]{#scale[0.76]{Preliminary}}")

h.Draw()
h_95.Draw("hist same")
h_68.Draw("hist same")
h.Draw("same")
h.Draw("same axis")
line.Draw("same")
title.Draw()
lumi.Draw()

c.SetLogy()
c.Print("dnll_expected.pdf")

fill = toy.Clone()
fill_data = toy.Clone()
x = toy.GetXaxis().GetBinUpEdge(toy.FindBin(medianq[0]))
line = ROOT.TLine(x, 0, x, toy.GetBinContent(toy.FindBin(medianq[0])))
line.SetLineColor(ROOT.TColor.GetColor("#3B8686"))
line.SetLineWidth(3)
x_data = toy.GetXaxis().GetBinUpEdge(toy.FindBin(q_data))
for i in range(1, fill.GetNbinsX() + 1):
    up = fill.GetXaxis().GetBinUpEdge(i)

    if up <= x:
        fill.SetBinContent(i, 0)
        fill.SetBinError(i, 0)

    if up <= x_data:
        fill_data.SetBinContent(i, 0)
        fill_data.SetBinError(i, 0)

fill.SetFillStyle(1001)
fill.SetFillColor(ROOT.TColor.GetColor("#79BD9A"))
fill.SetLineColor(fill.GetFillColor())

fill_data.SetFillStyle(3004)
fill_data.SetFillColor(ROOT.TColor.GetColor("#490A3D"))
fill_data.SetLineWidth(0)

toy.Draw()
fill.Draw("hist same")
toy.Draw("same")
line.Draw("same")

toy.Draw("axis same")

line_data = ROOT.TLine(x_data, 0, x_data, toy.GetBinContent(toy.FindBin(q_data)))
line_data.SetLineColor(ROOT.TColor.GetColor("#490A3D"))
line_data.SetLineStyle(ROOT.kDashed)
line_data.SetLineWidth(3)
line_data.Draw("same")
fill_data.Draw("hist same")
title.Draw()
lumi.Draw()

c.Print("dnll_toy.pdf")
