#! /bin/env python

import os, subprocess, datetime, tempfile, sys, math, json

import sys
sys.argv.append( '-b' )

from ROOT import TH1F, TCanvas, TColor, TBox, gROOT, TPaveText

gROOT.ProcessLine(".L tdrstyle.cc")
gROOT.ProcessLine("setTDRStyle();")

rates = {
    'w_rate': -0.6919109486304862,
    'ttbar_rate': 0.8137285421522549,
    'st_rate': (-1.068182582563093 + -1.2071352088945098) / 2,
    'z_rate': -0.5477218588756969,
    'lumi': -0.135092913449828,
#    'sat_rate': -1.2071352088945098,
    'diboson_rate': -0.26225554316630095,
    }

def getPrettyName(name):
    if name == "w_rate":
        return "W + jets"
    elif name == "ttbar_rate":
        return "tt"
    elif name == "st_rate" or name == "sat_rate":
        return "Single top"
    elif name == "lumi":
        return "Lumi"
    elif name == "z_rate":
        return "Z + jets"
    elif name == "diboson_rate":
        return "di-bosons"

h = TH1F("h", "h", 6, 0, 6)

index = 1
for name, rate in rates.items():
    h.GetXaxis().SetBinLabel(index, getPrettyName(name))
    h.SetBinContent(index, rate)

    index += 1

h.SetMaximum(3)
h.SetMinimum(-3)
h.GetYaxis().SetTitleOffset(1.2)
h.GetXaxis().SetTitleOffset(1.4)
h.GetYaxis().SetTitleSize(0.04)
h.GetXaxis().SetTitleSize(0.04)
h.GetYaxis().SetLabelSize(0.04)
h.GetYaxis().SetLabelOffset(0.02)
h.GetXaxis().SetLabelOffset(0.01)
h.GetXaxis().SetLabelSize(0.04)
#h.SetMarkerStyle(21)
#h.SetMarkerSize(1)
h.SetLineColor(TColor.GetColor('#031634'))

h.GetYaxis().SetTitle("Deviation in units of #sigma")
h.GetXaxis().SetTitle("Nuisance parameters")

c = TCanvas("c", "c", int(800 * 1.6180339887), 800)
c.SetLeftMargin(0.10)
c.SetTicky(1)
c.SetTickx(1)
h.Draw("L")

two_sigma = TBox(0, -2, 6, 2)
two_sigma.SetFillStyle(1001)
two_sigma.SetFillColor(TColor.GetColor("#A8DBA8"))
two_sigma.Draw()

one_sigma = TBox(0, -1, 6, 1)
one_sigma.SetFillStyle(1001)
one_sigma.SetFillColor(TColor.GetColor("#79BD9A"))
one_sigma.Draw()

h.Draw("L same")
h.Draw("same axis")

title = TPaveText(c.GetLeftMargin(), 1 - c.GetTopMargin() + 0.005, 1 - c.GetRightMargin(), 1, "brNDC")
title.SetFillStyle(0)
title.SetBorderSize(0)
title.SetMargin(0)
title.SetTextSize(0.03)
title.SetTextFont(42)
title.SetTextAlign(32)
title.AddText("CMS preliminary, L = 19.67 fb^{-1}, #sqrt{s} = 8 TeV")
title.Draw()

c.SaveAs("normalization_systs.pdf")
