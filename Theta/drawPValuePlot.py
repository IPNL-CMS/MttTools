#!/usr/bin/env python

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
    pvalues_all = json.load(json_data)

graph68  = ROOT.TGraphAsymmErrors()
#graph68.SetFillColor(ROOT.TColor.GetColor("#ECD078"))
graph68.SetFillColor(ROOT.TColor.GetColor("#79BD9A"))


graph95  = ROOT.TGraphAsymmErrors()
#graph95.SetFillColor(ROOT.TColor.GetColor("#CFF09E"))
graph95.SetFillColor(ROOT.TColor.GetColor("#A8DBA8"))


graphMed = ROOT.TGraphErrors()
graphMed.SetLineColor(ROOT.TColor.GetColor("#3B8686"))
graphMed.SetLineWidth(3)

graphObs = ROOT.TGraphErrors()
graphObs.SetLineColor(ROOT.TColor.GetColor("#0B486B"))
graphObs.SetLineWidth(3)

options.expectedOnly = True

index = 0
for mass in sorted(pvalues_all):
    regex = r'(\d+)'
    m = float(re.search(regex, mass).group(0))

    pvalues = pvalues_all[mass]

    graphMed.SetPoint(index, m, float(Z_to_p(pvalues["median"][0])))
    graph68.SetPoint(index, m, float(Z_to_p(pvalues["median"][0])))
    graph95.SetPoint(index, m, float(Z_to_p(pvalues["median"][0])))
    if pvalues["observed"][0] is not None:
        graphObs.SetPoint(index, m, float(Z_to_p(pvalues["observed"][0])))
        options.expectedOnly = False

    median = Z_to_p(pvalues["median"][0])
    median_up = Z_to_p(pvalues["median_up"][0])
    median_down = Z_to_p(pvalues["median_low"][0])

    median_up_2sigma = Z_to_p(pvalues["median_up_2sigma"][0])
    median_down_2sigma = Z_to_p(pvalues["median_low_2sigma"][0])

    median_s = pvalues["median"][0]
    median_up_s = pvalues["median_up"][0]
    median_down_s = pvalues["median_low"][0]

    # Significance and p-values are ordered in the opposite way
    error_up = median_down - median
    error_low = median - median_up

    error_up_95 = median_down_2sigma - median
    error_low_95 = median - median_up_2sigma

    graph68.SetPointError(index, 0, 0, error_low, error_up)
    graph95.SetPointError(index, 0, 0, error_low_95, error_up_95)
    #graph68.SetPointError(index, 0, 0, float(Z_to_p(pvalues["median_down"][0])), float(Z_to_p(pvalues["median_up"][0])))

    index += 1

MG = ROOT.TMultiGraph()
MG.Add(graph95)
MG.Add(graph68)
MG.Add(graphMed)
if not options.expectedOnly:
    MG.Add(graphObs)

FONTSIZE = 0.030
MINPV = 1 * 1e-5
MAXPV = 1.5
Lines = [1.,2.,3.,4.,5.]
MINMH = 400
MAXMH = 800

def MakePvalPlot(MG):

    c = ROOT.TCanvas("c","c", 800, 800)

    dhist = ROOT.TH1F("dh", "dh", 100, MINMH, MAXMH)
    dhist.GetYaxis().SetTitleOffset(1.5)
    dhist.GetXaxis().SetTitleOffset(1.2)
    dhist.GetYaxis().SetTitleSize(0.04)
    dhist.GetXaxis().SetTitleSize(0.04)
    dhist.GetYaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetLabelSize(0.04)
    dhist.GetXaxis().SetRangeUser(MINMH, MAXMH)
    dhist.GetYaxis().SetRangeUser(MINPV, MAXPV)
    dhist.GetXaxis().SetTitle("m_{t#bar{t}} (GeV)")
    dhist.GetYaxis().SetTitle("Local p-value")
    dhist.Draw("AXIS")

    # ------------------------------------------------------------------------
    # Additional Lines stored in --addline -----------------------------------
    #for lineF in options.addline:

        ## Parse the string, should be file.root:color:linestyle:legend entry    
        #vals = lineF.split(":")
        #ftmp = ROOT.TFile(vals[0])
        #grext = ftmp.Get("observed")
        #grext.SetLineColor(int(vals[1]))
        #grext.SetLineStyle(int(vals[2]))
        #grext.SetLineWidth(2)
        #legend.AddEntry(grext,vals[3],"L")
        #grext.Draw("same")
    # ------------------------------------------------------------------------

    text = ROOT.TLatex()
    text.SetTextColor(ROOT.TColor.GetColor("#53777A"))
    text.SetTextSize(FONTSIZE)
    text.SetTextFont(42)
        
    pvalues = [ROOT.RooStats.SignificanceToPValue(L) for L in Lines]
    TLines = [ROOT.TLine(MINMH, V, MAXMH, V) for V in pvalues]

    MG.Draw("L3")

    for j,TL in enumerate(TLines):
        TL.SetLineStyle(ROOT.kDashed)
        TL.SetLineColor(ROOT.TColor.GetColor("#53777A"))
        TL.SetLineWidth(2)
        TL.Draw("same")
        text.DrawLatex(MAXMH + 4, pvalues[j] * 0.95, "%d#sigma" % Lines[j])

    #c.SetGrid(True)
    c.SetLogy();

    dhist.Draw("AXIS SAME")

    legend=ROOT.TLegend(0.55,0.17+0.03,0.89,0.35+0.03)
    legend.SetFillColor(0)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetTextSize(FONTSIZE)
    if not options.expectedOnly:
        legend.AddEntry(graphObs, "Observed", "L")
    legend.AddEntry(graphMed, "Expected", "L")
    legend.AddEntry(graph68, "#pm 1#sigma expected", "f")
    legend.AddEntry(graph95, "#pm 2#sigma expected", "f")

    legend.Draw()
    c.Update()

    lumi = ROOT.TPaveText(c.GetLeftMargin(), 1 - 0.5 * c.GetTopMargin(), 1 - c.GetRightMargin(), 1, "brNDC")
    lumi.SetFillStyle(0)
    lumi.SetBorderSize(0)
    lumi.SetMargin(0)
    lumi.SetTextSize(0.6 * c.GetTopMargin())
    lumi.SetTextFont(42)
    lumi.SetTextAlign(33)
    lumi.AddText("19.67 fb^{-1} (8 TeV)")
    lumi.Draw()

    title = ROOT.TPaveText(c.GetLeftMargin(), 1 - 0.5 * c.GetTopMargin(), 1 - c.GetRightMargin(), 1, "brNDC")
    title.SetFillStyle(0)
    title.SetBorderSize(0)
    title.SetMargin(0)
    title.SetTextSize(0.75 * c.GetTopMargin())
    title.SetTextFont(62)
    title.SetTextAlign(13)
    title.AddText("CMS #font[52]{#scale[0.76]{Preliminary}}")
    title.Draw()

    c.Print("plot.pdf")

    #mytext= ROOT.TLatex()
    ##mytext.SetTextSize(FONTSIZE)
    #mytext.SetTextFont(42)
    #mytext.SetNDC()

    ##box = ROOT.TPave(0.19,0.17,0.42,0.3,2,"NDC")
    ##box.SetLineColor(1)
    ##box.SetFillColor(0)
    ##box.SetShadowColor(0)
    ##if not options.nobox: box.Draw()
    #mytext.DrawLatex(0.2,0.26,"CMS Preliminary")
    ##for t,lineT in enumerate(options.addtxt):
        ##mytext.DrawLatex(0.2,0.25-(t+1)*0.04,"%s"%(lineT))
    #legend.Draw()
    #ROOT.gPad.RedrawAxis();
    
    #raw_input("Looks Ok?")
    #c.SaveAs("pvaluesplot.pdf")
    #c.SaveAs("pvaluesplot.png")


MakePvalPlot(MG)
