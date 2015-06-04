#! /usr/bin/env python

from optparse import OptionParser
parser = OptionParser()

parser.add_option("-i", "--input", dest="input")
parser.add_option("-n", "--gen", dest="gen", type = int)
parser.add_option("", "--scale", action="store_true", dest="scale", default=False)

(options, args) = parser.parse_args()

from ROOT import *
import math

f = TFile.Open(options.input)

h = f.Get("mttSelected_btag_sel_reco_fullsel")
h2 = f.Get("mttSelected_btag_sel_reco_fullsel_no_gen_weight")

lumi = 1
xsection = 1
Ngen = 1

if options.scale:
    lumi = 19667
    if "scalar" in options.input:
        if "400" in options.input:
            xsection = 0.5289 * 5.340
            Ngen = 2074732
        elif "600" in options.input:
            xsection = 0.1891 * 4.886
            Ngen = 2354567
        elif "800" in options.input:
            xsection = 0.09436 * 2.888
            Ngen = 1999938
    elif "pseudoscalar" in options.input:
        if "400" in options.input:
            xsection = 1.169 * 5.039
            Ngen = 2276384
        elif "600" in options.input:
            xsection = 0.4414 * 2.709
            Ngen = 1999949
        elif "800" in options.input:
            xsection = 0.2325 * 1.697
            Ngen = 1999940


h.Scale(xsection * lumi / (Ngen))
entries_pos = 0
entries_neg = 0
for i in range(1, h.GetNbinsX() + 1):
    c = h.GetBinContent(i)
    if c > 0:
        entries_pos += c
    else:
        entries_neg += c

n = abs(entries_pos) + abs(entries_neg)

eff = n / options.gen
err = math.sqrt( (eff * (1 - eff)) / options.gen )

print eff * 100
print err * 100
print "\\num{%.3f \pm %.3f}" % (eff * 100, err * 100)

eff = float(h2.Integral()) / options.gen
err = math.sqrt( (eff * (1 - eff)) / options.gen )
print "\\num{%.3f \pm %.3f}" % (eff * 100, err * 100)
