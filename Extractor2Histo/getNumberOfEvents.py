#! /usr/bin/env python

from optparse import OptionParser
parser = OptionParser()

parser.add_option("-i", "--input", dest="input")
parser.add_option("-n", "--gen", dest="gen", type = int)

(options, args) = parser.parse_args()

from ROOT import *
import math

f = TFile.Open(options.input)

h = f.Get("mttSelected_btag_sel_reco_fullsel")
h2 = f.Get("mttSelected_btag_sel_reco_fullsel_no_gen_weight")

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
