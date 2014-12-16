#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()

parser.add_option("-i", "--input", dest="input")
parser.add_option("", "--sigma", dest="sigma", action="store_true")

(options, args) = parser.parse_args()

import scipy, scipy.stats
def Z_to_p(z_value):
    """
    Convert a Z-value to a p-value; inverse of :meth:`p_to_Z`
    """
    return scipy.stats.norm.sf(z_value)

def Z_to(z_value):
    return z_value if options.sigma else Z_to_p(z_value)

# Standard Imports and calculators
import json, re
import ROOT
import array, sys, numpy
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine(".L tdrstyle.cc")
ROOT.gROOT.ProcessLine("setTDRStyle();")

with open(options.input) as json_data:
    pvalues_all = json.load(json_data)

options.expectedOnly = True

index = 0
for mass in sorted(pvalues_all):
    regex = r'(\d+)'
    m = float(re.search(regex, mass).group(0))

    pvalues = pvalues_all[mass]

    sys.stdout.write("\sz (m = %d\,\GeV) & " % (m))

    sys.stdout.write(" %.4f $\pm$ %.4f & " % (float(Z_to(pvalues["median"][0])), float(Z_to(pvalues["median"][1]))))
    median = Z_to(pvalues["median"][0])
    median_up = Z_to(pvalues["median_up"][0])
    median_down = Z_to(pvalues["median_low"][0])

    median_up_2sigma = Z_to(pvalues["median_up_2sigma"][0])
    median_down_2sigma = Z_to(pvalues["median_low_2sigma"][0])

    # Significance and p-values are ordered in the opposite way
    error_up = median_down - median
    error_low = median - median_up

    error_up_95 = median_down_2sigma - median
    error_low_95 = median - median_up_2sigma

    sys.stdout.write(" %.4f -- %.4f & " % (median_up, median_down))
    sys.stdout.write(" %.4f -- %.4f & " % (median_up_2sigma, median_down_2sigma))

    sys.stdout.write(" %.4f \\\\" % (float(Z_to(pvalues["observed"][0]))))

    print("")

    index += 1
