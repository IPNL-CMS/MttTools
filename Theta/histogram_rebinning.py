#! /usr/bin/env python

import sys
sys.argv.append('-b')

from ROOT import TFile, TCanvas

import math
import array

def findMaximum(histogram):
    index = 1
    maximum = 0
    for i in range(1,histogram.GetNbinsX()+1):
        if histogram.GetBinContent(i) > maximum:
            maximum = histogram.GetBinContent(i)
            index = i
    return index, maximum


def findLowIndex(histogram, rerror):
    value = 0.0
    error = 0.0
    index, maximum = findMaximum(histogram)
    for i in range(1,histogram.GetNbinsX()+1):
        value = value + histogram.GetBinContent(i)
        error = math.sqrt(error**2+histogram.GetBinError(i)**2)
        if value > maximum and i > 1:
            return i-1
        elif value > maximum and i < 1:
            return 1
        ratio = 1.0
        if value != 0: ratio = error/value
        if ratio < rerror:
            return i


def findBinSize(histogram, highindexes, rerror, minvalue, maxbinsize, start, stop):
    #print "minvalue %f, maxbinsize %d, start %d, stop %d : " % (minvalue, maxbinsize, start, stop), highindexes
    value = 0.0
    error = 0.0
    for i in range(start, stop, -1):
        value = value + math.fabs(histogram.GetBinContent(i))
        error = math.sqrt(error**2+histogram.GetBinError(i)**2)
        ratio = 1.0
        binsize = i - 1
        if value != 0: ratio = error/math.fabs(value)
        #print "value %f ; error %f ; ratio : %f" % (value, error, ratio)
        # if ratio < rerror and value*(1+ratio) >= minvalue:
        if ratio < rerror:
            if binsize <= maxbinsize:
                highindexes.append(i)
                if not findBinSize(histogram, highindexes, rerror, value, binsize, i-1, stop):
                    highindexes.pop()
                    continue
                return True
            else:
                return False
    # highindexes.append(stop+1)
    return True        


def computeBinning(histogram, rerror):
    highindexes = []
    lowindex = findLowIndex(histogram, rerror)
    maxindex, maximum = findMaximum(histogram)
    findBinSize(histogram, highindexes, rerror, 0, histogram.GetNbinsX(), histogram.GetNbinsX(), lowindex)     
    highindexes = sorted(highindexes)
    binning = [histogram.GetBinLowEdge(1), histogram.GetBinLowEdge(lowindex)+histogram.GetBinWidth(lowindex)]
    print highindexes
    for i in highindexes[1:]:
        binning.append(histogram.GetBinLowEdge(i))
    binning.append(histogram.GetBinLowEdge(histogram.GetNbinsX())+histogram.GetBinWidth(histogram.GetNbinsX()))
    print binning
    return binning


f = TFile.Open("theta_histos_merged_BW5.root")
f_ = f.Get("mtt_mu_2btag__DATA");

binning = array.array('d', computeBinning(f_, 0.3))

c = TCanvas("t", "t", 800, 800)

#f_.DrawCopy()

f_ = f_.Rebin(len(binning) - 1, "test", binning)
f_.SetLineColor(1)

f_.Draw()
c.Print("test.pdf")
