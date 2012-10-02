#! /bin/env python

from __future__ import division
import argparse, os, array

# Batch mode
import sys
sys.argv.append('-b')

from ROOT import TChain, TH1F, TFile

from ROOT import gROOT
gDirectory = gROOT.GetGlobal("gDirectory")

parser = argparse.ArgumentParser(description='Submit jobs on the grid.')
parser.add_argument('dataset', nargs=1)
parser.add_argument('-b', action = "store_true")
args = parser.parse_args()

dataset = args.dataset[0]
inputList = "files_%s.list" % dataset

if not os.path.exists(inputList):
  print("Error: input list %s not found." % inputList)
  exit(1)

files = [line.strip() for line in open(inputList)]

chain = TChain("event", "event")
for file in files:
  print("Adding file '%s'" % file)
  chain.Add(file)

chain.SetBranchStatus("*", 0)
chain.SetBranchStatus("nTrueInteractions", 1);

chain.Draw("nTrueInteractions>>pileup(70, 0, 70)")

pu = gDirectory.Get("pileup");

scale = 1 / pu.Integral();
pu.Scale(scale);

filename = "summer12_computed_mc_%s_pu_truth_70bins.root" % dataset

output = TFile.Open(filename, "recreate")
pu.Write()
output.Close()
