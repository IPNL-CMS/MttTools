#! /bin/env python

from __future__ import division
import argparse, os, array
from ROOT import TChain, TH1F, TFile

parser = argparse.ArgumentParser(description='Submit jobs on the grid.')
parser.add_argument('dataset', nargs=1)
args = parser.parse_args()

dataset = args.dataset[0]
inputList = "files_%s.list" % dataset

if not os.path.exists(inputList):
  print("Error: input list %s not found." % inputList)
  exit(1)

files = [line.strip() for line in open(inputList)]

chain = TChain("event", "event")
for file in files:
  chain.AddFile(file)

#chain.SetBranchStatus("*", 0)
#chain.SetBranchStatus("nTrueInteractions", 1);
branch = chain.GetBranch("nTrueInteractions")

pu = TH1F("pileup", "MC Pileup truth", 70, 0, 70);

entries = chain.GetEntries()
for i in xrange(entries):
  branch.GetEntry(i)
  
  if i % 1000000 == 0:
     print("Iteration %d over %d; %f %%" % (i + 1, entries, (i + 1) / entries * 100))

  pu.Fill(chain.nTrueInteractions, 1)

scale = 1 / pu.Integral();
pu.Scale(scale);

filename = "summer12_computed_mc_%s_pu_truth_70bins.root" % dataset

output = TFile.Open(filename, "recreate")
pu.Write()
output.Close()
