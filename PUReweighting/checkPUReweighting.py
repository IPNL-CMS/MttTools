#! /bin/env python

from __future__ import division
import argparse, os, array

from ROOT import TChain, TH1F, TFile, TCanvas

from ROOT import gROOT
gDirectory = gROOT.GetGlobal("gDirectory")

parser = argparse.ArgumentParser(description='Check PU reweighting.')
parser.add_argument('data_pu_profile', nargs=1)
parser.add_argument('mc_input_file', nargs=1)
parser.add_argument('-b', action = "store_true")
args = parser.parse_args()

data_pu_profile_name = args.data_pu_profile[0]
mc_input_list = args.mc_input_list[0]

if not os.path.exists(mc_input_list):
  print("Error: input list %s not found." % mc_input_list)
  exit(1)

files = [line.strip() for line in open(mc_input_list)]

chain = TChain("event", "event")
for file in files:
  print("Adding file '%s'" % file)
  chain.Add(file)

chain.SetBranchStatus("*", 0)
chain.SetBranchStatus("nTrueInteractions", 1);

chain.Draw("nTrueInteractions>>pileup(70, 0, 70)")

pu = gDirectory.Get("pileup")

scale = 1 / pu.Integral()
pu.Scale(scale)

data_profile_file = TFile(data_pu_profile_name)
data_pu_profile = data_profile_file.Get("pileup")

scale = 1 / data_pu_profile.Integral()
data_pu_profile.Scale(scale)

data_pu_profile.SetMarkerStyle(20)
data_pu_profile.SetMarkerSize(0.8)
data_pu_profile.Draw("P")

pu.Draw("same")
