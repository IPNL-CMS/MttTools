#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os, math, subprocess, shutil, sys, tempfile, argparse, json

from ROOT import TH1F, TCanvas, TF1, TFile

if sys.version_info<(2,7,0):
  sys.stderr.write("You need python 2.7 or later to run this script\n")
  sys.exit(1)

parser = argparse.ArgumentParser(description='Compute JEC systematic using toys.')
parser.add_argument("--b-tag", dest="btag", required=True, type=int)
parser.add_argument("-m", dest="mass", required=True, type=int)
parser.add_argument("-i", dest="input", required=True, type=str)
parser.add_argument("-o", dest="output", required=True, type=str)
args = parser.parse_args()

systs = TH1F("systs", "systs", 200, 0, 5)

# Strategy: generate a toy dataset from data, and then compute JEC syst on it. Do that 100 times
for i in range(0, 10):
  print("Computing syst for toy #%d" % (i + 1))
  tmpfile = tempfile.NamedTemporaryFile(dir = '/tmp/', delete = False, suffix = '.root')
  tmpfile.close()

  print("Generating dataset ...")
  cmd = ["./generateToyDataset", "-i", args.input, "-o", tmpfile.name]
  with open(os.devnull, "w") as fnull: 
    if subprocess.call(cmd, stdout = fnull, stderr = fnull) != 0:
      raise SystemExit("Cannot generate toy dataset")

  tmpjson = tempfile.NamedTemporaryFile(dir = '/tmp/', delete = True, suffix = '.json')

  print("Fitting ...")
  cmd = ["./computeJECSystWithToys", "-m", str(args.mass), "--b-tag", str(args.btag), "-i", tmpfile.name, "-o", tmpjson.name]
  with open(os.devnull, "w") as fnull: 
    if subprocess.call(cmd, stdout = fnull, stderr = fnull) != 0:
      raise SystemExit("Cannot fit generated toy dataset")

  results = json.load(tmpjson)
  syst = results["result"]
  print("Done. Syst: %f pb" % syst)

  systs.Fill(syst)
  
  tmpjson.close()
  os.unlink(tmpfile.name)

#g = TF1("g", "gaus")
#g.SetParLimits(1, 0, 100) # Be sure out mean is always positive

#systs.Fit(g, "Q")
#print("Mean: %f" % g.GetParameter(1))

file = TFile.Open(args.output, "recreate")
systs.Write()
file.Close()
