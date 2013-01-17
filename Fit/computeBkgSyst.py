#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os, math, subprocess, shutil, sys, tempfile, argparse, json

from ROOT import TFile, RooWorkspace, RooRealVar, RooArgSet, TH1F

def runFitMtt(name, value):
  tmpjson = tempfile.NamedTemporaryFile(dir = '/tmp/', delete = True, suffix = '.json')
  cmds = ["./fitMtt", "-m", str(args.mass), "--b-tag", str(args.btag), "-i", args.input, "--bkg-parameter-name", name, "--bkg-parameter-value", "%.15f" % value, "--temp-output-file", tmpjson.name]
  with open(os.devnull, "w") as fnull:
    subprocess.call(cmds, stdout = fnull, stderr = fnull)

  results = json.load(tmpjson)
  sigma = results["result"]

  tmpjson.close()

  return math.fabs(sigma)
  
def getSigmaRef():
  jsonFile = open("%s/sigma_reference.json" % base_path)
  jsonValues = json.load(jsonFile)
  jsonFile.close()
  
  return jsonValues[analysisUUID][str(args.mass)][str(args.btag)]["sigma"]

def scanParameter(name, value, error):
  numPoints = 4

  #print("Value: %f, Error: %f" % (value, error))

  hist = TH1F("hist", "histo", numPoints + 1, value - 1.5 * error, value + 1.5 * error)

  sigma_ref = math.fabs(getSigmaRef())
 
  for i in range(1, int(((numPoints / 2)) + 1)):
    new_value = value - error / i
    sigma = runFitMtt(name, new_value)
    hist.Fill(new_value, math.fabs((sigma_ref - sigma)))
    #print("Value: %f; sigma: %f, syst: %f" % (new_value, sigma, math.fabs(sigma_ref - sigma)))

  for i in range(int(numPoints / 2), 0, -1):
    new_value = value + error / i
    sigma = runFitMtt(name, new_value)
    hist.Fill(new_value, math.fabs(sigma_ref - sigma))
    #print("Value: %f; sigma: %f, syst: %f" % (new_value, sigma, math.fabs(sigma_ref - sigma)))

  return hist.Integral() / (numPoints + 1);

  #f = TFile.Open("test.root", "recreate")
  #hist.Write();
  #f.Close()

if sys.version_info<(2,7,0):
  sys.stderr.write("You need python 2.7 or later to run this script\n")
  sys.exit(1)

parser = argparse.ArgumentParser(description='Compute background systematic.')
parser.add_argument("--b-tag", dest="btag", required=True, type=int)
parser.add_argument("-m", dest="mass", required=True, type=int)
parser.add_argument("-i", dest="input", required=True, type=str)
#parser.add_argument("-o", dest="output", required=True, type=str)
args = parser.parse_args()

f = open("analysis.json")
params = json.load(f)
f.close()

current_analysis = params["current_analysis"]
analysisUUID = params["analysis"][current_analysis].keys()[0]

base_path = "analysis/%s" % analysisUUID

sigma_ref = getSigmaRef()

f = TFile.Open("%s/background_parameters_%d_%d_btag.root" % (base_path, args.mass, args.btag))
w = f.Get("w")

vars = w.allVars()

iterator = vars.createIterator()

var = iterator.Next()

sigmas = []
systs = []

while var != None:
  name = var.GetName()
  value = var.getVal()
  error = var.getError()

  syst = scanParameter(name, value, error)
  print("For %s: %f" % (name, syst))
  systs.append(syst)
  var = iterator.Next()
  continue

  print("")

  sigma_up = runFitMtt(name, value + error)
  sigma_down = runFitMtt(name, value - error)

  syst_up = math.fabs(math.fabs(sigma_up) - math.fabs(sigma_ref))
  syst_down = math.fabs(math.fabs(sigma_down) - math.fabs(sigma_ref))

  print("For %s: %f (up: %f %f, down: %f %f)" % (name, (syst_up + syst_down) / 2., syst_up, sigma_up, syst_down, sigma_down))

  sigmas.append((syst_up + syst_down) / 2.)

print("")

syst = 0
for s in systs:
  syst = syst + (s * s)

print("Quadrature sum: %f" % (math.sqrt(syst)))

syst = 0
for s in systs:
  syst = syst + s

print("Mean: %f" % (syst / len(systs)))
#sigma = 0
#for s in sigmas:
  #sigma = sigma + s

#print("Systematic error: %f" % (sigma / len(sigmas)))

## Strategy: generate a toy dataset from data, and then compute JEC syst on it. Do that 100 times
#for i in range(0, 10):
  #print("Computing syst for toy #%d" % (i + 1))
  #tmpfile = tempfile.NamedTemporaryFile(dir = '/tmp/', delete = False, suffix = '.root')
  #tmpfile.close()

  #print("Generating dataset ...")
  #cmd = ["./generateToyDataset", "-i", args.input, "-o", tmpfile.name]
  #with open(os.devnull, "w") as fnull: 
    #if subprocess.call(cmd, stdout = fnull, stderr = fnull) != 0:
      #raise SystemExit("Cannot generate toy dataset")

  #tmpjson = tempfile.NamedTemporaryFile(dir = '/tmp/', delete = True, suffix = '.json')

  #print("Fitting ...")
  #cmd = ["./computeJECSyst", "-m", str(args.mass), "--b-tag", str(args.btag), "-i", tmpfile.name, "-o", tmpjson.name]
  #with open(os.devnull, "w") as fnull: 
    #if subprocess.call(cmd, stdout = fnull, stderr = fnull) != 0:
      #raise SystemExit("Cannot fit generated toy dataset")

  #results = json.load(tmpjson)
  #syst = results["result"]
  #print("Done. Syst: %f pb" % syst)

  #systs.Fill(syst)
  
  #tmpjson.close()
  #os.unlink(tmpfile.name)

##g = TF1("g", "gaus")
##g.SetParLimits(1, 0, 100) # Be sure out mean is always positive

##systs.Fit(g, "Q")
##print("Mean: %f" % g.GetParameter(1))

#file = TFile.Open(args.output, "recreate")
#systs.Write()
#file.Close()
