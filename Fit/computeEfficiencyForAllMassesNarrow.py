#!/bin/env python

import subprocess
import json, uuid, sys, os, datetime

if not os.path.exists("analysis.json"):
  print("No analysis currently defined. Please use startAnalysis first")
  sys.exit(1)

json_data = open("analysis.json")

try:
  data = json.load(json_data)
except:
  print("No analysis currently defined. Please use startAnalysis first")
  sys.exit(1)

json_data.close()

index = data["current_analysis"]
type = data["analysis"][index][data["analysis"][index].keys()[0]]["type"]

btag = [1, 2]

if type == "zprime":
    masses = {
            500: "5",
            750: "7p5",
            1000: "10",
            1250: "12p5",
            1500: "15",
            2000: "20"
            }
    for mass, width in masses.items():
        for b in btag:
            args = ["./computeEfficiencies", "-m", str(mass), "--b-tag", str(b), "--gen-file", "data/zprime_narrow_gen.json", "--trigger-eff", "data/zprime_trigger_efficiencies.json"]
            subprocess.call(args)

elif type == "higgs":
    masses = [400,500,600,700,800]
    for mass in masses:
        for b in btag:
            args = ["./computeEfficiencies", "-m", str(mass), "--b-tag", str(b), "--gen-file", "data/s0_pseudoscalar_gen.json"]
            subprocess.call(args)




