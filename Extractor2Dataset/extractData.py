#!/usr/bin/env python

import os, subprocess

inputs = {
    "MTT_Muons_2012.root": "input_data_semimu.list",
    "MTT_Electrons_2012.root": "input_data_semie.list"
    }

jobs = []

def launch(input, output):
  args = ["./extractor2Dataset", "--input-list", input, "-o", output, "--data"]
  p = subprocess.Popen(args)
  jobs.append(p)

print("Extracting dataset ...")
for output, input in inputs.items():
  launch(input, output)

for job in jobs:
  job.wait()

# All is done, merge

print("Merging ...")
args = ["hadd", "MTT_Data_merged_2012.root"]
for output in inputs.keys():
  args.append(output)

subprocess.call(args)

for output in inputs.keys():
  os.remove(output)
