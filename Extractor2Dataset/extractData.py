#!/usr/bin/env python

import os, subprocess

inputs = [
    ["MTT_Muons_2012.root", "input_data_semimu.list", "semimu"],
    ["MTT_Electrons_2012.root", "input_data_semie.list", "semie"]
    ]

jobs = []

def launch(input, output, type):
  args = ["./extractor2Dataset", "--input-list", input, "-o", output, "--data", "--type", type]
  p = subprocess.Popen(args)
  jobs.append(p)

print("Extracting dataset ...")
for input in inputs:
  launch(input[1], input[0], input[2])

for job in jobs:
  job.wait()

# All is done, merge

print("Merging ...")
args = ["hadd", "MTT_Data_merged_2012.root"]
for output in inputs:
  args.append(output[0])

subprocess.call(args)

for output in inputs.keys():
  os.remove(output)
