#!/usr/bin/env python

import os, subprocess, tempfile

files = [
    ["MTT_Zprime_750_Narrow_2012_dataset_%s.root", "MC/MTT_Zprime_750_Narrow_full_stat_%s.list"],
    ["MTT_Zprime_1000_Narrow_2012_dataset_%s.root", "MC/MTT_Zprime_1000_Narrow_full_stat_%s.list"],
    ["MTT_Zprime_1250_Narrow_2012_dataset_%s.root", "MC/MTT_Zprime_1250_Narrow_full_stat_%s.list"],
    ["MTT_Zprime_1500_Narrow_2012_dataset_%s.root", "MC/MTT_Zprime_1500_Narrow_full_stat_%s.list"]
  ]

def launch(input, output):
  args = ["./extractor2Dataset", "--input-list", input, "-o", output, "--mc"]
  if "semie" in input:
    args.append("--type semie")
  elif "semimu" in input:
    args.append("--type semimu")

  return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

print("Extracting dataset ...")
for file in files:
  for type in ["semie", "semimu"]:
    tmpfile.write(launch(file[1] % type, file[0] % type) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 
subprocess.call(args)

# Merge files
for file in files:
  merged_file = file[0].replace("_%s", "_merged")
  args = ["hadd", merged_file, file[0] % "semimu", file[0] % "semie"]
  subprocess.call(args)
