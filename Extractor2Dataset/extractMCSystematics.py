#!/usr/bin/env python

import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
    #["MTT_Signal_Zprime_500_Narrow_2012_dataset_%s_JECup.root",  "Systematics/JECup/MTT_Zprime_500_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_750_Narrow_2012_dataset_%s_JECup.root",  "Systematics/JECup/MTT_Zprime_750_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_1000_Narrow_2012_dataset_%s_JECup.root", "Systematics/JECup/MTT_Zprime_1000_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_1250_Narrow_2012_dataset_%s_JECup.root", "Systematics/JECup/MTT_Zprime_1250_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_1500_Narrow_2012_dataset_%s_JECup.root", "Systematics/JECup/MTT_Zprime_1500_Narrow_full_stat_%s.list"],

    #["MTT_Signal_Zprime_500_Narrow_2012_dataset_%s_JECdown.root",  "Systematics/JECdown/MTT_Zprime_500_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_750_Narrow_2012_dataset_%s_JECdown.root",  "Systematics/JECdown/MTT_Zprime_750_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_1000_Narrow_2012_dataset_%s_JECdown.root", "Systematics/JECdown/MTT_Zprime_1000_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_1250_Narrow_2012_dataset_%s_JECdown.root", "Systematics/JECdown/MTT_Zprime_1250_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_1500_Narrow_2012_dataset_%s_JECdown.root", "Systematics/JECdown/MTT_Zprime_1500_Narrow_full_stat_%s.list"]

    ["MTT_Signal_Zprime_500_Large_2012_dataset_%s_JECup.root",  "Systematics/JECup/MTT_Zprime_500_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_750_Large_2012_dataset_%s_JECup.root",  "Systematics/JECup/MTT_Zprime_750_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_1000_Large_2012_dataset_%s_JECup.root", "Systematics/JECup/MTT_Zprime_1000_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_1250_Large_2012_dataset_%s_JECup.root", "Systematics/JECup/MTT_Zprime_1250_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_1500_Large_2012_dataset_%s_JECup.root", "Systematics/JECup/MTT_Zprime_1500_Large_full_stat_%s.list"],

    ["MTT_Signal_Zprime_500_Large_2012_dataset_%s_JECdown.root",  "Systematics/JECdown/MTT_Zprime_500_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_750_Large_2012_dataset_%s_JECdown.root",  "Systematics/JECdown/MTT_Zprime_750_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_1000_Large_2012_dataset_%s_JECdown.root", "Systematics/JECdown/MTT_Zprime_1000_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_1250_Large_2012_dataset_%s_JECdown.root", "Systematics/JECdown/MTT_Zprime_1250_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_1500_Large_2012_dataset_%s_JECdown.root", "Systematics/JECdown/MTT_Zprime_1500_Large_full_stat_%s.list"]
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
  merged_file = file[0].replace("_%s", "_merged_%s" % d)
  args = ["hadd", merged_file, file[0] % "semimu", file[0] % "semie"]
  subprocess.call(args)

for file in files:
  os.remove(file[0] % "semimu")
  os.remove(file[0] % "semie")
