#!/usr/bin/env python

import os, subprocess, tempfile

files = [
    "MTT_Zprime_500_Narrow_2012_08Nov_%s.root",
    "MTT_Zprime_750_Narrow_2012_08Nov_%s.root",
    "MTT_Zprime_1000_Narrow_2012_08Nov_%s.root",
    "MTT_Zprime_1250_Narrow_2012_08Nov_%s.root",
    "MTT_Zprime_1500_Narrow_2012_08Nov_%s.root",
    "MTT_Zprime_2000_Narrow_2012_08Nov_%s.root"

#    "MTT_DYJetsToLL_M-50_2012_v1_%s.root",
#    "MTT_QCD_Pt_20_30_EMEnriched_2012_v1_semie.root",
#    "MTT_QCD_Pt_20_MuEnriched_2012_v1_semimu.root",
#    "MTT_QCD_Pt_30_80_EMEnriched_2012_v1_semie.root",
#    "MTT_QCD_Pt_80_170_EMEnriched_2012_v1_semie.root",
#    "MTT_Tbar_s-channel_2012_v1_%s.root",
#    "MTT_Tbar_t-channel_2012_v1_%s.root",
#    "MTT_Tbar_tW-channel_2012_v1_%s.root",
#    "MTT_TTJets_2012_v1_%s.root",
#    "MTT_T_tW-channel_2012_v1_%s.root",
#    "MTT_WJetsToLNu_2012_v1_%s.root",
  ]

def launch(input, output):
  args = ["./extractor2Dataset", "-i", "input/" + input, "-o", output, "--mc"]
  if "semie" in input:
    args.append("--type semie")
  elif "semimu" in input:
    args.append("--type semimu")

  return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

print("Extracting dataset ...")
for file in files:
  tmpfile.write(launch(file, file) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 
subprocess.call(args)

# Merge files
for file in files:
  merged_file = file.replace("_%s", "_merged")
  args = ["hadd", merged_file, file % "semimu", file % "semie"]
  subprocess.call(args)
