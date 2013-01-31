#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
    #["MTT_Signal_Zprime_500_Narrow_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_500_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_750_Narrow_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_750_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_1000_Narrow_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_1000_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_1250_Narrow_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_1250_Narrow_full_stat_%s.list"],
    #["MTT_Signal_Zprime_1500_Narrow_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_1500_Narrow_full_stat_%s.list"]

    ["MTT_Signal_Zprime_500_Large_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_500_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_750_Large_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_750_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_1000_Large_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_1000_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_1250_Large_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_1250_Large_full_stat_%s.list"],
    ["MTT_Signal_Zprime_1500_Large_2012_dataset_%s_nominal.root", "MC/MTT_Zprime_1500_Large_full_stat_%s.list"]

    #["MTT_MC_TTJets_powheg_18Jan13_dataset_%s_nominal.root", "MC/MTT_TT_powheg_18Jan13_%s.list", 20729745, 234.],
    #["MTT_MC_TTJets_MassiveBinDECAY_19Dec12_dataset_%s_nominal.root", "MC/MTT_TTJets_MassiveBinDECAY_19Dec12_%s.list", 6908904., 234.],
    #["MTT_MC_DYJetsToLL_M-50_19Dec12_dataset_%s_nominal.root", "MC/MTT_DYJetsToLL_M-50_19Dec12_%s.list", 30439508, 3503.71],
    #["MTT_MC_WJetsToLNu_19Dec12_dataset_%s_nominal.root", "MC/MTT_WJetsToLNu_19Dec12_%s.list", 57687669, 37509.0],

    #["MTT_MC_Tbar_s-channel_19Dec12_dataset_%s_nominal.root", "MC/MTT_Tbar_s-channel_19Dec12_%s.list", 139803, 1.75776],
    #["MTT_MC_Tbar_t-channel_19Dec12_dataset_%s_nominal.root", "MC/MTT_Tbar_t-channel_19Dec12_%s.list", 1932762, 30.7],
    #["MTT_MC_Tbar_tW-channel_19Dec12_dataset_%s_nominal.root", "MC/MTT_Tbar_tW-channel_19Dec12_%s.list", 492537, 11.1773],
    #["MTT_MC_T_s-channel_19Dec12_dataset_%s_nominal.root", "MC/MTT_T_s-channel_19Dec12_%s.list", 259571, 3.79],
    #["MTT_MC_T_t-channel_19Dec12_dataset_%s_nominal.root", "MC/MTT_T_t-channel_19Dec12_%s.list", 3722922, 56.4],
    #["MTT_MC_T_tW-channel_19Dec12_dataset_%s_nominal.root", "MC/MTT_T_tW-channel_19Dec12_%s.list", 496669, 11.1773],
  ]

def launch(input, output, weight):
  args = ["./extractor2Dataset", "--input-list", input, "-o", output, "--mc"]
  if "semie" in input:
    args.append("--type semie")
  elif "semimu" in input:
    args.append("--type semimu")

  args.append("--weight %.15f" % weight)

  return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

print("Extracting dataset ...")
for file in files:
  events = file[2] if len(file) >= 3 else -1
  xsection = file[3] if len(file) >= 4 else -1
  weight = 19580. * xsection / events if events > 0 else 1;
  print("Weight for %s: %.15f" % (file[1], weight))
  for type in ["semie", "semimu"]:
    tmpfile.write(launch(file[1] % type, file[0] % type, weight) + "\n");

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
