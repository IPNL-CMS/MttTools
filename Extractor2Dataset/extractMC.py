#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
        # Background + Signal
        ["Signal_S0_S_i_M500_cpl1_scalar_dataset_nominal.root", "lists/MTT_S0_S_i_M500_cpl1_scalar_%s.list"],
        ["Signal_S0_S_i_M700_cpl1_scalar_dataset_nominal.root", "lists/MTT_S0_S_i_M700_cpl1_scalar_%s.list"],
        ["Signal_S0_S_i_M500_cpl1_pseudoscalar_dataset_nominal.root", "lists/MTT_S0_S_i_M500_cpl1_pseudoscalar_%s.list"],
        ["Signal_S0_S_i_M700_cpl1_pseudoscalar_dataset_nominal.root", "lists/MTT_S0_S_i_M700_cpl1_pseudoscalar_%s.list"],

        # Background
        ["MC_TT_powheg_dataset_nominal.root", "lists/MTT_TT_powheg_%s.list", 21675970, 245.8],

        ["MC_T_tW-channel_dataset_nominal.root", "lists/MTT_T_tW-channel_%s.list", 497658, 11.1],
        ["MC_T_s-channel_dataset_nominal.root", "lists/MTT_T_s-channel_%s.list", 259961, 3.79],
        ["MC_T_t-channel_dataset_nominal.root", "lists/MTT_T_t-channel_%s.list", 3728227, 56.4],

        ["MC_Tbar_tW-channel_dataset_nominal.root", "lists/MTT_Tbar_tW-channel_%s.list", 493460, 11.1],
        ["MC_Tbar_s-channel_dataset_nominal.root", "lists/MTT_Tbar_s-channel_%s.list", 139974, 1.76],
        ["MC_Tbar_t-channel_dataset_nominal.root", "lists/MTT_Tbar_t-channel_%s.list", 1935072, 30.7],

        ["MC_DY1JetsToLL_M-50_dataset_nominal.root", "lists/MTT_DY1JetsToLL_M-50_%s.list", 24045248, 666.3],
        ["MC_DY2JetsToLL_M-50_dataset_nominal.root", "lists/MTT_DY2JetsToLL_M-50_%s.list", 21852156, 215.0],
        ["MC_DY3JetsToLL_M-50_dataset_nominal.root", "lists/MTT_DY3JetsToLL_M-50_%s.list", 11015445, 60.7],
        ["MC_DY4JetsToLL_M-50_dataset_nominal.root", "lists/MTT_DY4JetsToLL_M-50_%s.list", 6402827, 27.3],

        ["MC_W1JetsToLNu_dataset_nominal.root", "lists/MTT_W1JetsToLNu_%s.list", 23141598, 6662.8],
        ["MC_W2JetsToLNu_dataset_nominal.root", "lists/MTT_W2JetsToLNu_%s.list", 34044921, 2159.2],
        ["MC_W3JetsToLNu_dataset_nominal.root", "lists/MTT_W3JetsToLNu_%s.list", 15539503, 640.4],
        ["MC_W4JetsToLNu_dataset_nominal.root", "lists/MTT_W4JetsToLNu_%s.list", 13382803, 264.0],
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

# Build output tree structure
for type in ["semie", "semimu"]:
    path = "datasets/%s/%s" % (d, type)
    try:
        os.makedirs(path)
    except:
        pass

print("Extracting dataset ...")

for file in files:
    events = file[2] if len(file) >= 3 else -1
    xsection = file[3] if len(file) >= 4 else -1
    weight = 19667. * xsection / events if events > 0 else 1;
    print("Weight for %s: %.15f" % (file[1], weight))
    for type in ["semie", "semimu"]:
        path = "datasets/%s/%s" % (d, type)
        tmpfile.write(launch(file[1] % type, os.path.join(path, file[0]), weight) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 
subprocess.call(args)

# Merge files
for file in files:
    path = "datasets/%s/" % (d)
    merged_file = os.path.join(path, file[0].replace(".root", "_merged.root"))
    args = ["hadd", "-f", merged_file]
    for type in ["semie", "semimu"]:
        args.append(os.path.join(path, type, file[0]))

    subprocess.call(args)

# Merge background files
path = "datasets/%s/" % (d)
merged_file = os.path.join(path, "MC_background.root")
args = ["hadd", "-f", merged_file]
for file in files:
    if not "Signal" in file[0]:
        merged_file = os.path.join(path, file[0].replace(".root", "_merged.root"))
        args.append(merged_file)

    subprocess.call(args)
