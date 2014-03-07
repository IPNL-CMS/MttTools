#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
        # Background
        ["MC_TT_powheg_dataset_nominal.root", "skims/%s/MC_TT_powheg_skims_nominal.root", 21675970, 245.8],

        ["MC_T_tW-channel_dataset_nominal.root", "skims/%s/MC_T_tW-channel_skims_nominal.root", 497658, 11.1],
        ["MC_T_s-channel_dataset_nominal.root", "skims/%s/MC_T_s-channel_skims_nominal.root", 259961, 3.79],
        ["MC_T_t-channel_dataset_nominal.root", "skims/%s/MC_T_t-channel_skims_nominal.root", 3728227, 56.4],

        ["MC_Tbar_tW-channel_dataset_nominal.root", "skims/%s/MC_Tbar_tW-channel_skims_nominal.root", 493460, 11.1],
        ["MC_Tbar_s-channel_dataset_nominal.root", "skims/%s/MC_Tbar_s-channel_skims_nominal.root", 139974, 1.76],
        ["MC_Tbar_t-channel_dataset_nominal.root", "skims/%s/MC_Tbar_t-channel_skims_nominal.root", 1935072, 30.7],

        ["MC_DY1JetsToLL_M-50_dataset_nominal.root", "skims/%s/MC_DY1JetsToLL_M-50_skims_nominal.root", 24045248, 666.3],
        ["MC_DY2JetsToLL_M-50_dataset_nominal.root", "skims/%s/MC_DY2JetsToLL_M-50_skims_nominal.root", 21852156, 215.0],
        ["MC_DY3JetsToLL_M-50_dataset_nominal.root", "skims/%s/MC_DY3JetsToLL_M-50_skims_nominal.root", 11015445, 60.7],
        ["MC_DY4JetsToLL_M-50_dataset_nominal.root", "skims/%s/MC_DY4JetsToLL_M-50_skims_nominal.root", 6402827, 27.3],

        ["MC_W1JetsToLNu_dataset_nominal.root", "skims/%s/MC_W1JetsToLNu_skims_nominal.root", 23141598, 6662.8],
        ["MC_W2JetsToLNu_dataset_nominal.root", "skims/%s/MC_W2JetsToLNu_skims_nominal.root", 34044921, 2159.2],
        ["MC_W3JetsToLNu_dataset_nominal.root", "skims/%s/MC_W3JetsToLNu_skims_nominal.root", 15539503, 640.4],
        ["MC_W4JetsToLNu_dataset_nominal.root", "skims/%s/MC_W4JetsToLNu_skims_nominal.root", 13382803, 264.0],

        # QCD
        ["MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_dataset_nominal.root", "skims/semimu/MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_skims_nominal.root"],
        ["MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_dataset_nominal.root", "skims/semimu/MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_skims_nominal.root"],
        ["MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_dataset_nominal.root", "skims/semimu/MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_skims_nominal.root"],
        ["MC_QCD_pt150_bEnriched_MuEnrichedPt14_dataset_nominal.root", "skims/semimu/MC_QCD_pt150_bEnriched_MuEnrichedPt14_skims_nominal.root"],

        ["MC_QCD_Pt_20_30_EMEnriched_dataset_nominal.root", "skims/semie/MC_QCD_Pt_20_30_EMEnriched_skims_nominal.root"],
        ["MC_QCD_Pt_30_80_EMEnriched_dataset_nominal.root", "skims/semie/MC_QCD_Pt_30_80_EMEnriched_skims_nominal.root"],
        ["MC_QCD_Pt_80_170_EMEnriched_dataset_nominal.root", "skims/semie/MC_QCD_Pt_80_170_EMEnriched_skims_nominal.root"],
        ["MC_QCD_Pt_170_250_EMEnriched_dataset_nominal.root", "skims/semie/MC_QCD_Pt_170_250_EMEnriched_skims_nominal.root"],
        ["MC_QCD_Pt_250_350_EMEnriched_dataset_nominal.root", "skims/semie/MC_QCD_Pt_250_350_EMEnriched_skims_nominal.root"],
        ["MC_QCD_Pt_350_EMEnriched_dataset_nominal.root", "skims/semie/MC_QCD_Pt_350_EMEnriched_skims_nominal.root"],

        ## Background + Signal
        ["Signal_S0_S_i_M400_cpl1_scalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M400_cpl1_scalar_skims_nominal.root"],
        ["Signal_S0_S_i_M500_cpl1_scalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M500_cpl1_scalar_skims_nominal.root"],
        ["Signal_S0_S_i_M600_cpl1_scalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M600_cpl1_scalar_skims_nominal.root"],
        ["Signal_S0_S_i_M700_cpl1_scalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M700_cpl1_scalar_skims_nominal.root"],
        ["Signal_S0_S_i_M800_cpl1_scalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M800_cpl1_scalar_skims_nominal.root"],

        ["Signal_S0_S_i_M400_cpl1_pseudoscalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M400_cpl1_pseudoscalar_skims_nominal.root"],
        ["Signal_S0_S_i_M500_cpl1_pseudoscalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M500_cpl1_pseudoscalar_skims_nominal.root"],
        ["Signal_S0_S_i_M600_cpl1_pseudoscalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M600_cpl1_pseudoscalar_skims_nominal.root"],
        ["Signal_S0_S_i_M700_cpl1_pseudoscalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M700_cpl1_pseudoscalar_skims_nominal.root"],
        ["Signal_S0_S_i_M800_cpl1_pseudoscalar_dataset_nominal.root", "skims/%s/Signal_S0_S_i_M800_cpl1_pseudoscalar_skims_nominal.root"],

        ["Signal_ZPrimeToTTJets_M500GeV_W5GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M500GeV_W5GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M1000GeV_W10GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1000GeV_W10GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M1500GeV_W15GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1500GeV_W15GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M2000GeV_W20GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M2000GeV_W20GeV_skims_nominal_merged.root"],
	
        ["Signal_ZPrimeToTTJets_M500GeV_W50GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M500GeV_W50GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M750GeV_W75GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M750GeV_W75GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M1000GeV_W100GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1000GeV_W100GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M1250GeV_W125GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1250GeV_W125GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M1500GeV_W150GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1500GeV_W150GeV_skims_nominal_merged.root"],
        ["Signal_ZPrimeToTTJets_M2000GeV_W200GeV_dataset_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M2000GeV_W200GeV_skims_nominal_merged.root"],
        ]

def launch(input, output, weight):
  args = ["./extractor2Dataset", "-i", input, "-o", output, "--mc", "--skim"]
  if "semie" in input:
    args.append("--type semie")
  elif "semimu" in input:
    args.append("--type semimu")

  args.append("--weight %.15f" % weight)

  return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

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
        if not "%" in file[1] and not type in file[1]:
            continue
        path = "datasets/%s/%s" % (d, type)
        if not "%" in file[1]:
            tmpfile.write(launch(file[1], os.path.join(path, file[0]), weight) + "\n");
        else:
            tmpfile.write(launch(file[1] % type, os.path.join(path, file[0]), weight) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 

print args
subprocess.call(args)

# Merge files
for file in files:
    path = "datasets/%s/" % (d)
    merged_file = os.path.join(path, file[0].replace(".root", "_merged.root"))
    args = ["hadd", "-f", merged_file]
    for type in ["semie", "semimu"]:
        if not "%" in file[1] and not type in file[1]:
            continue
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
