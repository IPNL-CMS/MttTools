#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
        # Background
        #["MC_TT_powheg_histos_nominal.root", "skims/%s/MC_TT_powheg_skims_nominal.root"],

        #["MC_T_tW-channel_histos_nominal.root", "skims/%s/MC_T_tW-channel_skims_nominal.root"],
        #["MC_T_s-channel_histos_nominal.root", "skims/%s/MC_T_s-channel_skims_nominal.root"],
        #["MC_T_t-channel_histos_nominal.root", "skims/%s/MC_T_t-channel_skims_nominal.root"],

        #["MC_Tbar_tW-channel_histos_nominal.root", "skims/%s/MC_Tbar_tW-channel_skims_nominal.root"],
        #["MC_Tbar_s-channel_histos_nominal.root", "skims/%s/MC_Tbar_s-channel_skims_nominal.root"],
        #["MC_Tbar_t-channel_histos_nominal.root", "skims/%s/MC_Tbar_t-channel_skims_nominal.root"],

        #["MC_DY1JetsToLL_M-50_histos_nominal.root", "skims/%s/MC_DY1JetsToLL_M-50_skims_nominal.root"],
        #["MC_DY2JetsToLL_M-50_histos_nominal.root", "skims/%s/MC_DY2JetsToLL_M-50_skims_nominal.root"],
        #["MC_DY3JetsToLL_M-50_histos_nominal.root", "skims/%s/MC_DY3JetsToLL_M-50_skims_nominal.root"],
        #["MC_DY4JetsToLL_M-50_histos_nominal.root", "skims/%s/MC_DY4JetsToLL_M-50_skims_nominal.root"],

        #["MC_W1JetsToLNu_histos_nominal.root", "skims/%s/MC_W1JetsToLNu_skims_nominal.root"],
        #["MC_W2JetsToLNu_histos_nominal.root", "skims/%s/MC_W2JetsToLNu_skims_nominal.root"],
        #["MC_W3JetsToLNu_histos_nominal.root", "skims/%s/MC_W3JetsToLNu_skims_nominal.root"],
        ["MC_W4JetsToLNu_histos_nominal.root", "skims/%s/MC_W4JetsToLNu_skims_nominal.root"],

        ## Background + Signal
        #["Signal_S0_S_i_M500_cpl1_scalar_histos_nominal.root", "skims/%s/Signal_S0_S_i_M500_cpl1_scalar_skims_nominal.root"],
        #["Signal_S0_S_i_M700_cpl1_scalar_histos_nominal.root", "skims/%s/Signal_S0_S_i_M700_cpl1_scalar_skims_nominal.root"],
        #["Signal_S0_S_i_M500_cpl1_pseudoscalar_histos_nominal.root", "skims/%s/Signal_S0_S_i_M500_cpl1_pseudoscalar_skims_nominal.root"],
        #["Signal_S0_S_i_M700_cpl1_pseudoscalar_histos_nominal.root", "skims/%s/Signal_S0_S_i_M700_cpl1_pseudoscalar_skims_nominal.root"],
        ]

def launch(input, output, btag):
    args = ["./extractorToHisto", "-i", input, "-o", output, "--mc", "--skim", "--b-tag", str(btag)]
    if "semie" in input:
        args.append("--semie")
    elif "semimu" in input:
        args.append("--semimu")

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

# Build output tree structure
for btag in [1, 2]:
    for type in ["semie", "semimu"]:
        path = "plots/%s/%d-btag/%s" % (d, btag, type)
        try:
            os.makedirs(path)
        except:
            pass

print("Extracting datasets...")

for file in files:
    for btag in [1, 2]:
    #for btag in [2]:
        #for type in ["semie", "semimu"]:
        for type in ["semie", "semimu"]:
            path = "plots/%s/%d-btag/%s" % (d, btag, type)
            tmpfile.write(launch(file[1] % type, "%s/%s" % (path, file[0]), btag) + "\n");

tmpfile.flush()

print tmpfile.name
args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 
subprocess.call(args)
