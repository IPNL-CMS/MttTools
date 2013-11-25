#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
        # Background
        ["MC_TT_powheg_histos_nominal.root", "lists/MTT_TT_powheg_%s.list"],

        ["MC_T_tW-channel_histos_nominal.root", "lists/MTT_T_tW-channel_%s.list"],
        ["MC_T_s-channel_histos_nominal.root", "lists/MTT_T_s-channel_%s.list"],
        ["MC_T_t-channel_histos_nominal.root", "lists/MTT_T_t-channel_%s.list"],

        ["MC_Tbar_tW-channel_histos_nominal.root", "lists/MTT_Tbar_tW-channel_%s.list"],
        ["MC_Tbar_s-channel_histos_nominal.root", "lists/MTT_Tbar_s-channel_%s.list"],
        ["MC_Tbar_t-channel_histos_nominal.root", "lists/MTT_Tbar_t-channel_%s.list"],

        ["MC_DY1JetsToLL_M-50_histos_nominal.root", "lists/MTT_DY1JetsToLL_M-50_%s.list"],
        ["MC_DY2JetsToLL_M-50_histos_nominal.root", "lists/MTT_DY2JetsToLL_M-50_%s.list"],
        ["MC_DY3JetsToLL_M-50_histos_nominal.root", "lists/MTT_DY3JetsToLL_M-50_%s.list"],
        ["MC_DY4JetsToLL_M-50_histos_nominal.root", "lists/MTT_DY4JetsToLL_M-50_%s.list"],

        ["MC_W1JetsToLNu_histos_nominal.root", "lists/MTT_W1JetsToLNu_%s.list"],
        ["MC_W2JetsToLNu_histos_nominal.root", "lists/MTT_W2JetsToLNu_%s.list"],
        ["MC_W3JetsToLNu_histos_nominal.root", "lists/MTT_W3JetsToLNu_%s.list"],
        ["MC_W4JetsToLNu_histos_nominal.root", "lists/MTT_W4JetsToLNu_%s.list"],

        # Background + Signal
        ["Signal_S0_S_i_M500_cpl1_scalar_histos_nominal.root", "lists/MTT_S0_S_i_M500_cpl1_scalar_%s.list"],
        ["Signal_S0_S_i_M700_cpl1_scalar_histos_nominal.root", "lists/MTT_S0_S_i_M700_cpl1_scalar_%s.list"],
        ["Signal_S0_S_i_M500_cpl1_pseudoscalar_histos_nominal.root", "lists/MTT_S0_S_i_M500_cpl1_pseudoscalar_%s.list"],
        ["Signal_S0_S_i_M700_cpl1_pseudoscalar_histos_nominal.root", "lists/MTT_S0_S_i_M700_cpl1_pseudoscalar_%s.list"],
        ]

def launch(input, output, btag):
    args = ["./extractorToHisto", "--input-list", input, "-o", output, "--mc", "--b-tag", str(btag)]
    if "semie" in input:
        args.append("--semie")
    elif "semimu" in input:
        args.append("--semimu")

    #args.append("--weight %.15f" % weight)

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
    #for btag in [1, 2]:
    for btag in [1]:
        #for type in ["semie", "semimu"]:
        for type in ["semie", "semimu"]:
            path = "plots/%s/%d-btag/%s" % (d, btag, type)
            tmpfile.write(launch(file[1] % type, "%s/%s" % (path, file[0]), btag) + "\n");

tmpfile.flush()

print tmpfile.name
args = ["parallel", "-u", "-a", tmpfile.name, "-j", "32"] 
subprocess.call(args)
