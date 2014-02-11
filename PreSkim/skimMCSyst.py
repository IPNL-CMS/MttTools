#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
        # Background + Signal
        ["Signal_S0_S_i_M400_cpl1_scalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M400_cpl1_scalar_%s_%s.list"],
        ["Signal_S0_S_i_M500_cpl1_scalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M500_cpl1_scalar_%s_%s.list"],
        ["Signal_S0_S_i_M600_cpl1_scalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M600_cpl1_scalar_%s_%s.list"],
        ["Signal_S0_S_i_M700_cpl1_scalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M700_cpl1_scalar_%s_%s.list"],
        ["Signal_S0_S_i_M800_cpl1_scalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M800_cpl1_scalar_%s_%s.list"],

        ["Signal_S0_S_i_M400_cpl1_pseudoscalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M400_cpl1_pseudoscalar_%s_%s.list"],
        ["Signal_S0_S_i_M500_cpl1_pseudoscalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M500_cpl1_pseudoscalar_%s_%s.list"],
        ["Signal_S0_S_i_M600_cpl1_pseudoscalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M600_cpl1_pseudoscalar_%s_%s.list"],
        ["Signal_S0_S_i_M700_cpl1_pseudoscalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M700_cpl1_pseudoscalar_%s_%s.list"],
        ["Signal_S0_S_i_M800_cpl1_pseudoscalar_skims_%s.root", "lists/Systematics/MTT_S0_S_i_M800_cpl1_pseudoscalar_%s_%s.list"],

        # Signal Z'
        ["Signal_ZPrimeToTTJets_M500GeV_W50GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M500GeV_W50GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M750GeV_W75GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M750GeV_W75GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1000GeV_W100GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1000GeV_W100GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1250GeV_W125GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1250GeV_W125GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1500GeV_W150GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1500GeV_W150GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M2000GeV_W200GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M2000GeV_W200GeV_%s_%s.list"],

        ["Signal_ZPrimeToTTJets_M500GeV_W50GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M500GeV_W50GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M750GeV_W75GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M750GeV_W75GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1000GeV_W100GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1000GeV_W100GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1250GeV_W125GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1250GeV_W125GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1500GeV_W150GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1500GeV_W150GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M2000GeV_W200GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M2000GeV_W200GeV_ext_%s_%s.list"],

        ["Signal_ZPrimeToTTJets_M500GeV_W5GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M500GeV_W5GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M750GeV_W7p5GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1000GeV_W10GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1000GeV_W10GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1250GeV_W12p5GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1500GeV_W15GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1500GeV_W15GeV_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M2000GeV_W20GeV_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M2000GeV_W20GeV_%s_%s.list"],

        ["Signal_ZPrimeToTTJets_M500GeV_W5GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M500GeV_W5GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M750GeV_W7p5GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1000GeV_W10GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1000GeV_W10GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1250GeV_W12p5GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M1500GeV_W15GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M1500GeV_W15GeV_ext_%s_%s.list"],
        ["Signal_ZPrimeToTTJets_M2000GeV_W20GeV_ext_skims_%s.root", "lists/Systematics/MTT_ZPrimeToTTJets_M2000GeV_W20GeV_ext_%s_%s.list"],

        # Background
        ["MC_TT_powheg_skims_%s.root", "lists/Systematics/MTT_TT_powheg_%s_%s.list"],

        ["MC_T_tW-channel_skims_%s.root", "lists/Systematics/MTT_T_tW-channel_%s_%s.list"],
        ["MC_T_s-channel_skims_%s.root", "lists/Systematics/MTT_T_s-channel_%s_%s.list"],
        ["MC_T_t-channel_skims_%s.root", "lists/Systematics/MTT_T_t-channel_%s_%s.list"],

        ["MC_Tbar_tW-channel_skims_%s.root", "lists/Systematics/MTT_Tbar_tW-channel_%s_%s.list"],
        ["MC_Tbar_s-channel_skims_%s.root", "lists/Systematics/MTT_Tbar_s-channel_%s_%s.list"],
        ["MC_Tbar_t-channel_skims_%s.root", "lists/Systematics/MTT_Tbar_t-channel_%s_%s.list"],

        ["MC_DY1JetsToLL_M-50_skims_%s.root", "lists/Systematics/MTT_DY1JetsToLL_M-50_%s_%s.list"],
        ["MC_DY2JetsToLL_M-50_skims_%s.root", "lists/Systematics/MTT_DY2JetsToLL_M-50_%s_%s.list"],
        ["MC_DY3JetsToLL_M-50_skims_%s.root", "lists/Systematics/MTT_DY3JetsToLL_M-50_%s_%s.list"],
        ["MC_DY4JetsToLL_M-50_skims_%s.root", "lists/Systematics/MTT_DY4JetsToLL_M-50_%s_%s.list"],

        ["MC_W1JetsToLNu_skims_%s.root", "lists/Systematics/MTT_W1JetsToLNu_%s_%s.list"],
        ["MC_W2JetsToLNu_skims_%s.root", "lists/Systematics/MTT_W2JetsToLNu_%s_%s.list"],
        ["MC_W3JetsToLNu_skims_%s.root", "lists/Systematics/MTT_W3JetsToLNu_%s_%s.list"],
        ["MC_W4JetsToLNu_skims_%s.root", "lists/Systematics/MTT_W4JetsToLNu_%s_%s.list"],

        ["MC_QCD_Pt_20_30_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_20_30_EMEnriched_%s_semie.list"],
        ["MC_QCD_Pt_30_80_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_30_80_EMEnriched_%s_semie.list"],
        ["MC_QCD_Pt_80_170_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_80_170_EMEnriched_%s_semie.list"],
        ["MC_QCD_Pt_170_250_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_170_250_EMEnriched_%s_semie.list"],
        ["MC_QCD_Pt_250_350_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_250_350_EMEnriched_%s_semie.list"],
        ["MC_QCD_Pt_350_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_350_EMEnriched_%s_semie.list"],

        ["MC_QCD_Pt_30_80_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_30_80_BCtoE_%s_semie.list"],
        ["MC_QCD_Pt_80_170_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_80_170_BCtoE_%s_semie.list"],
        ["MC_QCD_Pt_170_250_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_170_250_BCtoE_%s_semie.list"],
        ["MC_QCD_Pt_250_350_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_250_350_BCtoE_%s_semie.list"],
        ["MC_QCD_Pt_350_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_350_BCtoE_%s_semie.list"],

        ["MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_skims_%s.root", "lists/Systematics/MTT_QCD_pt15to30_bEnriched_MuEnrichedPt14_%s_semimu.list"],
        ["MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_skims_%s.root", "lists/Systematics/MTT_QCD_pt30to50_bEnriched_MuEnrichedPt14_%s_semimu.list"],
        ["MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_skims_%s.root", "lists/Systematics/MTT_QCD_pt50to150_bEnriched_MuEnrichedPt14_%s_semimu.list"],
        ["MC_QCD_pt150_bEnriched_MuEnrichedPt14_skims_%s.root", "lists/Systematics/MTT_QCD_pt150_bEnriched_MuEnrichedPt14_%s_semimu.list"],
        ]

def launch(input, output):
    args = ["./preSkim", "--input-list", input, "-o", output, "--mc"]
    if "semie" in input:
        args.append("--semie")
    elif "semimu" in input:
        args.append("--semimu")

    #args.append("--weight %.15f" % weight)

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

systs = ["JECup", "JECdown", "JERup", "JERdown"]

# Build output tree structure
for type in ["semie", "semimu"]:
    path = "skims/%s/%s/Systematics" % (d, type)
    try:
        os.makedirs(path)
    except:
        pass

print("Skimming datasets...")

for file in files:
    for type in ["semie", "semimu"]:
        path = "skims/%s/%s/Systematics" % (d, type)
        for syst in systs:
            if "Enriched" in file[1] or "BCtoE" in file[1]:
                if ("MuEnriched" in file[1] and type == "semimu") or (("EMEnriched" in file[1] or "BCtoE" in file[1]) and type == "semie"):
                    tmpfile.write(launch(file[1] % syst, os.path.join(path, file[0] % syst)) + "\n");
            else:
                tmpfile.write(launch(file[1] % (syst, type), os.path.join(path, file[0] % syst)) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "8"] 
subprocess.call(args)