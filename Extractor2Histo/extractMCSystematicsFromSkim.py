#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
        # Background + Signal
        ["Signal_S0_S_i_M400_cpl1_scalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M400_cpl1_scalar_skims_%s.root"],
        ["Signal_S0_S_i_M500_cpl1_scalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M500_cpl1_scalar_skims_%s.root"],
        ["Signal_S0_S_i_M600_cpl1_scalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M600_cpl1_scalar_skims_%s.root"],
        ["Signal_S0_S_i_M700_cpl1_scalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M700_cpl1_scalar_skims_%s.root"],
        ["Signal_S0_S_i_M800_cpl1_scalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M800_cpl1_scalar_skims_%s.root"],

        ["Signal_S0_S_i_M400_cpl1_pseudoscalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_skims_%s.root"],
        ["Signal_S0_S_i_M500_cpl1_pseudoscalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_skims_%s.root"],
        ["Signal_S0_S_i_M600_cpl1_pseudoscalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_skims_%s.root"],
        ["Signal_S0_S_i_M700_cpl1_pseudoscalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_skims_%s.root"],
        ["Signal_S0_S_i_M800_cpl1_pseudoscalar_histos_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_skims_%s.root"],

        # Signal Z'
        ["Signal_ZPrimeToTTJets_M500GeV_W50GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M500GeV_W50GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M750GeV_W75GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M750GeV_W75GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M1000GeV_W100GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1000GeV_W100GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M1250GeV_W125GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1250GeV_W125GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M1500GeV_W150GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1500GeV_W150GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M2000GeV_W200GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M2000GeV_W200GeV_skims_%s_merged.root"],

        ["Signal_ZPrimeToTTJets_M500GeV_W5GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M500GeV_W5GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M1000GeV_W10GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1000GeV_W10GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M1500GeV_W15GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1500GeV_W15GeV_skims_%s_merged.root"],
        ["Signal_ZPrimeToTTJets_M2000GeV_W20GeV_histos_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M2000GeV_W20GeV_skims_%s_merged.root"],

        # Background
        ["MC_TT_powheg_histos_%s.root", "skims/%s/Systematics/MC_TT_powheg_skims_%s.root"],

        ["MC_T_tW-channel_histos_%s.root", "skims/%s/Systematics/MC_T_tW-channel_skims_%s.root"],
        ["MC_T_s-channel_histos_%s.root", "skims/%s/Systematics/MC_T_s-channel_skims_%s.root"],
        ["MC_T_t-channel_histos_%s.root", "skims/%s/Systematics/MC_T_t-channel_skims_%s.root"],

        ["MC_Tbar_tW-channel_histos_%s.root", "skims/%s/Systematics/MC_Tbar_tW-channel_skims_%s.root"],
        ["MC_Tbar_s-channel_histos_%s.root", "skims/%s/Systematics/MC_Tbar_s-channel_skims_%s.root"],
        ["MC_Tbar_t-channel_histos_%s.root", "skims/%s/Systematics/MC_Tbar_t-channel_skims_%s.root"],

        ["MC_DY1JetsToLL_M-50_histos_%s.root", "skims/%s/Systematics/MC_DY1JetsToLL_M-50_skims_%s.root"],
        ["MC_DY2JetsToLL_M-50_histos_%s.root", "skims/%s/Systematics/MC_DY2JetsToLL_M-50_skims_%s.root"],
        ["MC_DY3JetsToLL_M-50_histos_%s.root", "skims/%s/Systematics/MC_DY3JetsToLL_M-50_skims_%s.root"],
        ["MC_DY4JetsToLL_M-50_histos_%s.root", "skims/%s/Systematics/MC_DY4JetsToLL_M-50_skims_%s.root"],

        ["MC_W1JetsToLNu_histos_%s.root", "skims/%s/Systematics/MC_W1JetsToLNu_skims_%s.root"],
        ["MC_W2JetsToLNu_histos_%s.root", "skims/%s/Systematics/MC_W2JetsToLNu_skims_%s.root"],
        ["MC_W3JetsToLNu_histos_%s.root", "skims/%s/Systematics/MC_W3JetsToLNu_skims_%s.root"],
        ["MC_W4JetsToLNu_histos_%s.root", "skims/%s/Systematics/MC_W4JetsToLNu_skims_%s.root"],

        ["MC_QCD_Pt_20_30_EMEnriched_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_20_30_EMEnriched_skims_%s.root"],
        ["MC_QCD_Pt_30_80_EMEnriched_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_30_80_EMEnriched_skims_%s.root"],
        ["MC_QCD_Pt_80_170_EMEnriched_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_80_170_EMEnriched_skims_%s.root"],
        ["MC_QCD_Pt_170_250_EMEnriched_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_170_250_EMEnriched_skims_%s.root"],
        ["MC_QCD_Pt_250_350_EMEnriched_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_250_350_EMEnriched_skims_%s.root"],
        ["MC_QCD_Pt_350_EMEnriched_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_350_EMEnriched_skims_%s.root"],

        ["MC_QCD_Pt_30_80_BCtoE_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_30_80_BCtoE_skims_%s.root"],
        ["MC_QCD_Pt_80_170_BCtoE_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_80_170_BCtoE_skims_%s.root"],
        ["MC_QCD_Pt_170_250_BCtoE_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_170_250_BCtoE_skims_%s.root"],
        ["MC_QCD_Pt_250_350_BCtoE_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_250_350_BCtoE_skims_%s.root"],
        ["MC_QCD_Pt_350_BCtoE_histos_%s.root", "skims/semie/Systematics/MC_QCD_Pt_350_BCtoE_skims_%s.root"],

        ["MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_histos_%s.root", "skims/semimu/Systematics/MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_skims_%s.root"],
        ["MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_histos_%s.root", "skims/semimu/Systematics/MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_skims_%s.root"],
        ["MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_histos_%s.root", "skims/semimu/Systematics/MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_skims_%s.root"],
        ["MC_QCD_pt150_bEnriched_MuEnrichedPt14_histos_%s.root", "skims/semimu/Systematics/MC_QCD_pt150_bEnriched_MuEnrichedPt14_skims_%s.root"],
        ]

systs = {"JECup": ["JECup", "--jec-syst up"], "JECdown": ["JECdown", "--jec-syst down"], "JERup": ["JERup", ""], "JERdown": ["JERdown", ""], "puUp": ["nominal", "--pileup-syst up"], "puDown": ["nominal", "--pileup-syst down"]}

if True:
    systs["trigUp"] = ["nominal", "--trigger-syst up"]
    systs["trigDown"] = ["nominal", "--trigger-syst down"]

if False:
    systs["pdfUp"] = ["nominal", "--pdf-syst up"]
    systs["pdfDown"] = ["nominal", "--pdf-syst down"]

def launch(input, output, btag, extra):
    args = ["./extractorToHisto", "-i", input, "-o", output, "--mc", "--skim", "--b-tag", str(btag)]

    if len(extra) > 0:
        args.append(extra)

    if "semie" in input:
        args.append("--semie")
    elif "semimu" in input:
        args.append("--semimu")

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

# Build output tree structure
for btag in [1, 2]:
    for type in ["semie", "semimu"]:
        path = "plots/%s/%d-btag/%s/Systematics" % (d, btag, type)
        try:
            os.makedirs(path)
        except:
            pass

print("Extracting datasets...")

for file in files:
    for btag in [1, 2]:
        for type in ["semie", "semimu"]:
            for syst, extra in systs.items():
                if not "skims/%s" in file[1] and not type in file[1]:
                    continue

                path = "plots/%s/%d-btag/%s/Systematics" % (d, btag, type)
                if not "skims/%s" in file[1]:
                    inputFile = file[1] % extra[0]
                else:
                    inputFile = file[1] % (type, extra[0])

                if "nominal" in extra[0]:
                    # Remove the "Systematics/" part of the path
                    inputFile = inputFile.replace("Systematics/", "")

                tmpfile.write(launch(inputFile, os.path.join(path, file[0] % syst), btag, extra[1]) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "30"] 
subprocess.call(args)
