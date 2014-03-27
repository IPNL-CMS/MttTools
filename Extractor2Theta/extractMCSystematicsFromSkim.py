#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
        # Background + Signal
        #["Signal_S0_S_i_M400_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M400_cpl1_scalar_skims_%s.root"],
        #["Signal_S0_S_i_M500_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M500_cpl1_scalar_skims_%s.root"],
        #["Signal_S0_S_i_M600_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M600_cpl1_scalar_skims_%s.root"],
        #["Signal_S0_S_i_M700_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M700_cpl1_scalar_skims_%s.root"],
        #["Signal_S0_S_i_M800_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M800_cpl1_scalar_skims_%s.root"],

        #["Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_skims_%s.root"],
        #["Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_skims_%s.root"],
        #["Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_skims_%s.root"],
        #["Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_skims_%s.root"],
        #["Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_skims_%s.root"],

        # Signal Z'
        #["Signal_ZPrimeToTTJets_M500GeV_W50GeV_dataset_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M500GeV_W50GeV_skims_%s_merged.root"],
        #["Signal_ZPrimeToTTJets_M750GeV_W75GeV_dataset_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M750GeV_W75GeV_skims_%s_merged.root"],
        #["Signal_ZPrimeToTTJets_M1000GeV_W100GeV_dataset_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1000GeV_W100GeV_skims_%s_merged.root"],
        #["Signal_ZPrimeToTTJets_M1250GeV_W125GeV_dataset_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1250GeV_W125GeV_skims_%s_merged.root"],
        #["Signal_ZPrimeToTTJets_M1500GeV_W150GeV_dataset_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1500GeV_W150GeV_skims_%s_merged.root"],
        #["Signal_ZPrimeToTTJets_M2000GeV_W200GeV_dataset_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M2000GeV_W200GeV_skims_%s_merged.root"],

        #["signal_zprimetottjets_m500gev_w5gev_dataset_%s.root", "skims/%s/systematics/signal_zprimetottjets_m500gev_w5gev_skims_%s_merged.root"],
        #["signal_zprimetottjets_m750gev_w7p5gev_dataset_%s.root", "skims/%s/systematics/signal_zprimetottjets_m750gev_w7p5gev_skims_%s_merged.root"],
        #["signal_zprimetottjets_m1000gev_w10gev_dataset_%s.root", "skims/%s/systematics/signal_zprimetottjets_m1000gev_w10gev_skims_%s_merged.root"],
        #["signal_zprimetottjets_m1250gev_w12p5gev_dataset_%s.root", "skims/%s/systematics/signal_zprimetottjets_m1250gev_w12p5gev_skims_%s_merged.root"],
        #["signal_zprimetottjets_m1500gev_w15gev_dataset_%s.root", "skims/%s/systematics/signal_zprimetottjets_m1500gev_w15gev_skims_%s_merged.root"],
        #["signal_zprimetottjets_m2000gev_w20gev_dataset_%s.root", "skims/%s/systematics/signal_zprimetottjets_m2000gev_w20gev_skims_%s_merged.root"],


        ## Background
        ["MC_TT_powheg_theta_%s.root", "skims/%s/Systematics/MC_TT_powheg_skims_%s.root", 21675970, 245.8],

        ["MC_T_tW-channel_theta_%s.root", "skims/%s/Systematics/MC_T_tW-channel_skims_%s.root", 497658, 11.1],
        ["MC_T_s-channel_theta_%s.root", "skims/%s/Systematics/MC_T_s-channel_skims_%s.root", 259961, 3.79],
        ["MC_T_t-channel_theta_%s.root", "skims/%s/Systematics/MC_T_t-channel_skims_%s.root", 3728227, 56.4],

        ["MC_Tbar_tW-channel_theta_%s.root", "skims/%s/Systematics/MC_Tbar_tW-channel_skims_%s.root", 493460, 11.1],
        ["MC_Tbar_s-channel_theta_%s.root", "skims/%s/Systematics/MC_Tbar_s-channel_skims_%s.root", 139974, 1.76],
        ["MC_Tbar_t-channel_theta_%s.root", "skims/%s/Systematics/MC_Tbar_t-channel_skims_%s.root", 1935072, 30.7],

        ["MC_DY1JetsToLL_M-50_theta_%s.root", "skims/%s/Systematics/MC_DY1JetsToLL_M-50_skims_%s.root", 24045248, 666.3],
        ["MC_DY2JetsToLL_M-50_theta_%s.root", "skims/%s/Systematics/MC_DY2JetsToLL_M-50_skims_%s.root", 21852156, 215.0],
        ["MC_DY3JetsToLL_M-50_theta_%s.root", "skims/%s/Systematics/MC_DY3JetsToLL_M-50_skims_%s.root", 11015445, 60.7],
        ["MC_DY4JetsToLL_M-50_theta_%s.root", "skims/%s/Systematics/MC_DY4JetsToLL_M-50_skims_%s.root", 6402827, 27.3],

        ["MC_W1JetsToLNu_theta_%s.root", "skims/%s/Systematics/MC_W1JetsToLNu_skims_%s.root", 23141598, 6662.8],
        ["MC_W2JetsToLNu_theta_%s.root", "skims/%s/Systematics/MC_W2JetsToLNu_skims_%s.root", 34044921, 2159.2],
        ["MC_W3JetsToLNu_theta_%s.root", "skims/%s/Systematics/MC_W3JetsToLNu_skims_%s.root", 15539503, 640.4],
        ["MC_W4JetsToLNu_theta_%s.root", "skims/%s/Systematics/MC_W4JetsToLNu_skims_%s.root", 13382803, 264.0],

        #["MC_QCD_Pt_20_30_EMEnriched_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_20_30_EMEnriched_skims_%s.root"],
        #["MC_QCD_Pt_30_80_EMEnriched_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_30_80_EMEnriched_skims_%s.root"],
        #["MC_QCD_Pt_80_170_EMEnriched_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_80_170_EMEnriched_skims_%s.root"],
        #["MC_QCD_Pt_170_250_EMEnriched_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_170_250_EMEnriched_skims_%s.root"],
        #["MC_QCD_Pt_250_350_EMEnriched_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_250_350_EMEnriched_skims_%s.root"],
        #["MC_QCD_Pt_350_EMEnriched_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_350_EMEnriched_skims_%s.root"],

        #["MC_QCD_Pt_30_80_BCtoE_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_30_80_BCtoE_skims_%s.root"],
        #["MC_QCD_Pt_80_170_BCtoE_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_80_170_BCtoE_skims_%s.root"],
        #["MC_QCD_Pt_170_250_BCtoE_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_170_250_BCtoE_skims_%s.root"],
        #["MC_QCD_Pt_250_350_BCtoE_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_250_350_BCtoE_skims_%s.root"],
        #["MC_QCD_Pt_350_BCtoE_theta_%s.root", "skims/semie/Systematics/MC_QCD_Pt_350_BCtoE_skims_%s.root"],

        #["MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_theta_%s.root", "skims/semimu/Systematics/MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_skims_%s.root"],
        #["MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_theta_%s.root", "skims/semimu/Systematics/MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_skims_%s.root"],
        #["MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_theta_%s.root", "skims/semimu/Systematics/MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_skims_%s.root"],
        #["MC_QCD_pt150_bEnriched_MuEnrichedPt14_theta_%s.root", "skims/semimu/Systematics/MC_QCD_pt150_bEnriched_MuEnrichedPt14_skims_%s.root"],
        ]

systs = {}

if True:
    systs = {"JECup": ["JECup", "--jec-syst up"], "JECdown": ["JECdown", "--jec-syst down"], "JERup": ["JERup", ""], "JERdown": ["JERdown", ""], "puUp": ["nominal", "--pileup-syst up"], "puDown": ["nominal", "--pileup-syst down"]}

if True:
    systs["trigUp"] = ["nominal", "--trigger-syst up"]
    systs["trigDown"] = ["nominal", "--trigger-syst down"]

if True:
    systs["leptUp"] = ["nominal", "--lepton-syst up"]
    systs["leptDown"] = ["nominal", "--lepton-syst down"]

if True:
    systs["btagUp"] = ["nominal", "--btag-syst up"]
    systs["btagDown"] = ["nominal", "--btag-syst down"]

if False:
    systs["pdfUp"] = ["nominal", "--pdf-syst up"]
    systs["pdfDown"] = ["nominal", "--pdf-syst down"]

def launch(input, output, weight, extra):
    args = ["./extractor2Theta", "-i", input, "-o", output, "--mc", "--skim"]

    if len(extra) > 0:
        args.append(extra)

    if "semie" in input:
        args.append("--type semie")
    elif "semimu" in input:
        args.append("--type semimu")

    args.append("--weight %.15f" % weight)

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

# Build output tree structure
for type in ["semie", "semimu"]:
    path = "theta/%s/%s/Systematics" % (d, type)
    try:
        os.makedirs(path)
    except:
        pass

print("Extracting datasets...")

for file in files:
    for type in ["semie", "semimu"]:
        for syst, extra in systs.items():
            if not "skims/%s" in file[1] and not type in file[1]:
                continue

            path = "theta/%s/%s/Systematics" % (d, type)
            if not "skims/%s" in file[1]:
                inputFile = file[1] % extra[0]
            else:
                inputFile = file[1] % (type, extra[0])

            if "nominal" in extra[0]:
                # Remove the "Systematics/" part of the path
                inputFile = inputFile.replace("Systematics/", "")

            events = file[2] if len(file) >= 3 else -1
            xsection = file[3] if len(file) >= 4 else -1
            weight = 19667. * xsection / events if events > 0 else 1;
            print("Weight for %s: %.15f" % (file[1], weight))

            tmpfile.write(launch(inputFile, os.path.join(path, file[0] % syst), weight, extra[1]) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "20"] 
subprocess.call(args)
