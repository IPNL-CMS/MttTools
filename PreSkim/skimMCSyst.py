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

        # Signal RS gluons
        ["Signal_RSGluonToTT_M-700_skims_%s.root", "lists/Systematics/MTT_RSGluonToTT_M-700_%s_%s.list"],
        ["Signal_RSGluonToTT_M-1000_skims_%s.root", "lists/Systematics/MTT_RSGluonToTT_M-1000_%s_%s.list"],
        ["Signal_RSGluonToTT_M-1200_skims_%s.root", "lists/Systematics/MTT_RSGluonToTT_M-1200_%s_%s.list"],
        ["Signal_RSGluonToTT_M-1500_skims_%s.root", "lists/Systematics/MTT_RSGluonToTT_M-1500_%s_%s.list"],
        ["Signal_RSGluonToTT_M-2000_skims_%s.root", "lists/Systematics/MTT_RSGluonToTT_M-2000_%s_%s.list"],

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
        ["MC_TT_madgraph_skims_%s.root", "lists/Systematics/MTT_TT_madgraph_%s_%s.list"],
        ["MC_TT_mcatnlo_skims_%s.root", "lists/Systematics/MTT_TT_mcatnlo_%s_%s.list"],
        ["MC_TT_madgraph_semilept_skims_%s.root", "lists/Systematics/MTT_TT_madgraph_semilept_%s_%s.list"],
        ["MC_TT_madgraph_dilept_skims_%s.root", "lists/Systematics/MTT_TT_madgraph_dilept_%s_%s.list"],
        ["MC_TT_madgraph_hadronic_skims_%s.root", "lists/Systematics/MTT_TT_madgraph_hadronic_%s_%s.list"],

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

        ["MC_WW_skims_%s.root", "lists/Systematics/MTT_WW_%s_%s.list"],
        ["MC_WZ_skims_%s.root", "lists/Systematics/MTT_WZ_%s_%s.list"],
        ["MC_ZZ_skims_%s.root", "lists/Systematics/MTT_ZZ_%s_%s.list"],

        #["MC_QCD_Pt_20_30_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_20_30_EMEnriched_%s_semie.list"],
        #["MC_QCD_Pt_30_80_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_30_80_EMEnriched_%s_semie.list"],
        #["MC_QCD_Pt_80_170_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_80_170_EMEnriched_%s_semie.list"],
        #["MC_QCD_Pt_170_250_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_170_250_EMEnriched_%s_semie.list"],
        #["MC_QCD_Pt_250_350_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_250_350_EMEnriched_%s_semie.list"],
        #["MC_QCD_Pt_350_EMEnriched_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_350_EMEnriched_%s_semie.list"],

        #["MC_QCD_Pt_30_80_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_30_80_BCtoE_%s_semie.list"],
        #["MC_QCD_Pt_80_170_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_80_170_BCtoE_%s_semie.list"],
        #["MC_QCD_Pt_170_250_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_170_250_BCtoE_%s_semie.list"],
        #["MC_QCD_Pt_250_350_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_250_350_BCtoE_%s_semie.list"],
        #["MC_QCD_Pt_350_BCtoE_skims_%s.root", "lists/Systematics/MTT_QCD_Pt_350_BCtoE_%s_semie.list"],

        #["MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_skims_%s.root", "lists/Systematics/MTT_QCD_pt15to30_bEnriched_MuEnrichedPt14_%s_semimu.list"],
        #["MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_skims_%s.root", "lists/Systematics/MTT_QCD_pt30to50_bEnriched_MuEnrichedPt14_%s_semimu.list"],
        #["MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_skims_%s.root", "lists/Systematics/MTT_QCD_pt50to150_bEnriched_MuEnrichedPt14_%s_semimu.list"],
        #["MC_QCD_pt150_bEnriched_MuEnrichedPt14_skims_%s.root", "lists/Systematics/MTT_QCD_pt150_bEnriched_MuEnrichedPt14_%s_semimu.list"],

        ## Macthing / scale  up & down
        #["MC_TT_madgraph_skims_%s.root", "lists/Systematics/MTT_TTJets_%s_%s.list"],
        #["MC_DYJetsToLL_M-50_skims_%s.root", "lists/Systematics/MTT_DYJetsToLL_M-50_%s_%s.list"],
        #["MC_WJetsToLNu_skims_%s.root", "lists/Systematics/MTT_WJetsToLNu_%s_%s.list"]
        ]

def launch(input, output):
    if not os.path.exists(input):
        print("Warning: input file '%s' not found. Skipping this job." % input)
        return ""

    args = ["./preSkim", "--input-list", input, "-o", output, "--mc"]
    if "semie" in input:
        args.append("--semie")
    elif "semimu" in input:
        args.append("--semimu")

    #args.append("--weight %.15f" % weight)

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

systs = ["JECup", "JECdown", "JERup", "JERdown"]

if True:
    systs += ["matchingup", "matchingdown", "scaleup", "scaledown"]

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

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "15"]
subprocess.call(args)

masses_narrow = {
        500: "5",
        750: "7p5",
        1000: "10",
        1250: "12p5",
        1500: "15",
        2000: "20"
        }

masses_large = {
        500: "50",
        750: "75",
        1000: "100",
        1250: "125",
        1500: "150",
        2000: "200"
        }

masses = [masses_narrow, masses_large]

myFinalList = []
for type in ["semie", "semimu"]:
    path = "skims/%s/%s/Systematics" % (d, type)
    for syst in systs:
        for widths in masses:
            for mass, width in widths.items():
                aList = []
                for file in files:
                    skim = os.path.join(path, file[0] % syst)
                    if "ZPrimeToTTJets" in skim and ("M"+str(mass)+"G") in skim and ("W"+str(width)+"G") in skim:
                        aList.append(skim)
                if len(aList) != 0:
                    myFinalList.append(aList)

for pair in myFinalList:
    merged_file = pair[0].replace(".root", "_merged.root")
    args = ["hadd", "-f", merged_file]
    for input in pair:
        if os.path.exists(input):
            args.append(input)

    if len(args) > 3:
        subprocess.call(args)
