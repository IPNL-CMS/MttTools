#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime
from optparse import OptionParser

d = datetime.datetime.now().strftime("%d%b%y")

parser = OptionParser()
parser.add_option("", "--mva", action="store_true", dest="mva", default=False, help="Use MVA sorting algorithm")
parser.add_option("", "--chi2", action="store_true", dest="chi2", default=False, help="Use Chi2 sorting algorithm")
parser.add_option("", "--kf", action="store_true", dest="kf", default=False, help="Use KF sorting algorithm")
parser.add_option("", "--hybrid", action="store_true", dest="hybrid", default=False, help="Use hybrid sorting algorithm")
parser.add_option("", "--zprime", action="store_true", dest="zprime", default=False, help="Change mtt range for Zprime analysis")
parser.add_option("", "--higgs", action="store_true", dest="higgs", default=False, help="Change mtt range for Higgs analysis")

(option, args) = parser.parse_args()

files = [
        # Background + Signal
        ["Signal_S0_S_i_M400_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M400_cpl1_scalar_skims_%s.root", 2074732, 0.5289 * 5.340],
        ["Signal_S0_S_i_M500_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M500_cpl1_scalar_skims_%s.root", 1990686, 0.3023 * 6.088],
        ["Signal_S0_S_i_M600_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M600_cpl1_scalar_skims_%s.root", 2354567, 0.1891 * 4.886],
        ["Signal_S0_S_i_M700_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M700_cpl1_scalar_skims_%s.root", 1999935, 0.1295 * 3.654],
        ["Signal_S0_S_i_M800_cpl1_scalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M800_cpl1_scalar_skims_%s.root", 1999938, 0.09436 * 2.888],

        ["Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M400_cpl1_pseudoscalar_skims_%s.root", 2276384, 1.169 * 5.039],
        ["Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M500_cpl1_pseudoscalar_skims_%s.root", 1999940, 0.6756 * 3.803],
        ["Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M600_cpl1_pseudoscalar_skims_%s.root", 1999949, 0.4414 * 2.709],
        ["Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M700_cpl1_pseudoscalar_skims_%s.root", 1999932, 0.3121 * 2.057],
        ["Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_%s.root", "skims/%s/Systematics/Signal_S0_S_i_M800_cpl1_pseudoscalar_skims_%s.root", 1999940, 0.2325 * 1.697],

        ["Signal_ZPrimeToTTJets_M500GeV_W5GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M500GeV_W5GeV_skims_%s_merged.root", 118975 + 119053, 1],
        ["Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_skims_%s_merged.root", 108827 + 108802, 1],
        ["Signal_ZPrimeToTTJets_M1000GeV_W10GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1000GeV_W10GeV_skims_%s_merged.root", 103095 + 103751, 1],
        ["Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_skims_%s_merged.root", 97864 + 99558, 1],
        ["Signal_ZPrimeToTTJets_M1500GeV_W15GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1500GeV_W15GeV_skims_%s_merged.root", 97349 + 388760, 1],
        ["Signal_ZPrimeToTTJets_M2000GeV_W20GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M2000GeV_W20GeV_skims_%s_merged.root", 94817 + 94705, 1],

        ["Signal_ZPrimeToTTJets_M500GeV_W50GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M500GeV_W50GeV_skims_%s_merged.root", 110507 + 119260],
        ["Signal_ZPrimeToTTJets_M750GeV_W75GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M750GeV_W75GeV_skims_%s_merged.root",  106587 + 106587],
        ["Signal_ZPrimeToTTJets_M1000GeV_W100GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1000GeV_W100GeV_skims_%s_merged.root", 104043 + 103387],
        ["Signal_ZPrimeToTTJets_M1250GeV_W125GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1250GeV_W125GeV_skims_%s_merged.root", 100805 + 100615],
        ["Signal_ZPrimeToTTJets_M1500GeV_W150GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M1500GeV_W150GeV_skims_%s_merged.root", 98775 + 97976],
        ["Signal_ZPrimeToTTJets_M2000GeV_W200GeV_theta_%s.root", "skims/%s/Systematics/Signal_ZPrimeToTTJets_M2000GeV_W200GeV_skims_%s_merged.root", 97240 + 387897],

        # Background
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

        ["MC_WW_theta_%s.root", "skims/%s/Systematics/MC_WW_skims_%s.root", 10000431, 56.0],
        ["MC_ZZ_theta_%s.root", "skims/%s/Systematics/MC_ZZ_skims_%s.root", 9799908, 7.6],
        ["MC_WZ_theta_%s.root", "skims/%s/Systematics/MC_WZ_skims_%s.root", 10000283, 33.6],

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


generator_syst_files = [

        # Special systematics
        ["MC_TT_madgraph_theta_matchingup.root", "skims/%s/Systematics/MC_TT_madgraph_skims_matchingup.root", 5415010, 245.8],
        ["MC_TT_madgraph_theta_matchingdown.root", "skims/%s/Systematics/MC_TT_madgraph_skims_matchingdown.root", 5476728, 245.8],

        ["MC_TT_madgraph_theta_scaleup.root", "skims/%s/Systematics/MC_TT_madgraph_skims_scaleup.root", 5009488, 245.8],
        ["MC_TT_madgraph_theta_scaledown.root", "skims/%s/Systematics/MC_TT_madgraph_skims_scaledown.root", 5387181, 245.8],

        ["MC_WJetsToLNu_theta_scaleup.root", "skims/%s/Systematics/MC_WJetsToLNu_skims_scaleup.root", 20784770, 9726.4],
        ["MC_WJetsToLNu_theta_scaledown.root", "skims/%s/Systematics/MC_WJetsToLNu_skims_scaledown.root", 20760884, 9726.4],
        ["MC_WJetsToLNu_theta_matchingup.root", "skims/%s/Systematics/MC_WJetsToLNu_skims_matchingup.root", 20976082, 9726.4],
        ["MC_WJetsToLNu_theta_matchingdown.root", "skims/%s/Systematics/MC_WJetsToLNu_skims_matchingdown.root", 21364637, 9726.4],

        ["MC_DYJetsToLL_M-50_theta_scaleup.root", "skims/%s/Systematics/MC_DYJetsToLL_M-50_skims_scaleup.root", 2170270, 969.3],
        ["MC_DYJetsToLL_M-50_theta_scaledown.root", "skims/%s/Systematics/MC_DYJetsToLL_M-50_skims_scaledown.root", 1934901, 969.3],
        ["MC_DYJetsToLL_M-50_theta_matchingup.root", "skims/%s/Systematics/MC_DYJetsToLL_M-50_skims_matchingup.root", 1985529, 969.3],
        ["MC_DYJetsToLL_M-50_theta_matchingdown.root", "skims/%s/Systematics/MC_DYJetsToLL_M-50_skims_matchingdown.root", 2112387, 969.3],

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

if True:
    systs["pdfOnlyUp"] = ["nominal", "--pdf-syst up"]
    systs["pdfOnlyDown"] = ["nominal", "--pdf-syst down"]

if True:
    systs["alphasOnlyUp"] = ["nominal", "--alphas-syst up"]
    systs["alphasOnlyDown"] = ["nominal", "--alphas-syst down"]


sortingAlgoArg = ""
if option.mva:
    sortingAlgoArg = "--mva"
elif option.kf:
    sortingAlgoArg = "--kf"
elif option.chi2:
    sortingAlgoArg = "--chi2"
elif option.hybrid:
    sortingAlgoArg = "--hybrid"

typeArg = ""
if option.zprime:
    typeArg = "--zprime"
elif option.higgs:
    typeArg = "--higgs"

def launch(input, output, weight, extra):
    args = ["./extractor2Theta", "-i", input, "-o", output, "--mc", "--skim", sortingAlgoArg, typeArg]

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
            xsection = file[3] if len(file) >= 4 else 1
            weight = 19667. * xsection / events if events > 0 else 1;
            print("Weight for %s: %.15f" % (file[1], weight))

            tmpfile.write(launch(inputFile, os.path.join(path, file[0] % syst), weight, extra[1]) + "\n");

for file in generator_syst_files:
    for type in ["semie", "semimu"]:
        path = "theta/%s/%s/Systematics" % (d, type)
        inputFile = file[1] % type

        events = file[2] if len(file) >= 3 else -1
        xsection = file[3] if len(file) >= 4 else -1
        weight = 19667. * xsection / events if events > 0 else 1;
        print("Weight for %s: %.15f" % (file[1], weight))

        tmpfile.write(launch(inputFile, os.path.join(path, file[0]), weight, "") + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "20"] 
#print args
subprocess.call(args)
