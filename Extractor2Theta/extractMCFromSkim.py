#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
        # Background
        ["MC_TT_powheg_theta_nominal.root", "skims/%s/MC_TT_powheg_skims_nominal.root", 21675970, 245.8],

        ["MC_T_tW-channel_theta_nominal.root", "skims/%s/MC_T_tW-channel_skims_nominal.root", 497658, 11.1],
        ["MC_T_s-channel_theta_nominal.root", "skims/%s/MC_T_s-channel_skims_nominal.root", 259961, 3.79],
        ["MC_T_t-channel_theta_nominal.root", "skims/%s/MC_T_t-channel_skims_nominal.root", 3728227, 56.4],

        ["MC_Tbar_tW-channel_theta_nominal.root", "skims/%s/MC_Tbar_tW-channel_skims_nominal.root", 493460, 11.1],
        ["MC_Tbar_s-channel_theta_nominal.root", "skims/%s/MC_Tbar_s-channel_skims_nominal.root", 139974, 1.76],
        ["MC_Tbar_t-channel_theta_nominal.root", "skims/%s/MC_Tbar_t-channel_skims_nominal.root", 1935072, 30.7],

        ["MC_DY1JetsToLL_M-50_theta_nominal.root", "skims/%s/MC_DY1JetsToLL_M-50_skims_nominal.root", 24045248, 666.3],
        ["MC_DY2JetsToLL_M-50_theta_nominal.root", "skims/%s/MC_DY2JetsToLL_M-50_skims_nominal.root", 21852156, 215.0],
        ["MC_DY3JetsToLL_M-50_theta_nominal.root", "skims/%s/MC_DY3JetsToLL_M-50_skims_nominal.root", 11015445, 60.7],
        ["MC_DY4JetsToLL_M-50_theta_nominal.root", "skims/%s/MC_DY4JetsToLL_M-50_skims_nominal.root", 6402827, 27.3],

        ["MC_W1JetsToLNu_theta_nominal.root", "skims/%s/MC_W1JetsToLNu_skims_nominal.root", 23141598, 6662.8],
        ["MC_W2JetsToLNu_theta_nominal.root", "skims/%s/MC_W2JetsToLNu_skims_nominal.root", 34044921, 2159.2],
        ["MC_W3JetsToLNu_theta_nominal.root", "skims/%s/MC_W3JetsToLNu_skims_nominal.root", 15539503, 640.4],
        ["MC_W4JetsToLNu_theta_nominal.root", "skims/%s/MC_W4JetsToLNu_skims_nominal.root", 13382803, 264.0],

        ["MC_WW_theta_nominal.root", "skims/%s/MC_WW_skims_nominal.root", 10000431, 56.0],
        ["MC_ZZ_theta_nominal.root", "skims/%s/MC_ZZ_skims_nominal.root", 9799908, 7.6],
        ["MC_WZ_theta_nominal.root", "skims/%s/MC_WZ_skims_nominal.root", 10000283, 33.6],

        ## QCD
        #["MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_theta_nominal.root", "skims/semimu/MC_QCD_pt15to30_bEnriched_MuEnrichedPt14_skims_nominal.root"],
        #["MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_theta_nominal.root", "skims/semimu/MC_QCD_pt30to50_bEnriched_MuEnrichedPt14_skims_nominal.root"],
        #["MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_theta_nominal.root", "skims/semimu/MC_QCD_pt50to150_bEnriched_MuEnrichedPt14_skims_nominal.root"],
        #["MC_QCD_pt150_bEnriched_MuEnrichedPt14_theta_nominal.root", "skims/semimu/MC_QCD_pt150_bEnriched_MuEnrichedPt14_skims_nominal.root"],

        #["MC_QCD_Pt_20_30_EMEnriched_theta_nominal.root", "skims/semie/MC_QCD_Pt_20_30_EMEnriched_skims_nominal.root"],
        #["MC_QCD_Pt_30_80_EMEnriched_theta_nominal.root", "skims/semie/MC_QCD_Pt_30_80_EMEnriched_skims_nominal.root"],
        #["MC_QCD_Pt_80_170_EMEnriched_theta_nominal.root", "skims/semie/MC_QCD_Pt_80_170_EMEnriched_skims_nominal.root"],
        #["MC_QCD_Pt_170_250_EMEnriched_theta_nominal.root", "skims/semie/MC_QCD_Pt_170_250_EMEnriched_skims_nominal.root"],
        #["MC_QCD_Pt_250_350_EMEnriched_theta_nominal.root", "skims/semie/MC_QCD_Pt_250_350_EMEnriched_skims_nominal.root"],
        #["MC_QCD_Pt_350_EMEnriched_theta_nominal.root", "skims/semie/MC_QCD_Pt_350_EMEnriched_skims_nominal.root"],

        # Background + Signal
        ["Signal_S0_S_i_M400_cpl1_scalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M400_cpl1_scalar_skims_nominal.root", 2074732, 0.5289 * 5.340],
        ["Signal_S0_S_i_M500_cpl1_scalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M500_cpl1_scalar_skims_nominal.root", 1990686, 0.3023 * 6.088],
        ["Signal_S0_S_i_M600_cpl1_scalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M600_cpl1_scalar_skims_nominal.root", 2354567, 0.1891 * 4.886],
        ["Signal_S0_S_i_M700_cpl1_scalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M700_cpl1_scalar_skims_nominal.root", 1999935, 0.1295 * 3.654],
        ["Signal_S0_S_i_M800_cpl1_scalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M800_cpl1_scalar_skims_nominal.root", 1999938, 0.09436 * 2.888],

        ["Signal_S0_S_i_M400_cpl1_pseudoscalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M400_cpl1_pseudoscalar_skims_nominal.root", 2276384, 1.169 * 5.039],
        ["Signal_S0_S_i_M500_cpl1_pseudoscalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M500_cpl1_pseudoscalar_skims_nominal.root", 1999940, 0.6756 * 3.803],
        ["Signal_S0_S_i_M600_cpl1_pseudoscalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M600_cpl1_pseudoscalar_skims_nominal.root", 1999949, 0.4414 * 2.709],
        ["Signal_S0_S_i_M700_cpl1_pseudoscalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M700_cpl1_pseudoscalar_skims_nominal.root", 1999932, 0.3121 * 2.057],
        ["Signal_S0_S_i_M800_cpl1_pseudoscalar_theta_nominal.root", "skims/%s/Signal_S0_S_i_M800_cpl1_pseudoscalar_skims_nominal.root", 1999940, 0.2325 * 1.697],

        ["Signal_ZPrimeToTTJets_M500GeV_W5GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M500GeV_W5GeV_skims_nominal_merged.root", 118975 + 119053, 1],
        ["Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M750GeV_W7p5GeV_skims_nominal_merged.root", 108827 + 108802, 1],
        ["Signal_ZPrimeToTTJets_M1000GeV_W10GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1000GeV_W10GeV_skims_nominal_merged.root", 103095 + 103751, 1],
        ["Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1250GeV_W12p5GeV_skims_nominal_merged.root", 97864 + 99558, 1],
        ["Signal_ZPrimeToTTJets_M1500GeV_W15GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1500GeV_W15GeV_skims_nominal_merged.root", 97349 + 388760, 1],
        ["Signal_ZPrimeToTTJets_M2000GeV_W20GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M2000GeV_W20GeV_skims_nominal_merged.root", 94817 + 94705, 1],

        ["Signal_ZPrimeToTTJets_M500GeV_W50GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M500GeV_W50GeV_skims_nominal_merged.root", 110507 + 119260],
        ["Signal_ZPrimeToTTJets_M750GeV_W75GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M750GeV_W75GeV_skims_nominal_merged.root",  106587 + 106587],
        ["Signal_ZPrimeToTTJets_M1000GeV_W100GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1000GeV_W100GeV_skims_nominal_merged.root", 104043 + 103387],
        ["Signal_ZPrimeToTTJets_M1250GeV_W125GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1250GeV_W125GeV_skims_nominal_merged.root", 100805 + 100615],
        ["Signal_ZPrimeToTTJets_M1500GeV_W150GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M1500GeV_W150GeV_skims_nominal_merged.root", 98775 + 97976],
        ["Signal_ZPrimeToTTJets_M2000GeV_W200GeV_theta_nominal.root", "skims/%s/Signal_ZPrimeToTTJets_M2000GeV_W200GeV_skims_nominal_merged.root", 97240 + 387897],

        # For scale and matching systematics
        ["MC_TT_madgraph_theta_nominal.root", "skims/%s/MC_TT_madgraph_skims_nominal.root", 6901778, 245.8],
        ]

def launch(input, output, weight):
  args = ["./extractor2Theta", "-i", input, "-o", output, "--mc", "--skim"]
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
    path = "theta/%s/%s" % (d, type)
    try:
        os.makedirs(path)
    except:
        pass

print("Extracting dataset ...")

for file in files:
    events = file[2] if len(file) >= 3 else -1
    xsection = file[3] if len(file) >= 4 else 1
    weight = 19667. * xsection / events if events > 0 else 1;
    print("Weight for %s: %.15f" % (file[1], weight))

    for type in ["semie", "semimu"]:
        if not "%" in file[1] and not type in file[1]:
            continue
        path = "theta/%s/%s" % (d, type)
        if not "%" in file[1]:
            tmpfile.write(launch(file[1], os.path.join(path, file[0]), weight) + "\n");
        else:
            tmpfile.write(launch(file[1] % type, os.path.join(path, file[0]), weight) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "30"] 
subprocess.call(args)
