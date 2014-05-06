#! /usr/bin/env python

import subprocess
from optparse import OptionParser
parser = OptionParser()

parser.add_option("-d", "--dataset", dest="dataset")
parser.add_option("-o", "--output", dest="output")

datasets = {
        "/S0_S_i_M400_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M400_cpl1_pseudoscalar_01Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M400_cpl1_pseudoscalar",
        "/S0_S_i_M500_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M500_cpl1_pseudoscalar_03Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M500_cpl1_pseudoscalar",
        "/S0_S_i_M600_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M600_cpl1_pseudoscalar_02Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M600_cpl1_pseudoscalar",
        "/S0_S_i_M700_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M700_cpl1_pseudoscalar_02Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M700_cpl1_pseudoscalar",
        "/S0_S_i_M800_cpl1_pseudoscalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M800_cpl1_pseudoscalar_05Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M800_cpl1_pseudoscalar",

        "/S0_S_i_M400_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M400_cpl1_scalar_15Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M400_cpl1_scalar",
        "/S0_S_i_M500_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M500_cpl1_scalar_08Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M500_cpl1_scalar",
        "/S0_S_i_M600_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M600_cpl1_scalar_11Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M600_cpl1_scalar",
        "/S0_S_i_M700_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M700_cpl1_scalar_11Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M700_cpl1_scalar",
        "/S0_S_i_M800_cpl1_scalar_15Dec13_START53_V7C-GEN/sbrochet-S0_S_i_M800_cpl1_scalar_11Jan14_START53_V19-AOD-SIM-1af1161dfdea4f0d99d6c3f8450e776b/USER": "S0_S_i_M800_cpl1_scalar"
        }

for dataset, output in datasets.items():

    f = '%s_cross-section_ratio.json' % output

    args = ["./computeCrossSectionRatio.py", "-d", dataset, "-o", f]
    subprocess.call(args)
