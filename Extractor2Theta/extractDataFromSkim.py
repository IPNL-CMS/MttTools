#!/usr/bin/env python

import os, subprocess, datetime, tempfile
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

inputs = [
        ['MTT_MuHad_Run2012A-22Jan2013.root', 'skims/data/MTT_MuHad_Run2012A-22Jan2013.root', 'semimu'],
        ['MTT_SingleMu_Run2012B-TOPMuPlusJets-22Jan2013.root', 'skims/data/MTT_SingleMu_Run2012B-TOPMuPlusJets-22Jan2013.root', 'semimu'],
        ['MTT_SingleMu_Run2012C-TOPMuPlusJets-22Jan2013.root', 'skims/data/MTT_SingleMu_Run2012C-TOPMuPlusJets-22Jan2013.root', 'semimu'],
        ['MTT_SingleMu_Run2012D-TOPMuPlusJets-22Jan2013.root', 'skims/data/MTT_SingleMu_Run2012D-TOPMuPlusJets-22Jan2013.root', 'semimu'],

        ['MTT_ElectronHad_Run2012A-22Jan2013.root', 'skims/data/MTT_ElectronHad_Run2012A-22Jan2013.root', 'semie'],
        ['MTT_SingleElectron_Run2012B-TOPElePlusJets-22Jan2013.root', 'skims/data/MTT_SingleElectron_Run2012B-TOPElePlusJets-22Jan2013.root', 'semie'],
        ['MTT_SingleElectron_Run2012C-TOPElePlusJets-22Jan2013.root', 'skims/data/MTT_SingleElectron_Run2012C-TOPElePlusJets-22Jan2013.root', 'semie'],
        ['MTT_SingleElectron_Run2012D-TOPElePlusJets-22Jan2013.root', 'skims/data/MTT_SingleElectron_Run2012D-TOPElePlusJets-22Jan2013.root', 'semie'],
        ]

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

bdt_weights = "/gridgroup/cms/brochet/HTT/CMSSW_analysis/SL6/MttTools/SelectionMVA/BkgVsTT/bdt_trained/16Oct14/all-btag/weights/BDT_all-btag_BDT_boost_grad_0p2.weights.xml"

def launch(input, output, type):
  args = ["./extractor2Theta", "-i", input, "-o", output, "--data", "--type", type, "--skim", sortingAlgoArg, typeArg, "--bdt-weights", bdt_weights]
  
  return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

# Build output tree structure
path = "theta/%s/data" % (d)
try:
    os.makedirs(path)
except:
    pass

print("Extracting dataset ...")
for input in inputs:
  tmpfile.write(launch(input[1], os.path.join(path, input[0]), input[2]) + "\n")

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "20"] 
subprocess.call(args)

print("Merging ...")

for type in ["semie", "semimu"]:
    args = ["hadd", "-f", os.path.join(path, "MTT_Data_%s.root" % type)]
    args += ([os.path.join(path, output[0]) for output in inputs if type in output[2]])

    subprocess.call(args)


