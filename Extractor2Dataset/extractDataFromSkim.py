#!/usr/bin/env python

import os, subprocess, datetime, tempfile
from optparse import OptionParser

d = datetime.datetime.now().strftime("%d%b%y")

parser = OptionParser()
parser.add_option("", "--mva", action="store_true", dest="mva", default=False, help="Use MVA sorting algorithm")
parser.add_option("", "--chi2", action="store_true", dest="chi2", default=False, help="Use Chi2 sorting algorithm")
parser.add_option("", "--kf", action="store_true", dest="kf", default=False, help="Use KF sorting algorithm")
parser.add_option("", "--hybrid", action="store_true", dest="hybrid", default=False, help="Use hybrid sorting algorithm")
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

def launch(input, output, type):
  args = ["./extractor2Dataset", "-i", input, "-o", output, "--data", "--type", type, "--skim", sortingAlgoArg]
  
  return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

# Build output tree structure
path = "datasets/%s/data" % (d)
try:
    os.makedirs(path)
except:
    pass

print("Extracting dataset ...")
for input in inputs:
  tmpfile.write(launch(input[1], os.path.join(path, input[0]), input[2]) + "\n")

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "8"] 
subprocess.call(args)

# All is done, merge

print("Merging ...")
args = ["hadd", "-f", os.path.join(path, "MTT_Data_merged.root")]
for output in inputs:
  args.append(os.path.join(path, output[0]))

subprocess.call(args)
