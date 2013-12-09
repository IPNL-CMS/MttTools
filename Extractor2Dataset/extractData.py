#!/usr/bin/env python

import os, subprocess, datetime, tempfile

d = datetime.datetime.now().strftime("%d%b%y")

inputs = [
        ['MTT_MuHad_Run2012A-22Jan2013.root', 'lists/MTT_MuHad_Run2012A-22Jan2013_02Dec.list', 'semimu'],
        ['MTT_SingleMu_Run2012B-TOPMuPlusJets-22Jan2013.root', 'lists/MTT_SingleMu_Run2012B-TOPMuPlusJets-22Jan2013_02Dec.list', 'semimu'],
        ['MTT_SingleMu_Run2012C-TOPMuPlusJets-22Jan2013.root', 'lists/MTT_SingleMu_Run2012C-TOPMuPlusJets-22Jan2013_02Dec.list', 'semimu'],
        ['MTT_SingleMu_Run2012D-TOPMuPlusJets-22Jan2013.root', 'lists/MTT_SingleMu_Run2012D-TOPMuPlusJets-22Jan2013_02Dec.list', 'semimu'],

        ['MTT_ElectronHad_Run2012A-22Jan2013.root', 'lists/MTT_ElectronHad_Run2012A-22Jan2013_02Dec.list', 'semie'],
        ['MTT_SingleElectron_Run2012B-TOPElePlusJets-22Jan2013.root', 'lists/MTT_SingleElectron_Run2012B-TOPElePlusJets-22Jan2013_02Dec.list', 'semie'],
        ['MTT_SingleElectron_Run2012C-TOPElePlusJets-22Jan2013.root', 'lists/MTT_SingleElectron_Run2012C-TOPElePlusJets-22Jan2013_02Dec.list', 'semie'],
        ['MTT_SingleElectron_Run2012D-TOPElePlusJets-22Jan2013.root', 'lists/MTT_SingleElectron_Run2012D-TOPElePlusJets-22Jan2013_02Dec.list', 'semie'],
        ]

def launch(input, output, type):
  args = ["./extractor2Dataset", "--input-list", input, "-o", output, "--data", "--type", type]
  
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
args = ["hadd", "-f", os.path.join(path, "MTT_Data_merged_2012.root")]
for output in inputs:
  args.append(os.path.join(path, output[0]))

subprocess.call(args)
