#!/usr/bin/env python

from __future__ import division
import os, subprocess, tempfile, datetime

d = datetime.datetime.now().strftime("%d%b%y")

files = [
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
    args = ["./preSkim", "--input-list", input, "-o", output, "--data"]
    if "semie" in type:
        args.append("--semie")
    elif "semimu" in type:
        args.append("--semimu")

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

# Build output tree structure
path = "skims/%s/data" % (d)
try:
    os.makedirs(path)
except:
    pass

print("Skimming datasets...")

for file in files:
    tmpfile.write(launch(file[1], os.path.join(path, file[0]), file[2]) + "\n");

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "8"] 
subprocess.call(args)
