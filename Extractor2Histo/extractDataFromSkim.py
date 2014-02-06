#!/usr/bin/env python

import os, subprocess, datetime, tempfile

d = datetime.datetime.now().strftime("%d%b%y")

inputs = [
        ['Data_MuHad_Run2012A-22Jan2013.root', 'skims/data/MTT_MuHad_Run2012A-22Jan2013.root', 'semimu'],
        ['Data_SingleMu_Run2012B-TOPMuPlusJets-22Jan2013.root', 'skims/data/MTT_SingleMu_Run2012B-TOPMuPlusJets-22Jan2013.root', 'semimu'],
        ['Data_SingleMu_Run2012C-TOPMuPlusJets-22Jan2013.root', 'skims/data/MTT_SingleMu_Run2012C-TOPMuPlusJets-22Jan2013.root', 'semimu'],
        ['Data_SingleMu_Run2012D-TOPMuPlusJets-22Jan2013.root', 'skims/data/MTT_SingleMu_Run2012D-TOPMuPlusJets-22Jan2013.root', 'semimu'],

        ['Data_ElectronHad_Run2012A-22Jan2013.root', 'skims/data/MTT_ElectronHad_Run2012A-22Jan2013.root', 'semie'],
        ['Data_SingleElectron_Run2012B-TOPElePlusJets-22Jan2013.root', 'skims/data/MTT_SingleElectron_Run2012B-TOPElePlusJets-22Jan2013.root', 'semie'],
        ['Data_SingleElectron_Run2012C-TOPElePlusJets-22Jan2013.root', 'skims/data/MTT_SingleElectron_Run2012C-TOPElePlusJets-22Jan2013.root', 'semie'],
        ['Data_SingleElectron_Run2012D-TOPElePlusJets-22Jan2013.root', 'skims/data/MTT_SingleElectron_Run2012D-TOPElePlusJets-22Jan2013.root', 'semie'],
        ]

def launch(input, output, type, btag):
    args = ["./extractorToHisto", "-i", input, "-o", output, "--data", "--skim", "--%s" % type, "--b-tag", str(btag)]

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

# Build output tree structure
for btag in [1, 2]:
    path = "plots/%s/%d-btag/data" % (d, btag)
    try:
        os.makedirs(path)
    except:
        pass

print("Extracting dataset ...")
for input in inputs:
    for btag in [1, 2]:
        path = "plots/%s/%d-btag/data" % (d, btag)
        tmpfile.write(launch(input[1], os.path.join(path, input[0]), input[2], btag) + "\n")

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "4"] 
subprocess.call(args)

## All is done, merge

print("Merging ...")
for btag in [1, 2]:
    path = "plots/%s/%d-btag/data" % (d, btag)
    for type in ["semie", "semimu"]:
        args = ["hadd", "-f", os.path.join(path, "Data_SingleMu.root" if type == "semimu" else "Data_SingleElectron.root")]
        for output in inputs:
            if type == output[2]:
                args.append(os.path.join(path, output[0]))

        subprocess.call(args)
