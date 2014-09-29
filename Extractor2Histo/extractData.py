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
        #['MTT_SingleElectron_Run2012D-TOPElePlusJets-22Jan2013.root', 'lists/MTT_SingleElectron_Run2012D-TOPElePlusJets-22Jan2013_02Dec.list', 'semie'],
        ]

def launch(input, output, type, btag):
    args = ["./extractorToHisto", "--input-list", input, "-o", output, "--data", "--%s" % type, "--b-tag", str(btag)]

    return " ".join(args)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = False)

# Build output tree structure
for btag in [0, 1, 2]:
    path = "plots/%s/%d-btag/data" % (d, btag)
    try:
        os.makedirs(path)
    except:
        pass

print("Extracting dataset ...")
for input in inputs:
    for btag in [0, 1, 2]:
        path = "plots/%s/%d-btag/data" % (d, btag)
        tmpfile.write(launch(input[1], os.path.join(path, input[0]), input[2], btag) + "\n")

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "8"] 
subprocess.call(args)

## All is done, merge

#print("Merging ...")
#for btag in [1, 2]:
    #args = ["hadd", "MTT_Data_merged_2012_%d_btag_%s.root" % (btag, d)]
    #for output in inputs:
        #args.append(output[0])

#subprocess.call(args)

#for output in inputs:
#  os.remove(output[0])
