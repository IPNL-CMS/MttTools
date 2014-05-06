#! /usr/bin/env python

import subprocess
from optparse import OptionParser
parser = OptionParser()

parser.add_option("-d", "--dataset", dest="dataset")
parser.add_option("-o", "--output", dest="output")

(options, args) = parser.parse_args()

# Get list of input file
args = ["/gridgroup/cms/brochet/bin/das_client.py", "--limit", "30", "--format=json", "--query=instance=prod/phys02 file dataset=%s" % options.dataset]
output = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0]

import json

files = []
js = json.loads(output)
for f in js["data"]:
    files.append('root://lyogrid06.in2p3.fr//dpm/in2p3.fr/home/cms/data%s' % f["file"][0]["name"])

# Load CMSSW libraries
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable();
ROOT.gSystem.Load("libDataFormatsFWLite.so")
from DataFormats.FWLite import Events, Handle

genHandle = Handle('GenEventInfoProduct')
genLabel = ('generator')

genParticlesHandle = Handle('std::vector<reco::GenParticle>')
genParticlesLabel = ('genParticles')

import multiprocessing

def processFile(f):
    print "Opening '%s'" % f

    mtts = []
    weights = []

    events = Events(str(f))
    for event in events:
        event.getByLabel(genLabel, genHandle)
        event.getByLabel(genParticlesLabel, genParticlesHandle)

        gen = genHandle.product()
        genParticles = genParticlesHandle.product()

        weight = gen.weights()[0]

        mtt_p4 = ROOT.TLorentzVector(0, 0, 0, 0)
        for i in range(0, genParticles.size()):
            p = genParticles[i]
            if p.status() != 3:
                continue

            if abs(p.pdgId()) != 6:
                continue

            p4 = ROOT.TLorentzVector(p.px(), p.py(), p.pz(), p.energy())
            mtt_p4 += p4

        mtts.append(mtt_p4.M())
        weights.append(weight)

    return mtts, weights

pool = multiprocessing.Pool(30)
mtt_and_weight = pool.map(processFile, files)

mtts, weights = zip(*mtt_and_weight)

mtt = ROOT.TH1F("mtt", "mtt", 500, 250, 1250)
for i in range(0, len(mtts)):
    for j in range(0, len(mtts[i])):
        mtt.Fill(mtts[i][j], weights[i][j])

nEntries = 0
for bin in range(1, mtt.GetNbinsX() + 1):
    nEntries += abs(mtt.GetBinContent(bin))

js = {
        "entries": mtt.GetEntries(),
        "integral": mtt.Integral(),
        "effective_entries": nEntries
      }

with open(options.output, 'w') as f:
    json.dump(js, f)
    print "Output saved as '%s'" % options.output

f = ROOT.TFile.Open("test.root", "recreate")
mtt.Write()
f.Close()
