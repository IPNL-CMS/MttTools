#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os, math, subprocess, shutil, sys, tempfile, json

if sys.version_info<(2,7,0):
  sys.stderr.write("You need python 2.7 or later to run this script\n")
  sys.exit(1)

import argparse

parser = argparse.ArgumentParser(description='Compute JEC systematics for all masses with toys.')
parser.add_argument("--b-tag", dest="btag", required=True, type=int)
parser.add_argument("-i", dest="i", required=True, type=str)
args = parser.parse_args()

f = open("analysis.json")
params = json.load(f)
f.close()

current_analysis = params["current_analysis"]
analysisUUID = params["analysis"][current_analysis].keys()[0]

base_path = "analysis/%s" % analysisUUID

tmpfile = tempfile.NamedTemporaryFile(dir = '/tmp/', delete = False)

for m in range(750, 1501, 250):
  for i in range(0, 10):
    output = "%s/jec_systematics_%d_%d_btag_file_%d.root" % (base_path, m, args.btag, i)
    tmpfile.write("./computeJECSystWithToys.py -m %d --b-tag %d -i %s -o %s\n" % (m, args.btag, args.i, output))

tmpfile.flush()

cmds = ["parallel", "-u", "-a", tmpfile.name] 
subprocess.call(cmds)

tmpfile.close()

# Merge files
for m in range(750, 1501, 250):
  files_to_merge = []
  for i in range(0, 10):
    files_to_merge.append("%s/jec_systematics_%d_%d_btag_file_%d.root" % (base_path, m, args.btag, i))

  output = "%s/jec_systematics_%d_%d_btag.root" % (base_path, m, args.btag)
  cmds = ["hadd", "-f", output]
  cmds.extend(files_to_merge)
  subprocess.call(cmds)
