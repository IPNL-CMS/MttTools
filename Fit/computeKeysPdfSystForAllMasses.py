#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os, math, subprocess, shutil, sys, tempfile

if sys.version_info<(2,7,0):
  sys.stderr.write("You need python 2.7 or later to run this script\n")
  sys.exit(1)

import argparse

parser = argparse.ArgumentParser(description='Compute KeysPdf systematics for all masses.')
parser.add_argument("--fix-background", dest="fix_background", required=False, action='store_true')
parser.add_argument("--b-tag", dest="btag", required=True, type=int)
parser.add_argument("-i", dest="i", required=True, type=str)
parser.add_argument("--signal", dest="signal", required=True, type=str)
args = parser.parse_args()

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = True)

for m in range(750, 1501, 250):
  tmpfile.write("./computeKeysPdfSyst -m %d --b-tag %d -i %s --signal %s%s\n" % (m, args.btag, args.i, args.signal % m, " --fix-background" if args.fix_background else ""))

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 
subprocess.call(args)

tmpfile.close()
