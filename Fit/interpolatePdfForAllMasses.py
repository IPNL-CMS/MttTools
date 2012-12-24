#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os, math, subprocess, shutil, sys, tempfile

if sys.version_info<(2,7,0):
  sys.stderr.write("You need python 2.7 or later to run this script\n")
  sys.exit(1)

import argparse

parser = argparse.ArgumentParser(description='Generate PDFs interpolation for all masses.')
parser.add_argument("--b-tag", dest="btag", required=True, type=int)
args = parser.parse_args()

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = True)

for m in range(750, 1500, 50):
  if m == 750 or m == 1000 or m == 1250 or m == 1500:
    continue
  tmpfile.write("./interpolatePdf -m %d --b-tag %d\n" % (m, args.btag))

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 
subprocess.call(args)

tmpfile.close()
