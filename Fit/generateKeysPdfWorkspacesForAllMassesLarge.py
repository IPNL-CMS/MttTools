#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os, math, subprocess, shutil, sys, tempfile

if sys.version_info<(2,6,0):
  sys.stderr.write("You need python 2.6 or later to run this script\n")
  sys.exit(1)

from optparse import OptionParser

parser = OptionParser()
parser.add_option("", "--b-tag", dest="btag", type = int)
parser.add_option("-i", "--input", dest="i", type = str)
parser.add_option("-s", "--signal", dest="signal", type = str)

(args, pos) = parser.parse_args()

if not "%d" in args.signal or not "%s" in args.signal:
    print("Error: you must have a '%d' (mass) and a '%s' (width) in your signal filename")
    exit(1)

tmpfile = tempfile.NamedTemporaryFile(dir = '/scratch/', delete = True)

masses = {500: "50", 750: "75", 1000: "100", 1250: "125", 1500: "150", 2000: "200"}
for m, w in masses.items():
  tmpfile.write("./generateKeysPdfWorkspaces -m %d --b-tag %d -i %s --signal %s\n" % (m, args.btag, args.i, args.signal % (m, w)))

tmpfile.flush()

args = ["parallel", "-u", "-a", tmpfile.name, "-j", "6"] 
subprocess.call(args)

tmpfile.close()
