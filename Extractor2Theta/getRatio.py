#! /usr/bin/env python

import subprocess, json
from optparse import OptionParser
parser = OptionParser()

(options, args) = parser.parse_args()

js = {}
with open(args[0]) as f:
    js = json.load(f)

print "Ratio: %.3f" % (js["effective_entries"] / js["integral"])
print "Ratio (GetEntries): %.3f" % (js["entries"] / js["integral"])
