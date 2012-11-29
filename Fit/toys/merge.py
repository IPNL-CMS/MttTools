#! /usr/bin/env python

import os, subprocess, sys, json, glob
from subprocess import Popen, PIPE, STDOUT

def test_root_file(root_file):

  root_code="""{
    TFile *f = TFile::Open("%s");
    if (! f) {
      return 1;
    } else {
      f->Close();
      delete f;
      return 0;
    }
  }""" % root_file

  f = open("root_file_test.C", "w")
  f.write(root_code)
  f.close()

  with open(os.devnull, "w") as fnull:
    returncode = subprocess.call(['root', '-b', '-q', '-l', '-n', 'root_file_test.C'], stdout = fnull, stderr = fnull)

  os.remove("root_file_test.C")

  print(returncode)

working_dir = os.getcwd() + "/works/"
input_dir = os.getcwd() + "/inputs/"
output_dir = sys.argv[1]

if not os.path.exists("../analysis.json"):
  print("No analysis currently defined. Please use startAnalysis first")
  sys.exit(1)

json_data = open("../analysis.json")

try:
  data = json.load(json_data)
except:
  print("No analysis currently defined. Please use startAnalysis first")
  sys.exit(1)

json_data.close()

current_analysis = data["current_analysis"]
uuid = data["analysis"][current_analysis].keys()[0]
analysisTuple = data["analysis"][current_analysis][uuid]
analysisName = analysisTuple["name"]

print("Merging...")
masses = [
    750,
    1000,
    1250,
    1500
    ]

for mass in masses:
  filenames = glob.glob("%s/data_2012_nominal_%d_toylimit_%s_*.root" % (output_dir, mass, analysisName))
  merged_filename = "%s/data_2012_nominal_%d_toylimit_%s.root" % (output_dir, mass, analysisName)

  if len(filenames) == 0:
    print("No files found for %d, continuing." % mass)
    continue

  filenames = ' '.join(x)

  with open(os.devnull, "w") as fnull:
    errorcode = subprocess.call(['hadd', merged_filename, filenames], stdout = fnull, stderr = fnull)

  if errorcode == 0 and os.path.exists(merged_filename):
    if test_root_file(merged_filename) == 0:
      # We can safely remove the files
      print("Would have deleted: %s" % filenames)
      print("Done.")

print("All done!")
