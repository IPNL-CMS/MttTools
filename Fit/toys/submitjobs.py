#!/usr/bin/env python

from __future__ import division
import os, math, subprocess, shutil, sys, glob, stat

if sys.version_info<(2,7,0):
  sys.stderr.write("You need python 2.7 or later to run this script\n")
  sys.exit(1)

import json, argparse

masses = [
    750,
    1000,
    1250,
    1500
    ]

parser = argparse.ArgumentParser(description='Submit jobs on the grid.')
parser.add_argument('--python', dest='print_python', action='store_true')
parser.add_argument('--b-tag', dest='btag', type=int, required=True)
args = parser.parse_args()

printPython = args.print_python

num_toys = 1000
num_jobs = 100

num_toys_per_job = int(math.ceil(num_toys / num_jobs))
num_toys = num_jobs * num_toys_per_job

if printPython:
  print "num_toys = %d" % num_toys
  print "num_jobs = %d" % num_jobs
  print "num_toys_per_job = %d" % num_toys_per_job
  exit(0)

f = open("../parameters.json")
params = json.load(f)
pdfSignalName = params["parameters"]["pdf"]["signal"]
f.close()

data = [
  'data/ds_data_2011-nominal_*.txt']

fit_configuration_files = [
    'fit_configuration/*'
    ]

input_files = [
  'efficiencies.json',
  'sigma_reference.json',
  'systematics.json',
  'parameters.json',
  'pdf_parameters.json',
  'fit_pdf.json',
  '*.script']

input_files_mass_dependant = [
    "data_2011_nominal_%%(mass)d_%(signalPdf)s_%(btag)d_btag/data_2011_nominal_%%(mass)d_fitRes_%(signalPdf)s.root" % {'signalPdf': pdfSignalName, 'btag': args.btag}
    ]

frit_files = [
    "frit/nominal-Zprime%%(mass)d_%(signalPdf)s_*_workspace.root" % {'signalPdf': pdfSignalName}
    ]


works_files = [
    "../fitMtt",
    "templates/wrapper.sh"
    ]

print "Submitting %d jobs with %d toys per job: Total: %d toys for each mass" % (num_jobs, num_toys_per_job, num_toys)

working_dir = os.getcwd() + "/works/"
input_dir = os.getcwd() + "/inputs/"
output_dir = os.getcwd() + "/results/%d-btag/" % args.btag
wrapper_dir = os.getcwd() + "/templates/wrapper.sh"

os.path.exists(os.path.join(input_dir, "data")) or os.mkdir(os.path.join(input_dir, "data"))
os.path.exists(os.path.join(input_dir, "frit")) or os.mkdir(os.path.join(input_dir, "frit"))
os.path.exists(os.path.join(input_dir, "fit_configuration")) or os.mkdir(os.path.join(input_dir, "fit_configuration"))
os.path.exists(output_dir) or os.mkdir(output_dir)

os.chmod(output_dir, 0777)

def copy_file(files, dest):
  for file in glob.glob(files):
    if not os.path.isfile(file):
      print "ERROR: '%s' not found." % file
      sys.exit(1)
    print("Copy %s to %s" % (file, dest))
    shutil.copy(file, dest)

def copy_data():
  global input_dir
  for input_file in data:
    correct_file = "../%s" % input_file
    copy_file(correct_file, os.path.join(input_dir, "data"))

def copy_fit_configuration_files():
  global input_dir
  for input_file in fit_configuration_files:
    correct_file = "../%s" % input_file
    copy_file(correct_file, os.path.join(input_dir, "fit_configuration"))

def copy_deps():
  global input_files, input_dir
  for input_file in input_files:
    correct_file = "../%s" % input_file
    copy_file(correct_file, input_dir)

def copy_mass_deps(mass):
  global input_files, input_dir
  for input_file in input_files_mass_dependant:
    correct_file = "../%s" % (input_file % {'mass': mass})
    copy_file(correct_file, input_dir)

  for input_file in frit_files:
    correct_file = "../%s" % (input_file % {'mass': mass})
    copy_file(correct_file, os.path.join(input_dir, "frit"))

def copy_works_files():
  global works_files, working_dir
  for file in works_files:
    copy_file(file, working_dir)

def build_parameter(mass, index):
  global num_toys_per_job, input_dir, working_dir, output_dir, args
  return [mass, num_toys_per_job, index, working_dir, input_dir, output_dir, args.btag]

def create_ipnl_job():
  global working_dir
  j = Job()
  j.application = Executable(exe=File(wrapper_dir))
  j.backend = 'CREAM'
  j.backend.CE = 'lyogrid07.in2p3.fr:8443/cream-pbs-cms'
  return j

# Remove all jobs
# jobs.remove()

copy_data()
copy_fit_configuration_files()
copy_deps()
copy_works_files()

for mass in masses:

  print "Submitting for M=%d" % mass
  copy_mass_deps(mass)

  # Create arguments
  parameters = []
  for i in range(0, num_jobs):
  #for i in range(0, 2):
    parameters.append(build_parameter(mass, i))

  j = create_ipnl_job()
  j.name = "Zprime_%d" % mass

  s = ArgSplitter(args = parameters)
  j.splitter = s

  j.submit()

print "All done!"
