#!/usr/bin/env python

from __future__ import division
import os, math, subprocess, shutil, sys, glob, stat, json
from optparse import OptionParser

if sys.version_info<(2,7,0):
  sys.stderr.write("You need python 2.7 or later to run this script\n")
  sys.exit(1)

masses = [
    750,
    1000,
    1250,
    1500
    ]

#parser = argparse.ArgumentParser(description='Submit jobs on the grid.')
#parser.add_argument('--python', dest='print_python', action='store_true')
#parser.add_argument('--b-tag', dest='btag', type=int, required=True)
#parser.add_argument('-i', dest='inputfile', required=True)
#args = parser.parse_args()

parser = OptionParser()
parser.add_option("--python", dest="print_python", action="store_true")
parser.add_option("--b-tag", dest="btag", type="int")
parser.add_option("-i", dest="inputfile")

(args, other) = parser.parse_args()

printPython = args.print_python

num_toys = 1000
num_jobs = 100

num_toys_per_job = int(math.ceil(num_toys / num_jobs))
num_toys = num_jobs * num_toys_per_job

import socket
isLXPLUS = False
if "lxplus" in socket.gethostname():
  isLXPLUS = True

if printPython:
  print "num_toys = %d" % num_toys
  print "num_jobs = %d" % num_jobs
  print "num_toys_per_job = %d" % num_toys_per_job
  exit(0)

if not args.btag:
  parser.error("No b-tag given")

if not args.inputfile:
  parser.error("No input file given")


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

data = args.inputfile

fit_configuration_files = [
    'fit_configuration/*'
    ]

input_files = [
  'analysis.json',
  'efficiencies.json',
  'sigma_reference.json',
  'systematics.json',
  'parameters.json',
  'pdf_parameters.json',
  'fit_pdf.json',
  '*.script']

input_files_mass_dependant = [
    "data_2012_nominal_%%(mass)d_%(name)s_%(btag)d_btag/data_2012_nominal_%%(mass)d_fitRes_%(name)s.root" % {'name': analysisName, 'btag': args.btag}
    ]

frit_files = [
    "frit/nominal-Zprime%%(mass)d_%(name)s_*_workspace.root" % {'name': analysisName}
    ]


works_files = [
    "../../../fitMtt",
    "../analysis.json"
    ]

print "Submitting %d jobs with %d toys per job: Total: %d toys for each mass" % (num_jobs, num_toys_per_job, num_toys)

working_dir = os.getcwd() + "/works/"
input_dir = os.getcwd() + "/inputs/"
output_dir = os.getcwd() + "/results/%d-btag/" % args.btag

wrapper = os.getcwd() + "/templates/wrapper.sh"
if isLXPLUS:
  wrapper = os.getcwd() + "/templates/wrapper_lxplus.sh"

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
  copy_file(data, os.path.join(input_dir, "data"))

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

def build_parameter(input_file, mass, index):
  global num_toys_per_job, input_dir, working_dir, output_dir, args
  return [input_file, mass, num_toys_per_job, index, working_dir, input_dir, output_dir, args.btag]

def create_ipnl_job():
  global working_dir
  j = Job(do_auto_resubmit=True)
  j.do_auto_resubmit = True
  j.application = Executable(exe=File(wrapper))

  if not isLXPLUS:
    j.backend = 'CREAM'
    j.backend.CE = 'lyogrid07.in2p3.fr:8443/cream-pbs-cms'
  else:
    j.backend = LSF()
    j.backend.queue = '8nh'

  return j

# Remove all jobs
# jobs.remove()

copy_data()
copy_fit_configuration_files()
copy_deps()
copy_works_files()

data_file_name = os.path.basename(data)
remote_data = os.path.join(input_dir, os.path.join("data", data_file_name))

for mass in masses:

  print "Submitting for M=%d" % mass
  copy_mass_deps(mass)

  # Create arguments
  parameters = []
  for i in range(0, num_jobs):
  #for i in range(0, 2):
    parameters.append(build_parameter(remote_data, mass, i))

  j = create_ipnl_job()
  j.name = "Zprime_%d" % mass

  s = ArgSplitter(args = parameters)
  j.splitter = s

  #j.submit(keep_going = True)
  j.submit()

print "All done!"
