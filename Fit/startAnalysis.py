#! /usr/bin/env python

import json, uuid, sys, os, datetime, shutil, glob
from pprint import pprint

import sys

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

def query(question, answsers):
    """Ask a question via raw_input() and return their answer.

    "question" is a string that is presented to the user.

    The "answer" return value is one of "answers".
    """
    prompt = " [%s] " % ("/".join(answsers))

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if choice in map(str.lower, answsers):
            return choice
        else:
            sys.stdout.write("Please respond with %s\n" % (", ".join(answsers)))

open("analysis.json", "a").close() # "touch" file

analysisName = raw_input("Enter analysis name: ")
analysisDescription = raw_input("Enter analysis description: ")
analysisType = query("What analysis do you want to do?", ["Higgs", "Zprime"])

useAsCurrentAnalysis = query_yes_no("Set this new analysis as the current one?")

analysisUUID = str(uuid.uuid1())

now = datetime.datetime.now()
date = now.strftime("%d%b%Y")

json_data = open("analysis.json")

try:
  data = json.load(json_data)
except:
  data = {}
  data["analysis"] = []
  data["current_analysis"] = -1

json_data.close()

new_analysis = {
    analysisUUID: {
      "name": analysisName,
      "description": analysisDescription,
      "date": date,
      "type": analysisType
      }
    }

data["analysis"].append(new_analysis)

if useAsCurrentAnalysis:
  data["current_analysis"] = len(data["analysis"]) - 1

with open('analysis.json', 'w') as f:
  json.dump(data, f, indent=2)

# Create folder structure
os.makedirs("analysis/%s" % analysisUUID)
os.mkdir("analysis/%s/frit" % analysisUUID)
shutil.copy("analysis.json", "analysis/%s" % analysisUUID)

shutil.copytree("fit_configuration", "analysis/%s/fit_configuration" % analysisUUID)
shutil.copytree("toys", "analysis/%s/toys" % analysisUUID, ignore = shutil.ignore_patterns('*.root'))

for filename in glob.glob('*.script'):
  shutil.copy(filename, "analysis/%s" % analysisUUID)

print("Added new analysis:")
print("\tID: %s" % analysisUUID)
print("\tName: %s" % analysisName)
print("\tDescription: %s" % analysisDescription)
print("\tDate: %s" % date)
print("\tCurrent analysis? %s" % useAsCurrentAnalysis)
