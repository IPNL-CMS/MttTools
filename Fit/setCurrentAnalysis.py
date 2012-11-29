#! /usr/bin/env python

import json, uuid, sys, os, datetime
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


if not os.path.exists("analysis.json"):
  print("No analysis currently defined. Please use startAnalysis first")
  sys.exit(1)

json_data = open("analysis.json")

try:
  data = json.load(json_data)
except:
  print("No analysis currently defined. Please use startAnalysis first")
  sys.exit(1)

json_data.close()

index = 0
for analysis in data["analysis"]:
  uuid = analysis.keys()[0]
  analysisTuple = analysis[uuid]
  current = " - current" if index == data["current_analysis"] else ""
  print("[%d%s] %s - %s (%s)" % (index, current, analysisTuple["name"], analysisTuple["date"], uuid))
  print("\t%s" % analysisTuple["description"])
  index += 1

while True:
  i = raw_input("Enter new current analysis index [0-%d]: " % (index - 1))
  if not i.isdigit():
    continue

  i = int(i)
  if i >= index:
    print("Please use an index in the range 0 - %d" % (index -1))
    continue

  data["current_analysis"] = i
  break

with open('analysis.json', 'w') as f:
  json.dump(data, f, indent=2)
