#! /usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import division
import os, math, subprocess, shutil, sys, tempfile

if sys.version_info<(2,7,0):
  sys.stderr.write("You need python 2.7 or later to run this script\n")
  sys.exit(1)

import json, datetime, argparse
from string import Template

parser = argparse.ArgumentParser(description='Bookmark analysis.')
parser.add_argument('--intro', dest='do_intro', action='store_true')
parser.add_argument("--b-tag", dest="btag", required=True, type=int)
args = parser.parse_args()

doIntro = args.do_intro

f = open("analysis.json")
params = json.load(f)
f.close()

current_analysis = params["current_analysis"]
analysisUUID = params["analysis"][current_analysis].keys()[0]
analysisTuple = params["analysis"][current_analysis][analysisUUID]

analysisName = analysisTuple["name"]
analysisDescription = analysisTuple["description"]
analysisDate = analysisTuple["date"]

useSystematics = analysisTuple["systematics"]
useJECSyst = analysisTuple["jec_syst"]
useSignalSyst = analysisTuple["signal_syst"]
useBkgSyst = analysisTuple["background_syst"]

masses = [750, 1000, 1250, 1500]
jecs = ["nominal", "JECup", "JECdown"]
pwd = os.getcwd() + "/analysis/" + analysisUUID

btag = str(args.btag)

tmp = tempfile.mkdtemp(dir="/scratch")

PARAMS = {
    'pwd': pwd,
    'analysis_name': analysisName,
    'btag': btag
}

def run(program, *args):
  pid = os.fork()
  if not pid:
    os.execvp(program, (program,) + args)
  return os.wait()[1]

#if doIntro:
  #intro_file = tempfile.mkstemp(dir="/scratch/")
  #run("vim", intro_file[1])
  #shutil.copy(intro_file[1], tmp + "/intro.tex")
#else:
  #open(tmp + "/intro.tex", "w").close()

intro = open(tmp + "/intro.tex", "w")
template = Template(r"""\begin{itemize}
    \item Analysis UUID: ${uuid}
    \item Analysis name: ${name}
    \item Analysis description: ${description}
    \item Analysis date: {$date}
\end{itemize}

Analysis done with ${btag} b-tag.""")

btag_str = ""
minBTag = maxBTag = 0
if args.btag <= 2:
  btag_str = str(btag)
  minBTag = maxBTag = int(btag)
else:
  if args.btag == 3:
    minBTag = 1
    maxBTag = 2
    btag_str = "1 + 2"

intro.write(template.substitute(uuid = analysisUUID, name = analysisName.replace("_", "\_"), description = analysisDescription, date = analysisDate, btag = btag_str))

intro.close()

if not useSystematics:
  f = open(tmp + "/intro.tex", "w+")
  f.write(r"\begin{center}\textcolor{red}{WARNING: Analysis ran without systematics!}\end{center}")
  f.close()

# First, frit

for i in range(minBTag, maxBTag + 1):
  template_full = Template(r"""\subsubsection{${btag} b-tag}
  \begin{minipage}{0.33\textwidth} \centering
  \includegraphics[width=0.99\textwidth]{${pwd}/frit/nominal-Zprime${mass}_${analysis_name}_${btag}_btag_fit.pdf}\\
  Nominal\\
  $$\chi^2 = ${chi2_nominal}$$
  \end{minipage}%
  \begin{minipage}{0.33\textwidth} \centering
  \includegraphics[width=0.99\textwidth]{${pwd}/frit/JECup-Zprime${mass}_${analysis_name}_${btag}_btag_fit.pdf}\\
  JEC Up\\
  $$\chi^2 = ${chi2_JECup}$$
  \end{minipage}%
  \begin{minipage}{0.33\textwidth} \centering
  \includegraphics[width=0.99\textwidth]{${pwd}/frit/JECdown-Zprime${mass}_${analysis_name}_${btag}_btag_fit.pdf}\\
  JEC Down\\
  $$\chi^2 = ${chi2_JECdown}$$
  \end{minipage}""")

  template_reduced = Template(r"""\subsubsection{${btag} b-tag}
  \begin{center}\begin{minipage}{0.50\textwidth} \centering
  \includegraphics[width=0.99\textwidth]{${pwd}/frit/nominal-Zprime${mass}_${analysis_name}_${btag}_btag_fit.pdf}\\
  Nominal\\
  $$\chi^2 = ${chi2_nominal}$$
  \end{minipage}\end{center}""")
  
  for mass in masses:
    jsonFile = open(pwd + "/frit_efficiencies.json")
    chi2 = {}
    jsonValues = json.load(jsonFile)
    jsonFile.close()

    for jec in jecs:
      if jec in jsonValues[analysisUUID][str(mass)][str(i)]:
        chi2["chi2_" + jec] = jsonValues[analysisUUID][str(mass)][str(i)][jec]["chi2"]
      else:
        chi2["chi2_" + jec] = "Not computed"

    chi2.update(PARAMS)
    f = open(tmp + "/frit_%d.tex" % mass, "a+")

    if useJECSyst:
      f.write(template_full.substitute(chi2, btag = i, mass = mass))
    else:
      f.write(template_reduced.substitute(chi2, btag = i, mass = mass))

    f.close()

  shutil.copy(pwd + ("/efficiencies_table_%s_%s_btag.tex" % (analysisName, str(i))), tmp + "/efficiencies_table_%dbtag.tex" % i)

  # Compute final efficiency

  template = Template(r"""\subsubsection{${btag} b-tag}
  \begin{tabular}{|c|c|c|c|c|}
  \hline
  \mtt & 750 GeV & 1000 GeV & 1250 GeV & 1500 GeV\\
  \hline
  $$\epsilon(Z^{\prime})$$, semi-mu & $eff_mu_750 & $eff_mu_1000 & $eff_mu_1250 & $eff_mu_1500\\
  $$\epsilon(Z^{\prime})$$, semi-e & $eff_e_750 & $eff_e_1000 & $eff_e_1250 & $eff_e_1500\\
  \hline
  \end{tabular}""")

  selEff_mu = {}
  selEff_e = {}
  hltEff_mu = {}
  hltEff_e = {}
  eff = {}
  for mass in masses:
    jsonFile = open(pwd + "/efficiencies.json")
    jsonValues = json.load(jsonFile)
    jsonFile.close()
  
    strMass = str(mass)
    selEff_mu[strMass] = jsonValues[analysisUUID][strMass][str(i)]["nominal"][0]
    selEff_e[strMass] = jsonValues[analysisUUID][strMass][str(i)]["nominal"][1]
    hltEff_mu[strMass] = jsonValues[analysisUUID][strMass][str(i)]["nominal"][2]
    hltEff_e[strMass] = jsonValues[analysisUUID][strMass][str(i)]["nominal"][3]
  
  from ctypes import cdll, c_double
  lib = cdll.LoadLibrary("./libUtils.so")
  
  lib.computeEfficiencyMuons_1btag.restype = c_double
  lib.computeEfficiencyMuons_2btag.restype = c_double
  lib.computeEfficiencyElectrons_1btag.restype = c_double
  lib.computeEfficiencyElectrons_2btag.restype = c_double
  for mass in masses:
    strMass = str(mass)
  
    if i == 1:
      eff["eff_mu_" + strMass] = round(lib.computeEfficiencyMuons_1btag(c_double(selEff_mu[strMass]), c_double(hltEff_mu[strMass])) * 100, 2) 
      eff["eff_e_" + strMass] = round(lib.computeEfficiencyElectrons_1btag(c_double(selEff_e[strMass]), c_double(hltEff_e[strMass])) * 100, 2)
    else:
      eff["eff_mu_" + strMass] = round(lib.computeEfficiencyMuons_2btag(c_double(selEff_mu[strMass]), c_double(hltEff_mu[strMass])) * 100, 2) 
      eff["eff_e_" + strMass] = round(lib.computeEfficiencyElectrons_2btag(c_double(selEff_e[strMass]), c_double(hltEff_e[strMass])) * 100, 2)
  
  f = open(tmp + "/total_eff.tex", "a+")
  f.write(template.substitute(eff, btag = i))
  f.close()

# data_2012_nominal_1000_crystalball_faltB_2_btag/
# Second, sigma ref
template = Template(ur"""\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${analysis_name}_${btag}_btag/data_2012_nominal_${mass}_fitRes_${analysis_name}.pdf}\\
Nominal\\
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${analysis_name}_${btag}_btag/data_2012_nominal_${mass}_fitRes_${analysis_name}_log.pdf}\\
Nominal, échelle log\\
\end{minipage}

\begin{itemize}
\item $$\chi^2 = ${chi2}$$
\item Statut du fit : ${fit}
\end{itemize}""")
for mass in masses:
  jsonFile = open(pwd + "/sigma_reference.json")
  jsonValues = json.load(jsonFile)
  jsonFile.close()

  chi2 = jsonValues[analysisUUID][str(mass)][btag]["chi2"]
  fit = "OK" if (jsonValues[analysisUUID][str(mass)][btag]["fit_status"] == 0 and jsonValues[analysisUUID][str(mass)][btag]["fit_covQual"] == 3) else "Echec"

  f = open(tmp + "/sigma_ref_%d.tex" % mass, "w")
  f.write(template.substitute(PARAMS, mass = mass, chi2 = chi2, fit = fit).encode("utf-8"))
  f.close()

template = Template(r"""\begin{tabular}{|c|c|c|c|c|}
\hline
\mtt & 750 GeV & 1000 GeV & 1250 GeV & 1500 GeV\\
\hline
$$\sigma$$ (pb) & $sigma_750 & $sigma_1000 & $sigma_1250 & $sigma_1500\\
\hline
\end{tabular}""")

datas = {}
jsonFile = open(pwd + "/sigma_reference.json")
jsonValues = json.load(jsonFile)
jsonFile.close()
for mass in masses:
  datas["chi2_" + str(mass)] = jsonValues[analysisUUID][str(mass)][btag]["chi2"]
  datas["sigma_" + str(mass)] = jsonValues[analysisUUID][str(mass)][btag]["sigma"]

f = open(tmp + "/sigma_ref.tex", "w")
f.write(template.substitute(datas))
f.close()

# Third, systematics
datas = {}
if useSystematics:
  jsonFile = open(pwd + "/systematics_parameters.json")
  jsonValues = json.load(jsonFile)
  jsonFile.close()
  jsonFile = open(pwd + "/systematics.json")
  jsonSyst = json.load(jsonFile)
  jsonFile.close()

reducedJEC = ["up", "down"]

if useJECSyst:
  # JEC
  template = Template(r"""
  \begin{adjustwidth}{-2cm}{-2cm}
  \begin{center}
  \begin{tabular}{|c|c|c|c|c|c|c|c|c|}
  \hline
  \mtt & \multicolumn{2}{|c|}{750 GeV} & \multicolumn{2}{|c|}{1000 GeV} & \multicolumn{2}{|c|}{1250 GeV} & \multicolumn{2}{|c|}{1500 GeV}\\
  \hline
   & JEC up & JEC down & JEC up & JEC down & JEC up & JEC down & JEC up & JEC down\\
  \hline
  \hline
  $$\chi^2$$ & $c_750_up & $c_750_down & $c_1000_up & $c_1000_down & $c_1250_up & $c_1250_down & $c_1500_up & $c_1500_down\\
  Fit & $f_750_up & $f_750_down & $f_1000_up & $f_1000_down & $f_1250_up & $f_1250_down & $f_1500_up & $f_1500_down\\
  $$\sigma$$ (pb) & $s_750_up & $s_750_down & $s_1000_up & $s_1000_down & $s_1250_up & $s_1250_down & $s_1500_up & $s_1500_down\\
  \hline
  $$\sigma_{syst}$$ (pb) & \multicolumn{2}{|c|}{$s_750} & \multicolumn{2}{|c|}{$s_1000} & \multicolumn{2}{|c|}{$s_1250} & \multicolumn{2}{|c|}{$s_1500}\\
  \hline
  \end{tabular}
  \end{center}
  \end{adjustwidth}""")
  
  for mass in masses:
    if "jec" in jsonValues[analysisUUID][str(mass)][btag]:
      for jec in reducedJEC:
        if "JEC" + jec in jsonValues[analysisUUID][str(mass)][btag]["jec"]:
          datas["c_" + str(mass) + "_" + jec] = round(jsonValues[analysisUUID][str(mass)][btag]["jec"]["JEC" + jec]["chi2"], 4)
          datas["s_" + str(mass) + "_" + jec] = round(jsonValues[analysisUUID][str(mass)][btag]["jec"]["JEC" + jec]["sigma"], 4)
          datas["f_" + str(mass) + "_" + jec] = "OK" if (jsonValues[analysisUUID][str(mass)][btag]["jec"]["JEC" + jec]["fit_status"] == 0 and jsonValues[analysisUUID][str(mass)][btag]["jec"]["JEC" + jec]["fit_covQual"] == 3) else "Echec"
        else:
          datas["c_" + str(mass) + "_" + jec] = "N/A"
          datas["s_" + str(mass) + "_" + jec] = "N/A"
          datas["f_" + str(mass) + "_" + jec] = "N/A"


      datas["s_" + str(mass)] = round(jsonSyst[analysisUUID][str(mass)][btag]["jec"], 4)
    else:
      for jec in reducedJEC:
        datas["c_" + str(mass) + "_" + jec] = "N/A"
        datas["s_" + str(mass) + "_" + jec] = "N/A"
        datas["f_" + str(mass) + "_" + jec] = "N/A"

      datas["s_" + str(mass)] = "N/A"

  f = open(tmp + "/syst_jec.tex", "w")
  f.write(template.substitute(datas))
  f.close()
else:
  f = open(tmp + "/syst_jec.tex", "w")
  f.write(r"\textbf{No JEC systematics for this analysis}")
  f.close()
  
if useSignalSyst:
  # Signal
  template = Template(ur"""
  \begin{center}
  \begin{longtable}{|c|c|c|c|c|}
  \hline Paramètre & Variation & $$\chi^2$$ & $$\sigma$$ (pb) & Statut du fit \\ \hline \hline
  \endfirsthead
  
  \hline Paramètre & Variation & $$\chi^2$$ & $$\sigma$$ (pb) & Statut du fit \\ \hline
  \endhead
  
  \hline \multicolumn{5}{|r|}{La suite page suivante} \\ \hline
  \endfoot
  
  \endlastfoot
  
  $content
  \end{longtable}
  \end{center}""")
  
  arrayLine = ""
  for mass in masses:
    arrayLine = arrayLine + "\multicolumn{5}{|l|}{$m = %d$ GeV}\\\\\n\\hline" % mass
    if "signal" in jsonValues[analysisUUID][str(mass)][btag]:
      for param, values in jsonValues[analysisUUID][str(mass)][btag]["signal"].items():
        arrayLine = arrayLine + "\multirow{2}{*}{%s}" % param.replace("_", r"\_")
        for var in reducedJEC:
          chi2 = round(values["chi2"][var], 4)
          sigma = round(values["sigma"][var], 4)
          fit = "OK" if (values["fit_status"][var] == 0 and values["fit_covQual"][var] == 3) else "Echec"
          arrayLine = arrayLine + " & %s & %.04f & %.04f & %s\\\\\n" % (var, chi2, sigma, fit)
          if var == "up":
            arrayLine = arrayLine + "\\cline{2-5}"
        
        arrayLine = arrayLine + "\\hline\\hline\n"
      arrayLine = arrayLine + "\multicolumn{5}{|l|}{$\sigma_{syst} = %.04f$ pb}\\\\ \\hline\\hline" % jsonSyst[analysisUUID][str(mass)][btag]["signal_pdf"]  
    elif "signal_pdf" in jsonSyst[analysisUUID][str(mass)][btag]:
      arrayLine = "\multicolumn{5}{|l|}{$\sigma_{syst} = %.04f$ pb}\\\\ \\hline\\hline" % jsonSyst[analysisUUID][str(mass)][btag]["signal_pdf"]
  
  f = open(tmp + "/syst_signal.tex", "w")
  f.write(template.substitute(content = arrayLine).encode("utf-8"))
  f.close()
else:
  f = open(tmp + "/syst_signal.tex", "w")
  f.write(r"\textbf{No signal systematics for this analysis}")
  f.close()
  
if useBkgSyst:
  # Background
  template = Template(ur"""
  \begin{center}
  \begin{longtable}{|c|c|c|c|}
  \hline Fonction de bkg & $$\chi^2$$ & $$\sigma$$ (pb) & Statut du fit \\ \hline \hline
  \endfirsthead
  
  \hline Fonction de bkg & $$\chi^2$$ & $$\sigma$$ (pb) & Statut du fit \\ \hline
  \endhead
  
  \hline \multicolumn{4}{|r|}{La suite page suivante} \\ \hline
  \endfoot
  
  \endlastfoot
  
  $content
  \end{longtable}
  \end{center}""")
  
  arrayLine = ""
  for mass in masses:
    arrayLine = arrayLine + "\multicolumn{4}{|l|}{$m = %d$ GeV}\\\\\n\\hline\n" % mass
    for param, values in jsonValues[analysisUUID][str(mass)][btag]["background"].items():
      arrayLine = arrayLine + param.replace("_", r"\_")
      chi2 = round(values["chi2"], 4)
      sigma = round(values["sigma"], 4)
      fit = "OK" if (values["fit_status"] == 0 and values["fit_covQual"] == 3) else "Echec"
      arrayLine = arrayLine + " & %.04f & %.04f & %s\\\\\n" % (chi2, sigma, fit)
  
      arrayLine = arrayLine + "\\hline\\hline\n"
    arrayLine = arrayLine + "\multicolumn{4}{|l|}{$\sigma_{syst} = %.04f$ pb}\\\\ \\hline\\hline\n" % jsonSyst[analysisUUID][str(mass)][btag]["background_pdf"]
  
  f = open(tmp + "/syst_bkg.tex", "w")
  f.write(template.substitute(content = arrayLine).encode("utf-8"))
  f.close()
else:
  f = open(tmp + "/syst_bkg.tex", "w")
  f.write(r"\textbf{No background systematics for this analysis}")
  f.close()
  
##################################################################################
##################################################################################
##################################################################################
# Likelihood scan

template = Template(ur"""\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${analysis_name}_${btag}_btag/data_2012_nominal_${mass}_likscan_${analysis_name}.pdf}\\
Likelihood scan
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${analysis_name}_${btag}_btag/data_2012_nominal_${mass}_pdfscan_${analysis_name}.pdf}\\
PDF scan
\end{minipage}

\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${analysis_name}_${btag}_btag/data_2012_nominal_${mass}_pdfscan_wsyst_${analysis_name}.pdf}\\
PDF scan + systématiques
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${analysis_name}_${btag}_btag/data_2012_nominal_${mass}_pdfscan_wsyst_cut_${analysis_name}.pdf}\\
PDF scan + systématiques pour $$N_{sig} > 0$$
\end{minipage}

""")

jsonFile = open(pwd + "/likelihood_scan.json")
jsonValues = json.load(jsonFile)
jsonFile.close()
values = {}

for mass in masses:
  f = open(tmp + "/likscan_%d.tex" % mass, "w")
  f.write(template.substitute(PARAMS, mass = mass).encode("utf-8"))
  f.close()

  values["limit_" + str(mass)] = round(jsonValues[analysisUUID][str(mass)][btag]["scan_wsyst_cut_limit"], 4)
  
template = Template(ur"""\begin{tabular}{|c|c|c|c|c|}
\hline
\mtt & 750 GeV & 1000 GeV & 1250 GeV & 1500 GeV\\
\hline
Limite observée (pb) & $limit_750 & $limit_1000 & $limit_1250 & $limit_1500\\
\hline
\end{tabular}""")

f = open(tmp + "/limites_obs.tex", "w")
f.write(template.substitute(values).encode("utf-8"))
f.close()

observed_limits = values

############################################################################
############################################################################
## Toy MC

# This store num_jobs, num_toys_per_job and num_toys
p = subprocess.Popen(["python", "toys/submitjobs.py", "--python", "--b-tag", btag], stdout=subprocess.PIPE)
for line in p.stdout.readlines():
  exec(line)

template = Template(r"""\begin{itemize}
\item Nombre de toys par masse : $num_toys
\item Nombre de jobs par masse : $num_jobs
\item Nombre de toys par jobs : $num_toys_per_job
\end{itemize}""")

f = open(tmp + "/toys.tex", "w")
f.write(template.substitute(num_jobs=num_jobs, num_toys = num_toys, num_toys_per_job = num_toys_per_job).encode("utf-8"))
f.close()

template = Template(ur"""\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${analysis_name}_LimitNLLToyExp.pdf}\\
Nll Toy exp
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${analysis_name}_pull.pdf}\\
Pull
\end{minipage}
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${analysis_name}_LimitPlotZ.pdf}\\
Limite Z'
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${analysis_name}_LimitErrors.pdf}\\
Erreur sur limite
\end{minipage}
\begin{center}
  \includegraphics[width=0.50\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${analysis_name}_LimitPlotsEMu.pdf}\\
  Plots $$e$$ $$\mu$$
\end{center}
""")

jsonFile = open(pwd + "/expected_limits.json")
jsonValues = json.load(jsonFile)
jsonFile.close()
values = {}

for mass in masses:
  f = open(tmp + "/toys_%d.tex" % mass, "w")
  f.write(template.substitute(PARAMS, mass = mass).encode("utf-8"))
  f.close()

  values["elimit_" + str(mass)] = round(jsonValues[analysisUUID][str(mass)][btag]["median"], 4)
  values["m68_" + str(mass)] = round(jsonValues[analysisUUID][str(mass)][btag]["widthM68"], 4)
  values["p68_" + str(mass)] = round(jsonValues[analysisUUID][str(mass)][btag]["widthP68"], 4)
  values["m95_" + str(mass)] = round(jsonValues[analysisUUID][str(mass)][btag]["widthM95"], 4)
  values["p95_" + str(mass)] = round(jsonValues[analysisUUID][str(mass)][btag]["widthP95"], 4)

template = Template(ur"""{
\renewcommand{\arraystretch}{2}
\begin{tabular}{|c|c|c|c|c|}
\hline
\mtt & 750 GeV & 1000 GeV & 1250 GeV & 1500 GeV\\
\hline
Limite attendue (pb) & $elimit_750 & $elimit_1000 & $elimit_1250 & $elimit_1500\\
\hline \hline
Bande d'exclusion (68\%) (pb) & $$^{+${p68_750}}_{-${m68_750}}$$ & $$^{+${p68_1000}}_{-${m68_1000}}$$ & $$^{+${p68_1250}}_{-${m68_1250}}$$ & $$^{+${p68_1500}}_{-${m68_1500}}$$\\
Bande d'exclusion (95\%) (pb) & $$^{+${p95_750}}_{-${m95_750}}$$ & $$^{+${p95_1000}}_{-${m95_1000}}$$ & $$^{+${p95_1250}}_{-${m95_1250}}$$ & $$^{+${p95_1500}}_{-${m95_1500}}$$\\
\hline
\end{tabular}
}""")

f = open(tmp + "/limites_exp.tex", "w")
f.write(template.substitute(values).encode("utf-8"))
f.close()

expected_limits = values

#########################################################################
#########################################################################
#### Limites

template = Template(ur"""{
\renewcommand{\arraystretch}{2}
\begin{tabular}{|c|c|c|c|c|}
\hline
\mtt & 750 GeV & 1000 GeV & 1250 GeV & 1500 GeV\\
\hline
Limite observée (pb) & $limit_750 & $limit_1000 & $limit_1250 & $limit_1500\\
\hline
Limite attendue (pb) & $elimit_750 & $elimit_1000 & $elimit_1250 & $elimit_1500\\
\hline \hline
Bande d'exclusion (68\%) (pb) & $$^{+${p68_750}}_{-${m68_750}}$$ & $$^{+${p68_1000}}_{-${m68_1000}}$$ & $$^{+${p68_1250}}_{-${m68_1250}}$$ & $$^{+${p68_1500}}_{-${m68_1500}}$$\\
Bande d'exclusion (95\%) (pb) & $$^{+${p95_750}}_{-${m95_750}}$$ & $$^{+${p95_1000}}_{-${m95_1000}}$$ & $$^{+${p95_1250}}_{-${m95_1250}}$$ & $$^{+${p95_1500}}_{-${m95_1500}}$$\\
\hline
\end{tabular}

\begin{center}
  \includegraphics[width=0.70\textwidth]{${pwd}/limitCurve_2012_${analysis_name}_${btag}btag.pdf}
\end{center}
}""")

values = {}
values.update(observed_limits)
values.update(expected_limits)
values.update(PARAMS)

f = open(tmp + "/limites.tex", "w")
f.write(template.substitute(values).encode("utf-8"))
f.close()

shutil.copy(os.getcwd() + "/template/analysis_summary.tex", tmp)
os.chdir(tmp)
os.system("pdflatex analysis_summary.tex")
os.system("pdflatex analysis_summary.tex")
os.chdir(pwd)

bookmark = pwd + "/bookmark_%s_%s.pdf" % (analysisName, datetime.date.today().isoformat())
shutil.copy(tmp + "/analysis_summary.pdf", bookmark)

print ""
print "Bookmark saved as %s" % bookmark
