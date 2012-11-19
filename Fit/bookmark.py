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

pdfBackgroundName = "faltB"

f = open("parameters.json")
params = json.load(f)
pdfSignalName = params["parameters"]["pdf"]["signal"]
f.close()

masses = [750, 1000, 1250, 1500]
jecs = ["nominal", "JECup", "JECdown"]
pwd = os.getcwd()

btag = str(args.btag)

tmp = tempfile.mkdtemp(dir="/scratch")

PARAMS = {
    'pwd': pwd,
    'sig_pdf': pdfSignalName,
    'bkg_pdf': pdfBackgroundName,
    'btag': btag
}

# Check if we are using systematics or not
if subprocess.call(["grep", "^#define NO_SYST", os.path.join(pwd, "fitMtt.cc")]) == 0:
  NO_SYST = True
else:
  NO_SYST = False

def run(program, *args):
  pid = os.fork()
  if not pid:
    os.execvp(program, (program,) + args)
  return os.wait()[1]

if doIntro:
  intro_file = tempfile.mkstemp(dir="/scratch/")
  run("vim", intro_file[1])
  shutil.copy(intro_file[1], tmp + "/intro.tex")
else:
  open(tmp + "/intro.tex", "w").close()

if NO_SYST:
  f = open(tmp + "/intro.tex", "w+")
  f.write(r"\begin{center}\textcolor{red}{WARNING: Analysis ran without systematics!}\end{center}")
  f.close()

# First, frit

if btag != "3":
  template_full = Template(r"""\begin{minipage}{0.33\textwidth} \centering
  \includegraphics[width=0.99\textwidth]{${pwd}/frit/nominal-Zprime${mass}_${sig_pdf}_${btag}_btag_fitCB.pdf}\\
  Nominal\\
  $$\chi^2 = ${chi2_nominal}$$
  \end{minipage}%
  \begin{minipage}{0.33\textwidth} \centering
  \includegraphics[width=0.99\textwidth]{${pwd}/frit/JECup-Zprime${mass}_${sig_pdf}_${btag}_btag_fitCB.pdf}\\
  JEC Up\\
  $$\chi^2 = ${chi2_JECup}$$
  \end{minipage}%
  \begin{minipage}{0.33\textwidth} \centering
  \includegraphics[width=0.99\textwidth]{${pwd}/frit/JECdown-Zprime${mass}_${sig_pdf}_${btag}_btag_fitCB.pdf}\\
  JEC Down\\
  $$\chi^2 = ${chi2_JECdown}$$
  \end{minipage}""")

  template_reduced = Template(r"""\begin{center}\begin{minipage}{0.50\textwidth} \centering
  \includegraphics[width=0.99\textwidth]{${pwd}/frit/nominal-Zprime${mass}_${sig_pdf}_${btag}_btag_fitCB.pdf}\\
  Nominal\\
  $$\chi^2 = ${chi2_nominal}$$
  \end{minipage}\end{center}""")
  
  for mass in masses:
    jsonFile = open("frit_efficiencies.json")
    chi2 = {}
    jsonValues = json.load(jsonFile)
    jsonFile.close()
  
    for jec in jecs:
      chi2["chi2_" + jec] = jsonValues[str(mass)][btag][jec]["chi2"]
  
    chi2.update(PARAMS)
    f = open(tmp + "/frit_%d.tex" % mass, "w")

    if (not NO_SYST) and (os.path.exists("${pwd}/frit/JECup-Zprime${mass}_${sig_pdf}_${btag}_btag_fitCB.pdf" % (chi2))):
      f.write(template_full.substitute(chi2, mass = mass))
    else:
      f.write(template_reduced.substitute(chi2, mass = mass))

    f.close()

  shutil.copy(pwd + ("/efficiencies_table_%s_btag.tex" % btag), tmp + "/efficiencies_table.tex")

  # Compute final efficiency

  template = Template(r"""\begin{tabular}{|c|c|c|c|c|}
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
    jsonFile = open("efficiencies.json")
    jsonValues = json.load(jsonFile)
    jsonFile.close()
  
    strMass = str(mass)
    selEff_mu[strMass] = jsonValues[strMass][btag]["nominal"][0]
    selEff_e[strMass] = jsonValues[strMass][btag]["nominal"][1]
    hltEff_mu[strMass] = jsonValues[strMass][btag]["nominal"][2]
    hltEff_e[strMass] = jsonValues[strMass][btag]["nominal"][3]
  
  from ctypes import cdll, c_double
  lib = cdll.LoadLibrary("./libUtils.so")
  
  lib.computeEfficiency.restype = c_double
  for mass in masses:
    strMass = str(mass)
  
    eff["eff_mu_" + strMass] = round(lib.computeEfficiency(c_double(selEff_mu[strMass]), c_double(hltEff_mu[strMass])) * 100, 2) 
    eff["eff_e_" + strMass] = round(lib.computeEfficiency(c_double(selEff_e[strMass]), c_double(hltEff_e[strMass])) * 100, 2)
  
  f = open(tmp + "/total_eff.tex", "w")
  f.write(template.substitute(eff))
  f.close()

# data_2012_nominal_1000_crystalball_faltB_2_btag/
# Second, sigma ref
template = Template(ur"""\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${sig_pdf}_${btag}_btag/data_2012_nominal_${mass}_fitRes_${sig_pdf}.pdf}\\
Nominal\\
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${sig_pdf}_${btag}_btag/data_2012_nominal_${mass}_fitRes_${sig_pdf}_log.pdf}\\
Nominal, échelle log\\
\end{minipage}

\begin{itemize}
\item $$\chi^2 = ${chi2}$$
\item Statut du fit : ${fit}
\end{itemize}""")
for mass in masses:
  jsonFile = open("sigma_reference.json")
  jsonValues = json.load(jsonFile)
  jsonFile.close()

  chi2 = jsonValues[str(mass)][btag]["chi2"]
  fit = "OK" if (jsonValues[str(mass)][btag]["fit_status"] == 0 and jsonValues[str(mass)][btag]["fit_covQual"] == 3) else "Echec"

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
jsonFile = open("sigma_reference.json")
jsonValues = json.load(jsonFile)
jsonFile.close()
for mass in masses:
  datas["chi2_" + str(mass)] = jsonValues[str(mass)][btag]["chi2"]
  datas["sigma_" + str(mass)] = jsonValues[str(mass)][btag]["sigma"]

f = open(tmp + "/sigma_ref.tex", "w")
f.write(template.substitute(datas))
f.close()

# Third, systematics

if (not NO_SYST) and (btag != "3"):

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
  
  datas = {}
  jsonFile = open("systematics_parameters.json")
  jsonValues = json.load(jsonFile)
  jsonFile.close()
  jsonFile = open("systematics.json")
  jsonSyst = json.load(jsonFile)
  jsonFile.close()
  
  reducedJEC = ["up", "down"]
  for mass in masses:
    for jec in reducedJEC:
      datas["c_" + str(mass) + "_" + jec] = round(jsonValues[str(mass)][btag]["jec"]["JEC" + jec]["chi2"], 4)
      datas["s_" + str(mass) + "_" + jec] = round(jsonValues[str(mass)][btag]["jec"]["JEC" + jec]["sigma"], 4)
      datas["f_" + str(mass) + "_" + jec] = "OK" if (jsonValues[str(mass)][btag]["jec"]["JEC" + jec]["fit_status"] == 0 and jsonValues[str(mass)][btag]["jec"]["JEC" + jec]["fit_covQual"] == 3) else "Echec"
      
    datas["s_" + str(mass)] = round(jsonSyst[str(mass)][btag]["jec"], 4)

  f = open(tmp + "/syst_jec.tex", "w")
  f.write(template.substitute(datas))
  f.close()
  
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
    for param, values in jsonValues[str(mass)][btag]["signal"].items():
      arrayLine = arrayLine + "\multirow{2}{*}{%s}" % param.replace("_", r"\_")
      for var in reducedJEC:
        chi2 = round(values["chi2"][var], 4)
        sigma = round(values["sigma"][var], 4)
        fit = "OK" if (values["fit_status"][var] == 0 and values["fit_covQual"][var] == 3) else "Echec"
        arrayLine = arrayLine + " & %s & %.04f & %.04f & %s\\\\\n" % (var, chi2, sigma, fit)
        if var == "up":
          arrayLine = arrayLine + "\\cline{2-5}"
      
      arrayLine = arrayLine + "\\hline\\hline\n"
    arrayLine = arrayLine + "\multicolumn{5}{|l|}{$\sigma_{syst} = %.04f$ pb}\\\\ \\hline\\hline" % jsonSyst[str(mass)][btag]["signal_pdf"]
  
  
  f = open(tmp + "/syst_signal.tex", "w")
  f.write(template.substitute(content = arrayLine).encode("utf-8"))
  f.close()
  
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
    for param, values in jsonValues[str(mass)][btag]["background"].items():
      arrayLine = arrayLine + param.replace("_", r"\_")
      chi2 = round(values["chi2"], 4)
      sigma = round(values["sigma"], 4)
      fit = "OK" if (values["fit_status"] == 0 and values["fit_covQual"] == 3) else "Echec"
      arrayLine = arrayLine + " & %.04f & %.04f & %s\\\\\n" % (chi2, sigma, fit)
  
      arrayLine = arrayLine + "\\hline\\hline\n"
    arrayLine = arrayLine + "\multicolumn{4}{|l|}{$\sigma_{syst} = %.04f$ pb}\\\\ \\hline\\hline\n" % jsonSyst[str(mass)][btag]["background_pdf"]
  
  f = open(tmp + "/syst_bkg.tex", "w")
  f.write(template.substitute(content = arrayLine).encode("utf-8"))
  f.close()
  
##################################################################################
##################################################################################
##################################################################################
# Likelihood scan

template = Template(ur"""\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${sig_pdf}_${btag}_btag/data_2012_nominal_${mass}_likscan_${sig_pdf}.pdf}\\
Likelihood scan
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${sig_pdf}_${btag}_btag/data_2012_nominal_${mass}_pdfscan_${sig_pdf}.pdf}\\
PDF scan
\end{minipage}

\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${sig_pdf}_${btag}_btag/data_2012_nominal_${mass}_pdfscan_wsyst_${sig_pdf}.pdf}\\
PDF scan + systématiques
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/data_2012_nominal_${mass}_${sig_pdf}_${btag}_btag/data_2012_nominal_${mass}_pdfscan_wsyst_cut_${sig_pdf}.pdf}\\
PDF scan + systématiques pour $$N_{sig} > 0$$
\end{minipage}

""")

jsonFile = open("likelihood_scan.json")
jsonValues = json.load(jsonFile)
jsonFile.close()
values = {}

for mass in masses:
  f = open(tmp + "/likscan_%d.tex" % mass, "w")
  f.write(template.substitute(PARAMS, mass = mass).encode("utf-8"))
  f.close()

  values["limit_" + str(mass)] = round(jsonValues[str(mass)][btag]["scan_wsyst_cut_limit"], 4)
  
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
\includegraphics[width=0.99\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${sig_pdf}_LimitNLLToyExp.pdf}\\
Nll Toy exp
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${sig_pdf}_pull.pdf}\\
Pull
\end{minipage}
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${sig_pdf}_LimitPlotZ.pdf}\\
Limite Z'
\end{minipage}%
\begin{minipage}{0.49\textwidth} \centering
\includegraphics[width=0.99\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${sig_pdf}_LimitErrors.pdf}\\
Erreur sur limite
\end{minipage}
\begin{center}
  \includegraphics[width=0.50\textwidth]{${pwd}/toys/plots/${btag}-btag/data_2012_Zprime${mass}_${sig_pdf}_LimitPlotsEMu.pdf}\\
  Plots $$e$$ $$\mu$$
\end{center}
""")

jsonFile = open("expected_limits.json")
jsonValues = json.load(jsonFile)
jsonFile.close()
values = {}

for mass in masses:
  f = open(tmp + "/toys_%d.tex" % mass, "w")
  f.write(template.substitute(PARAMS, mass = mass).encode("utf-8"))
  f.close()

  values["elimit_" + str(mass)] = round(jsonValues[str(mass)][btag]["median"], 4)
  values["m68_" + str(mass)] = round(jsonValues[str(mass)][btag]["widthM68"], 4)
  values["p68_" + str(mass)] = round(jsonValues[str(mass)][btag]["widthP68"], 4)
  values["m95_" + str(mass)] = round(jsonValues[str(mass)][btag]["widthM95"], 4)
  values["p95_" + str(mass)] = round(jsonValues[str(mass)][btag]["widthP95"], 4)

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
  \includegraphics[width=0.70\textwidth]{${pwd}/limitCurve_2012_${sig_pdf}_${bkg_pdf}_${btag}btag.pdf}
\end{center}
}""")

values = {}
values.update(observed_limits)
values.update(expected_limits)
values.update(PARAMS)

f = open(tmp + "/limites.tex", "w")
f.write(template.substitute(values).encode("utf-8"))
f.close()

shutil.copy(pwd + "/template/analysis_summary.tex", tmp)
os.chdir(tmp)
os.system("pdflatex analysis_summary.tex")
os.system("pdflatex analysis_summary.tex")
os.chdir(pwd)

bookmark = "bookmark_%s_%s_%s.pdf" % (pdfSignalName, pdfBackgroundName, datetime.date.today().isoformat())
shutil.copy(tmp + "/analysis_summary.pdf", bookmark)

print ""
print "Bookmark saved as %s" % bookmark
