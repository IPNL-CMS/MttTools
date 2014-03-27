# Model

The file `model.py` contains the theta definition of the model, along with all the systematic uncertainties on the normalization

Naming schema is defined below :

## Observables

We have four observables, based on number of b-tagged jets and lepton flavor :

 - `mtt_mu_1btag`: muon, 1btag
 - `mtt_e_1btag`: electron, 1btag
 - `mtt_mu_2btag`: muon, 2btag
 - `mtt_e_2btag`: electron, 2btag

## Processus

### MC

 - `TT`: TTbar powheg
 - `T[bar]_s`: single [anti-]top, s-channel
 - `T[bar]_t`: single [anti-]top, t-channel
 - `T[bar]_tW`: single [anti-]top, tW-channel
 - `Z*Jets`: Z + * jets (DY)
 - `W*Jets`: W + * jets

### Signal
 - `H*_scalar`: Massive Higgs of mass \*, scalar resonance
 - `H*_pseudoscalar`: Massive Higgs of mass \*, pseudoscalar resonance

## Systematics

 - `jec`: jet energy correction, +/- 1 sigma
 - `jer`: jet energy resolution, +/- 1 sigma
 - `trig`: trigger scale factor, +/- 1 sigma
 - `lept`: lepton scale factor, +/- 1 sigma
 - `btag`: b-tagging scale factor, +/- 1 sigma
 - `pu`: PU reweighting, +/- 1 sigma

# Normalization systematic uncertainties

TODO

# How-to run

## Pre-processing

Save analysis summary:
```shell
theta_driver preprocess --model 'model: type = scalar' --workdir summary_scalar --analysis summary file.root
```

Perform a maximum likehood fit of all the parameters :
```shell
theta_driver preprocess --model 'model: type = scalar' --workdir mle_scalar --analysis mle file.root
```

Compute expected and observed limit. This step only create necessary files to launch crab
```shell
theta_driver preprocess --model 'model: type = scalar' --workdir bayesian_scalar --analysis 'bayesian: iterations = 160000' file.root
```
