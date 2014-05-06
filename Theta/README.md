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
 - `single_top`: single top
 - `single_antitop`: single anti-top
 - `dibosons`: di-bosons
 - `zjets`: Z + * jets (DY)
 - `wjets`: W + * jets

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
 - `matching`: matching generator uncertainties
 - `scale`: renormalization scale uncertainties

# Normalization systematic uncertainties

 - `lumi`: 2.6%
 - `TT`: 15%
 - `single-top`: 30%
 - `w jets / z jets`: 50%
 - `dibosons`: 50%

# How-to run

Possible model type:
 - `scalar`: Higgs scalar
 - `pseudoscalar`: Higgs pseudo-scalar
 - `zprime_narrow`: Z' narrow
 - `zprime_large`: Z' large 

Possible options (all booleans):
 - `stat_only`: run without any systematics
 - `lumi_only`: run only with luminosity uncertainty
 - `no_data`: run only expected, no observed
 - `no_shape_syst`: run without any shape uncertainties, only normalisation uncertainties

## Steps

Save analysis summary:
```shell
theta_driver preprocess --model 'model: type = pseudoscalar' --workdir summary_scalar --analysis summary file.root
```

Perform a maximum likehood fit of all the parameters :
```shell
theta_driver preprocess --model 'model: type = pseudoscalar' --workdir mle_scalar --analysis mle file.root
```

Compute expected and observed limit.
```shell
theta_driver preprocess --model 'model: type = pseudoscalar' --workdir bayesian_scalar --analysis 'bayesian: iterations = 160000' file.root
```

Compute pvalues
```shell
theta_driver preprocess --model 'model: type = pseudoscalar' --workdir pvalues2_scalar --analysis pvalues2 file.root
```
