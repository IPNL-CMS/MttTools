systs = ['jec', 'jer', 'pu', 'lept']

def select_scalar_signal_hist(name):

    # Keep only jec syst
    #if 'up' in name or 'down' in name:
        #ret = False
        #for syst in systs:
            #ret = ret or syst in name

        #if not ret:
            #return False

    # Accept MC
    #if '__T_tW' in name:
        #return False

    #if '__Tbar_' in name:
        #return False

    #if '__W' in name:
        #return False

    #if '__Z' in name:
        #return False

    if not "scalar" in name:
        return True

    if '_pseudoscalar' in name:
        return False

    return True

def select_pseudoscalar_signal_hist(name):

    # Accept MC
    if not "scalar" in name:
        return True

    if '_scalar' in name:
        return False

    return True

def build_theta_model(file, filter, signal):

    model = build_model_from_rootfile(file, filter, include_mc_uncertainties = True)
    model.fill_histogram_zerobins()
    model.set_signal_processes(signal)

    for p in model.processes:
        model.add_lognormal_uncertainty('lumi', math.log(1.025), p)

    model.add_lognormal_uncertainty('ttbar_rate', math.log(1.15), 'TT')
    model.add_lognormal_uncertainty('st_rate', math.log(1.50), 'singletop')

    model.add_lognormal_uncertainty('w_rate', math.log(1.50), 'wjets')
    model.add_lognormal_uncertainty('z_rate', math.log(2.00), 'zjets')

    #for p in model.distribution.get_parameters():
        #d = model.distribution.get_distribution(p)
        #if d['typ'] == 'gauss' and d['mean'] == 0.0 and d['width'] == 1.0:
            #model.distribution.set_distribution_parameters(p, range = [-10.0, 10.0])

    return model

def build_model(type):

    if type == "scalar":
        return build_theta_model("higgs_analysis_lepton_merged.root", select_scalar_signal_hist, "H*_scalar")
    elif type == "pseudoscalar":
        return build_theta_model("higgs_analysis_lepton_merged.root", select_pseudoscalar_signal_hist, "H*_pseudoscalar")

    return None
