STAT_ONLY = False
LUMI_ONLY = False
NO_DATA = False
NO_SHAPE_SYST = False

def getSystematic(name):
    p = name.split("__")
    if len(p) > 2:
        return p[2]
    else:
        return None

def getProcess(name):
    p = name.split("__")
    return p[1]

def default_filter(name):
    global STAT_ONLY
    global LUMI_ONLY
    global NO_DATA
    global NO_SHAPE_SYST

    if (STAT_ONLY or LUMI_ONLY or NO_SHAPE_SYST) and getSystematic(name) is not None:
        return False

    if NO_DATA and "DATA" in name:
        return False

    # Ignore scale & matching uncertainties, except for tt
    syst = getSystematic(name)
    if syst is not None:
        if 'scale' in syst or 'matching' in syst:
            process = getProcess(name)
            if 'TT' not in process:
                return False

    return True


def select_scalar_signal_hist(name):

    if not default_filter(name):
        return False

    if "zp" in name:
        return False

    if '_pseudoscalar' in name:
        return False

    return True

def select_pseudoscalar_signal_hist(name):

    if not default_filter(name):
        return False

    # Exclude Z'
    if "zp" in name:
        return False

    if '_scalar' in name:
        return False

    #if "pseudoscalar" in name and (not '400' in name):
        #return False

    return True

def select_zp_narrow_signal_hist(name):

    if not default_filter(name):
        return False

    # Exclude Higgs
    if "scalar" in name:
        return False

    if "large" in name:
        return False

    if 'zp' in name:
        p = getProcess(name)
        m = p.split('_')[0]
        m = int(m[2:])

        if m > 1000:
            return False

    return True

def build_theta_model(file, filter, signal):
    global LUMI_ONLY
    global STAT_ONLY

    model = build_model_from_rootfile(file, filter, include_mc_uncertainties = True)
    model.fill_histogram_zerobins()
    model.set_signal_processes(signal)

    if not STAT_ONLY:
        for p in model.processes:
            model.add_lognormal_uncertainty('lumi', math.log(1.026), p)

    if not STAT_ONLY and not LUMI_ONLY:
        model.add_lognormal_uncertainty('ttbar_rate', math.log(1.15), 'TT')

        model.add_lognormal_uncertainty('st_rate', math.log(1.30), 'single_top')
        model.add_lognormal_uncertainty('sat_rate', math.log(1.30), 'single_antitop')

        model.add_lognormal_uncertainty('w_rate', math.log(1.50), 'wjets')
        model.add_lognormal_uncertainty('z_rate', math.log(1.50), 'zjets')

        model.add_lognormal_uncertainty('diboson_rate', math.log(1.50), 'dibosons')

    for p in model.distribution.get_parameters():
        d = model.distribution.get_distribution(p)
        if d['typ'] == 'gauss' and d['mean'] == 0.0 and d['width'] == 1.0:
            model.distribution.set_distribution_parameters(p, range = [-30.0, 30.0])

    return model

def build_model(type, file, stat_only = False, lumi_only = False, no_data = False, no_shape_syst = False):
    global STAT_ONLY
    global LUMI_ONLY
    global NO_DATA
    global NO_SHAPE_SYST

    STAT_ONLY = stat_only
    LUMI_ONLY = lumi_only
    NO_DATA = no_data
    NO_SHAPE_SYST = no_shape_syst

    options = Options()
    options.set('global', 'debug', 'True')

    if type == "scalar":
        return build_theta_model(file, select_scalar_signal_hist, "H*_scalar")
    elif type == "pseudoscalar":
        return build_theta_model(file, select_pseudoscalar_signal_hist, "H*_pseudoscalar")
    elif type == "zprime_narrow":
        return build_theta_model(file, select_zp_narrow_signal_hist, "zp*_narrow")
    elif type == "zprime_large":
        return build_theta_model(file, select_zp_large_signal_hist, "zp*_large")

    return None

