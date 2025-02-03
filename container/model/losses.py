# based on https://github.com/runopti/stg/blob/master/python/stg/losses.py 

import torch 
from lifelines.utils import concordance_index


def NegPartialLogLikelihood(logits, fail_indicator, fail_time):
    '''
    fail_indicator: 1 if the sample fails, 0 if the sample is censored.
    fail_time: time until failure or censorship
    logits: raw output from model
    '''
    _sort_idx = torch.argsort(fail_time, descending=True)
    _logits = logits[_sort_idx]
    _fail_indicator = fail_indicator[_sort_idx]

    _hazard_ratio = torch.exp(_logits)

    ### contains now sum of probabilities after time Ti
    ### Hense, need to sort in descending order
    _cumsum_hazard_ratio = torch.cumsum(_hazard_ratio, 0)
    
    _log_at_risk = torch.log(_cumsum_hazard_ratio)
    _likelihood = _logits - _log_at_risk

    # sum over events 
    _uncensored_likelihood = _likelihood * _fail_indicator
    _logL = torch.sum(_uncensored_likelihood)

    # negative average log-likelihood
    _events = torch.sum(_fail_indicator)
    return -1.0*_logL / _events


def calc_concordance_index(logits, fail_indicator, fail_time):
    """
    Compute the concordance-index value.
    Parameters:
        fail_indicator -  event happened or censored
        fail_time - event/censor time
        logits: - predictive proportional risk of network.
    Returns:
        concordance index.
    """
    # need to convert hazard risk to the same order as survival times. 
    # so, higher hazard -> less time
    hr_pred = -logits
    #
    ci = concordance_index(fail_time,
                           hr_pred,
                           fail_indicator)
    return ci
