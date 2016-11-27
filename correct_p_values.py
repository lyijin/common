#!/usr/bin/env python

"""
> correct_p_values.py <

Has a def to correct P values (to produce q values), which by default
is Benjamini-Hochberg (1995). Also incorporates the FDR adjustment from 
Yekuteli and Benjamini (1999), which ensures monotonicity (i.e. rankings are
identical pre- and post-correction). This was a problem in B-H because all
P values are multiplied by different values (sample population / rank), which 
might cause rankings to shift due to the different multiplications.

See:
http://stats.stackexchange.com/questions/870/multiple-hypothesis-testing-correction-with-benjamini-hochberg-p-values-or-q-va
R's code in implementing B-H
    BH = {
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(n/i * p[o]))[ro]
    }
"""
import collections

def correct_p_values(p_values_dict, method='BH'):
    """
    p_values_dict is defined as:
      p_values_dict[unique_key] = p_value (float)
    
    By default, Benjamini-Hochberg is used to correct the p values;
    alternatively, use 'Holm' for Bonferroni step-down correction,
    'Bonferroni' for Bonferroni correction.
    
    FDR adjustment simply implies that the corrected P value can never be
    smaller than the one that is ranked above it pre-correction.
    """

    ordered_p_values = collections.OrderedDict(
        sorted(p_values_dict.items(), key=lambda x: x[1]))
    
    # correct P values
    n = len(ordered_p_values)
    for rank, p in enumerate(ordered_p_values):
        # note that rank is [0, n)!
        if method == 'BH':
            ordered_p_values[p] *= n / (rank + 1)
        elif method == 'Holm':
            ordered_p_values[p] *= (n - rank)
        elif method == 'Bonferroni':
            # a bit dumb but...
            ordered_p_values[p] *= n
        elif method == 'None':
            # dummy function
            pass

        # cap P values at 1
        ordered_p_values[p] = min(ordered_p_values[p], 1)
    
    # preserve monotonicity
    if method == 'Holm':
        # Holm: if rank r has higher p value than r+1, r+1 takes higher p value
        prev_p_value = 0
        for p in ordered_p_values:
            ordered_p_values[p] = max(ordered_p_values[p], prev_p_value)
            prev_p_value = ordered_p_values[p]
    elif method == 'BH':
        # BH (1995) + YB (1999) -- do not confuse this with BY (2001):
        #   if rank r has higher p value than r+1, r takes smaller p value
        prev_p_value = 1
        for p in reversed(ordered_p_values):
            ordered_p_values[p] = min(ordered_p_values[p], prev_p_value)
            prev_p_value = ordered_p_values[p]
    
    return ordered_p_values
    
if __name__ == '__main__':
    test = [0.01, 0.011, 0.02, 0.021, 0.022, 0.022, 0.023, 0.024, 0.025, 0.1]
    test_dict = {n:x for n,x in enumerate(test)}
    print (test_dict)
    print (correct_p_values(test_dict, 'BH'))
    print (correct_p_values(test_dict, 'Bonferroni'))
    print (correct_p_values(test_dict, 'Holm'))