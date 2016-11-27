#!/usr/bin/env python3

"""
> calc_fishers_exact.py <

Does what it says on the tin - provides easier-to-remember defs that help 
calculate P values via Fisher's exact test.

Requires scipy.
"""
import scipy.stats

def fishers_exact_cdf(k, m, n, N):
    """
    Used to calculate p-values (e.g. tissue is 'brain', AffyCall is 'Down').
    k = genes in brain of AffyCall 'Down' targeted by that certain miR
    m = total genes in brain with AffyCall 'Down'
    n = all genes targeted by that miR
    N = all genes
    Take note that calculations should only consider genes that are common
    across both genesets - inaccuracies occur when N is larger than it should!
    
    Note: cdf = cumulative distribution function
    """
    # Fisher's exact test: table layout
    #   k     (m-k)        m
    #   (n-k) (N-n-m+k)    N-m
    #   n     (N-n)        N
    w = k
    x = m - k
    y = n - k
    z = N - n - m + k
    
    return scipy.stats.fisher_exact([[w, x], [y, z]])[-1]
