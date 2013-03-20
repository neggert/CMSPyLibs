import matplotlib.pyplot as plt
import numpy as np
import inspect

def hist_errorbars( data, xerrs=True, color="k", *args, **kwargs) :
    """Plot a histogram with error bars. Accepts any kwarg accepted by either numpy.histogram or pyplot.errorbar"""
    # pop off normed kwarg, since we want to handle it specially
    norm = False
    if 'normed' in kwargs.keys() :
        norm = kwargs.pop('normed')

    # retrieve the kwargs for numpy.histogram
    histkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(np.histogram).args :
            histkwargs[key] = value

    histvals, binedges = np.histogram( data, **histkwargs )

    weighted=False
    if "weights" in histkwargs.keys():
        hist_weights = np.asarray(histkwargs.pop('weights'))
        weighted = True

    if weighted:
        sumw2, binedges = np.histogram(data, weights=hist_weights**2, **histkwargs)
        yerrs = np.sqrt(sumw2)
    else:
        yerrs = np.sqrt(histvals.tolist()) # no effing idea why tolist is necessary


    if norm :
        nevents = float(sum(histvals))
        binwidth = (binedges[1]-binedges[0])
        histvals = histvals/nevents/binwidth
        yerrs = yerrs/nevents/binwidth

    bincenters = (binedges[1:]+binedges[:-1])/2

    if xerrs :
        xerrs = (binedges[1]-binedges[0])/2
    else :
        xerrs = None

    # retrieve the kwargs for errorbar
    ebkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(plt.errorbar).args :
            histkwargs[key] = value
    out = plt.errorbar(bincenters, histvals, yerrs, xerrs, fmt=".", color=color, **ebkwargs)

    if 'log' in kwargs.keys() :
        if kwargs['log'] :
            plt.yscale('log')

    if 'range' in kwargs.keys() :
        plt.xlim(*kwargs['range'])

    return (histvals, yerrs, out)

def hist_ratio( data_num, data_denom, weight_denom, weight_num=None, xerrs=True, color="k", *args, **kwargs) :
    """Plot a histogram with error bars. Accepts any kwarg accepted by either numpy.histogram or pyplot.errorbar"""
    # pop off normed kwarg, since we want to handle it specially
    norm = False
    if 'normed' in kwargs.keys() :
        norm = kwargs.pop('normed')

    # retrieve the kwargs for numpy.histogram
    histkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(np.histogram).args :
            histkwargs[key] = value

    histvals1, binedges = np.histogram( data_num, weights=weight_num, **histkwargs )
    yerrs1 = np.sqrt(histvals1.tolist()) # no effing idea why tolist is necessary
    if weight_num is not None:
        sumw2, binedges = np.histogram( data_num, weights=np.asarray(weight_num)**2, **histkwargs)
        yerrs1 = np.sqrt(sumw2) # no effing idea why tolist is necessary

    histvals2, binedges = np.histogram( data_denom, weights=weight_denom, **histkwargs )
    yerrs2 = np.sqrt(histvals2.tolist()) # no effing idea why tolist is necessary
    if weight_denom is not None:
        sumw2, binedges = np.histogram( data_denom, weights=np.asarray(weight_denom)**2, **histkwargs)
        yerrs2 = np.sqrt(sumw2)

    if norm:
        yerrs1 = yerrs1/sum(histvals1)
        yerrs2 = yerrs2/sum(histvals2)
        histvals1 = 1.*histvals1/sum(histvals1)
        histvals2 = 1.*histvals2/sum(histvals2)

    ratio = histvals1/histvals2
    ratio_err = np.sqrt((yerrs1/histvals2)**2+(histvals1/histvals2**2*yerrs2)**2)

    bincenters = (binedges[1:]+binedges[:-1])/2

    if xerrs :
        xerrs = (binedges[1]-binedges[0])/2
    else :
        xerrs = None

    ratio_err = ratio_err[ratio>0.]
    bincenters = bincenters[ratio>0.]
    ratio = ratio[ratio>0.]
    

    # retrieve the kwargs for errorbar
    ebkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(plt.errorbar).args :
            histkwargs[key] = value
    out = plt.errorbar(bincenters, ratio, ratio_err, xerrs, fmt=".", color=color, **ebkwargs)

    if 'log' in kwargs.keys() :
        if kwargs['log'] :
            plt.yscale('log')

    if 'range' in kwargs.keys() :
        plt.xlim(*kwargs['range'])


    return out
