import matplotlib.pyplot as plt
import numpy as np
import inspect

def hist_errorbars( data, *args, **kwargs) :
    """Plot a histogram with error bars. Accepts any kwarg accepted by either numpy.histogram or pyplot.errorbar"""
    # retrieve the kwargs for numpy.histogram
    histkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(np.histogram) :
            histkwargs[key] = value

    histvals, binedges = np.histogram( data, **histkwargs )

    bincenters = (binedges[1:]+binedges[:-1])/2

    yerrs = np.sqrt(histvals)
    xerrs = (binedges[1]-binedges[0])/2

    # retrieve the kwargs for errorbar
    ebkwargs = {}
    for key, value in kwargs.iteritems() :
        if key in inspect.getargspec(plt.errorbar) :
            histkwargs[key] = value
    plt.errorbar(bincenters, histvals, yerrs, xerrs, fmt=".", **ebkwargs)

    if 'log' in kwargs.keys() :
        if kwargs['log'] :
            plt.yscale('log')