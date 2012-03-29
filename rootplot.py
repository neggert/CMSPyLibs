import ROOT
from array import array

"""
Wrapper functions to produce ROOT histograms and graphs
using an interface similar to matplotlib
"""

def plot_hist_ROOT( data, bins=20, range=(0., 100.), name="hist") :
    """Plot a TH1D of data"""
    hist = ROOT.TH1D(name, name, bins, range[0], range[1])
    for d in data :
        hist.Fill(d)
    hist.SetMinimum(0.)
    hist.Draw()

    return hist

def plot_hist2d_ROOT( xdata, ydata, xbins=20, xrange=(0., 100.), ybins=20, yrange=(0., 100.), name="hist2d") :
    """Plot a TH2D of data"""
    hist = ROOT.TH2D(name, name, xbins, xrange[0], xrange[1], ybins, yrange[0], yrange[1])
    for x,y in zip(xdata,ydata) :
        hist.Fill(x,y)
    hist.Draw("COLZ")

    return hist

def plot_graph_ROOT( x, y ) :
    """Plot a TGraph of x vs y"""
    g = ROOT.TGraph(len(x), array("f", x), array("f", y))
    g.SetMarkerStyle(7)
    # g.SetMarkerSize(4)
    g.Draw("pa")

    return g