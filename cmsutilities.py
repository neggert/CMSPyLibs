#! /usr/bin/env python

from DataFormats.FWLite import Handle
from ROOT import TLorentzVector
from pylab import pi, sin, cos, sqrt
import ROOT
from copy import copy

def angle_0_2pi (phi) :
    while phi > 2*pi : phi -= 2*pi
    while phi < 0 : phi += 2*pi
    return phi

def get_list_from_handle( event, handle, collection_name) :
    event.getByLabel(collection_name, handle)
    return [item for item in handle.product()] # convert to a python list

# these probably cause memory leaks
def get_electrons_from_event(event, collection_name):
    """Return electrons in the collection named collection_name from event"""
    handle = Handle("std::vector<pat::Electron>")
    event.getByLabel(collection_name, handle)
    return [item for item in handle.product()] # convert to a python list

def get_muons_from_event(event, collection_name):
    """Return muons in the collection named collection_name from event"""
    handle = Handle("std::vector<pat::Muon>")
    event.getByLabel(collection_name, handle)
    return [item for item in handle.product()]

def get_jets_from_event(event, collection_name):
    """Return jets in the collection named collection_name from event"""
    handle = Handle("std::vector<pat::Jet>")
    event.getByLabel(collection_name, handle)
    return [item for item in handle.product()]

def get_met_from_event(event, collection_name):
    """Return met in the collection named collection_name from event"""
    handle = Handle("std::vector<pat::MET>")
    event.getByLabel(collection_name, handle)
    return handle.product()[0]

def get_TLorentzVector(cms_object):
    """Return a TLorentzVector made from the input object"""
    v = TLorentzVector(cms_object.px(), cms_object.py(), cms_object.pz(), cms_object.energy())
    return v

def get_effective_area(eta):
    """Return effective area for 2011 data"""
    aeta = abs(eta)
    if aeta < 1.0:
        return 0.13
    elif aeta < 1.479:
        return 0.14
    elif aeta < 2.0:
        return 0.07
    elif aeta < 2.2:
        return 0.09
    elif aeta < 2.3:
        return 0.11
    elif aeta < 2.4:
        return 0.11
    else:
        return 0.14

def get_PF_isolation(lepton):
    """Return the particle flow isolation for the input object"""
    iso = lepton.chargedHadronIso()+lepton.neutralHadronIso()+lepton.photonIso()
    return iso/lepton.pt()

def get_EA_PF_isolation(lepton, rho):
    """Return the particle flow isolation with the effective area correction for the input object"""
    ea = get_effective_area(lepton.eta())
    iso = lepton.chargedHadronIso()+max([lepton.neutralHadronIso()+lepton.photonIso()-ea*rho, 0.])
    return iso/lepton.pt()

def get_delta_beta_PF_isolation(lepton):
      iso_struct = lepton.pfIsolationR03()
      iso = iso_struct.sumChargedHadronPt+max([0.,iso_struct.sumNeutralHadronEt+iso_struct.sumPhotonEt-0.5*iso_struct.sumPUPt]);
      return iso/lepton.pt()

def get_5_X_isolation(lepton, rho):
    """Check particle type and get the appropriate isolation"""
    if abs(lepton.pdgId()) == 11:
        return get_EA_PF_isolation(lepton, rho)
    elif abs(lepton.pdgId()) == 13:
        return get_delta_beta_PF_isolation(lepton)
    else:
        raise RuntimeError("Trying to calculate isolation for non-lepton")

def get_upstream_phi_res(p4Up, sigMatrix) :
    """Get the upstream phi resoluiton from it's four-vector and significance matrix"""
    # rotate so that upstream pt is pointing in x direction
    phi  = -p4Up.Phi()
    R = ROOT.TMatrixD(2,2)
    R[0][0] = R[1][1]= cos(phi)
    R[0][1] = sin(phi)
    R[1][0] = -sin(phi)
    RI = copy(R) # R inverse
    RI.Invert()
    rotSig = RI*sigMatrix*R

    return sqrt(rotSig(1,1)/p4Up.Pt()**2)

def get_cov_matrix(p4, ptRes, phiRes):
    """Make a covariance matrix for an object"""
    # rotation matrix
    phi  = p4.Phi()
    R = ROOT.TMatrixD(2,2)
    R[0][0] = R[1][1]= cos(phi)
    R[0][1] = sin(phi)
    R[1][0] = -sin(phi)
    RI = copy(R) # R inverse
    RI.Invert()

    # significance matrix in frame where pt points in x direction
    s = ROOT.TMatrixD(2,2)
    s[0][0] = ptRes**2
    s[0][1] = s[1][0] = 0.
    s[1][1] = p4.Et()**2*phiRes**2

    # print s[0][0], s[0][1]
    # print s[1][0], s[1][1]

    # rotate matrix into x-y frame
    return RI*s*R
