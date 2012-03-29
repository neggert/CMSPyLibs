#! /usr/bin/env python

from DataFormats.FWLite import Events
from DataFormats.FWLite import Handle

from ..cmsutilities import get_PF_isolation, get_list_from_handle

ele_iso_cut = 0.17
muon_iso_cut = 0.2
jet_btag_cut = 1.7
jet_btag = "trackCountingHighEffBJetTags"

class CMSEvent(object):
    """General class to provide a nice python interface to a CMS event in FWLite"""
    def __init__(self, electrons, muons, jets, met):
        """Initialize with collections of objects, except MET, which is a single object"""
        self.electrons = electrons
        self.muons = muons
        self.jets = jets
        self.met = met
    
    def __repr__(self):
        """Print information about the event"""
        out = ""
        
        out += str(len(self.get_electrons()))+" electrons:\n"
        for i, electron in enumerate(self.get_electrons()):
            out+="Electron "+str(i)+"\tPt: "+str(electron.pt())+"\tEta: "+str(electron.eta())+"\tPhi: "+str(electron.phi())+"\n"
        
        out += "\n"+str(len(self.get_muons()))+" muons:\n"
        for i, muon in enumerate(self.get_muons()):
            out+="Muon "+str(i)+"\tPt: "+str(muon.pt())+"\tEta: "+str(muon.eta())+"\tPhi: "+str(muon.phi())+"\n"
        
        out += "\n"+str(len(self.get_jets()))+" jets:\n"
        for i, jet in enumerate(self.get_jets()):
            out+="Jet "+str(i)+"\tPt: "+str(jet.pt())+"\tEta: "+str(jet.eta())+"\tPhi: "+str(jet.phi())+"\n"
        
        met = self.get_met()
        out+="\nMET \tPt: "+str(met.pt())+"\tEta: "+str(met.eta())+"\tPhi: "+str(met.phi())+"\n"

        return out
        
    def get_electrons(self):
        """Return a list of the electrons in the event"""
        return self.electrons
    
    def get_muons(self):
        """Return a list of the muons in the event"""
        return self.muons
    
    def get_jets(self):
        """Return a list of the jets in the event"""
        return self.jets
    
    def get_met(self):
        """Return the MET object from the event"""
        return self.met
    
    def get_leptons(self):
        """Return a list containing the electrons and muons in the Event"""
        return self.get_electrons()+self.get_muons()
    
    def set_electrons(self, electrons):
        """Set the electrons for the event"""
        self.electrons = electrons
    
    def set_muons(self, muons):
        """Set the muons for the event"""
        self.muons = muons
    
    def set_jets(self, jets):
        """Set the jets for the event"""
        self.jets = jets
    
    def set_met(self, met):
        """Set the MET for the event"""
        self.met = met
    
    def get_isolated_electrons(self):
        """Return the electrons with relative PF isolation < ele_iso_cut"""
        return [e for e in self.get_electrons() if get_PF_isolation(e) < ele_iso_cut]
    
    def get_nonisolated_electrons(self):
        """Return the electrons with relative PF isolation > ele_iso_cut"""
        return [e for e in self.get_electrons() if get_PF_isolation(e) > ele_iso_cut]
    
    def get_isolated_muons(self):
        """Return the muons with relative PF isolation < muons_iso_cut"""
        return [m for m in self.get_muons() if get_PF_isolation(m) < muon_iso_cut]
    
    def get_nonisolated_muons(self):
        """Return the muons with relative PF isolation > muons_iso_cut"""
        return [m for m in self.get_muons() if get_PF_isolation(m) > muon_iso_cut]
    
    def get_isolated_leptons(self):
        """Return the leptons that pass the relative isolation cuts"""
        return self.get_isolated_electrons()+self.get_isolated_muons()
    
    def get_nonisolated_leptons(self):
        """Return the leptons that fail the relative isolation cuts"""
        return self.get_nonisolated_electrons()+self.get_nonisolated_muons()
    
    def get_tagged_jets(self):
        """Return the jets that pass the b-tagging cuts defined in the module"""
        return [j for j in self.get_jets() if j.bDiscriminator(jet_btag) > jet_btag_cut]
    
    def get_untagged_jets(self):
        """Return the jets that fail the b-tagging cuts defined in the module"""
        return [j for j in self.get_jets() if j.bDiscriminator(jet_btag) < jet_btag_cut]

class CMSEventGetter(object):
    """The main purpose of this class is to supply a generator (events) which feeds CMSEvents from a file.
    The rest of the class is to provide an interface to change parameters"""
    _parent = CMSEvent

    def __init__(self, files):
        """Initialize using files and default object collection names"""
        self.files = files
        self.electron_collection = "selectedPatElectronsPFlow"
        self.muon_collection = "selectedPatMuonsPFlow"
        self.jet_collection = "selectedPatJetsPFlow"
        self.met_collection = "patMETsPFlow"

    
    def set_electron_collection(self, collection_name):
        """Set the name of the electron collection"""
        self.electron_collection = collection_name
    
    def set_muon_collection(self, collection_name):
        """Set the name of the muon collection"""
        self.muon_collection = collection_name
    
    def set_jet_collection(self, collection_name):
        """Set the name of the jet collection"""
        self.jet_collection = collection_name
    
    def set_met_collection(self, collection_name):
        """Set the name of the MET collection"""
        self.met_collection = collection_name
    
    def make_event(self, fwlite_event, handles):
        """Create a CMSEvent from the input FWLite event"""

        ele_handle, mu_handle, jet_handle, met_handle = handles

        electrons = get_list_from_handle(fwlite_event, ele_handle, self.electron_collection)
        muons     = get_list_from_handle(fwlite_event, mu_handle, self.muon_collection)
        jets      = get_list_from_handle(fwlite_event, jet_handle, self.jet_collection)
        met       = get_list_from_handle(fwlite_event, met_handle, self.met_collection)[0]

        
        event = self._parent(electrons, muons, jets, met)
        del electrons, muons, jets, met

        return event
    
    def passes_cuts(self, cms_event):
        """Returns whether or not the event passes the cuts. Mostly intended to be implemented
        in inherited classes."""
        return True
    
    def events(self):
        """Generator for events"""
        fwlite_events = Events(self.files)

        ele_handle = Handle("std::vector<pat::Electron>")
        mu_handle = Handle("std::vector<pat::Muon>")
        jet_handle = Handle("std::vector<pat::Jet>")
        met_handle = Handle("std::vector<pat::MET>")
        handles = (ele_handle, mu_handle, jet_handle, met_handle)
        for fwlite_event in fwlite_events:
            event = self.make_event(fwlite_event, handles)
            if self.passes_cuts(event):
                yield event
            del event

def test():
    """Run some tests to make sure everything works"""
    files = ["/Users/nic/cms/August11MC/Signal/BsmMassesSkim_Summer11_Sync.root"]
    getter = CMSEventGetter(files)
    iterator = getter.events()
    for event in iterator :
        print event

def test_dumb():
    """Generator for events"""
    files = ["/Users/nic/cms/August11MC/Signal/BsmMassesSkim_Summer11_Sync.root"]

    fwlite_events = Events(files)
    handle = Handle("std::vector<pat::MET>")

    for fwlite_event in fwlite_events:
        fwlite_event.getByLabel("patMETsPFlow", handle)
        met = handle.product()[0]
        print met.pt()
        