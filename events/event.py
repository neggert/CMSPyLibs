#! /usr/bin/env python

from DataFormats.FWLite import Events
from DataFormats.FWLite import Handle

from ..cmsutilities import get_PF_isolation, get_list_from_handle

from fnmatch import fnmatch

class EventID(object) :
    """Class to contain run, luminosity block, and event number which uniquely identify an event"""
    def run_number():
        doc = "The run_number property."
        def fget(self):
            return self._run_number
        def fset(self, value):
            self._run_number = value
        def fdel(self):
            del self._run_number
        return locals()
    run_number = property(**run_number())
    def luminosity_block():
        doc = "The luminosity_block property."
        def fget(self):
            return self._luminosity_block
        def fset(self, value):
            self._luminosity_block = value
        def fdel(self):
            del self._luminosity_block
        return locals()
    luminosity_block = property(**luminosity_block())

    def event_number():
        doc = "The event_number property."
        def fget(self):
            return self._event_number
        def fset(self, value):
            self._event_number = value
        def fdel(self):
            del self._event_number
        return locals()
    event_number = property(**event_number())

    def __init__(self, event) :
        """Initialize from an fwlite event"""
        self.run_number = event.eventAuxiliary().run()
        self.luminosity_block = event.eventAuxiliary().luminosityBlock()
        self.event_number = event.eventAuxiliary().event()

    def __repr__(self) :
        return "Run {0}, Lumi {1}, Event {2}".format(self.run_number, self.luminosity_block, self.event_number)

    def __hash__(self) :
        return hash(str(self.run_number)+str(self.luminosity_block)+str(self.event_number))

    def __eq__(self, other) :
        return (self.run_number == other.run_number) and (self.luminosity_block == other.luminosity_block)\
                and (self.event_number == other.event_number)


class CMSEvent(object):
    """General class to provide a nice python interface to a CMS event in FWLite"""
    ele_iso_cut = 0.17
    muon_iso_cut = 0.2
    jet_btag_cut = 1.7
    jet_btag = "trackCountingHighEffBJetTags"

    def __init__(self, eventID, vertices, genparticles, electrons, muons, jets, met, triggerresults, triggernames, metadata={}):
        """Initialize with collections of objects, except MET, which is a single object"""
        self.eventID = eventID
        self.vertices = vertices
        self.genparticles = genparticles
        self.electrons = electrons
        self.muons = muons
        self.jets = jets
        self.met = met
        self.triggerresults = triggerresults
        self.triggernames = triggernames
        self.metadata = metadata

    def __repr__(self):
        """Print information about the event"""
        out = repr(self.eventID)

        out += "\n\n"

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

    def set_genparticles(self, genparticles):
        """Set the MET for the event"""
        self.genparticles = genparticles

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
        return [j for j in self.get_jets() if j.bDiscriminator(self.jet_btag) > self.jet_btag_cut]

    def get_untagged_jets(self):
        """Return the jets that fail the b-tagging cuts defined in the module"""
        return [j for j in self.get_jets() if j.bDiscriminator(self.jet_btag) < self.jet_btag_cut]

    def get_vertices(self):
        """Return the vertex collection"""
        return self.vertices

    def get_genparticles(self):
        """Return the genparticle collection"""
        return self.genparticles

    def passes_HLT(self, hltpath):
        """Check to see if the event passes the given HLT path. Unix filename style wildcards are allowed"""
        # first find triggers which pass the wildcard
        triggers = [i for i in range(len(self.triggernames)) if fnmatch(self.triggernames[i], hltpath)] # might be a nicer way
        for t in triggers:
            if self.triggerresults.accept(t):
                return True
        return False


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
        self.vertex_collection = "goodOfflinePrimaryVertices"
        self.genparticle_collection = "genParticles"
        self.gen_handle = Handle("std::vector<reco::GenParticle>")
        self.ele_handle = Handle("std::vector<pat::Electron>")
        self.mu_handle = Handle("std::vector<pat::Muon>")
        self.jet_handle = Handle("std::vector<pat::Jet>")
        self.met_handle = Handle("std::vector<pat::MET>")
        self.vertex_handle = Handle("std::vector<reco::Vertex>")
        self.pu_handle = Handle("std::vector<PileupSummaryInfo>")
        self.vdouble_handle = Handle("std::vector<double>")
        self.double_handle = Handle("double")
        self.triggerresults_handle = Handle("edm::TriggerResults")
        self.do_genparticles = True
        self.do_PU = False
        self.do_SMS = False


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

    def get_num_pu_vertices(self, fwlite_event):
        """ Get the number of PU vertices"""
        try :
            puInfos = get_list_from_handle(fwlite_event, self.pu_handle, "addPileupInfo")
        except Exception :
            return None
        return sum((pvi.getPU_NumInteractions() for pvi in puInfos if pvi.getBunchCrossing()==0))

    def get_true_pu_vertices(self, fwlite_event):
        try :
            puInfos = get_list_from_handle(fwlite_event, self.pu_handle, "addPileupInfo")
        except Exception :
            return None
        return puInfos[0].getTrueNumInteractions()

    def get_trigger_results(self, fwlite_event):
        """Get the trigger names and results"""
        fwlite_event.getByLabel("TriggerResults", "", "HLT", self.triggerresults_handle)
        tr = self.triggerresults_handle.product()
        tn = fwlite_event.object().triggerNames(tr)
        tn_l = []
        for i in range(tn.triggerNames().size()):
            tn_l.append(tn.triggerName(i))
        return tr, tn_l

    def get_sms_params(self, fwlite_event):
        """ Get the SMS model parameters"""
        sms_params = get_list_from_handle( fwlite_event, self.double_handle,
                       ['susyScanP1', ''])
        return sms_params

    def get_rho(self, fwlite_event):
        """Get the rho parameter"""
        fwlite_event.getByLabel( ["kt6PFJetsForIsolation", "rho"], self.double_handle )
        return self.double_handle.product()[0]

    def make_event(self, fwlite_event):
        """Create a CMSEvent from the input FWLite event"""

        if self.do_genparticles:
            gen_particles = get_list_from_handle(fwlite_event, self.gen_handle, self.genparticle_collection)
        else :
            gen_particles = []
        electrons = get_list_from_handle(fwlite_event, self.ele_handle, self.electron_collection)
        muons     = get_list_from_handle(fwlite_event, self.mu_handle, self.muon_collection)
        jets      = get_list_from_handle(fwlite_event, self.jet_handle, self.jet_collection)
        met       = get_list_from_handle(fwlite_event, self.met_handle, self.met_collection)[0]
        vertices  = get_list_from_handle(fwlite_event, self.vertex_handle, self.vertex_collection)

        eventID = EventID(fwlite_event)

        trigger_results, trigger_names = self.get_trigger_results(fwlite_event)


        event = self._parent(eventID, vertices, gen_particles, electrons, muons, jets, met, trigger_results, trigger_names)
        del gen_particles, electrons, muons, jets, met, vertices, trigger_results

        if self.do_PU :
            event.metadata['num_pu_vertices'] = self.get_num_pu_vertices( fwlite_event )
            event.metadata['num_true_pu_vertices'] = self.get_true_pu_vertices( fwlite_event )
        if self.do_SMS :
            event.metadata['modelParams'] = self.get_sms_params(fwlite_event)

        event.metadata['rho'] = self.get_rho(fwlite_event)

        return event

    def passes_cuts(self, cms_event):
        """Returns whether or not the event passes the cuts. Mostly intended to be implemented
        in inherited classes."""
        return True

    def events(self):
        """Generator for events"""
        fwlite_events = Events(self.files)

        for fwlite_event in fwlite_events:
            event = self.make_event(fwlite_event)
            if self.passes_cuts(event):
                yield event
            del event

def test():
    """Run some tests to make sure everything works"""
    files = ["root://osg-se.cac.cornell.edu//xrootd/path/cms/store/user/neggert/WWZNoGstarJets_8TeV-madgraph/OSDil_MCT_HCP2012_WWZNoGstarJets/47c3305489ab060ff201ab9227b9eec2/OSDil_MCT_HCP2012_1_1_o6a.root"]
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

