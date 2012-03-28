#! /usr/bin/env python

from .event import CMSEvent, CMSEventGetter
from  ..cmsutilities import get_TLorentzVector

class CMSDileptonEvent(CMSEvent):
    """Event class which guarantees the event has two leptons"""
    def __init__(self, electrons, muons, jets, met):
        """Make sure that event has at least two leptons"""
        if (len(electrons)+len(muons) < 2):
            raise ValueError('DileptonEvent must have at least two leptons')
        super(CMSDileptonEvent, self).__init__(electrons, muons, jets, met)
    
    def upstream_of_leptons(self):
        """Get 4-momentum upstream of dilepton pair and MET"""
        p4l1 = get_TLorentzVector(self.get_leptons()[0])
        p4l2 = get_TLorentzVector(self.get_leptons()[1])
        p4met = get_TLorentzVector(self.get_met())
        return -p4met-p4l1-p4l2

class CMSDileptonEventGetter(CMSEventGetter):
    """Provides a generator which feeds events from file(s). Skips events that don't have at least two leptons"""
    _parent = CMSDileptonEvent

    def passes_cuts(self, cms_event):
        if len(cms_event.get_leptons()) < 2 :
            return False
        else :
            return True

def test():
    """Run some tests to make sure everything works"""
    files = ["/Users/nic/cms/August11MC/Signal/BsmMassesSkim_Summer11_Sync.root"]
    getter = CMSDileptonEventGetter(files)
    iterator = getter.events()
    event = iterator.next()
    print event
    print "Upstream Phi:", event.upstream_of_leptons().Phi()
    return event


if __name__ == '__main__':
    test()