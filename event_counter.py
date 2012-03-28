class EventCounter :
    """Class that counts events"""
    maxEvents = -1
    eventNum = 1
    printOut = True
    printRate = 1
    def __init__(self) :
        self.maxEvents = -1
        self.eventNum = 1
        self.printOut = True
        self.constPrintRate = False
        self.printRate = 1
    def SetPrintRate(self, rate) :
        self.printRate = rate
    def SetPrinting(self, onoff) :
        self.printOut = onoff
    def SetConstPrintRate(self, onoff) :
        self.constPrintRate = onoff
    def SetMaxEvents(self, maxEvents) :
        self.maxEvents = maxEvents
    def Reset(self) :
        self.eventNum = 0
    def Increment(self) :
        """Increments event number. Will also print event number
        at regular intervals if printOut is True. Returns True if
        eventNum < maxEvents and False otherwise"""
        if (self.eventNum >= 10*self.printRate and not self.constPrintRate) :
            self.printRate *= 10
        if (self.eventNum%self.printRate == 0 and self.printOut) :
            print "Event number", self.eventNum
        self.eventNum += 1
        if (self.eventNum <= self.maxEvents or self.maxEvents < 0) :
            return True
        else :
            raise StopIteration
            return False
            

