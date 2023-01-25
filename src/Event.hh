#ifndef FDSI_EVENT_HH
#define FDSI_EVENT_HH

#include "DEGAi.hh"
#include "LaBr.hh"

#include "libpixie/event.hh"
#include "libpixie/reader.hh"

namespace FDSi {
  class Event {
  public:
    DEGAi degai;  //DEGAi array
    LaBrArray labr;  //DEGAi array
    unsigned long long time; //first time
    
    Event(DEGAi d, LaBrArray l) { degai = d; labr = l; }
    const PIXIE::Event& Set(const PIXIE::Reader &reader, int eventInd); 
  };
}

#endif
