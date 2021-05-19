#ifndef GENPROTONS_EVENT_H 
#define GENPROTONS_EVENT_H

//// system include files
#include <memory>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//GenParticle collection
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//
// class declaration
//

class GenProtonsEvent{
  public:
      GenProtonsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<std::vector< reco::GenParticle>>&);
      ~GenProtonsEvent();
      std::vector<const reco::GenParticle*> GetGenProtons();
                
  private:
      std::vector<const reco::GenParticle*> genparticleslist;

};

#endif

