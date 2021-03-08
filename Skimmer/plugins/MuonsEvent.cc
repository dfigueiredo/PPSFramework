// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/MuonsEvent
// Class:      MuonsEvent
// 
/**\class MuonsEvent MuonsEvent.cc CTPPSAnalysisCode/MuonsEvent/plugins/MuonsEvent.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

#include "MuonsEvent.h"

// Retrieve Muons
std::vector<const pat::Muon*> MuonsEvent::GetMuons(){
  return muonslist;
}

// Class Definition, loop over all muons events and fill vector.
MuonsEvent::MuonsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::View<pat::Muon>>& muonsToken_){

  muonslist.clear();

  try{
    // Get the Jet collection from the event
    edm::Handle<edm::View<pat::Muon> > muonColl; // PAT
    iEvent.getByToken( muonsToken_, muonColl );

    for (unsigned int i = 0; i < muonColl->size(); ++i ) {
      const pat::Muon* mu = &((*muonColl)[i]);
      muonslist.push_back(mu);
    }
    LogDebug( "Muons Event" ) << "Passed Loop on muons";
  }catch(...){
    LogDebug( "Muons Event" ) << "collection not present in the sample";
  }

}
