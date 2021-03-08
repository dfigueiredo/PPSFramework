// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/PhotonsEvent
// Class:      PhotonsEvent
// 
/**\class PhotonsEvent PhotonsEvent.cc CTPPSAnalysisCode/PhotonsEvent/plugins/PhotonsEvent.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

#include "PhotonsEvent.h"

// Retrieve Photons
std::vector<const pat::Photon*> PhotonsEvent::GetPhotons(){
  return photonslist;
}

// Class Definition, loop over all photons events and fill vector.
PhotonsEvent::PhotonsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::View<pat::Photon>>& photonsToken_){

  photonslist.clear();

  try{
    // Get the Jet collection from the event
    edm::Handle<edm::View<pat::Photon> > photonColl; // PAT
    iEvent.getByToken( photonsToken_, photonColl );

    for (unsigned int i = 0; i < photonColl->size(); ++i ) {
      const pat::Photon* photon = &((*photonColl)[i]);
      photonslist.push_back(photon);
    }
    LogDebug( "Photons Event" ) << "Passed Loop on photons";
  }catch(...){
    LogDebug( "Photons Event" ) << "collection not present in the sample";
  }

}
