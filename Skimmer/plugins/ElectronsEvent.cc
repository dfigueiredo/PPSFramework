// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/ElectronsEvent
// Class:      ElectronsEvent
// 
/**\class ElectronsEvent ElectronsEvent.cc CTPPSAnalysisCode/ElectronsEvent/plugins/ElectronsEvent.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

#include "ElectronsEvent.h"

// Retrieve Electrons
std::vector<const pat::Electron*> ElectronsEvent::GetElectrons(){
  return electronslist;
}

// Class Definition, loop over all electrons events and fill vector.
ElectronsEvent::ElectronsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::View<pat::Electron>>& electronsToken_){

  electronslist.clear();

  try{

    // Get the Jet collection from the event
    edm::Handle<edm::View<pat::Electron> > elecColl; // PAT
    iEvent.getByToken( electronsToken_, elecColl );

    for (unsigned int i = 0; i < elecColl->size(); ++i ) {
      const pat::Electron* ele = &((*elecColl)[i]);
      electronslist.push_back(ele);
    }
    LogDebug( "Electrons Event" ) << "Passed Loop on electrons";
  }catch(...){
    LogDebug( "Electrons Event" ) << "collection not present in the sample";
  }

}
