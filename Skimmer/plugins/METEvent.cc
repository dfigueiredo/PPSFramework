// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/METEvent
// Class:      METEvent
// 
/**\class METEvent METEvent.cc CTPPSAnalysisCode/METEvent/plugins/METEvent.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

#include "METEvent.h"

// Retrieve MET
edm::View<pat::MET>::const_iterator METEvent::GetMET(){
  return metlist;
}

// Class Definition, loop over all muons events and fill vector.
METEvent::METEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::View<pat::MET>>& metToken_){

  try{

    edm::Handle<edm::View<pat::MET> > MET;
    iEvent.getByToken( metToken_, MET);

    const edm::View<pat::MET>* metColl = MET.product();
    edm::View<pat::MET>::const_iterator met = metColl->begin();
    metlist = met;

    LogDebug( "MET Event" ) << "Passed Loop on muons";
  }catch(...){
    LogDebug( "MET Event" ) << "collection not present in the sample";
  }

}
