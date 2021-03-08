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

#include "GenMETEvent.h"

// Retrieve GenMET
edm::View<reco::GenMET>::const_iterator GenMETEvent::GetGenMET(){
  return genmetlist;
}

// Class Definition, loop over all muons events and fill vector.
GenMETEvent::GenMETEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::View<reco::GenMET>>& GenMETToken_){

  try{

    edm::Handle<edm::View<reco::GenMET> > GenMET;
    iEvent.getByToken( GenMETToken_, GenMET);

    const edm::View<reco::GenMET>* genmetColl = GenMET.product();
    edm::View<reco::GenMET>::const_iterator genmet = genmetColl->begin();
    genmetlist = genmet;

    LogDebug( "GenMET Event" ) << "Passed Loop on muons";
  }catch(...){
    LogDebug( "GenMET Event" ) << "collection not present in the sample";
  }

}
