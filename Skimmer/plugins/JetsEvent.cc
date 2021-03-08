// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/JetsEvent
// Class:      JetsEvent
// 
/**\class JetsEvent JetsEvent.cc CTPPSAnalysisCode/JetsEvent/plugins/JetsEvent.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

#include "JetsEvent.h"

// Retrieve Jets
std::vector<const pat::Jet*> JetsEvent::GetJets(){
  return patjetslist;
}

// Retrieve Jets
std::vector<const reco::PFJet*> JetsEvent::GetJetsPF(){
  return pfjetslist;
}

// Class Definition, loop over all jets events and fill vector.
JetsEvent::JetsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::View<pat::Jet>>& jetsToken_){

  patjetslist.clear();

  try{

    // Get the Jet collection from the event
    edm::Handle<edm::View<pat::Jet> > jetColl; // PAT
    iEvent.getByToken( jetsToken_, jetColl );

    for (unsigned int i = 0; i < jetColl->size(); ++i ) {
      const pat::Jet* jet = &((*jetColl)[i]);
      patjetslist.push_back(jet);
    }
    LogDebug( "Jets Event" ) << "Passed Loop on jets";
  }catch(...){
    LogDebug( "Jets Event" ) << "collection not present in the sample";
  }

}

// Class Definition, loop over all jets events and fill vector.
JetsEvent::JetsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<std::vector<reco::PFJet>>& jetsToken_){

  pfjetslist.clear();

  try{

    // Get the Jet collection from the event
    edm::Handle<reco::PFJetCollection > jetColl; // PAT
    iEvent.getByToken( jetsToken_, jetColl );

    for (unsigned int i = 0; i < jetColl->size(); ++i ) {
      const reco::PFJet* jet = &((*jetColl)[i]);
      pfjetslist.push_back(jet);
    }
    LogDebug( "Jets Event" ) << "Passed Loop on jets";
  }catch(...){
    LogDebug( "Jets Event" ) << "collection not present in the sample";
  }

}
