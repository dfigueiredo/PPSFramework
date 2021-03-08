// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/ParticleFlowEvent
// Class:      ParticleFlowEvent
// 
/**\class ParticleFlowEvent ParticleFlowEvent.cc CTPPSAnalysisCode/ParticleFlowEvent/plugins/ParticleFlowEvent.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

#include "ParticleFlowEvent.h"

// Retrieve PF
std::vector<const reco::PFCandidate*> ParticleFlowEvent::GetParticleFlow(){
  return pflist;
}

std::vector<const pat::PackedCandidate*> ParticleFlowEvent::GetPackedFlow(){
  return packedlist;
}

// Class Definition, loop over all particle flow events and fill vector.
ParticleFlowEvent::ParticleFlowEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<std::vector<reco::PFCandidate>>& pfToken_){

  pflist.clear();

  try{
    // Get the PF collection from the event
    edm::Handle<std::vector<reco::PFCandidate>> pfColl; // PAT
    iEvent.getByToken( pfToken_, pfColl );

    for (unsigned int i = 0; i < pfColl->size(); ++i ) {
      const reco::PFCandidate* pf = &((*pfColl)[i]);
      pflist.push_back(pf);
    }
    LogDebug( "Particle Flow Event" ) << "Passed Loop on PF";
  }catch(...){
    LogDebug( "Particle Flow Event" ) << "collection not present in the sample";
  }

}

// Class Definition, loop over all particle flow events and fill vector.
ParticleFlowEvent::ParticleFlowEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<std::vector<pat::PackedCandidate>>& packedToken_){

  packedlist.clear();

  try{
    // Get the PF collection from the event
    edm::Handle<std::vector<pat::PackedCandidate>> packedColl; // PAT
    iEvent.getByToken( packedToken_, packedColl );

    for (unsigned int i = 0; i < packedColl->size(); ++i ) {
      const pat::PackedCandidate* packed = &((*packedColl)[i]);
      packedlist.push_back(packed);
    }
    LogDebug( "Packed Flow Event" ) << "Passed Loop on Packed Candidate";
  }catch(...){
    LogDebug( "Packed Flow Event" ) << "collection not presented in the sample";
  }

}

