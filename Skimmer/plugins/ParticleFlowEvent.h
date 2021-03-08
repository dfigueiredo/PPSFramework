#ifndef PARTICLEFLOW_EVENT_H
#define PARTICLEFLOW_EVENT_H

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


// system include files
#include <memory>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// PF collection
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"

//
// class declaration
//

class ParticleFlowEvent{
  public:
    ParticleFlowEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<std::vector<pat::PackedCandidate>>&);
    ParticleFlowEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<std::vector<reco::PFCandidate>>&);
    ~ParticleFlowEvent();
    std::vector<const reco::PFCandidate*> GetParticleFlow();
    std::vector<const pat::PackedCandidate*> GetPackedFlow();

  private:
    std::vector<const reco::PFCandidate*> pflist;
    std::vector<const pat::PackedCandidate*> packedlist;

};

#endif
