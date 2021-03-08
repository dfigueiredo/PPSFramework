#ifndef GENJETS_EVENT_H
#define GENJETS_EVENT_H

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

// system include files
#include <memory>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Jets/MET collection
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

//
// class declaration
//

class GenJetsEvent{
  public:
    GenJetsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::View<reco::GenJet>>&);
    //GenJetsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<std::vector<reco::PFJet>>&);
    ~GenJetsEvent();
    std::vector<const reco::GenJet*> GetGenJets();
    //std::vector<const reco::PFJet*> GetJetsPF();

  private:
    std::vector<const reco::GenJet*> genjetslist;
    //std::vector<const reco::PFJet*> pfjetslist;

};

#endif
