#ifndef JETS_EVENT_H
#define JETS_EVENT_H

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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

//
// class declaration
//

class JetsEvent{
  public:
    JetsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::View<pat::Jet>>&);
    JetsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<std::vector<reco::PFJet>>&);
    ~JetsEvent();
    std::vector<const pat::Jet*> GetJets();
    std::vector<const reco::PFJet*> GetJetsPF();

  private:
    std::vector<const pat::Jet*> patjetslist;
    std::vector<const reco::PFJet*> pfjetslist;

};

#endif
