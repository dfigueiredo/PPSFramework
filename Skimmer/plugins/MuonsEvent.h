#ifndef MUONS_EVENT_H
#define MUONS_EVENT_H

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

// system include files
#include <memory>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Muons collection
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//
// class declaration
//

class MuonsEvent{
  public:
    MuonsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::View<pat::Muon>>&);
    ~MuonsEvent();
    std::vector<const pat::Muon*> GetMuons();

  private:
    std::vector<const pat::Muon*> muonslist;

};

#endif
