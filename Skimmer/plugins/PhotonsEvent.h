#ifndef PHOTONS_EVENT_H
#define PHOTONS_EVENT_H

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

// system include files
#include <memory>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Photons collection
#include "DataFormats/PatCandidates/interface/Photon.h"

//
// class declaration
//

class PhotonsEvent{
  public:
    PhotonsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::View<pat::Photon>>&);
    ~PhotonsEvent();
    std::vector<const pat::Photon*> GetPhotons();

  private:
    std::vector<const pat::Photon*> photonslist;

};

#endif
