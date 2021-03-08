#ifndef PROTONS_EVENT_H
#define PROTONS_EVENT_H

// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/ProtonsEvent
// Class:      ProtonsEvent
// 
/**\class ProtonsEvent ProtonsEvent.cc CTPPSAnalysisCode/ProtonsEvent/plugins/ProtonsEvent.cc

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

// CT-PPS objects
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "DataFormats/CTPPSReco/interface/CTPPSPixelLocalTrack.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"
#include "DataFormats/ProtonReco/interface/ForwardProtonFwd.h"


//
// class declaration
//

class ProtonsEvent{
  public:
    ProtonsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite>>&);
    ProtonsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack>>&);
    ProtonsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondLocalTrack>>&);
    ProtonsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack>>&);
    ProtonsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<reco::ForwardProtonCollection>&);

    ~ProtonsEvent();
    std::vector<const CTPPSLocalTrackLite*> GetProtonsLite();
    std::vector<const reco::ForwardProton*> GetProtonsReco();
    std::vector<std::pair<const CTPPSPixelLocalTrack*, const CTPPSDetId>> GetPixels();
    std::vector<std::pair<const CTPPSDiamondLocalTrack*, const CTPPSDetId>> GetTiming();
    std::vector<std::pair<const TotemRPLocalTrack*, const CTPPSDetId>> GetStrips();

  private:
    std::vector<const CTPPSLocalTrackLite*> protonslitelist;
    std::vector<const reco::ForwardProton*> protonsrecolist;
    std::vector<std::pair<const CTPPSPixelLocalTrack*, const CTPPSDetId>> pixelsppslist;
    std::vector<std::pair<const CTPPSDiamondLocalTrack*, const CTPPSDetId>> timingppslist;
    std::vector<std::pair<const TotemRPLocalTrack*, const CTPPSDetId>> stripstotemlist;

};

#endif
