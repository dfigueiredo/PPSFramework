#ifndef ELECTRONS_EVENT_H
#define ELECTRONS_EVENT_H

// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/ElectronsEvent
// Class:      ElectronsEvent
// 
/**\class ElectronsEvent ElectronsEvent.cc CTPPSAnalysisCode/ElectronsEvent/plugins/ElectronsEvent.cc

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

// Electrons collection
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//
// class declaration
//

class ElectronsEvent{
  public:
    ElectronsEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::View<pat::Electron>>&);
    ~ElectronsEvent();
    std::vector<const pat::Electron*> GetElectrons();

  private:
    std::vector<const pat::Electron*> electronslist;

};

#endif
