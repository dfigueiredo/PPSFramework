#ifndef GENMET_EVENT_H
#define GENMET_EVENT_H

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

// system include files
#include <memory>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// GenMET collection
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
//
// class declaration
//

class GenMETEvent{
  public:
    GenMETEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::View<reco::GenMET>>&);
    ~GenMETEvent();
    edm::View<reco::GenMET>::const_iterator GetGenMET();

  private:
    edm::View<reco::GenMET>::const_iterator genmetlist;

};

#endif
