#ifndef MET_EVENT_H
#define MET_EVENT_H

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

// MET collection
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

//
// class declaration
//

class METEvent{
  public:
    METEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::View<pat::MET>>&);
    ~METEvent();
    edm::View<pat::MET>::const_iterator GetMET();

  private:
    edm::View<pat::MET>::const_iterator metlist;

};

#endif
