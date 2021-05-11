#ifndef TRIGGER_EVENT_H
#define TRIGGER_EVENT_H

// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/TriggerEvent
// Class:      TriggerEvent
// 
/**\class TriggerEvent TriggerEvent.cc CTPPSAnalysisCode/TriggerEvent/plugins/TriggerEvent.cc

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

// L1 collections
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEmCand.h"

// HLT information
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"


//
// class declaration
//

class TriggerEvent{
  public:
    TriggerEvent(const edm::Event&, const edm::EventSetup&, HLTPrescaleProvider &hltPrescaleProvider_, edm::EDGetTokenT<edm::TriggerResults>&, std::vector<std::string>&, bool);
    ~TriggerEvent();
    std::vector<int> GetTrigger();
    std::vector<int> GetPrescalesL1();
    std::vector<int> GetPrescalesHLT();

  private:
    std::vector<int> triggerlist;
    std::vector<int> prescalesL1;
    std::vector<int> prescalesHLT;

};

#endif
