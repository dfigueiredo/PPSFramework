#ifndef VERTICES_EVENT_H
#define VERTICES_EVENT_H

// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/VerticesEvent
// Class:      VerticesEvent
// 
/**\class VerticesEvent VerticesEvent.cc CTPPSAnalysisCode/VerticesEvent/plugins/VerticesEvent.cc

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

// Vertices collection
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

//
// class declaration
//

class VerticesEvent{
  public:
    VerticesEvent(const edm::Event&, const edm::EventSetup&, edm::EDGetTokenT<edm::View<reco::Vertex>>&);
    ~VerticesEvent();
    std::vector<const reco::Vertex*> GetVertices();

  private:
    std::vector<const reco::Vertex*> verticeslist;

};

#endif
