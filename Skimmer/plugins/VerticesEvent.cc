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

#include "VerticesEvent.h"

// Retrieve Vertices
std::vector<const reco::Vertex*> VerticesEvent::GetVertices(){
  return verticeslist;
}

// Class Definition, loop over all vertices events and fill vector.
VerticesEvent::VerticesEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::View<reco::Vertex>>& verticesToken_){

  verticeslist.clear();

  try{
    // Get the Jet collection from the event
    edm::Handle<edm::View<reco::Vertex> > vertexColl; // PAT
    iEvent.getByToken( verticesToken_, vertexColl );

    for (unsigned int i = 0; i < vertexColl->size(); ++i ) {
      const reco::Vertex* vtx = &((*vertexColl)[i]);
      verticeslist.push_back(vtx);
    }
    LogDebug( "Vertices Event" ) << "Passed Loop on vertices";
  }catch(...){
    LogDebug( "Vertices Event" ) << "collection not present in the sample";
  }

}
