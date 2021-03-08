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

#include "ProtonsEvent.h"

// Retrieve Protons
std::vector<const CTPPSLocalTrackLite*> ProtonsEvent::GetProtonsLite(){
  return protonslitelist;
}

// Retrieve Protons Reco
std::vector<const reco::ForwardProton*> ProtonsEvent::GetProtonsReco(){
  return protonsrecolist;
}

// Retrieve Pixels
std::vector<std::pair<const CTPPSPixelLocalTrack*, const CTPPSDetId>> ProtonsEvent::GetPixels(){
  return pixelsppslist;
}

// Retrieve Timing
std::vector<std::pair<const CTPPSDiamondLocalTrack*, const CTPPSDetId>> ProtonsEvent::GetTiming(){
  return timingppslist;
}

// Retrieve Strips
std::vector<std::pair<const TotemRPLocalTrack*, const CTPPSDetId>> ProtonsEvent::GetStrips(){
  return stripstotemlist;
}

// Class Definition, loop over all pps pixels events and fill vector.
ProtonsEvent::ProtonsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack>>& pixelsppsToken_){

  pixelsppslist.clear();

  try{
    edm::Handle< edm::DetSetVector<CTPPSPixelLocalTrack> > rppixellocaltracks;
    iEvent.getByToken(pixelsppsToken_, rppixellocaltracks);

    for ( const auto& dsv : *rppixellocaltracks ) {
      const CTPPSDetId det_id( dsv.detId() );
      for ( const auto& trk : dsv ) {
	if ( !trk.isValid() ) continue;
	pixelsppslist.push_back( std::pair<const CTPPSPixelLocalTrack*, const CTPPSDetId>(&trk, det_id) );
      }
    }
    LogDebug( "Pixel Event" ) << "Passed Loop on Proton Candidate";
  }catch(...){
    LogDebug( "Pixel Event" ) << "collection not presented in the sample.";
  }

}

// Class Definition, loop over all pps timing events and fill vector.
ProtonsEvent::ProtonsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondLocalTrack>>& timingppsToken_){

  timingppslist.clear();

  try{
    edm::Handle< edm::DetSetVector<CTPPSDiamondLocalTrack> > rptiminglocaltracks;
    iEvent.getByToken(timingppsToken_, rptiminglocaltracks);

    for ( const auto& dsv : *rptiminglocaltracks ) {
      const CTPPSDetId det_id( dsv.detId() );
      for ( const auto& timing : dsv ) {
	if ( !timing.isValid() ) continue;
	timingppslist.push_back( std::pair<const CTPPSDiamondLocalTrack*, const CTPPSDetId>(&timing, det_id) );
      }
    }
    LogDebug( "Timing Event" ) << "Passed Loop on Proton Candidate";
  }catch(...){
    LogDebug( "Timing Event" ) << "collection not presented in the sample.";
  }

}

// Class Definition, loop over all pps timing events and fill vector.
ProtonsEvent::ProtonsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack>>& stripstotemToken_){

  stripstotemlist.clear();

  try{
    edm::Handle< edm::DetSetVector<TotemRPLocalTrack> > rpstripslocaltracks;
    iEvent.getByToken(stripstotemToken_, rpstripslocaltracks);

    for ( const auto& dsv : *rpstripslocaltracks ) {
      const CTPPSDetId det_id( dsv.detId() );
      for ( const auto& trk : dsv ) {
	if ( !trk.isValid() ) continue;
	stripstotemlist.push_back( std::pair<const TotemRPLocalTrack*, const CTPPSDetId>(&trk, det_id) );
      }
    }
    LogDebug( "Strips Event" ) << "Passed Loop on Proton Candidate";
  }catch(...){
    LogDebug( "Strips Event" ) << "collection not presented in the sample.";
  }

}

// Class Definition, loop over all protons events and fill vector.
ProtonsEvent::ProtonsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite>>& protonsliteToken_){

  protonslitelist.clear();

  try{
    // Get the Proton collection from the event
    edm::Handle<std::vector<CTPPSLocalTrackLite>> protonColl;
    iEvent.getByToken(protonsliteToken_, protonColl);

    for (unsigned int i = 0; i < protonColl->size(); ++i ) {
      const CTPPSLocalTrackLite* proton = &((*protonColl)[i]);
      protonslitelist.push_back(proton);
    }
    LogDebug( "Proton Event" ) << "Passed Loop on Proton Candidate";
  }catch(...){
    LogDebug( "Proton Event" ) << "collection not presented in the sample.";
  }

}

// Class Definition, loop over all reconstructed protons events and fill vector.
ProtonsEvent::ProtonsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<reco::ForwardProtonCollection>& protonsrecoToken_){

  protonslitelist.clear();

  try{
    // Get the Proton reconstructed from the event
    edm::Handle<reco::ForwardProtonCollection> protonReco;
    iEvent.getByToken(protonsrecoToken_, protonReco);

    for (unsigned int i = 0; i < protonReco->size(); ++i ) {
      const reco::ForwardProton* proton = &((*protonReco)[i]);
      protonsrecolist.push_back(proton);
    }
    LogDebug( "Proton Reco Event" ) << "Passed Loop on Proton Candidate";
  }catch(...){
    LogDebug( "Proton Reco Event" ) << "collection not presented in the sample.";
  }

}
