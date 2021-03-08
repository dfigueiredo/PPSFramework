// -*- C++ -*-
//
// Package:    Filter/ProtonTagFilter
// Class:      ProtonTagFilter
// 
/**\class ProtonTagFilter ProtonTagFilter.cc Filter/ProtonTagFilter/plugins/ProtonTagFilter.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Wed, 10 Oct 2018 00:10:06 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "ProtonsEvent.h"

//
// class declaration
//

using namespace edm;

class ProtonTagFilter : public edm::stream::EDFilter<> {
  public:
    explicit ProtonTagFilter(const edm::ParameterSet&);
    ~ProtonTagFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------

    // Switches
    bool debug_;
    bool singlePotMode_;
      bool includeProtonsReco_;

    // PPS/TOTEM Detectors
    edm::EDGetTokenT<reco::ForwardProtonCollection > protonSingleToken_;
    edm::EDGetTokenT<reco::ForwardProtonCollection > protonMultiToken_;
    edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack> > pixelsppsToken_;
    edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondLocalTrack> > timingppsToken_;
    edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack> > stripstotemToken_;

    // Filling vectors with CMS physics objects
    std::vector<const reco::ForwardProton*> protonSingleVec;
    std::vector<std::pair<const CTPPSPixelLocalTrack*, const CTPPSDetId>> pixelsVec;
    std::vector<std::pair<const CTPPSDiamondLocalTrack*, const CTPPSDetId>> timingVec;
    std::vector<std::pair<const TotemRPLocalTrack*, const CTPPSDetId>> stripsVec;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

