// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/MissingMassEfficiency
// Class:      MissingMassEfficiency
// 
/**\class MissingMassEfficiency MissingMassEfficiency.cc CTPPSAnalysisCode/MissingMassEfficiency/plugins/MissingMassEfficiency.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo and Nicola Turini
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//


// system include files
#include <memory>
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include <iostream>  
#include <typeinfo>
#include "TLorentzVector.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

// HLT
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

// CMS Calorimeters IDs
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

// Event
#include "TriggerEvent.h"
#include "JetsEvent.h"
#include "ElectronsEvent.h"
#include "MuonsEvent.h"
#include "PhotonsEvent.h"
#include "METEvent.h"
#include "ParticleFlowEvent.h"
#include "VerticesEvent.h"
#include "ProtonsEvent.h"
#include "GenProtonsEvent.h"
#include "GenParticlesEvent.h"
#include "GenMETEvent.h"
#include "GenJetsEvent.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#define ECM 13000.0 // GeV

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MissingMassEfficiency : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {

  public:
    explicit MissingMassEfficiency(const edm::ParameterSet&);
    ~MissingMassEfficiency();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    void endRun(const edm::Run&, const edm::EventSetup&);
    void beginRun(const edm::Run&, const edm::EventSetup&);

    // Creating CMS Physics Object
    void fetchGenParticles(const edm::Event&, const edm::EventSetup&);
    void fetchGenMET(const edm::Event&, const edm::EventSetup&);
    void fetchGenJets(const edm::Event&, const edm::EventSetup&);
    void fetchMuons(const edm::Event&, const edm::EventSetup&);
    void fetchElectrons(const edm::Event&, const edm::EventSetup&);
    void fetchJets(const edm::Event&, const edm::EventSetup&);
    void fetchPhotons(const edm::Event&, const edm::EventSetup&);
    void fetchMET(const edm::Event&, const edm::EventSetup&);
    void fetchVertices(const edm::Event&, const edm::EventSetup&);
    void fetchProtonsReco(const edm::Event&, const edm::EventSetup&);
    bool fetchTrigger(const edm::Event&, const edm::EventSetup&);
    void cleaning();
    void selection(const edm::Event&, const edm::EventSetup&);

    template <class T, class W>
      math::XYZTLorentzVector DiSystem(T obj1, W obj2);

    template <class T, class W>
      double TransverseMass(T lepton1, W lepton2);

    template <class T, class W>
      W Unmatching(T lepton, W jet, double radius);

    template <class T, class W>
      bool DiSystemFilter(T obj1, W obj2);

    template<typename T>
      constexpr T reduceRange(T x);

    // ----------member data ---------------------------

    bool debugging_;
    int verbosity_;

    // Interface to the HLT information
    HLTPrescaleProvider hltPrescaleProvider_;

    // Switches
    std::string physics_;
    bool unmatching_;
    bool orderDiscriminator_;
    bool enableMC_;
    bool enableTrigger_;
    bool enablePrescales_;

    // Trigger
    std::vector<std::string> triggersList_;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;

    //GenParticles
    edm::EDGetTokenT<std::vector<reco::GenParticle>> GenPartToken_;

    //GenMET
    edm::EDGetTokenT<edm::View<reco::GenMET>> GenMETToken_;

    //GenJets
    edm::EDGetTokenT<edm::View<reco::GenJet>> GenJetTokenA_;
    edm::EDGetTokenT<edm::View<reco::GenJet>> GenJetTokenB_;

    // Jets
    edm::EDGetTokenT<edm::View<pat::Jet> > patjetTokenA_;
    edm::EDGetTokenT<edm::View<pat::Jet> > patjetTokenB_;

    // Electrons
    edm::EDGetTokenT<edm::View<pat::Electron> > eleToken_;

    // Muons
    edm::EDGetTokenT<edm::View<pat::Muon> > muonToken_;

    // Photons
    edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;

    // MET
    edm::EDGetTokenT<edm::View<pat::MET> > patmetToken_;
    edm::EDGetTokenT<std::vector<reco::PFMET> > metToken_;

    // Vertex
    edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;

    // Protons Reco
    edm::EDGetTokenT<reco::ForwardProtonCollection > protonSingleToken_;
    edm::EDGetTokenT<reco::ForwardProtonCollection > protonMultiToken_;

    // Vectors with Collections
    std::vector<int> triggerVec;
    std::vector<int> prescalesL1Vec;
    std::vector<int> prescalesHLTVec;

    std::vector<const reco::GenParticle*> genpartVec;
    std::vector<const reco::GenParticle*> genprotonsVec;
    edm::View<reco::GenMET>::const_iterator genmet;
    std::vector<const reco::GenJet*> genjetsVecA;
    std::vector<const reco::GenJet*> genjetsVecB;

    std::vector<const pat::Jet*> jetsVecA;
    std::vector<const pat::Jet*> jetsVecB;

    std::vector<const pat::Electron*> electronsVec;
    std::vector<const pat::Muon*> muonsVec;
    std::vector<const pat::Photon*> photonsVec;
    std::vector<const reco::Vertex*> vtxVec;
    std::vector<const reco::ForwardProton*> protonSingleVec;
    std::vector<const reco::ForwardProton*> protonMultiVec;

    std::vector<const pat::Jet*> jetsak4Cand;
    std::vector<const pat::Jet*> jetsak8Cand;

    edm::View<pat::MET>::const_iterator met;

    // Cross-check histograms
    TH1F *selectionHisto_;

    TH1F *diMassHisto_;
    TH1F *diPhiHisto_;
    TH1F *diEtaHisto_;
    TH1F *diPtHisto_;
    TH1F *xMassHisto_;
    TH1F *xPhiHisto_;
    TH1F *xEtaHisto_;
    TH1F *xPtHisto_;
    TH1F *xiArm0Histo_;
    TH1F *xiArm1Histo_;
    TH1F *diffMassHisto_;

    TH1F *deltaPtHisto_;
    TH1F *deltaPhiHisto_;

    TH1F *xdeltaPtHisto_;
    TH1F *xdeltaPhiHisto_;


    TLorentzVector p1Multi, p2Multi, disystem, xsystem, cand1, cand2;

    // To be used with i.e. std::sort(v.begin(), v.end(), orderPT()), vector will be organized by pt.
    struct orderPt
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (vec1->pt() > vec2->pt());
	}
    };

    struct orderEta
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (vec1->eta() > vec2->eta());
	}
    };

    struct orderAbsolutPz
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (fabs(vec1->pz()) > fabs(vec2->pz()));
	}
    };

    struct orderVz
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (fabs(vec1->z()) > fabs(vec2->z()));
	}
    };

    struct orderDiscriminator
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return ((vec1->bDiscriminator("pfDeepCSVJetTags:probb") + vec1->bDiscriminator("pfDeepCSVJetTags:probbb")) > vec2->bDiscriminator("pfDeepCSVJetTags:probb") + vec2->bDiscriminator("pfDeepCSVJetTags:probbb"));
	}
    };


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

