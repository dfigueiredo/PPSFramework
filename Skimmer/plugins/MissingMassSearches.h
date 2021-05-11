// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/MissingMassSearches
// Class:      MissingMassSearches
// 
/**\class MissingMassSearches MissingMassSearches.cc CTPPSAnalysisCode/MissingMassSearches/plugins/MissingMassSearches.cc

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
#include "GenParticlesEvent.h"
#include "GenMETEvent.h"
#include "GenJetsEvent.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MissingMassSearches : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {

  public:
    explicit MissingMassSearches(const edm::ParameterSet&);
    ~MissingMassSearches();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    virtual void cleaning();
    virtual void eventClear();
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
    void fetchPF(const edm::Event&, const edm::EventSetup&);
    void fetchPPSInfo(const edm::Event&, const edm::EventSetup&);
    void fetchLHCInfo(const edm::EventSetup&);
    void fetchProtonsReco(const edm::Event&, const edm::EventSetup&);
    bool fetchTrigger(const edm::Event&, const edm::EventSetup&);
    void fetchDetector(const edm::Event&, const edm::EventSetup&);
    void fetchEventTagger(const edm::Event&);
    //void fetchPileUpInfo(const edm::EventSetup&, const edm::EventSetup&);
    void Debug();

    TTree * tree_;

    template <class T, class W>
      math::XYZTLorentzVector DiSystem(T obj1, W obj2);

    template <class T, class W>
      double TransverseMass(T lepton1, W lepton2);

    template <class T, class W>
      W Unmatching(T lepton, W jet, double radius);

    template <class T, class W>
      void fillingJets(T jet1, W jet2);


    // Only for PF Jets
    //template <class T, class W>
    //void fillingExtraInfoJets(T jet1, W jet2);

    // ----------member data ---------------------------

    // Interface to the HLT information
    HLTPrescaleProvider hltPrescaleProvider_;

    // Switches
    bool debug_;
    std::string tier_;
    std::string physics_;
    std::string year_;
    bool ppslite_;
    bool includeMuons_, includeElectrons_;
    bool includeJets_;
    bool includePhotons_;
    bool includeMET_;
    bool includeVertices_;
    bool includePF_;     
    bool includeProtonsReco_;
    bool includePPSInfo_;
    bool enableMC_;
    bool enableTrigger_;
    bool enablePrescales_;
    bool unmatching_;

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

    // PF
    edm::EDGetTokenT<std::vector< reco::PFCandidate > > pfToken_;
    edm::EDGetTokenT<std::vector< pat::PackedCandidate > > packedToken_;

    // Photons
    edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;

    // MET
    edm::EDGetTokenT<edm::View<pat::MET> > patmetToken_;
    edm::EDGetTokenT<std::vector<reco::PFMET> > metToken_;

    // Vertex
    edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;

    // Protons Lite
    edm::EDGetTokenT<std::vector<CTPPSLocalTrackLite> > protonToken_;

    // Protons Reco
    edm::EDGetTokenT<reco::ForwardProtonCollection > protonSingleToken_;
    edm::EDGetTokenT<reco::ForwardProtonCollection > protonMultiToken_;

    // PPS/TOTEM Detectors
    edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack> > pixelsppsToken_;
    edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondLocalTrack> > timingppsToken_;
    edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack> > stripstotemToken_;

    // Calo Towers
    edm::EDGetTokenT<CaloTowerCollection> calotowersToken_;
  
    //PileUpSummaryInfo 
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;

    // Energy Thresholds
    double energyThresholdHB_;
    double energyThresholdHE_;
    double energyThresholdHF_;
    double energyThresholdEB_;
    double energyThresholdEE_;

    // LHC Info Label
    std::string lhcInfoLabel_;
    
    // Vectors with Collections
    std::vector<int> triggerVec;
    std::vector<int> prescalesL1Vec;
    std::vector<int> prescalesHLTVec;

    std::vector<const reco::GenParticle*> genpartVec;
    edm::View<reco::GenMET>::const_iterator genmet;
    std::vector<const reco::GenJet*> genjetsVecA;
    std::vector<const reco::GenJet*> genjetsVecB;

    std::vector<const pat::Jet*> jetsVecA;
    std::vector<const pat::Jet*> jetsVecB;

    std::vector<const pat::Electron*> electronsVec;
    std::vector<const pat::Muon*> muonsVec;
    std::vector<const reco::PFCandidate*> pfVec;
    std::vector<const pat::PackedCandidate*> packedVec;
    std::vector<const pat::Photon*> photonsVec;
    std::vector<const reco::Vertex*> vtxVec;
    std::vector<const CTPPSLocalTrackLite*> ppsVec;
    std::vector<const reco::ForwardProton*> protonSingleVec;
    std::vector<const reco::ForwardProton*> protonMultiVec;
    std::vector<std::pair<const CTPPSPixelLocalTrack*, const CTPPSDetId>> pixelsVec;
    std::vector<std::pair<const CTPPSDiamondLocalTrack*, const CTPPSDetId>> timingVec;
    std::vector<std::pair<const TotemRPLocalTrack*, const CTPPSDetId>> stripsVec;

    std::vector<const pat::Jet*> jetsak4Cand;
    std::vector<const pat::Jet*> jetsak8Cand;

    edm::View<pat::MET>::const_iterator met;

    std::vector<std::pair<const pat::Muon*, const pat::Jet*>> MuonsJetsCand;
    std::vector<std::pair<const pat::Electron*, const pat::Jet*>> ElectronsJetsCand;

    // TTree Vectors

    long int *ev_;
    int *run_;
    int *lumiblock_;

    std::vector<int> *trigger_;
    std::vector<int> *prescalesL1_;
    std::vector<int> *prescalesHLT_;

    int *xangle_;
    double *betastar_;
    int *bunchesb1_;
    int *bunchesb2_;
    double *intensityb1_;
    double *intensityb2_;
    double *instlumi_;

    int *PUinter_;
    int *PUtrueinter_;

    std::vector<double> *vertex_x_;
    std::vector<double> *vertex_y_;
    std::vector<double> *vertex_z_;
    std::vector<int> *vertex_ntrack_;
    std::vector<double> *vertex_chi2_;
    std::vector<double> *vertex_ndof_;

    std::vector<int> *genleptons_pdgid_;
    std::vector<double> *genleptons_energy_;
    std::vector<double> *genleptons_pt_;
    std::vector<double> *genleptons_eta_;
    std::vector<double> *genleptons_phi_;
    std::vector<double> *genleptons_px_;
    std::vector<double> *genleptons_py_;
    std::vector<double> *genleptons_pz_;
    std::vector<double> *genleptons_charge_;
    std::vector<double> *genleptons_vx_;
    std::vector<double> *genleptons_vy_;
    std::vector<double> *genleptons_vz_;
 
    std::vector<double> *genprotons_energy_;
    std::vector<double> *genprotons_pt_;
    std::vector<double> *genprotons_eta_;
    std::vector<double> *genprotons_phi_;
    std::vector<double> *genprotons_px_;
    std::vector<double> *genprotons_py_;
    std::vector<double> *genprotons_pz_;
    std::vector<double> *genprotons_xi_;
  
    std::vector<double> *genphotons_pt_;
    std::vector<double> *genphotons_eta_;
    std::vector<double> *genphotons_phi_;
    std::vector<double> *genphotons_energy_;

    double *genmisset_;
    double *genmisset_phi_;

    std::vector<double> *genjetsak4_px_;
    std::vector<double> *genjetsak4_py_;
    std::vector<double> *genjetsak4_pz_;
    std::vector<double> *genjetsak4_pt_;
    std::vector<double> *genjetsak4_energy_;
    std::vector<double> *genjetsak4_phi_;
    std::vector<double> *genjetsak4_eta_;
    std::vector<double> *genjetsak4_vz_;

    std::vector<double> *genjetsak8_px_;
    std::vector<double> *genjetsak8_py_;
    std::vector<double> *genjetsak8_pz_;
    std::vector<double> *genjetsak8_pt_;
    std::vector<double> *genjetsak8_energy_;
    std::vector<double> *genjetsak8_phi_;
    std::vector<double> *genjetsak8_eta_;
    std::vector<double> *genjetsak8_vz_;

    std::vector<double> *jetsak4_px_;
    std::vector<double> *jetsak4_py_;
    std::vector<double> *jetsak4_pz_;
    std::vector<double> *jetsak4_pt_;
    std::vector<double> *jetsak4_energy_;
    std::vector<double> *jetsak4_phi_;
    std::vector<double> *jetsak4_eta_;
    std::vector<double> *jetsak4_vx_;
    std::vector<double> *jetsak4_vy_;
    std::vector<double> *jetsak4_vz_;
    std::vector<double> *jetsak4_bdis_;
    std::vector<double> *jetsak4_qgdis_;
    std::vector<double> *jetsak4_neutralEmFraction_;
    std::vector<double> *jetsak4_neutralHadFraction_;
    std::vector<double> *jetsak4_chargedEmFraction_;
    std::vector<double> *jetsak4_chargedHadFraction_;
    std::vector<double> *jetsak4_muonFraction_;
    std::vector<int> *jetsak4_neutralMultiplicity_;
    std::vector<int> *jetsak4_chargedMultiplicity_;
    std::vector<bool> *jetsak4_looseJetId_;
    std::vector<bool> *jetsak4_tightJetId_;
    std::vector<bool> *jetsak4_lepVetoJetId_;
    std::vector<double> *jetsak4_puIdfdisc_;
    std::vector<int> *jetsak4_puIdcbased_;
    std::vector<int> *jetsak4_puIdfid_;

    std::vector<double> *jetsak8_px_;
    std::vector<double> *jetsak8_py_;
    std::vector<double> *jetsak8_pz_;
    std::vector<double> *jetsak8_pt_;
    std::vector<double> *jetsak8_energy_;
    std::vector<double> *jetsak8_phi_;
    std::vector<double> *jetsak8_eta_;
    std::vector<double> *jetsak8_vx_;
    std::vector<double> *jetsak8_vy_;
    std::vector<double> *jetsak8_vz_;
    std::vector<double> *jetsak8_bdis_;
    std::vector<double> *jetsak8_neutralEmFraction_;
    std::vector<double> *jetsak8_neutralHadFraction_;
    std::vector<double> *jetsak8_chargedEmFraction_;
    std::vector<double> *jetsak8_chargedHadFraction_;
    std::vector<double> *jetsak8_muonFraction_;
    std::vector<int> *jetsak8_neutralMultiplicity_;
    std::vector<int> *jetsak8_chargedMultiplicity_;
    std::vector<bool> *jetsak8_looseJetId_;
    std::vector<bool> *jetsak8_tightJetId_;
    std::vector<bool> *jetsak8_lepVetoJetId_;

    std::vector<int> *leptons_pdgid_;
    std::vector<double> *leptons_energy_;
    std::vector<double> *leptons_pt_;
    std::vector<double> *leptons_eta_;
    std::vector<double> *leptons_phi_;
    std::vector<double> *leptons_px_;
    std::vector<double> *leptons_py_;
    std::vector<double> *leptons_pz_;
    std::vector<double> *leptons_charge_;
    std::vector<double> *leptons_vx_;
    std::vector<double> *leptons_vy_;
    std::vector<double> *leptons_vz_;
    std::vector<bool> *leptons_looseId_;
    std::vector<bool> *leptons_mediumId_;
    std::vector<bool> *leptons_tightId_;
    std::vector<bool> *leptons_pfIsoMedium_;
    std::vector<bool> *leptons_miniIsoTight_;
    std::vector<bool> *leptons_pfIsoVeryTight_;
    std::vector<double> *leptons_pfIso_;
    std::vector<double> *leptons_tkIso_;

    std::vector<double> *photons_pt_;
    std::vector<double> *photons_eta_;
    std::vector<double> *photons_phi_;
    std::vector<double> *photons_energy_;
    std::vector<double> *photons_r9_;
    std::vector<bool> *photons_looseId_;
    std::vector<bool> *photons_mediumId_;
    std::vector<bool> *photons_tightId_;

    double *misset_;
    double *misset_phi_;

    int *PUweight_;
    int *nGenParticles_;
    int *nGenMuons_;    
    int *nGenElectrons_;
    int *nGenProtons_;
    int *nPF_;
    int *nMuons_;
    int *nMuons_looseId_;
    int *nMuons_mediumId_;
    int *nMuons_tightId_;
    int *nMuons_pfIsoMedium_;
    int *nMuons_miniIsoTight_;
    int *nMuons_pfIsoVeryTight_;
    int *nElectrons_;
    int *nElectrons_looseId_;
    int *nElectrons_mediumId_;
    int *nElectrons_tightId_;
    int *nPhotons_;
    double *SumPF_energy_;

    double *pfHFHadPlus_energy_;
    double *pfHFHadPlus_pt_;
    double *pfHFHadMinus_energy_;
    double *pfHFHadMinus_pt_;
    double *pfHFEmPlus_energy_;
    double *pfHFEmPlus_pt_;
    double *pfHFEmMinus_energy_;
    double *pfHFEmMinus_pt_;
    double *pfEcalMinus_energy_;
    double *pfEcalPlus_energy_;
    double *pfHcalMinus_energy_;
    double *pfHcalPlus_energy_;
    double *pfHadNeutralMinus_energy_;
    double *pfHadNeutralPlus_energy_;
    double *SumChargedPFMultiPV_pt_Loose_;
    double *SumChargedPFMultiPV_pt_Tight_;
    double *SumChargedPFMultiPV_pt_UsedInFit_;
    double *SumChargedPFMultiPV_pt_Tight_Fit_;

    int *nChargedPFMultiPV_Loose_;
    int *nChargedPFMultiPV_Tight_;
    int *nChargedPFMultiPV_UsedInFit_;
    int *nChargedPFMultiPV_Tight_Fit_;
    int *nPFPlus_;
    int *nPFMinus_;
    int *nPFHFHadPlus_;
    int *nPFHFHadMinus_;
    int *nPFHFEmPlus_;
    int *nPFHFEmMinus_;

    double *pfSlice0P_energy_;
    double *pfSlice1P_energy_;
    double *pfSlice2P_energy_;
    double *pfSlice3P_energy_;
    double *pfSlice4P_energy_;
    double *pfSlice0N_energy_;
    double *pfSlice1N_energy_;
    double *pfSlice2N_energy_;
    double *pfSlice3N_energy_;
    double *pfSlice4N_energy_;

    double *pfSlice0P_pt_;
    double *pfSlice1P_pt_;
    double *pfSlice2P_pt_;
    double *pfSlice3P_pt_;
    double *pfSlice4P_pt_;
    double *pfSlice0N_pt_;
    double *pfSlice1N_pt_;
    double *pfSlice2N_pt_;
    double *pfSlice3N_pt_;
    double *pfSlice4N_pt_;

    int *npfSlice0P_;
    int *npfSlice1P_;
    int *npfSlice2P_;
    int *npfSlice3P_;
    int *npfSlice4P_;
    int *npfSlice0N_;
    int *npfSlice1N_;
    int *npfSlice2N_;
    int *npfSlice3N_;
    int *npfSlice4N_;

    int *nmuonSlice0P_;
    int *nmuonSlice1P_;
    int *nmuonSlice2P_;
    int *nmuonSlice3P_;
    int *nmuonSlice4P_;
    int *nmuonSlice0N_;
    int *nmuonSlice1N_;
    int *nmuonSlice2N_;
    int *nmuonSlice3N_;
    int *nmuonSlice4N_;

    int *nelectronSlice0P_;
    int *nelectronSlice1P_;
    int *nelectronSlice2P_;
    int *nelectronSlice3P_;
    int *nelectronSlice4P_;
    int *nelectronSlice0N_;
    int *nelectronSlice1N_;
    int *nelectronSlice2N_;
    int *nelectronSlice3N_;
    int *nelectronSlice4N_;

    std::vector<int> *singleProtonArm_;
    std::vector<int> *singleProtonStation_;
    std::vector<int> *singleProtonPot_;
    std::vector<double> *singleProtonXi_;
    std::vector<double> *singleProtonT_;
    std::vector<double> *singleProtonThetaX_;
    std::vector<double> *singleProtonThetaY_;

    std::vector<int> *multiProtonArm_;
    std::vector<double> *multiProtonXi_;
    std::vector<double> *multiProtonT_;
    std::vector<double> *multiProtonTime_;
    std::vector<double> *multiProtonTimeError_;
    std::vector<double> *multiProtonVz_;
    std::vector<double> *multiProtonThetaX_;
    std::vector<double> *multiProtonThetaY_;

    std::vector<bool> *protonsIsValid_;
    std::vector<int> *protonsArm_;
    std::vector<int> *protonsStation_;
    std::vector<int> *protonsRP_;
    std::vector<double> *protonsX_;
    std::vector<double> *protonsXUnc_;
    std::vector<double> *protonsY_;
    std::vector<double> *protonsYUnc_;
    std::vector<double> *protonsTx_;
    std::vector<double> *protonsTxUnc_;
    std::vector<double> *protonsTy_;
    std::vector<double> *protonsTyUnc_;
    std::vector<double> *protonsTime_;
    std::vector<double> *protonsTimeUnc_;

    // Cross-check histograms
    TH1F *hltTriggerNamesHisto_;

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

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

