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
// Original Authors:  Diego Figueiredo, Lorenzo Pagliai and Nicola Turini
//          Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

// user include files
#include "MissingMassSearches.h"

//
// constructors and destructor
//
MissingMassSearches::MissingMassSearches(const edm::ParameterSet& iConfig):
  hltPrescaleProvider_ ( iConfig, consumesCollector(), *this),
  debug_               ( iConfig.getParameter<bool>( "debugging" ) ),
  tier_                ( iConfig.getParameter<std::string>( "tier" ) ),
  physics_                ( iConfig.getParameter<std::string>( "physics" ) ),
  year_                ( iConfig.getParameter<std::string>( "year" ) ),
  ppslite_             ( iConfig.getParameter<bool>( "ppslite" ) ),
  includeMuons_        ( iConfig.getParameter<bool>( "includeMuons" ) ),
  includeElectrons_    ( iConfig.getParameter<bool>( "includeElectrons" ) ),
  includeJets_         ( iConfig.getParameter<bool>( "includeJets" ) ),
  includePhotons_      ( iConfig.getParameter<bool>( "includePhotons" ) ),
  includeMET_          ( iConfig.getParameter<bool>( "includeMET" ) ),
  includeVertices_     ( iConfig.getParameter<bool>( "includeVertices" ) ),
  includePF_           ( iConfig.getParameter<bool>( "includeParticleFlow" ) ),
  includeProtonsReco_      ( iConfig.getParameter<bool>( "includeProtonsReco" ) ),
  includePPSInfo_      ( iConfig.getParameter<bool>( "includePPSInfo" ) ),
  enableMC_            ( iConfig.getParameter<bool>( "enableMC" ) ),
  enableTrigger_            ( iConfig.getParameter<bool>( "enableTrigger" ) ),
  enablePrescales_            ( iConfig.getParameter<bool>( "enablePrescales" ) ),
  unmatching_            ( iConfig.getParameter<bool>( "unmatching" ) ),
  triggersList_        ( iConfig.getParameter<std::vector<std::string>>              ( "triggersList" ) ),
  triggerResultsToken_ ( consumes<edm::TriggerResults>                               ( iConfig.getParameter<edm::InputTag>( "triggerResults" ) ) ),
  GenPartToken_        ( consumes<std::vector<reco::GenParticle>>                    ( iConfig.getParameter<edm::InputTag>( "GenPartTag" ) ) ),
  GenMETToken_         ( consumes<edm::View<reco::GenMET>>                           ( iConfig.getParameter<edm::InputTag>( "GenMETTag" ) ) ),
  GenJetTokenA_        ( consumes<edm::View<reco::GenJet>>                           ( iConfig.getParameter<edm::InputTag>( "GenJetAlgoA" ) ) ),
  GenJetTokenB_        ( consumes<edm::View<reco::GenJet>>                           ( iConfig.getParameter<edm::InputTag>( "GenJetAlgoB" ) ) ),
  patjetTokenA_        ( consumes<edm::View<pat::Jet>>                               ( iConfig.getParameter<edm::InputTag>( "patJetAlgoA" ) ) ),
  patjetTokenB_        ( consumes<edm::View<pat::Jet>>                               ( iConfig.getParameter<edm::InputTag>( "patJetAlgoB" ) ) ),
  eleToken_            ( consumes<edm::View<pat::Electron>>                          ( iConfig.getParameter<edm::InputTag>( "electronTag" ) ) ),
  muonToken_           ( consumes<edm::View<pat::Muon>>                              ( iConfig.getParameter<edm::InputTag>( "muonTag" ) ) ),
  pfToken_             ( consumes<std::vector< reco::PFCandidate>>                   ( iConfig.getParameter<edm::InputTag>( "pfTag" ) ) ),
  packedToken_         ( consumes<std::vector< pat::PackedCandidate>>                ( iConfig.getParameter<edm::InputTag>( "packedTag" ) ) ),
  photonToken_         ( consumes<edm::View<pat::Photon>>                            ( iConfig.getParameter<edm::InputTag>( "photonTag" ) ) ),
  patmetToken_         ( consumes<edm::View<pat::MET>>                               ( iConfig.getParameter<edm::InputTag>( "patmetTag" ) ) ),
  vertexToken_         ( consumes<edm::View<reco::Vertex>>                           ( iConfig.getParameter<edm::InputTag>( "vertexTag" ) ) ),
  protonToken_         ( consumes<std::vector<CTPPSLocalTrackLite>>                  ( iConfig.getParameter<edm::InputTag>( "protonTag" ) ) ),
  protonSingleToken_   ( consumes<reco::ForwardProtonCollection>                     ( iConfig.getParameter<edm::InputTag>( "protonSingleTag" ) ) ),
  protonMultiToken_    ( consumes<reco::ForwardProtonCollection>                     ( iConfig.getParameter<edm::InputTag>( "protonMultiTag" ) ) ),
  pixelsppsToken_      ( consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>           ( iConfig.getParameter<edm::InputTag>( "pixelsppsTag" ) ) ),
  timingppsToken_      ( consumes<edm::DetSetVector<CTPPSDiamondLocalTrack>>         ( iConfig.getParameter<edm::InputTag>( "timingppsTag" ) ) ),
  stripstotemToken_    ( consumes<edm::DetSetVector<TotemRPLocalTrack>>              ( iConfig.getParameter<edm::InputTag>( "stripstotemTag" ) ) ),
  calotowersToken_     ( consumes<CaloTowerCollection>                               ( iConfig.getParameter<edm::InputTag>("caloTowerTag"))),
  puToken_             ( consumes<std::vector<PileupSummaryInfo>>                    (edm::InputTag("slimmedAddPileupInfo"))),
  energyThresholdHB_   ( iConfig.getParameter<double>("energyThresholdHB")),
  energyThresholdHE_   ( iConfig.getParameter<double>("energyThresholdHE")),
  energyThresholdHF_   ( iConfig.getParameter<double>("energyThresholdHF")),
  energyThresholdEB_   ( iConfig.getParameter<double>("energyThresholdEB")),
  energyThresholdEE_   ( iConfig.getParameter<double>("energyThresholdEE")),
  lhcInfoLabel_        ( iConfig.getParameter<std::string>("lhcInfoLabel"))         
{

  //now do what ever initialization is needed
  usesResource("TFileService");

  if(physics_ == "Muon" || physics_ == "MUON" || physics_ == "muon"){
    includeElectrons_ = false;
    includeMuons_ = true;
  }
  else if(physics_ == "Electron" || physics_ == "ELECTRON" || physics_ == "electron"){
    includeElectrons_ = true;
    includeMuons_ = false;
  }
  else if(physics_ == "DisplacedJet" || physics_ == "DISPLACEDJET" || physics_ == "displacedjet" || physics_ =="bjet" || physics_=="BJET" || physics_=="BJet"){
    includeElectrons_ = true;
    includeMuons_ = true;
  }
  else if(physics_ == "Electron" || physics_ == "ELECTRON" || physics_ == "electron"){
    includeElectrons_ = true;
    includeMuons_ = false;
  }
  else if(physics_.find("ZeroBias") != std::string::npos || physics_.find("ZEROBIAS") != std::string::npos || physics_.find("zerobias") != std::string::npos){
    includeElectrons_ = false;
    includeMuons_ = true;
  }
  else if(physics_ == "EMu" || physics_ == "emu" || physics_ == "EMU"){
    includeElectrons_ = true;
    includeMuons_ = true;
  }

  else throw cms::Exception( "MissingMassSearches" ) << "'physics' parameter should either be:\n"
    << "   * 'electron' or 'muon' or 'emu'";

  std::cout<<"\n<CMSSW Plugin> Missing Mass Skimmer..." << std::endl;
  std::cout<<"\t-->Enable MC: " << enableMC_ <<  std::endl;
  std::cout<<"\t-->Enable Trigger: " << enableTrigger_ <<  std::endl;
  std::cout<<"\t-->Enable Prescale: " << enablePrescales_ <<  std::endl;
  std::cout<<"\t-->Debug: " << debug_ <<  std::endl;              
  std::cout<<"\t-->Forward Proton Reco: " << includeProtonsReco_ <<  std::endl;              
  std::cout<<"\t-->PPS Info: " << includePPSInfo_ <<  std::endl;              
  std::cout<<"\t-->PPS Lite: " << ppslite_ <<  std::endl;
  std::cout<<"\t-->Physics Option: " << physics_ <<  std::endl;
  std::cout<<"\t-->Year: " << year_ <<  std::endl;
  std::cout<<"\t-->Electrons: " << includeElectrons_ <<  std::endl;
  std::cout<<"\t-->Muons: " << includeMuons_ << "\n" << std::endl;

}

MissingMassSearches::~MissingMassSearches()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//


// ------------ Cleaning  ------------
void MissingMassSearches::cleaning(){

  triggerVec.clear();
  genpartVec.clear();
  genprotonsVec.clear();
  genjetsVecA.clear();
  genjetsVecB.clear();
  jetsVecA.clear(); // AK4
  jetsVecB.clear(); // AK8
  electronsVec.clear();
  muonsVec.clear();
  packedVec.clear();
  pfVec.clear();
  photonsVec.clear();
  vtxVec.clear();
  ppsVec.clear();
  protonSingleVec.clear();
  protonMultiVec.clear();
  pixelsVec.clear();
  timingVec.clear();
  stripsVec.clear();
  jetsak4Cand.clear();
  jetsak8Cand.clear();
  MuonsJetsCand.clear();
  ElectronsJetsCand.clear();

}

// ------------ Cleaning  ------------
void MissingMassSearches::eventClear(){

  *run_=-999;
  *ev_=-999;
  *lumiblock_=-999;
  *misset_=-999.;
  *misset_phi_=-999.;
  *xangle_ = -999;
  *betastar_ = -999.;
  *bunchesb1_ = -999;
  *bunchesb2_ = -999;
  *intensityb1_ = -999.;
  *intensityb2_ = -999.;
  *instlumi_ = -999.;

  *PUinter_=-999;
  *PUtrueinter_=-999;

  (*trigger_).clear();
  (*prescalesL1_).clear();
  (*prescalesHLT_).clear();

  if(includeVertices_){
    (*vertex_x_).clear();
    (*vertex_y_).clear();
    (*vertex_z_).clear();
    (*vertex_ntrack_).clear();
    (*vertex_chi2_).clear();
    (*vertex_ndof_).clear();
  }

  if(enableMC_){

    *genmisset_=-999.;
    *genmisset_phi_=-999.;

    if(includeElectrons_ || includeMuons_){
      (*genleptons_pdgid_).clear();
      (*genleptons_status_).clear();
      (*genleptons_energy_).clear();
      (*genleptons_pt_).clear();
      (*genleptons_eta_).clear();
      (*genleptons_phi_).clear();
      (*genleptons_px_).clear();
      (*genleptons_py_).clear();
      (*genleptons_pz_).clear();
      (*genleptons_charge_).clear();
      (*genleptons_vx_).clear();
      (*genleptons_vy_).clear();
      (*genleptons_vz_).clear();
    }

    if(includeProtonsReco_){
      (*genprotons_status_).clear();
      (*genprotons_energy_).clear();
      (*genprotons_pt_).clear();
      (*genprotons_eta_).clear();
      (*genprotons_phi_).clear();
      (*genprotons_px_).clear();
      (*genprotons_py_).clear();
      (*genprotons_pz_).clear();
      (*genprotons_xi_).clear();
    }

    if(includePhotons_){
      (*genphotons_pt_).clear();
      (*genphotons_eta_).clear();
      (*genphotons_phi_).clear();
      (*genphotons_energy_).clear();
    }

    if(includeJets_){   
      (*genjetsak4_px_).clear();
      (*genjetsak4_py_).clear();
      (*genjetsak4_pz_).clear();
      (*genjetsak4_pt_).clear();
      (*genjetsak4_energy_).clear();
      (*genjetsak4_phi_).clear();
      (*genjetsak4_eta_).clear();
      (*genjetsak4_vz_).clear();

      (*genjetsak8_px_).clear();
      (*genjetsak8_py_).clear();
      (*genjetsak8_pz_).clear();
      (*genjetsak8_pt_).clear();
      (*genjetsak8_energy_).clear();
      (*genjetsak8_phi_).clear();
      (*genjetsak8_eta_).clear();
      (*genjetsak8_vz_).clear();
    }
  }

  if(includeJets_){
    (*jetsak4_px_).clear();
    (*jetsak4_py_).clear();
    (*jetsak4_pz_).clear();
    (*jetsak4_pt_).clear();
    (*jetsak4_energy_).clear();
    (*jetsak4_phi_).clear();
    (*jetsak4_eta_).clear();
    (*jetsak4_vx_).clear();
    (*jetsak4_vy_).clear();
    (*jetsak4_vz_).clear();
    (*jetsak4_bdis_).clear();
    (*jetsak4_qgdis_).clear();
    (*jetsak4_neutralEmFraction_).clear();
    (*jetsak4_neutralHadFraction_).clear();
    (*jetsak4_chargedEmFraction_).clear();
    (*jetsak4_chargedHadFraction_).clear();
    (*jetsak4_muonFraction_).clear();
    (*jetsak4_neutralMultiplicity_).clear();
    (*jetsak4_chargedMultiplicity_).clear();
    (*jetsak4_looseJetId_).clear();
    (*jetsak4_tightJetId_).clear();
    (*jetsak4_lepVetoJetId_).clear();
    (*jetsak4_puIdfdisc_).clear();
    (*jetsak4_puIdcbased_).clear();
    (*jetsak4_puIdfid_).clear();
    (*jetsak8_px_).clear();
    (*jetsak8_py_).clear();
    (*jetsak8_pz_).clear();
    (*jetsak8_pt_).clear();
    (*jetsak8_energy_).clear();
    (*jetsak8_phi_).clear();
    (*jetsak8_eta_).clear();
    (*jetsak8_vx_).clear();
    (*jetsak8_vy_).clear();
    (*jetsak8_vz_).clear();
    (*jetsak8_bdis_).clear();
    (*jetsak8_neutralEmFraction_).clear();
    (*jetsak8_neutralHadFraction_).clear();
    (*jetsak8_chargedEmFraction_).clear();
    (*jetsak8_chargedHadFraction_).clear();
    (*jetsak8_muonFraction_).clear();
    (*jetsak8_neutralMultiplicity_).clear();
    (*jetsak8_chargedMultiplicity_).clear();
    (*jetsak8_looseJetId_).clear();
    (*jetsak8_tightJetId_).clear();
    (*jetsak8_lepVetoJetId_).clear();
  }

  if(includeElectrons_ || includeMuons_){
    (*leptons_pdgid_).clear();
    (*leptons_energy_).clear();
    (*leptons_pt_).clear();
    (*leptons_eta_).clear();
    (*leptons_phi_).clear();
    (*leptons_px_).clear();
    (*leptons_py_).clear();
    (*leptons_pz_).clear();
    (*leptons_charge_).clear();
    (*leptons_vx_).clear();
    (*leptons_vy_).clear();
    (*leptons_vz_).clear();
    (*leptons_looseId_).clear();
    (*leptons_mediumId_).clear();
    (*leptons_tightId_).clear();
    (*leptons_pfIsoMedium_).clear();
    (*leptons_miniIsoTight_).clear();
    (*leptons_pfIsoVeryTight_).clear();
    (*leptons_pfIso_).clear();
    (*leptons_tkIso_).clear();
  }

  if(includePhotons_){
    (*photons_pt_).clear();
    (*photons_eta_).clear();
    (*photons_phi_).clear();
    (*photons_energy_).clear();
    (*photons_r9_).clear();
    (*photons_looseId_).clear();
    (*photons_mediumId_).clear();
    (*photons_tightId_).clear();
  }

  if(includePPSInfo_){
    (*protonsArm_).clear();
    (*protonsStation_).clear();
    (*protonsRP_).clear();
    (*protonsX_).clear();
    (*protonsXUnc_).clear();
    (*protonsY_).clear();
    (*protonsYUnc_).clear();
    if(!includeProtonsReco_){
      (*protonsIsValid_).clear();
      (*protonsTx_).clear();
      (*protonsTxUnc_).clear();
      (*protonsTy_).clear();
      (*protonsTyUnc_).clear();
      (*protonsTime_).clear();
      (*protonsTimeUnc_).clear();
    }
  }

  if(includeProtonsReco_){
    (*singleProtonArm_).clear();
    (*singleProtonStation_).clear();
    (*singleProtonPot_).clear();
    (*singleProtonXi_).clear();
    (*singleProtonT_).clear();
    (*singleProtonThetaX_).clear();
    (*singleProtonThetaY_).clear();
    (*multiProtonArm_).clear();
    (*multiProtonXi_).clear();
    (*multiProtonT_).clear();
    (*multiProtonTime_).clear();
    (*multiProtonTimeError_).clear();
    (*multiProtonVz_).clear();
    (*multiProtonThetaX_).clear();
    (*multiProtonThetaY_).clear();
  }

}

// ------------ running over genparticles -----
void MissingMassSearches::fetchGenParticles(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  GenParticlesEvent *genpart;
  genpart = new GenParticlesEvent(iEvent, iSetup, GenPartToken_);
  genpartVec = genpart->GetGenParticles();

  GenProtonsEvent *genproton;
  genproton = new GenProtonsEvent(iEvent, iSetup, GenPartToken_);
  genprotonsVec = genproton->GetGenProtons();

  int nGenMuons_cnt=0; 
  int nGenElectrons_cnt=0;

  for ( const auto& genpartEvt: genpartVec) {
    if(fabs(genpartEvt->pdgId()) == 11) ++nGenElectrons_cnt;
    if(fabs(genpartEvt->pdgId()) == 13) ++nGenMuons_cnt;
  }

  *nGenMuons_ = nGenMuons_cnt;
  *nGenElectrons_ = nGenElectrons_cnt;
  *nGenParticles_=genpartVec.size();
  *nGenProtons_=genprotonsVec.size();

}

// ------------- running over genmet ----------
void MissingMassSearches::fetchGenMET(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  GenMETEvent *genmissinget;
  genmissinget = new GenMETEvent(iEvent, iSetup, GenMETToken_);
  genmet = genmissinget->GetGenMET();

  *genmisset_ = genmet->et();
  *genmisset_phi_ = genmet->phi();

}

// ------------- running over genjets --------
void MissingMassSearches::fetchGenJets(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  GenJetsEvent *genjetAlgoA;
  genjetAlgoA = new GenJetsEvent(iEvent, iSetup, GenJetTokenA_);
  genjetsVecA = genjetAlgoA->GetGenJets();

  GenJetsEvent *genjetAlgoB;
  genjetAlgoB = new GenJetsEvent(iEvent, iSetup, GenJetTokenB_);
  genjetsVecB = genjetAlgoB->GetGenJets();

}

// ------------ running over muons  ------------
void MissingMassSearches::fetchMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  MuonsEvent *mu;
  mu = new MuonsEvent(iEvent, iSetup, muonToken_);
  muonsVec = mu->GetMuons();
  *nMuons_=muonsVec.size();

  int nMuons_looseId_cnt = 0;
  int nMuons_mediumId_cnt = 0;
  int nMuons_tightId_cnt = 0;
  int nMuons_pfIsoMedium_cnt = 0;
  int nMuons_miniIsoTight_cnt = 0;
  int nMuons_pfIsoVeryTight_cnt = 0;
  for ( const auto& leptonEvt: muonsVec) {
    if(leptonEvt->passed(reco::Muon::CutBasedIdLoose)) ++nMuons_looseId_cnt;
    if(leptonEvt->passed(reco::Muon::CutBasedIdMedium)) ++nMuons_mediumId_cnt;
    if(leptonEvt->passed(reco::Muon::CutBasedIdTight)) ++nMuons_tightId_cnt;
    if(leptonEvt->passed(reco::Muon::PFIsoMedium)) ++nMuons_pfIsoMedium_cnt;
    if(leptonEvt->passed(reco::Muon::MiniIsoTight)) ++nMuons_miniIsoTight_cnt;
    if(leptonEvt->passed(reco::Muon::PFIsoVeryTight)) ++nMuons_pfIsoVeryTight_cnt;
  }

  *nMuons_looseId_ = nMuons_looseId_cnt;
  *nMuons_mediumId_ = nMuons_mediumId_cnt;
  *nMuons_tightId_ = nMuons_tightId_cnt;
  *nMuons_pfIsoMedium_ = nMuons_pfIsoMedium_cnt;
  *nMuons_miniIsoTight_ = nMuons_miniIsoTight_cnt;
  *nMuons_pfIsoVeryTight_ = nMuons_pfIsoVeryTight_cnt;

}

// ------------ running over electrons  ------------
void MissingMassSearches::fetchElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  ElectronsEvent *ele;
  ele = new ElectronsEvent(iEvent, iSetup, eleToken_);
  electronsVec = ele->GetElectrons();

  *nElectrons_=electronsVec.size();

  int nElectrons_looseId_cnt = 0;
  int nElectrons_mediumId_cnt = 0;
  int nElectrons_tightId_cnt = 0;

  for ( const auto& leptonEvt: electronsVec) {
    if(leptonEvt->electronID("cutBasedElectronID-Fall17-94X-V1-loose")) ++nElectrons_looseId_cnt;
    if(leptonEvt->electronID("cutBasedElectronID-Fall17-94X-V1-medium")) ++nElectrons_mediumId_cnt;
    if(leptonEvt->electronID("cutBasedElectronID-Fall17-94X-V1-tight")) ++nElectrons_tightId_cnt;
  }

  *nElectrons_looseId_ = nElectrons_looseId_cnt;
  *nElectrons_mediumId_ = nElectrons_mediumId_cnt;
  *nElectrons_tightId_ = nElectrons_tightId_cnt;

}

// ------------ running over jets  ------------
void MissingMassSearches::fetchJets(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  JetsEvent *jetAlgoA;
  jetAlgoA = new JetsEvent(iEvent, iSetup, patjetTokenA_);
  jetsVecA = jetAlgoA->GetJets();

  JetsEvent *jetAlgoB;
  jetAlgoB = new JetsEvent(iEvent, iSetup, patjetTokenB_);
  jetsVecB = jetAlgoB->GetJets();

}

// ------------ running over photons  ------------
void MissingMassSearches::fetchPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  PhotonsEvent *photon;
  photon = new PhotonsEvent(iEvent, iSetup, photonToken_);
  photonsVec = photon->GetPhotons();

}

// ------------ running over MET  ------------
void MissingMassSearches::fetchMET(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  METEvent *missinget;
  missinget = new METEvent(iEvent, iSetup, patmetToken_);
  met = missinget->GetMET();

  *misset_ = met->et();
  *misset_phi_ = met->phi();

}

// ------------ running over vertices  ------------
void MissingMassSearches::fetchVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  VerticesEvent *vtx;
  vtx = new VerticesEvent(iEvent, iSetup, vertexToken_);
  vtxVec = vtx->GetVertices();

  if(includeVertices_){
    for ( const auto& vtxEvt: vtxVec){
      (*vertex_x_).push_back(vtxEvt->x());
      (*vertex_y_).push_back(vtxEvt->y());
      (*vertex_z_).push_back(vtxEvt->z());
      (*vertex_ntrack_).push_back(vtxEvt->nTracks());
      (*vertex_chi2_).push_back(vtxEvt->chi2());
      (*vertex_ndof_).push_back(vtxEvt->ndof());
    }
  }

}

// ------------ running over PF  ------------
void MissingMassSearches::fetchPF(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool debug = false;

  ParticleFlowEvent *packed;
  packed = new ParticleFlowEvent(iEvent, iSetup, packedToken_);
  packedVec = packed->GetPackedFlow();

  ParticleFlowEvent *pf;
  pf = new ParticleFlowEvent(iEvent, iSetup, pfToken_);
  pfVec = pf->GetParticleFlow();

  double sumEnergyPF = 0;
  double pfSlice0P_energy = 0;
  double pfSlice0P_pt = 0;
  double pfSlice1P_energy = 0;
  double pfSlice1P_pt = 0;
  double pfSlice2P_energy = 0;
  double pfSlice2P_pt = 0;
  double pfSlice3P_energy = 0;
  double pfSlice3P_pt = 0;
  double pfSlice4P_energy = 0;
  double pfSlice4P_pt = 0;
  double pfSlice0N_energy = 0;
  double pfSlice0N_pt = 0;
  double pfSlice1N_energy = 0;
  double pfSlice1N_pt = 0;
  double pfSlice2N_energy = 0;
  double pfSlice2N_pt = 0;
  double pfSlice3N_energy = 0;
  double pfSlice3N_pt = 0;
  double pfSlice4N_energy = 0;
  double pfSlice4N_pt = 0;
  int npfSlice0P = 0;
  int npfSlice1P = 0;
  int npfSlice2P = 0;
  int npfSlice3P = 0;
  int npfSlice4P = 0;
  int npfSlice0N = 0;
  int npfSlice1N = 0;
  int npfSlice2N = 0;
  int npfSlice3N = 0;
  int npfSlice4N = 0;
  int nmuonSlice0P = 0;
  int nmuonSlice1P = 0;
  int nmuonSlice2P = 0;
  int nmuonSlice3P = 0;
  int nmuonSlice4P = 0;
  int nmuonSlice0N = 0;
  int nmuonSlice1N = 0;
  int nmuonSlice2N = 0;
  int nmuonSlice3N = 0;
  int nmuonSlice4N = 0;
  int nelectronSlice0P = 0;
  int nelectronSlice1P = 0;
  int nelectronSlice2P = 0;
  int nelectronSlice3P = 0;
  int nelectronSlice4P = 0;
  int nelectronSlice0N = 0;
  int nelectronSlice1N = 0;
  int nelectronSlice2N = 0;
  int nelectronSlice3N = 0;
  int nelectronSlice4N = 0;
  double pfHFHadPlus_energy = 0;
  double pfHFHadMinus_energy = 0;
  double pfHFEmPlus_energy = 0;
  double pfHFEmMinus_energy = 0;
  double pfHFHadPlus_pt = 0;
  double pfHFHadMinus_pt = 0;
  double pfHFEmPlus_pt = 0;
  double pfHFEmMinus_pt = 0;
  double pfEcalMinus_energy = 0;
  double pfHcalMinus_energy = 0;
  double pfHadNeutralMinus_energy = 0;
  double pfEcalPlus_energy = 0;
  double pfHcalPlus_energy = 0;
  double pfHadNeutralPlus_energy = 0;
  double SumChargedPFMultiPV_pt_Loose = 0;
  double SumChargedPFMultiPV_pt_Tight = 0;
  double SumChargedPFMultiPV_pt_UsedInFit = 0;
  double SumChargedPFMultiPV_pt_Tight_Fit = 0;
  int nChargedPFMultiPV_Loose = 0;
  int nChargedPFMultiPV_Tight = 0;
  int nChargedPFMultiPV_UsedInFit = 0;
  int nChargedPFMultiPV_Tight_Fit = 0;
  int nPFPlus = 0;
  int nPFMinus = 0;
  int nPFHFHadPlus = 0;
  int nPFHFHadMinus = 0;
  int nPFHFEmPlus = 0;
  int nPFHFEmMinus = 0;
  math::XYZTLorentzVector Vertex_pt_Loose(0.,0.,0.,0.);
  math::XYZTLorentzVector Vertex_pt_Tight(0.,0.,0.,0.);
  math::XYZTLorentzVector Vertex_pt_UsedInFit(0.,0.,0.,0.);
  math::XYZTLorentzVector Vertex_pt_Tight_Fit(0.,0.,0.,0.);

  if(!(tier_.find("mini") != std::string::npos || tier_.find("MINI") != std::string::npos)){

    for ( const auto& pfEvt: pfVec){
      sumEnergyPF+=pfEvt->energy();

      if(pfEvt->eta()>0){
	pfEcalPlus_energy+=pfEvt->ecalEnergy();
	pfHcalPlus_energy+=pfEvt->hcalEnergy();
	pfHadNeutralPlus_energy+=pfEvt->hoEnergy();
	++nPFPlus;
      }
      else{
	pfEcalMinus_energy+=pfEvt->ecalEnergy();
	pfHcalMinus_energy+=pfEvt->hcalEnergy();
	pfHadNeutralMinus_energy+=pfEvt->hoEnergy();
	++nPFMinus;
      }
      // Hadronic HF energy
      if(pfEvt->particleId()==reco::PFCandidate::h_HF){
	if(pfEvt->eta()>0){
	  pfHFHadPlus_energy += pfEvt->energy();
	  pfHFHadPlus_pt += pfEvt->pt();
	  ++nPFHFHadPlus;
	}else{
	  pfHFHadMinus_energy += pfEvt->energy();
	  pfHFHadMinus_pt += pfEvt->pt();
	  ++nPFHFHadMinus;
	}
      }
      // Eletromagnetic HF energy
      if(pfEvt->particleId()==reco::PFCandidate::egamma_HF){
	if(pfEvt->eta()>0){
	  pfHFEmPlus_energy += pfEvt->energy();
	  pfHFEmPlus_pt += pfEvt->pt();
	  ++nPFHFEmPlus;
	}else{
	  pfHFEmMinus_energy += pfEvt->energy();
	  pfHFEmMinus_pt += pfEvt->pt();
	  ++nPFHFEmMinus;
	}
      }

      if(pfEvt->eta()>=0. && pfEvt->eta()<1.){
	pfSlice0P_energy += pfEvt->energy();
	pfSlice0P_pt += pfEvt->pt();
	++npfSlice0P;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice0P;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice0P;
	}
      }

      if(pfEvt->eta()>=1. && pfEvt->eta()<2.){
	pfSlice1P_energy += pfEvt->energy();
	pfSlice1P_pt += pfEvt->pt();
	++npfSlice1P;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice1P;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice1P;
	}
      }

      if(pfEvt->eta()>=2. && pfEvt->eta()<3.){
	pfSlice2P_energy += pfEvt->energy();
	pfSlice2P_pt += pfEvt->pt();
	++npfSlice2P;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice2P;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice2P;
	}
      }

      if(pfEvt->eta()>=3. && pfEvt->eta()<4.){
	pfSlice3P_energy += pfEvt->energy();
	pfSlice3P_pt += pfEvt->pt();
	++npfSlice3P;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice3P;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice3P;
	}
      }

      if(pfEvt->eta()>=4.){
	pfSlice4P_energy += pfEvt->energy();
	pfSlice4P_pt += pfEvt->pt();
	++npfSlice4P;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice4P;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice4P;
	}
      }

      if(pfEvt->eta()<0. && pfEvt->eta()>-1.){
	pfSlice0N_energy += pfEvt->energy();
	pfSlice0N_pt += pfEvt->pt();
	++npfSlice0N;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice0N;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice0N;
	}
      }

      if(pfEvt->eta()<=-1. && pfEvt->eta()>-2.){
	pfSlice1N_energy += pfEvt->energy();
	pfSlice1N_pt += pfEvt->pt();
	++npfSlice1N;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice1N;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice1N;
	}
      }

      if(pfEvt->eta()<=-2. && pfEvt->eta()>-3.){
	pfSlice2N_energy += pfEvt->energy();
	pfSlice2N_pt += pfEvt->pt();
	++npfSlice2N;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice2N;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice2N;
	}
      }

      if(pfEvt->eta()<=-3. && pfEvt->eta()>-4.){
	pfSlice3N_energy += pfEvt->energy();
	pfSlice3N_pt += pfEvt->pt();
	++npfSlice3N;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice3N;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice3N;
	}
      }

      if(pfEvt->eta()<=-4.){
	pfSlice4N_energy += pfEvt->energy();
	pfSlice4N_pt += pfEvt->pt();
	++npfSlice4N;
	if(pfEvt->particleId()==reco::PFCandidate::mu){
	  ++nmuonSlice4N;
	}
	if(pfEvt->particleId()==reco::PFCandidate::e){
	  ++nelectronSlice4N;
	}
      }

    } // end of loop

    *nPF_=pfVec.size();
    *SumPF_energy_=sumEnergyPF;
    *pfHFHadPlus_energy_=pfHFHadPlus_energy;
    *pfHFHadPlus_pt_=pfHFHadPlus_pt;
    *pfHFHadMinus_energy_=pfHFHadMinus_energy;
    *pfHFHadMinus_pt_=pfHFHadMinus_pt;
    *pfHFEmPlus_energy_=pfHFEmPlus_energy;
    *pfHFEmPlus_pt_=pfHFEmPlus_pt;
    *pfHFEmMinus_energy_=pfHFEmMinus_energy;
    *pfHFEmMinus_pt_=pfHFEmMinus_pt;
    *pfEcalMinus_energy_=pfEcalMinus_energy;
    *pfEcalPlus_energy_=pfEcalPlus_energy;
    *pfHcalMinus_energy_=pfHcalMinus_energy;
    *pfHcalPlus_energy_=pfHcalPlus_energy;
    *pfHadNeutralMinus_energy_=pfHadNeutralMinus_energy;
    *pfHadNeutralPlus_energy_=pfHadNeutralPlus_energy;
    *nPFPlus_=nPFPlus;
    *nPFMinus_=nPFMinus;
    *nPFHFHadPlus_=nPFHFHadPlus;
    *nPFHFHadMinus_=nPFHFHadMinus;
    *nPFHFEmPlus_=nPFHFEmPlus;
    *nPFHFEmMinus_=nPFHFEmMinus;
    *pfSlice0P_energy_ = pfSlice0P_energy;
    *pfSlice0P_pt_ = pfSlice0P_pt;
    *pfSlice1P_energy_ = pfSlice1P_energy;
    *pfSlice1P_pt_ = pfSlice1P_pt;
    *pfSlice2P_energy_ = pfSlice2P_energy;
    *pfSlice2P_pt_ = pfSlice2P_pt;
    *pfSlice3P_energy_ = pfSlice3P_energy;
    *pfSlice3P_pt_ = pfSlice3P_pt;
    *pfSlice4P_energy_ = pfSlice4P_energy;
    *pfSlice4P_pt_ = pfSlice4P_pt;
    *pfSlice0N_energy_ = pfSlice0N_energy;
    *pfSlice0N_pt_ = pfSlice0N_pt;
    *pfSlice1N_energy_ = pfSlice1N_energy;
    *pfSlice1N_pt_ = pfSlice1N_pt;
    *pfSlice2N_energy_ = pfSlice2N_energy;
    *pfSlice2N_pt_ = pfSlice2N_pt;
    *pfSlice3N_energy_ = pfSlice3N_energy;
    *pfSlice3N_pt_ = pfSlice3N_pt;
    *pfSlice4N_energy_ = pfSlice4N_energy;
    *pfSlice4N_pt_ = pfSlice4N_pt;
    *npfSlice0P_=npfSlice0P;
    *npfSlice1P_=npfSlice1P;
    *npfSlice2P_=npfSlice2P;
    *npfSlice3P_=npfSlice3P;
    *npfSlice4P_=npfSlice4P;
    *npfSlice0N_=npfSlice0N;
    *npfSlice1N_=npfSlice1N;
    *npfSlice2N_=npfSlice2N;
    *npfSlice3N_=npfSlice3N;
    *npfSlice4N_=npfSlice4N;
    *nmuonSlice0P_=nmuonSlice0P;
    *nmuonSlice1P_=nmuonSlice1P;
    *nmuonSlice2P_=nmuonSlice2P;
    *nmuonSlice3P_=nmuonSlice3P;
    *nmuonSlice4P_=nmuonSlice4P;
    *nmuonSlice0N_=nmuonSlice0N;
    *nmuonSlice1N_=nmuonSlice1N;
    *nmuonSlice2N_=nmuonSlice2N;
    *nmuonSlice3N_=nmuonSlice3N;
    *nmuonSlice4N_=nmuonSlice4N;
    *nelectronSlice0P_=nelectronSlice0P;
    *nelectronSlice1P_=nelectronSlice1P;
    *nelectronSlice2P_=nelectronSlice2P;
    *nelectronSlice3P_=nelectronSlice3P;
    *nelectronSlice4P_=nelectronSlice4P;
    *nelectronSlice0N_=nelectronSlice0N;
    *nelectronSlice1N_=nelectronSlice1N;
    *nelectronSlice2N_=nelectronSlice2N;
    *nelectronSlice3N_=nelectronSlice3N;
    *nelectronSlice4N_=nelectronSlice4N;

  }else{

    for (const auto& packedEvt: packedVec){

      sumEnergyPF+=packedEvt->energy();
      if(fabs(packedEvt->eta())<4.7 && packedEvt->charge()!=0){
	bool passChargeSelPack(packedEvt->pt()>0.9 && fabs(packedEvt->eta())<2.5);
	const pat::PackedCandidate::PVAssoc pvassocpack=packedEvt->fromPV();
	if(passChargeSelPack && pvassocpack==pat::PackedCandidate::PVLoose){
	  ++nChargedPFMultiPV_Loose;
	  Vertex_pt_Loose += packedEvt->p4();
	}
	if(passChargeSelPack && pvassocpack==pat::PackedCandidate::PVTight){
	  ++nChargedPFMultiPV_Tight;
	  Vertex_pt_Tight += packedEvt->p4();
	}
	if(passChargeSelPack && pvassocpack==pat::PackedCandidate::PVUsedInFit){
	  ++nChargedPFMultiPV_UsedInFit;
	  Vertex_pt_UsedInFit += packedEvt->p4();
	}
	if(passChargeSelPack && pvassocpack>=pat::PackedCandidate::PVTight){
	  ++nChargedPFMultiPV_Tight_Fit;
	  Vertex_pt_Tight_Fit += packedEvt->p4();
	}	
      }

      // EM Energy
      if(abs(packedEvt->pdgId()) == 11 || abs(packedEvt->pdgId()) == 22 || abs(packedEvt->pdgId()) == 2){
	// ECAL
	if(packedEvt->eta()>0. && packedEvt->eta()<3.){
	  pfEcalPlus_energy+=packedEvt->energy();
	  ++nPFPlus;
	}
	if(packedEvt->eta()<=0. && packedEvt->eta()>-3.){
	  pfEcalMinus_energy+=packedEvt->energy();
	  ++nPFMinus;
	}
	// HF
	if(packedEvt->eta()<=-3.){
	  pfHFEmMinus_energy+=packedEvt->energy();
	  pfHFEmMinus_pt+=packedEvt->pt();
	  ++nPFHFEmMinus;
	  ++nPFMinus;
	}
	if(packedEvt->eta()>=3.){
	  pfHFEmPlus_energy+=packedEvt->energy();
	  pfHFEmPlus_pt+=packedEvt->pt();
	  ++nPFHFEmPlus;
	  ++nPFPlus;
	}
      }

      // HADRONIC Energy, charged
      if(abs(packedEvt->pdgId()) == 211 || abs(packedEvt->pdgId()) == 321){
	// HCAL
	if(packedEvt->eta()>0. && packedEvt->eta()<3.){
	  pfHcalPlus_energy+=packedEvt->energy();
	  ++nPFPlus;
	}
	if(packedEvt->eta()<=0. && packedEvt->eta()>-3.){
	  pfHcalMinus_energy+=packedEvt->energy();
	  ++nPFMinus;
	}
	// HF
	if(packedEvt->eta()<=-3.){
	  pfHFHadMinus_energy+=packedEvt->energy();
	  pfHFHadMinus_pt+=packedEvt->pt();
	  ++nPFHFHadMinus;
	  ++nPFMinus;
	}
	if(packedEvt->eta()>=3.){
	  pfHFHadPlus_energy+=packedEvt->energy();
	  pfHFHadPlus_pt+=packedEvt->pt();
	  ++nPFHFHadPlus;
	  ++nPFPlus;
	}
      }

      // HADRONIC Energy, neutral
      if(packedEvt->pdgId() == 111 || packedEvt->pdgId() == 311){
	if(packedEvt->eta()>0.){
	  pfHadNeutralPlus_energy+=packedEvt->energy();
	  ++nPFPlus;
	}
	if(packedEvt->eta()<=0.){
	  pfHadNeutralMinus_energy+=packedEvt->energy();
	  ++nPFMinus;
	}
	// HF
	if(packedEvt->eta()<=-3.){
	  pfHFHadMinus_energy+=packedEvt->energy();
	  pfHFHadMinus_pt+=packedEvt->pt();
	  ++nPFHFHadMinus;
	  ++nPFMinus;
	}
	if(packedEvt->eta()>=3.){
	  pfHFHadPlus_energy+=packedEvt->energy();
	  pfHFHadPlus_pt+=packedEvt->pt();
	  ++nPFHFHadPlus;
	  ++nPFPlus;
	}
      }

      if(packedEvt->eta()>=0. && packedEvt->eta()<1.){
	pfSlice0P_energy += packedEvt->energy();
	pfSlice0P_pt += packedEvt->pt();
	++npfSlice0P;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice0P;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice0P;
	}
      }

      if(packedEvt->eta()>=1. && packedEvt->eta()<2.){
	pfSlice1P_energy += packedEvt->energy();
	pfSlice1P_pt += packedEvt->pt();
	++npfSlice1P;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice1P;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice1P;
	}
      }

      if(packedEvt->eta()>=2. && packedEvt->eta()<3.){
	pfSlice2P_energy += packedEvt->energy();
	pfSlice2P_pt += packedEvt->pt();
	++npfSlice2P;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice2P;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice2P;
	}
      }

      if(packedEvt->eta()>=3. && packedEvt->eta()<4.){
	pfSlice3P_energy += packedEvt->energy();
	pfSlice3P_pt += packedEvt->pt();
	++npfSlice3P;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice3P;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice3P;
	}
      }

      if(packedEvt->eta()>=4.){
	pfSlice4P_energy += packedEvt->energy();
	pfSlice4P_pt += packedEvt->pt();
	++npfSlice4P;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice4P;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice4P;
	}
      }

      if(packedEvt->eta()<0. && packedEvt->eta()>-1.){
	pfSlice0N_energy += packedEvt->energy();
	pfSlice0N_pt += packedEvt->pt();
	++npfSlice0N;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice0N;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice0N;
	}
      }

      if(packedEvt->eta()<=-1. && packedEvt->eta()>-2.){
	pfSlice1N_energy += packedEvt->energy();
	pfSlice1N_pt += packedEvt->pt();
	++npfSlice1N;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice1N;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice1N;
	}
      }

      if(packedEvt->eta()<=-2. && packedEvt->eta()>-3.){
	pfSlice2N_energy += packedEvt->energy();
	pfSlice2N_pt += packedEvt->pt();
	++npfSlice2N;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice2N;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice2N;
	}
      }

      if(packedEvt->eta()<=-3. && packedEvt->eta()>-4.){
	pfSlice3N_energy += packedEvt->energy();
	pfSlice3N_pt += packedEvt->pt();
	++npfSlice3N;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice3N;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice3N;
	}
      }

      if(packedEvt->eta()<=-4.){
	pfSlice4N_energy += packedEvt->energy();
	pfSlice4N_pt += packedEvt->pt();
	++npfSlice4N;
	if(fabs(packedEvt->pdgId())==13){
	  ++nmuonSlice4N;
	}
	if(fabs(packedEvt->pdgId())==11){
	  ++nelectronSlice4N;
	}
      }

    } // end PF loop

    SumChargedPFMultiPV_pt_Loose = sqrt(Vertex_pt_Loose.Px()*Vertex_pt_Loose.Px() + Vertex_pt_Loose.Py()*Vertex_pt_Loose.Py());
    SumChargedPFMultiPV_pt_Tight = sqrt(Vertex_pt_Tight.Px()*Vertex_pt_Tight.Px() + Vertex_pt_Tight.Py()*Vertex_pt_Tight.Py());
    SumChargedPFMultiPV_pt_UsedInFit = sqrt(Vertex_pt_UsedInFit.Px()*Vertex_pt_UsedInFit.Px() + Vertex_pt_UsedInFit.Py()*Vertex_pt_UsedInFit.Py());
    SumChargedPFMultiPV_pt_Tight_Fit = sqrt(Vertex_pt_Tight_Fit.Px()*Vertex_pt_Tight_Fit.Px() + Vertex_pt_Tight_Fit.Py()*Vertex_pt_Tight_Fit.Py());
    *nPF_=packedVec.size();
    *SumPF_energy_=sumEnergyPF;
    *pfHFHadPlus_energy_=pfHFHadPlus_energy;
    *pfHFHadPlus_pt_=pfHFHadPlus_pt;
    *pfHFHadMinus_energy_=pfHFHadMinus_energy;
    *pfHFHadMinus_pt_=pfHFHadMinus_pt;
    *pfHFEmPlus_energy_=pfHFEmPlus_energy;
    *pfHFEmPlus_pt_=pfHFEmPlus_pt;
    *pfHFEmMinus_energy_=pfHFEmMinus_energy;
    *pfHFEmMinus_pt_=pfHFEmMinus_pt;
    *pfEcalMinus_energy_=pfEcalMinus_energy;
    *pfEcalPlus_energy_=pfEcalPlus_energy;
    *pfHcalMinus_energy_=pfHcalMinus_energy;
    *pfHcalPlus_energy_=pfHcalPlus_energy;
    *pfHadNeutralMinus_energy_=pfHadNeutralMinus_energy;
    *pfHadNeutralPlus_energy_=pfHadNeutralPlus_energy;
    *SumChargedPFMultiPV_pt_Loose_=SumChargedPFMultiPV_pt_Loose;
    *SumChargedPFMultiPV_pt_Tight_=SumChargedPFMultiPV_pt_Tight;
    *SumChargedPFMultiPV_pt_UsedInFit_=SumChargedPFMultiPV_pt_UsedInFit;
    *SumChargedPFMultiPV_pt_Tight_Fit_=SumChargedPFMultiPV_pt_Tight_Fit;
    *nChargedPFMultiPV_Loose_=nChargedPFMultiPV_Loose;
    *nChargedPFMultiPV_Tight_=nChargedPFMultiPV_Tight;
    *nChargedPFMultiPV_UsedInFit_=nChargedPFMultiPV_UsedInFit;
    *nChargedPFMultiPV_Tight_Fit_=nChargedPFMultiPV_Tight_Fit;
    *nPFPlus_=nPFPlus;
    *nPFMinus_=nPFMinus;
    *nPFHFHadPlus_=nPFHFHadPlus;
    *nPFHFHadMinus_=nPFHFHadMinus;
    *nPFHFEmPlus_=nPFHFEmPlus;
    *nPFHFEmMinus_=nPFHFEmMinus;
    *pfSlice0P_energy_ = pfSlice0P_energy;
    *pfSlice0P_pt_ = pfSlice0P_pt;
    *pfSlice1P_energy_ = pfSlice1P_energy;
    *pfSlice1P_pt_ = pfSlice1P_pt;
    *pfSlice2P_energy_ = pfSlice2P_energy;
    *pfSlice2P_pt_ = pfSlice2P_pt;
    *pfSlice3P_energy_ = pfSlice3P_energy;
    *pfSlice3P_pt_ = pfSlice3P_pt;
    *pfSlice4P_energy_ = pfSlice4P_energy;
    *pfSlice4P_pt_ = pfSlice4P_pt;
    *pfSlice0N_energy_ = pfSlice0N_energy;
    *pfSlice0N_pt_ = pfSlice0N_pt;
    *pfSlice1N_energy_ = pfSlice1N_energy;
    *pfSlice1N_pt_ = pfSlice1N_pt;
    *pfSlice2N_energy_ = pfSlice2N_energy;
    *pfSlice2N_pt_ = pfSlice2N_pt;
    *pfSlice3N_energy_ = pfSlice3N_energy;
    *pfSlice3N_pt_ = pfSlice3N_pt;
    *pfSlice4N_energy_ = pfSlice4N_energy;
    *pfSlice4N_pt_ = pfSlice4N_pt;
    *npfSlice0P_=npfSlice0P;
    *npfSlice1P_=npfSlice1P;
    *npfSlice2P_=npfSlice2P;
    *npfSlice3P_=npfSlice3P;
    *npfSlice4P_=npfSlice4P;
    *npfSlice0N_=npfSlice0N;
    *npfSlice1N_=npfSlice1N;
    *npfSlice2N_=npfSlice2N;
    *npfSlice3N_=npfSlice3N;
    *npfSlice4N_=npfSlice4N;
    *nmuonSlice0P_=nmuonSlice0P;
    *nmuonSlice1P_=nmuonSlice1P;
    *nmuonSlice2P_=nmuonSlice2P;
    *nmuonSlice3P_=nmuonSlice3P;
    *nmuonSlice4P_=nmuonSlice4P;
    *nmuonSlice0N_=nmuonSlice0N;
    *nmuonSlice1N_=nmuonSlice1N;
    *nmuonSlice2N_=nmuonSlice2N;
    *nmuonSlice3N_=nmuonSlice3N;
    *nmuonSlice4N_=nmuonSlice4N;
    *nelectronSlice0P_=nelectronSlice0P;
    *nelectronSlice1P_=nelectronSlice1P;
    *nelectronSlice2P_=nelectronSlice2P;
    *nelectronSlice3P_=nelectronSlice3P;
    *nelectronSlice4P_=nelectronSlice4P;
    *nelectronSlice0N_=nelectronSlice0N;
    *nelectronSlice1N_=nelectronSlice1N;
    *nelectronSlice2N_=nelectronSlice2N;
    *nelectronSlice3N_=nelectronSlice3N;
    *nelectronSlice4N_=nelectronSlice4N;
  }

  if(debug){
    std::cout << "PF HF CMS+ Had Energy [GeV]: " << pfHFHadPlus_energy << std::endl;
    std::cout << "PF HF CMS+ Had pT [GeV]: " << pfHFHadPlus_pt << std::endl;
    std::cout << "PF HF CMS- Had Energy [GeV]: " << pfHFHadMinus_energy << std::endl;
    std::cout << "PF HF CMS- Had pT [GeV]: " << pfHFHadMinus_pt << std::endl;
    std::cout << "PF HF CMS+ Em Energy [GeV]: " << pfHFEmPlus_energy << std::endl;
    std::cout << "PF HF CMS+ Em pT [GeV]: " << pfHFEmPlus_pt << std::endl;
    std::cout << "PF HF CMS- Em Energy [GeV]: " << pfHFEmMinus_energy << std::endl;
    std::cout << "PF HF CMS- Em pT [GeV]: " << pfHFEmMinus_pt << std::endl;
    std::cout << "Ecal CMS- Energy [GeV]: " << pfEcalMinus_energy << std::endl;
    std::cout << "Ecal CMS+ Energy [GeV]: " << pfEcalPlus_energy << std::endl;
    std::cout << "Hcal CMS- Energy [GeV]: " << pfHcalMinus_energy << std::endl;
    std::cout << "Hcal CMS+ Energy [GeV]: " << pfHcalPlus_energy << std::endl;
    std::cout << "Hneutral CMS- Energy [GeV]: " << pfHadNeutralMinus_energy << std::endl;
    std::cout << "Hneutral CMS+ Energy [GeV]: " << pfHadNeutralPlus_energy << std::endl;
    std::cout << "nPF CMS+: " << nPFPlus << std::endl;
    std::cout << "nPF CMS-: " << nPFMinus << std::endl;
    std::cout << "nPF HF Had CMS+: " << nPFHFHadPlus << std::endl;
    std::cout << "nPF HF Had CMS-: " << nPFHFHadMinus << std::endl;
    std::cout << "nPF HF Em CMS+: " << nPFHFEmPlus << std::endl;
    std::cout << "nPF HF Em CMS-: " << nPFHFEmMinus << std::endl;
    std::cout << "Energy [GeV], Slice CMS+: (" << pfSlice0P_energy << ", " << pfSlice1P_energy << ", " << pfSlice2P_energy << ", " << pfSlice3P_energy << ", " << pfSlice4P_energy << ")" << std::endl;
    std::cout << "pT [GeV], Slice CMS+: (" << pfSlice0P_pt << ", " << pfSlice1P_pt << ", " << pfSlice2P_pt << ", " << pfSlice3P_pt << ", " << pfSlice4P_pt << ")" << std::endl;
    std::cout << "Multiplicity, Slice CMS+: (" << npfSlice0P << ", " << npfSlice1P << ", " << npfSlice2P << ", " << npfSlice3P << ", " << npfSlice4P << ")" << std::endl;
    std::cout << "Energy [GeV], Slice CMS-: (" << pfSlice0N_energy << ", " << pfSlice1N_energy << ", " << pfSlice2N_energy << ", " << pfSlice3N_energy << ", " << pfSlice4N_energy << ")" << std::endl;
    std::cout << "pT [GeV], Slice CMS-: (" << pfSlice0N_pt << ", " << pfSlice1N_pt << ", " << pfSlice2N_pt << ", " << pfSlice3N_pt << ", " << pfSlice4N_pt << ")" << std::endl;
    std::cout << "Multiplicity, Slice CMS-: (" << npfSlice0N << ", " << npfSlice1N << ", " << npfSlice2N << ", " << npfSlice3N << ", " << npfSlice4N << ")" << std::endl;
    std::cout << "Multiplicity, Muon CMS+: (" << nmuonSlice0P << ", " << nmuonSlice1P << ", " << nmuonSlice2P << ", " << nmuonSlice3P << ", " << nmuonSlice4P << ")" << std::endl;
    std::cout << "Multiplicity, Muon CMS-: (" << nmuonSlice0N << ", " << nmuonSlice1N << ", " << nmuonSlice2N << ", " << nmuonSlice3N << ", " << nmuonSlice4N << ")" << std::endl;
    std::cout << "Multiplicity, Electron CMS+: (" << nelectronSlice0P << ", " << nelectronSlice1P << ", " << nelectronSlice2P << ", " << nelectronSlice3P << ", " << nelectronSlice4P << ")" << std::endl;
    std::cout << "Multiplicity, Electron CMS-: (" << nelectronSlice0N << ", " << nelectronSlice1N << ", " << nelectronSlice2N << ", " << nelectronSlice3N << ", " << nelectronSlice4N << ")" << std::endl;
  }

}

// ------------ fetching over LHC info ----------------
void MissingMassSearches::fetchLHCInfo(const edm::EventSetup& iSetup){

  bool debug = false;

  edm::ESHandle<LHCInfo> hLHCInfo;
  iSetup.get<LHCInfoRcd>().get(lhcInfoLabel_, hLHCInfo);

  *xangle_ = hLHCInfo->crossingAngle();
  *betastar_ = hLHCInfo->betaStar();
  *bunchesb1_ = hLHCInfo->bunchesInBeam1();
  *bunchesb2_ = hLHCInfo->bunchesInBeam2();
  *intensityb1_ = hLHCInfo->intensityForBeam1();
  *intensityb2_ = hLHCInfo->intensityForBeam2();
  *instlumi_ = hLHCInfo->instLumi();

  if(debug){
    std::cout << "Colliding bunches: " << hLHCInfo->collidingBunches() << std::endl;
    std::cout << "# bunches Beam1: " << hLHCInfo->bunchesInBeam1() << std::endl;
    std::cout << "# bunches Beam2: " << hLHCInfo->bunchesInBeam2() << std::endl;
    std::cout << "Intensity beam1: " << hLHCInfo->intensityForBeam1() << std::endl;
    std::cout << "Intensity beam2: " << hLHCInfo->intensityForBeam2() << std::endl;
    std::cout << "Inst. Luminosity: " << hLHCInfo->instLumi() << std::endl;
    std::cout << "XAngle (LHC) [murad]: " << hLHCInfo->crossingAngle() << std::endl;
    std::cout << "Beta Star: " << hLHCInfo->betaStar() << std::endl;
  }

}

// ------------ running over reco protons  ------------
void MissingMassSearches::fetchProtonsReco(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  ProtonsEvent *protonsSingle;
  protonsSingle = new ProtonsEvent(iEvent, iSetup, protonSingleToken_);
  protonSingleVec = protonsSingle->GetProtonsReco();

  ProtonsEvent *protonsMulti;
  protonsMulti = new ProtonsEvent(iEvent, iSetup, protonMultiToken_);
  protonMultiVec = protonsMulti->GetProtonsReco();

  for ( const auto& protonEvt: protonSingleVec) {
    const CTPPSDetId det_id_single((*protonEvt->contributingLocalTracks().begin())->getRPId());
    (*singleProtonArm_).push_back(det_id_single.arm());
    (*singleProtonStation_).push_back(det_id_single.station());
    (*singleProtonPot_).push_back(det_id_single.rp());
    (*singleProtonXi_).push_back(protonEvt->xi());
    (*singleProtonT_).push_back(protonEvt->t());
    (*singleProtonThetaX_).push_back(protonEvt->thetaX());
    (*singleProtonThetaY_).push_back(protonEvt->thetaY());
  }

  for ( const auto& protonEvt: protonMultiVec) {
    const CTPPSDetId det_id_multi((*protonEvt->contributingLocalTracks().begin())->getRPId());
    (*multiProtonArm_).push_back(det_id_multi.arm());
    (*multiProtonXi_).push_back(protonEvt->xi());
    (*multiProtonT_).push_back(protonEvt->t());
    (*multiProtonTime_).push_back(protonEvt->time());
    (*multiProtonTimeError_).push_back(protonEvt->timeError());
    (*multiProtonVz_).push_back(protonEvt->vz());
    (*multiProtonThetaX_).push_back(protonEvt->thetaX());
    (*multiProtonThetaY_).push_back(protonEvt->thetaY());
  }

}

// ------------ running over protons  ------------
void MissingMassSearches::fetchPPSInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  if(ppslite_){
    ProtonsEvent *protons;
    protons = new ProtonsEvent(iEvent, iSetup, protonToken_);
    ppsVec = protons->GetProtonsLite();
    for ( const auto& ppsEvt: ppsVec ) {
      const CTPPSDetId det_id( ppsEvt->getRPId() );
      (*protonsArm_).push_back(det_id.arm());
      (*protonsStation_).push_back(det_id.station());
      (*protonsRP_).push_back(det_id.rp());
      (*protonsX_).push_back(ppsEvt->getX());
      (*protonsXUnc_).push_back(ppsEvt->getXUnc());
      (*protonsY_).push_back(ppsEvt->getY());
      (*protonsYUnc_).push_back(ppsEvt->getYUnc());
      if(!includeProtonsReco_){
	(*protonsIsValid_).push_back(true);
	(*protonsTx_).push_back(-1);
	(*protonsTxUnc_).push_back(-1);
	(*protonsTy_).push_back(-1);
	(*protonsTyUnc_).push_back(-1);
	(*protonsTime_).push_back(ppsEvt->getTime());
	(*protonsTimeUnc_).push_back(ppsEvt->getTimeUnc());
      }
    }
  }else{
    // Not in Miniaod!
    ProtonsEvent *pixels;
    pixels = new ProtonsEvent(iEvent, iSetup, pixelsppsToken_);
    pixelsVec = pixels->GetPixels();

    // Not in Miniaod!
    ProtonsEvent *timing;
    timing = new ProtonsEvent(iEvent, iSetup, timingppsToken_);
    timingVec = timing->GetTiming();

    // Not in Miniaod!
    ProtonsEvent *strips;
    strips = new ProtonsEvent(iEvent, iSetup, stripstotemToken_);
    stripsVec = strips->GetStrips();

    for ( const auto& pixelsEvt: pixelsVec ) {
      (*protonsArm_).push_back(pixelsEvt.second.arm());
      (*protonsStation_).push_back(pixelsEvt.second.station());
      (*protonsRP_).push_back(pixelsEvt.second.rp());
      (*protonsX_).push_back(pixelsEvt.first->getX0());
      (*protonsXUnc_).push_back(pixelsEvt.first->getX0Sigma());
      (*protonsY_).push_back(pixelsEvt.first->getY0());
      (*protonsYUnc_).push_back(pixelsEvt.first->getY0Sigma());
      if(!includeProtonsReco_){
	(*protonsIsValid_).push_back(pixelsEvt.first->isValid());
	(*protonsTx_).push_back(pixelsEvt.first->getTx());
	(*protonsTxUnc_).push_back(pixelsEvt.first->getTxSigma());
	(*protonsTy_).push_back(pixelsEvt.first->getTy());
	(*protonsTyUnc_).push_back(pixelsEvt.first->getTySigma());
	(*protonsTime_).push_back(-1.);
	(*protonsTimeUnc_).push_back(-1.);
      }
    }
    for ( const auto& stripsEvt: stripsVec ) {
      (*protonsArm_).push_back(stripsEvt.second.arm());
      (*protonsStation_).push_back(stripsEvt.second.station());
      (*protonsRP_).push_back(stripsEvt.second.rp());
      (*protonsX_).push_back(stripsEvt.first->getX0());
      (*protonsXUnc_).push_back(stripsEvt.first->getX0Sigma());
      (*protonsY_).push_back(stripsEvt.first->getY0());
      (*protonsYUnc_).push_back(stripsEvt.first->getY0Sigma());
      if(!includeProtonsReco_){
	(*protonsIsValid_).push_back(stripsEvt.first->isValid());
	(*protonsTx_).push_back(stripsEvt.first->getTx());
	(*protonsTxUnc_).push_back(stripsEvt.first->getTxSigma());
	(*protonsTy_).push_back(stripsEvt.first->getTy());
	(*protonsTyUnc_).push_back(stripsEvt.first->getTySigma());
	(*protonsTime_).push_back(-1.);
	(*protonsTimeUnc_).push_back(-1.);
      }
    }
    for ( const auto& timingEvt: timingVec ) {
      (*protonsArm_).push_back(timingEvt.second.arm());
      (*protonsStation_).push_back(timingEvt.second.station());
      (*protonsRP_).push_back(timingEvt.second.rp());
      (*protonsX_).push_back(timingEvt.first->getX0());
      (*protonsXUnc_).push_back(timingEvt.first->getX0Sigma());
      (*protonsY_).push_back(timingEvt.first->getY0());
      (*protonsYUnc_).push_back(timingEvt.first->getY0Sigma());
      if(!includeProtonsReco_){
	(*protonsIsValid_).push_back(timingEvt.first->isValid());
	(*protonsTx_).push_back(-1);
	(*protonsTxUnc_).push_back(-1);
	(*protonsTy_).push_back(-1);
	(*protonsTyUnc_).push_back(-1);
	(*protonsTime_).push_back(timingEvt.first->getT());
	(*protonsTimeUnc_).push_back(timingEvt.first->getTSigma());
      }
    }
  }

}

// ------------ saving trigger  ------------
bool MissingMassSearches::fetchTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool trigger_fired = false;

  TriggerEvent *trigger;
  trigger = new TriggerEvent(iEvent, iSetup, hltPrescaleProvider_, triggerResultsToken_, triggersList_, enablePrescales_);
  triggerVec = trigger->GetTrigger();

  if(enablePrescales_){
    prescalesL1Vec = trigger->GetPrescalesL1();
    prescalesHLTVec = trigger->GetPrescalesHLT();
  }

  for ( const auto& triggerEvt: triggerVec) {
    if(triggerEvt==1) trigger_fired = true;
  }

  if(!enableTrigger_) trigger_fired = true;

  return trigger_fired;

}

// ------------ saving trigger  ------------
void MissingMassSearches::fetchDetector(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool debug = false;

  std::vector<double> energy_tower;
  std::vector<double> eta_tower;

  energy_tower.clear();
  eta_tower.clear();

  double Epz_plus=0.;  
  double Epz_minus=0.;  
  double Et_eta_plus=0.;
  double Et_eta_minus=0.;

  int nTowersHF_plus = 0;
  int nTowersHF_minus = 0;
  int nTowersHE_plus = 0;
  int nTowersHE_minus = 0;
  int nTowersHB_plus = 0;
  int nTowersHB_minus = 0;
  int nTowersEE_plus = 0;
  int nTowersEE_minus = 0;
  int nTowersEB_plus = 0;
  int nTowersEB_minus = 0;

  //Sum(E)
  double sumEHF_plus = 0.;
  double sumEHF_minus = 0.;
  double sumEHF_L_plus = 0.;
  double sumEHF_L_minus = 0.;
  double sumEHF_S_plus = 0.;
  double sumEHF_S_minus = 0.;
  double sumEHE_plus = 0.;
  double sumEHE_minus = 0.;
  double sumEHB_plus = 0.;
  double sumEHB_minus = 0.;
  double sumEEE_plus = 0.;
  double sumEEE_minus = 0.;
  double sumEEB_plus = 0.;
  double sumEEB_minus = 0.;

  // Sum(ET)
  double sumETHF_plus = 0.;
  double sumETHF_minus = 0.;
  double sumETHE_plus = 0.;
  double sumETHE_minus = 0.;
  double sumETHB_plus = 0.;
  double sumETHB_minus = 0.;
  double sumETEB_plus = 0.;
  double sumETEB_minus = 0.;
  double sumETEE_plus = 0.;
  double sumETEE_minus = 0.;

  edm::Handle<CaloTowerCollection> towersColl;
  iEvent.getByToken(calotowersToken_,towersColl);

  int towersize = towersColl->size();
  int itTower;
  int counter_tower=0;

  for(itTower=0; itTower < towersize; ++itTower){
    const CaloTower* calotower = &((*towersColl)[itTower]);

    //bool hasHCAL = false;
    bool hasHF = false;
    bool hasHE = false;
    bool hasHB = false;
    //bool hasHO = false;
    //bool hasECAL = false;
    bool hasEE = false;
    bool hasEB = false;     

    for(size_t iconst = 0; iconst < calotower->constituentsSize(); iconst++){
      DetId adetId = calotower->constituent(iconst);
      if(adetId.det()==DetId::Hcal){
	//hasHCAL = true;
	HcalDetId hcalDetId(adetId);
	if(hcalDetId.subdet()==HcalForward) hasHF = true;
	else if(hcalDetId.subdet()==HcalEndcap) hasHE = true;
	else if(hcalDetId.subdet()==HcalBarrel) hasHB = true;
	//else if(hcalDetId.subdet()==HcalOuter) hasHO = true;  
      } 
      else if(adetId.det()==DetId::Ecal){
	//hasECAL = true;
	EcalSubdetector ecalSubDet = (EcalSubdetector)adetId.subdetId();
	if(ecalSubDet == EcalEndcap) hasEE = true;
	else if(ecalSubDet == EcalBarrel) hasEB = true;
      }
    }

    int zside = calotower->zside();
    double caloTowerEnergy = calotower->energy();
    double caloTowerEmEnergy = calotower->emEnergy();
    double caloTowerHadEnergy = calotower->hadEnergy();
    double caloTowerPz = calotower->pz();
    double caloTowerEt = calotower->et();
    double caloTowerEmEt = calotower->emEt();
    double caloTowerHadEt = calotower->hadEt();
    double EHF_S = 0;
    double EHF_L = 0;

    bool CalAboveTh = false;

    if( hasHF && !hasHE )
    {
      if( caloTowerEnergy > energyThresholdHF_)
      {
	CalAboveTh = true;

	if (debug) std::cout << "HF>> " <<  calotower->id() << " HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl; 

	// Added Long and short fibers
	// emc=L-S
	// hac=2*S
	// Tot = L+S
	// S = hac/2
	// L = Tot - S

	EHF_S = caloTowerHadEnergy/2;
	EHF_L = caloTowerEnergy - caloTowerHadEnergy/2;

	if(zside >= 0)
	{
	  ++nTowersHF_plus;
	  sumEHF_plus += caloTowerEnergy;
	  sumEHF_S_plus += EHF_S;
	  sumEHF_L_plus += EHF_L;
	  sumETHF_plus += caloTowerEt;
	}
	else
	{
	  ++nTowersHF_minus;
	  sumEHF_minus += caloTowerEnergy;
	  sumEHF_S_minus += EHF_S;
	  sumEHF_L_minus += EHF_L;
	  sumETHF_minus += caloTowerEt;
	}
      }
    }
    else if( hasHE && !hasHF && !hasHB )
    {
      if( caloTowerHadEnergy > energyThresholdHE_)
      {
	CalAboveTh = true;

	if (debug) std::cout << "HE>> " <<  calotower->id() << "  HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;

	if(zside >= 0)
	{
	  ++nTowersHE_plus;
	  sumEHE_plus += caloTowerHadEnergy;
	  sumETHE_plus += caloTowerHadEt;
	}
	else
	{
	  ++nTowersHE_minus;
	  sumEHE_minus += caloTowerHadEnergy;
	  sumETHE_minus += caloTowerHadEt;
	}
      }
    }
    else if( hasHB && !hasHE )
    {
      if( caloTowerHadEnergy > energyThresholdHB_)
      {
	CalAboveTh = true;

	if (debug) std::cout << "HB>> " <<  calotower->id() << "  HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;

	if(zside >= 0)
	{
	  ++nTowersHB_plus;
	  sumEHB_plus += caloTowerHadEnergy;
	  sumETHB_plus += caloTowerHadEt;
	}
	else
	{
	  ++nTowersHB_minus;
	  sumEHB_minus += caloTowerHadEnergy;
	  sumETHB_minus += caloTowerHadEt;
	}
      }
    }

    if( hasEE && !hasEB )
    {
      if( caloTowerEmEnergy >= energyThresholdEE_)
      {

	CalAboveTh = true;

	if (debug) std::cout << "EE>> " <<  calotower->id() << "  HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl;

	if(zside >= 0)
	{
	  ++nTowersEE_plus;
	  sumEEE_plus += caloTowerEmEnergy;
	  sumETEE_plus += caloTowerEmEt;
	}
	else
	{
	  ++nTowersEE_minus;
	  sumEEE_minus += caloTowerEmEnergy;
	  sumETEE_minus += caloTowerEmEt;
	}
      }
    }
    else if( hasEB && !hasEE )
    {
      if( caloTowerEmEnergy >= energyThresholdEB_)
      {
	CalAboveTh = true;

	if (debug) std::cout << "EB>> " <<  calotower->id() << " HAD energy "     << caloTowerHadEnergy << " EM energy " << caloTowerEmEnergy << " energy " << caloTowerEnergy << std::endl; 

	if(zside >= 0)
	{
	  ++nTowersEB_plus;
	  sumEEB_plus += caloTowerEmEnergy;
	  sumETEB_plus += caloTowerEmEt;
	}
	else
	{
	  ++nTowersEB_minus;
	  sumEEB_minus += caloTowerEmEnergy;
	  sumETEB_minus += caloTowerEmEt;
	}
      }
    }

    if(CalAboveTh)
    {
      ++counter_tower;
      energy_tower.push_back(calotower->energy());
      eta_tower.push_back(calotower->eta());
      Et_eta_minus += caloTowerEt * pow(2.71,-calotower->eta());
      Et_eta_plus += caloTowerEt * pow(2.71,calotower->eta());
      Epz_plus  += caloTowerEnergy + caloTowerPz;
      Epz_minus += caloTowerEnergy - caloTowerPz;
    }

  }  ////has to close calotower loop

  if(debug){
    std::cout << "Calorimeter, Et*exp(eta): " << Et_eta_plus << std::endl;
    std::cout << "Calorimeter, Et*exp(-eta): " << Et_eta_minus << std::endl;
    std::cout << "Calorimeter, E + pz: " << Epz_plus << std::endl;
    std::cout << "Calorimeter, E - pz: " << Epz_minus << std::endl;
    std::cout << "Calorimeter, Sum HF+: " << sumEHF_plus << " [GeV]" << std::endl;
    std::cout << "Calorimeter, Sum HF-: " << sumEHF_minus << " [GeV]"<< std::endl;
  }
}

// ------------ event tagger ------------
void MissingMassSearches::fetchEventTagger(const edm::Event& iEvent){

  // Force pT sorting... (CMSSW does the job)
  std::sort(muonsVec.begin(), muonsVec.end(), orderPt());
  std::sort(electronsVec.begin(), electronsVec.end(), orderPt());
  std::sort(jetsVecA.begin(), jetsVecA.end(), orderPt());
  std::sort(jetsVecB.begin(), jetsVecB.end(), orderPt());
  std::sort(photonsVec.begin(), photonsVec.end(), orderPt());

  if(enableMC_){

    std::sort(genprotonsVec.begin(), genprotonsVec.end(), orderAbsolutPz());
    std::sort(genpartVec.begin(), genpartVec.end(), orderPt());
    std::sort(genjetsVecA.begin(), genjetsVecA.end(), orderPt());
    std::sort(genjetsVecB.begin(), genjetsVecB.end(), orderPt());

    for (const auto& genpartEvt: genpartVec) {
      if(genpartEvt->pdgId() == 22 && includePhotons_){
	(*genphotons_pt_).push_back(genpartEvt->pt());
	(*genphotons_eta_).push_back(genpartEvt->eta());
	(*genphotons_phi_).push_back(genpartEvt->phi());
	(*genphotons_energy_).push_back(genpartEvt->energy());
      }
      if(fabs(genpartEvt->pdgId()) == 11 && includeElectrons_){
	(*genleptons_pdgid_).push_back(11);
	(*genleptons_status_).push_back(genpartEvt->status());
	(*genleptons_energy_).push_back(genpartEvt->energy());
	(*genleptons_pt_).push_back(genpartEvt->pt());
	(*genleptons_eta_).push_back(genpartEvt->eta());
	(*genleptons_phi_).push_back(genpartEvt->phi());
	(*genleptons_px_).push_back(genpartEvt->px());
	(*genleptons_py_).push_back(genpartEvt->py());
	(*genleptons_pz_).push_back(genpartEvt->pz());
	(*genleptons_charge_).push_back(genpartEvt->charge());
	(*genleptons_vx_).push_back(genpartEvt->vertex().x());
	(*genleptons_vy_).push_back(genpartEvt->vertex().y());
	(*genleptons_vz_).push_back(genpartEvt->vertex().z());
      }
      if(fabs(genpartEvt->pdgId()) == 13 && includeMuons_){
	(*genleptons_pdgid_).push_back(13);
	(*genleptons_status_).push_back(genpartEvt->status());
	(*genleptons_energy_).push_back(genpartEvt->energy());
	(*genleptons_pt_).push_back(genpartEvt->pt());
	(*genleptons_eta_).push_back(genpartEvt->eta());
	(*genleptons_phi_).push_back(genpartEvt->phi());
	(*genleptons_px_).push_back(genpartEvt->px());
	(*genleptons_py_).push_back(genpartEvt->py());
	(*genleptons_pz_).push_back(genpartEvt->pz());
	(*genleptons_charge_).push_back(genpartEvt->charge());
	(*genleptons_vx_).push_back(genpartEvt->vertex().x());
	(*genleptons_vy_).push_back(genpartEvt->vertex().y());
	(*genleptons_vz_).push_back(genpartEvt->vertex().z());
      } 
    }

    if(includeProtonsReco_){
      for (const auto& genprotonsEvt: genprotonsVec){
	(*genprotons_status_).push_back(genprotonsEvt->status());
	(*genprotons_energy_).push_back(genprotonsEvt->energy());
	(*genprotons_pt_).push_back(genprotonsEvt->pt());
	(*genprotons_eta_).push_back(genprotonsEvt->eta());
	(*genprotons_phi_).push_back(genprotonsEvt->phi());
	(*genprotons_px_).push_back(genprotonsEvt->px());
	(*genprotons_py_).push_back(genprotonsEvt->py());
	(*genprotons_pz_).push_back(genprotonsEvt->pz());
	(*genprotons_xi_).push_back(1 - abs(genprotonsEvt->pz()) / 6500.);
      }
    }

    if (includeJets_){
      for ( const auto& genjetsak4Evt: genjetsVecA) {
	(*genjetsak4_px_).push_back(genjetsak4Evt->px());
	(*genjetsak4_py_).push_back(genjetsak4Evt->py());
	(*genjetsak4_pz_).push_back(genjetsak4Evt->pz());
	(*genjetsak4_pt_).push_back(genjetsak4Evt->pt());
	(*genjetsak4_energy_).push_back(genjetsak4Evt->energy());
	(*genjetsak4_phi_).push_back(genjetsak4Evt->phi());
	(*genjetsak4_eta_).push_back(genjetsak4Evt->eta());
	(*genjetsak4_vz_).push_back(genjetsak4Evt->vz());
      }
      for ( const auto& genjetsak8Evt: genjetsVecB) {
	(*genjetsak8_px_).push_back(genjetsak8Evt->px());
	(*genjetsak8_py_).push_back(genjetsak8Evt->py());
	(*genjetsak8_pz_).push_back(genjetsak8Evt->pz());
	(*genjetsak8_pt_).push_back(genjetsak8Evt->pt());
	(*genjetsak8_energy_).push_back(genjetsak8Evt->energy());
	(*genjetsak8_phi_).push_back(genjetsak8Evt->phi());
	(*genjetsak8_eta_).push_back(genjetsak8Evt->eta());
	(*genjetsak8_vz_).push_back(genjetsak8Evt->vz());
      }
    } 
  }

  if(includePhotons_){
    for ( const auto& photonsEvt: photonsVec) {
      (*photons_pt_).push_back(photonsEvt->pt());
      (*photons_eta_).push_back(photonsEvt->eta());
      (*photons_phi_).push_back(photonsEvt->phi());
      (*photons_energy_).push_back(photonsEvt->energy());
      (*photons_r9_).push_back(photonsEvt->r9());
      (*photons_looseId_).push_back(photonsEvt->photonID("cutBasedPhotonID-Fall17-94X-V1-loose"));
      (*photons_mediumId_).push_back(photonsEvt->photonID("cutBasedPhotonID-Fall17-94X-V1-medium"));
      (*photons_tightId_).push_back(photonsEvt->photonID("cutBasedPhotonID-Fall17-94X-V1-tight"));
    }
  }

  if(includeElectrons_){
    for ( const auto& leptonEvt: electronsVec) {
      (*leptons_pdgid_).push_back(11);
      (*leptons_energy_).push_back(leptonEvt->energy());
      (*leptons_pt_).push_back(leptonEvt->pt());
      (*leptons_eta_).push_back(leptonEvt->eta());
      (*leptons_phi_).push_back(leptonEvt->phi());
      (*leptons_px_).push_back(leptonEvt->px());
      (*leptons_py_).push_back(leptonEvt->py());
      (*leptons_pz_).push_back(leptonEvt->pz());
      (*leptons_charge_).push_back(leptonEvt->charge());
      (*leptons_vx_).push_back(leptonEvt->vertex().x());
      (*leptons_vy_).push_back(leptonEvt->vertex().y());
      (*leptons_vz_).push_back(leptonEvt->vertex().z());
      (*leptons_looseId_).push_back(leptonEvt->electronID("cutBasedElectronID-Fall17-94X-V1-loose"));
      (*leptons_mediumId_).push_back(leptonEvt->electronID("cutBasedElectronID-Fall17-94X-V1-medium"));
      (*leptons_tightId_).push_back(leptonEvt->electronID("cutBasedElectronID-Fall17-94X-V1-tight"));
      (*leptons_pfIsoMedium_).push_back(0);
      (*leptons_miniIsoTight_).push_back(0);
      (*leptons_pfIsoVeryTight_).push_back(0);
      (*leptons_pfIso_).push_back(0);
      (*leptons_tkIso_).push_back(0);
    }
    try{

      if(unmatching_){
	if(debug_) std::cout << "Running unmatching..." << std::endl; 
	// storing jets AK4 candidates which are not matching each lepton candidate (highest and second highest pt)
	jetsak4Cand = Unmatching(electronsVec, jetsVecA, 0.5);
	// storing jets AK8 candidates which are not matching each lepton candidate (highest and second highest pt)
	jetsak8Cand = Unmatching(electronsVec, jetsVecB, 0.9);
      }else{
	if(debug_) std::cout << "_NOT_ running unmatching..." << std::endl; 
	jetsak4Cand = jetsVecA;
	jetsak8Cand = jetsVecB;
      }


    }catch(...){}
  }

  if(includeMuons_){
    for ( const auto& leptonEvt: muonsVec) {
      (*leptons_pdgid_).push_back(13);
      (*leptons_energy_).push_back(leptonEvt->energy());
      (*leptons_pt_).push_back(leptonEvt->pt());
      (*leptons_eta_).push_back(leptonEvt->eta());
      (*leptons_phi_).push_back(leptonEvt->phi());
      (*leptons_px_).push_back(leptonEvt->px());
      (*leptons_py_).push_back(leptonEvt->py());
      (*leptons_pz_).push_back(leptonEvt->pz());
      (*leptons_charge_).push_back(leptonEvt->charge());
      (*leptons_vx_).push_back(leptonEvt->vertex().x());
      (*leptons_vy_).push_back(leptonEvt->vertex().y());
      (*leptons_vz_).push_back(leptonEvt->vertex().z());
      (*leptons_looseId_).push_back(leptonEvt->passed(reco::Muon::CutBasedIdLoose));
      (*leptons_mediumId_).push_back(leptonEvt->passed(reco::Muon::CutBasedIdMedium));
      (*leptons_tightId_).push_back(leptonEvt->passed(reco::Muon::CutBasedIdTight));

      double pfIso = (leptonEvt->pfIsolationR04().sumChargedHadronPt + std::max(0., leptonEvt->pfIsolationR04().sumNeutralHadronEt + leptonEvt->pfIsolationR04().sumPhotonEt - 0.5*leptonEvt->pfIsolationR04().sumPUPt))/leptonEvt->pt();
      double tkIso = (leptonEvt->isolationR03().sumPt)/(leptonEvt->pt());

      (*leptons_pfIsoMedium_).push_back(leptonEvt->passed(reco::Muon::PFIsoMedium));
      (*leptons_miniIsoTight_).push_back(leptonEvt->passed(reco::Muon::MiniIsoTight));
      (*leptons_pfIsoVeryTight_).push_back(leptonEvt->passed(reco::Muon::PFIsoVeryTight));
      (*leptons_pfIso_).push_back(pfIso);
      (*leptons_tkIso_).push_back(tkIso);
    }
    try{

      if(unmatching_){
	if(debug_) std::cout << "Running unmatching..." << std::endl; 
	// storing jets AK4 candidates which are not matching each lepton candidate (highest and second highest pt)
	jetsak4Cand = Unmatching(muonsVec, jetsVecA, 0.5);
	// storing jets AK8 candidates which are not matching each lepton candidate (highest and second highest pt)
	jetsak8Cand = Unmatching(muonsVec, jetsVecB, 0.9);
      }else{
	if(debug_) std::cout << "_NOT_ running unmatching..." << std::endl; 
	jetsak4Cand = jetsVecA;
	jetsak8Cand = jetsVecB;
      }

    }catch(...){}
  }

  // Force pT sorting again...
  std::sort(jetsak4Cand.begin(), jetsak4Cand.end(), orderPt());
  std::sort(jetsak8Cand.begin(), jetsak8Cand.end(), orderPt());

  // ak4, ak8 input, flag for slimmed Jets.
  fillingJets(jetsak4Cand, jetsak8Cand);

}

// ------------ debugging  ------------
void MissingMassSearches::Debug(){

  bool debug_deep = false;
  bool debug_protons = false;

  for ( const auto& triggerEvt: triggerVec ) {
    std::cout << "Trigger: " << triggerEvt << std::endl;
  }

  if(enableMC_){

    for ( const auto& genpartEvt: genpartVec) {
      std::cout << "GenParticle pT: " << genpartEvt->pt() << " [GeV], eta: " << genpartEvt->eta() << ", phi: "<< genpartEvt->phi() << ", ID: " << genpartEvt->pdgId() << ", Status: "<< genpartEvt->status() <<std::endl;
    }

    for ( const auto& genprotonEvt: genprotonsVec) {
      std::cout << "GenProtons pZ: " << genprotonEvt->pz() << "[GeV], pT: "<< genprotonEvt->pt() << " [GeV], eta: " << genprotonEvt->eta() << ", phi: "<< genprotonEvt->phi() << ", ID: " << genprotonEvt->pdgId() << ", Status: "<< genprotonEvt->status() <<std::endl;
    }

    std::cout << "GenMET eT: " << genmet->et() << " [GeV], phi: " << genmet->phi() << ", significance: " << genmet->significance() << std::endl;

    for ( const auto& genjetsak4Evt: genjetsVecA) {
      std::cout << "GenJets (AK4) pT: " << genjetsak4Evt->pt() << " [GeV], eta: " << genjetsak4Evt->eta() << ", phi: "<< genjetsak4Evt->phi() << std::endl;
    }

    for ( const auto& genjetsak8Evt: genjetsVecB) {
      std::cout << "GenJets (AK8) pT: " << genjetsak8Evt->pt() << " [GeV], eta: " << genjetsak8Evt->eta() << ", phi: "<< genjetsak8Evt->phi() << std::endl;
    }
  }

  for ( const auto& jetsak4Evt: jetsVecA) {
    std::cout << "Jets (AK4) pT: " << jetsak4Evt->pt() << " [GeV], eta: " << jetsak4Evt->eta() << ", phi: "<< jetsak4Evt->phi() << std::endl;
  }

  for ( const auto& jetsak8Evt: jetsVecB) {
    std::cout << "Jets (AK8) pT: " << jetsak8Evt->pt() << " [GeV], eta: " << jetsak8Evt->eta() << ", phi: "<< jetsak8Evt->phi() << std::endl;
  }

  for ( const auto& electronsEvt: electronsVec ) {
    std::cout << "Electrons pT: " << electronsEvt->pt() << " [GeV], eta: " << electronsEvt->eta() << ", phi: "<< electronsEvt->phi() << std::endl;
  }

  for ( const auto& muonsEvt: muonsVec ) {
    std::cout << "Muons pT: " << muonsEvt->pt() << " [GeV], eta: " << muonsEvt->eta() << ", phi: "<< muonsEvt->phi() << std::endl;
  }

  for ( const auto& photonsEvt: photonsVec ) {
    std::cout << "Photons pT: " << photonsEvt->pt() << " [GeV], eta: " << photonsEvt->eta() << ", phi: "<< photonsEvt->phi() << std::endl;
  }

  std::cout << "MET eT: " << met->et() << " [GeV], phi: " << met->phi() << ", significance: " << met->significance() << std::endl; 

  if(debug_deep){
    for ( const auto& pfEvt: pfVec ) {
      std::cout << "PF pT: " << pfEvt->pt() << " [GeV], eta: " << pfEvt->eta() << ", phi: "<< pfEvt->phi() << std::endl;
    }

    for ( const auto& packedEvt: packedVec ) {
      std::cout << "Packed PF pT: " << packedEvt->pt() << " [GeV], eta: " << packedEvt->eta() << ", phi: "<< packedEvt->phi() << ", vertex (x,y,z) [cm]: " << packedEvt->vertex() << std::endl;
    }

    for ( const auto& vtxEvt: vtxVec ) {
      std::cout << "Vertices position (x,y,z) [cm]: " << vtxEvt->x() << ", " << vtxEvt->y() << ", " << vtxEvt->z() << std::endl;
    }
  } //debug deep

  if(debug_protons){
    for ( const auto& protonEvt: protonSingleVec ) {
      const CTPPSDetId det_id((*protonEvt->contributingLocalTracks().begin())->getRPId());
      std::cout << "Proton Reco SingleRP: xi [a.u]: " << protonEvt->xi() << ", time [ns]: " << protonEvt->time() << ", Theta_y: " <<  protonEvt->thetaY() << " Theta_x: " <<  protonEvt->thetaX() << ", arm: " << det_id.arm() << ", station: " << det_id.station() << ", pot: " << det_id.rp() << std::endl;
    }

    for ( const auto& protonEvt: protonMultiVec ) {
      const CTPPSDetId det_id((*protonEvt->contributingLocalTracks().begin())->getRPId());
      std::cout << "Proton Reco MultiRP: xi [a.u]: " << protonEvt->xi() << ", time [ns]: " << protonEvt->time() << ", Theta_y: " <<  protonEvt->thetaY() << " Theta_x: "<<  protonEvt->thetaX() << ", arm: " << det_id.arm() << ", station: " << det_id.station() << ", pot: " << det_id.rp() << std::endl;
    }

    for ( const auto& ppsEvt: ppsVec ) {
      const CTPPSDetId det_id( ppsEvt->getRPId() );
      std::cout << "PPS Hit (x,y)[mm]: (" << ppsEvt->getX() << ", " << ppsEvt->getY() << "), arm: " << det_id.arm() << ", station: " << det_id.station() << ", pot: " << det_id.rp() << std::endl;
    }

    for ( const auto& pixelsEvt: pixelsVec ) {
      std::cout << "Pixels (x,y)[mm]: (" << pixelsEvt.first->getX0() << ", " << pixelsEvt.first->getY0() << "), Arm: " << pixelsEvt.second.arm() << ", Pot: " << pixelsEvt.second.rp() << ", Station: " << pixelsEvt.second.station() << std::endl;
    }

    for ( const auto& timingEvt: timingVec ) {
      std::cout << "Timing (x,y)[mm]: (" << timingEvt.first->getX0() << ", " << timingEvt.first->getY0() << "), Arm: " << timingEvt.second.arm() << ", Pot: " << timingEvt.second.rp() << ", Station: " << timingEvt.second.station() << std::endl;
    }

    for ( const auto& stripsEvt: stripsVec ) {
      std::cout << "Strips (x,y)[mm]: (" << stripsEvt.first->getX0() << ", " << stripsEvt.first->getY0() << "), Arm: " << stripsEvt.second.arm() << ", Pot: " << stripsEvt.second.rp() << ", Station: " << stripsEvt.second.station() << std::endl;
    }

  } //debug protons

}

// ------------ method called for each event  ------------
  void
MissingMassSearches::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // Cleaning vectors per event
  cleaning();

  // Trigger fired!
  if(fetchTrigger(iEvent,iSetup)){

    for (std::size_t i = 0; i != triggerVec.size(); ++i) {
      if(enableTrigger_){
	hltTriggerNamesHisto_->Fill( triggersList_[i].c_str(),triggerVec[i]);
	(*trigger_).push_back(triggerVec[i]);
      }else{
	(*trigger_).push_back(-1);
      }
      if(enablePrescales_){
	(*prescalesL1_).push_back(prescalesL1Vec[i]);
	(*prescalesHLT_).push_back(prescalesHLTVec[i]);
      }else{
	(*prescalesL1_).push_back(-1);
	(*prescalesHLT_).push_back(-1);

      }
    }

    *run_ = iEvent.id().run();
    *ev_ = iEvent.id().event();
    *lumiblock_ = iEvent.luminosityBlock();

    fetchLHCInfo(iSetup);

    // Filling Vectors
    if(enableMC_) {
      fetchGenParticles(iEvent, iSetup);
      fetchGenMET(iEvent, iSetup);
      fetchGenJets(iEvent, iSetup);
      edm::Handle<std::vector<PileupSummaryInfo>> PupInfo;
      iEvent.getByToken(puToken_,PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator ipu;
      for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu){
	if (ipu->getBunchCrossing() != 0) continue; // storing detailed PU info only for BX=0
	*PUinter_ = ipu->getPU_NumInteractions();
	*PUtrueinter_ =ipu->getTrueNumInteractions();
      }

    }

    fetchMuons(iEvent, iSetup);
    fetchElectrons(iEvent, iSetup);
    if(includeJets_) fetchJets(iEvent, iSetup);
    if(includePhotons_) fetchPhotons(iEvent, iSetup);
    if(includeMET_) fetchMET(iEvent, iSetup);
    if(includeVertices_) fetchVertices(iEvent, iSetup);
    if(includePF_) fetchPF(iEvent, iSetup);
    if(includeProtonsReco_) fetchProtonsReco(iEvent, iSetup); 
    if(includePPSInfo_) fetchPPSInfo(iEvent, iSetup); 
    fetchEventTagger(iEvent);
    tree_->Fill();

    // Printout for Debugging
    if(debug_) Debug();

  } //end trigger

  eventClear();

}


// ------------ method called when a run finishes ------------
void MissingMassSearches::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}

// ------------ method called when a new run number begins  ------------
void MissingMassSearches::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){

  bool debug = false;

  using namespace std;
  using namespace edm;

  bool changed(true);

  if(enableTrigger_ && enablePrescales_){
    if (hltPrescaleProvider_.init(iRun,iSetup,"HLT", changed)) {
      hltPrescaleProvider_.hltConfigProvider();
      if(debug){
	HLTConfigProvider const&  hltConfig = hltPrescaleProvider_.hltConfigProvider();
	if (changed) {
	  for(vector<string>::const_iterator triggerName_ = triggersList_.begin(); triggerName_ != triggersList_.end(); ++triggerName_) {
	    const unsigned int n(hltConfig.size());
	    const unsigned int triggerIndex(hltConfig.triggerIndex(*triggerName_));
	    if (triggerIndex>=n) {
	      hltConfig.dump("Triggers");
	    }
	  }
	  hltConfig.dump("ProcessName");
	  hltConfig.dump("GlobalTag");
	  hltConfig.dump("TableName");
	  hltConfig.dump("Streams");
	  hltConfig.dump("Datasets");
	  hltConfig.dump("PrescaleTable");
	  hltConfig.dump("ProcessPSet");
	}
      }
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void MissingMassSearches::beginJob(){

  edm::Service<TFileService> fs;
  tree_=fs->make<TTree>("analyzer","analyzer");

  // Control Histograms
  if(enableTrigger_){
    TFileDirectory triggerDir = fs->mkdir("TriggerInfo");
    hltTriggerNamesHisto_ = triggerDir.make<TH1F>("HLTTriggerNames","HLTTriggerNames",1,0,1);
    hltTriggerNamesHisto_->SetCanExtend(TH1::kXaxis);
  }

  run_ = new int;
  ev_ = new long int;
  lumiblock_ = new int;

  trigger_ = new std::vector<int>;
  prescalesL1_ = new std::vector<int>;
  prescalesHLT_ = new std::vector<int>;

  xangle_ = new int;
  betastar_ = new double;
  bunchesb1_ = new int;
  bunchesb2_ = new int;
  intensityb1_ = new double;
  intensityb2_ = new double;
  instlumi_ = new double;

  PUinter_ = new int;
  PUtrueinter_ = new int;

  if(includeVertices_){
    vertex_x_ = new std::vector<double>;
    vertex_y_ = new std::vector<double>;
    vertex_z_ = new std::vector<double>;
    vertex_ntrack_ = new std::vector<int>;
    vertex_chi2_ = new std::vector<double>;
    vertex_ndof_ = new std::vector<double>;
  }

  if(enableMC_){

    nGenParticles_ = new int;

    if(includeElectrons_ || includeMuons_){
      nGenMuons_ = new int;
      nGenElectrons_ = new int;  
      genleptons_pdgid_ = new std::vector<int>;
      genleptons_status_ = new std::vector<int>;
      genleptons_energy_ = new std::vector<double>;
      genleptons_pt_ = new std::vector<double>;
      genleptons_eta_ = new std::vector<double>;
      genleptons_phi_ = new std::vector<double>;
      genleptons_px_ = new std::vector<double>;
      genleptons_py_ = new std::vector<double>;
      genleptons_pz_ = new std::vector<double>;
      genleptons_charge_ = new std::vector<double>;
      genleptons_vx_ = new std::vector<double>;
      genleptons_vy_ = new std::vector<double>;
      genleptons_vz_ = new std::vector<double>;
    }
    if(includePhotons_){
      genphotons_pt_ = new std::vector<double>;
      genphotons_eta_ = new std::vector<double>;
      genphotons_phi_ = new std::vector<double>;
      genphotons_energy_ = new std::vector<double>;
    }
    if(includeProtonsReco_){
      nGenProtons_ = new int;
      genprotons_status_ = new std::vector<int>;
      genprotons_energy_ = new std::vector<double>;
      genprotons_pt_ = new std::vector<double>;
      genprotons_eta_ = new std::vector<double>;
      genprotons_phi_ = new std::vector<double>;
      genprotons_px_ = new std::vector<double>;
      genprotons_py_ = new std::vector<double>;
      genprotons_pz_ = new std::vector<double>;
      genprotons_xi_ = new std::vector<double>;
    }
    if(includeMET_){
      genmisset_ = new double;
      genmisset_phi_ = new double;
    }   
    if(includeJets_){
      genjetsak4_px_ = new std::vector<double>;
      genjetsak4_py_ = new std::vector<double>;
      genjetsak4_pz_ = new std::vector<double>;
      genjetsak4_pt_ = new std::vector<double>;
      genjetsak4_energy_ = new std::vector<double>;
      genjetsak4_phi_ = new std::vector<double>;
      genjetsak4_eta_ = new std::vector<double>;
      genjetsak4_vz_ = new std::vector<double>;
      genjetsak8_px_ = new std::vector<double>;
      genjetsak8_py_ = new std::vector<double>;
      genjetsak8_pz_ = new std::vector<double>;
      genjetsak8_pt_ = new std::vector<double>;
      genjetsak8_energy_ = new std::vector<double>;
      genjetsak8_phi_ = new std::vector<double>;
      genjetsak8_eta_ = new std::vector<double>;
      genjetsak8_vz_ = new std::vector<double>;
    }
  }

  if(includeJets_){
    jetsak4_px_ = new std::vector<double>;
    jetsak4_py_ = new std::vector<double>;
    jetsak4_pz_ = new std::vector<double>;
    jetsak4_pt_ = new std::vector<double>;
    jetsak4_energy_ = new std::vector<double>;
    jetsak4_phi_ = new std::vector<double>;
    jetsak4_eta_ = new std::vector<double>;
    jetsak4_vx_ = new std::vector<double>;
    jetsak4_vy_ = new std::vector<double>;
    jetsak4_vz_ = new std::vector<double>;
    jetsak4_bdis_ = new std::vector<double>;
    jetsak4_qgdis_ = new std::vector<double>;
    jetsak4_neutralEmFraction_ = new std::vector<double>;
    jetsak4_neutralHadFraction_ = new std::vector<double>;
    jetsak4_chargedEmFraction_ = new std::vector<double>;
    jetsak4_chargedHadFraction_ = new std::vector<double>;
    jetsak4_muonFraction_ = new std::vector<double>;
    jetsak4_neutralMultiplicity_ = new std::vector<int>;
    jetsak4_chargedMultiplicity_ = new std::vector<int>;
    jetsak4_looseJetId_ = new std::vector<bool>;
    jetsak4_tightJetId_ = new std::vector<bool>;
    jetsak4_lepVetoJetId_ = new std::vector<bool>;
    jetsak4_puIdfdisc_ = new std::vector<double>;
    jetsak4_puIdcbased_ = new std::vector<int>;
    jetsak4_puIdfid_ = new std::vector<int>;
    jetsak8_px_ = new std::vector<double>;
    jetsak8_py_ = new std::vector<double>;
    jetsak8_pz_ = new std::vector<double>;
    jetsak8_pt_ = new std::vector<double>;
    jetsak8_energy_ = new std::vector<double>;
    jetsak8_phi_ = new std::vector<double>;
    jetsak8_eta_ = new std::vector<double>;
    jetsak8_vx_ = new std::vector<double>;
    jetsak8_vy_ = new std::vector<double>;
    jetsak8_vz_ = new std::vector<double>;
    jetsak8_bdis_ = new std::vector<double>;
    jetsak8_neutralEmFraction_ = new std::vector<double>;
    jetsak8_neutralHadFraction_ = new std::vector<double>;
    jetsak8_chargedEmFraction_ = new std::vector<double>;
    jetsak8_chargedHadFraction_ = new std::vector<double>;
    jetsak8_muonFraction_ = new std::vector<double>;
    jetsak8_neutralMultiplicity_ = new std::vector<int>;
    jetsak8_chargedMultiplicity_ = new std::vector<int>;
    jetsak8_looseJetId_ = new std::vector<bool>;
    jetsak8_tightJetId_ = new std::vector<bool>;
    jetsak8_lepVetoJetId_ = new std::vector<bool>;
  }

  if(includeElectrons_ || includeMuons_){
    nMuons_ = new int;
    nMuons_looseId_ = new int;
    nMuons_mediumId_ = new int;
    nMuons_tightId_ = new int;
    nMuons_pfIsoMedium_ = new int;
    nMuons_miniIsoTight_ = new int;
    nMuons_pfIsoVeryTight_ = new int;
    nElectrons_ = new int;
    nElectrons_looseId_ = new int;
    nElectrons_mediumId_ = new int;
    nElectrons_tightId_ = new int;
    leptons_pdgid_ = new std::vector<int>;
    leptons_energy_ = new std::vector<double>;
    leptons_pt_ = new std::vector<double>;
    leptons_eta_ = new std::vector<double>;
    leptons_phi_ = new std::vector<double>;
    leptons_px_ = new std::vector<double>;
    leptons_py_ = new std::vector<double>;
    leptons_pz_ = new std::vector<double>;
    leptons_charge_ = new std::vector<double>;
    leptons_vx_ = new std::vector<double>;
    leptons_vy_ = new std::vector<double>;
    leptons_vz_ = new std::vector<double>;
    leptons_looseId_ = new std::vector<bool>;
    leptons_mediumId_ = new std::vector<bool>;
    leptons_tightId_ = new std::vector<bool>;
    leptons_pfIsoMedium_ = new std::vector<bool>;
    leptons_miniIsoTight_ = new std::vector<bool>;
    leptons_pfIsoVeryTight_ = new std::vector<bool>;
    leptons_pfIso_ = new std::vector<double>;
    leptons_tkIso_ = new std::vector<double>;
  }

  if(includePhotons_){
    photons_pt_ = new std::vector<double>;
    photons_eta_ = new std::vector<double>;
    photons_phi_ = new std::vector<double>;
    photons_energy_ = new std::vector<double>;
    photons_r9_ = new std::vector<double>;
    photons_looseId_ = new std::vector<bool>;
    photons_mediumId_ = new std::vector<bool>;
    photons_tightId_ = new std::vector<bool>;
  }

  if(includeMET_){
    misset_ = new double;
    misset_phi_ = new double;
  }

  if(includePF_){
    nPF_ = new int;
    SumPF_energy_ = new double;
    pfHFHadPlus_energy_ = new double;
    pfHFHadPlus_pt_ = new double;
    pfHFHadMinus_energy_ = new double;
    pfHFHadMinus_pt_ = new double;
    pfHFEmPlus_energy_ = new double;
    pfHFEmPlus_pt_ = new double;
    pfHFEmMinus_energy_ = new double;
    pfHFEmMinus_pt_ = new double;
    pfEcalMinus_energy_ = new double;
    pfEcalPlus_energy_ = new double;
    pfHcalMinus_energy_ = new double;
    pfHcalPlus_energy_ = new double;
    pfHadNeutralMinus_energy_ = new double;
    pfHadNeutralPlus_energy_ = new double;
    SumChargedPFMultiPV_pt_Loose_ = new double;
    SumChargedPFMultiPV_pt_Tight_ = new double;
    SumChargedPFMultiPV_pt_UsedInFit_ = new double;
    SumChargedPFMultiPV_pt_Tight_Fit_ = new double;
    nChargedPFMultiPV_Loose_ = new int;
    nChargedPFMultiPV_Tight_ = new int;
    nChargedPFMultiPV_UsedInFit_ = new int;
    nChargedPFMultiPV_Tight_Fit_ = new int;
    nPFPlus_ = new int;
    nPFMinus_ = new int;
    nPFHFHadPlus_ = new int;
    nPFHFHadMinus_ = new int;
    nPFHFEmPlus_ = new int;
    nPFHFEmMinus_ = new int;
    pfSlice0P_energy_ = new double;
    pfSlice1P_energy_ = new double;
    pfSlice2P_energy_ = new double;
    pfSlice3P_energy_ = new double;
    pfSlice4P_energy_ = new double;
    pfSlice0N_energy_ = new double;
    pfSlice1N_energy_ = new double;
    pfSlice2N_energy_ = new double;
    pfSlice3N_energy_ = new double;
    pfSlice4N_energy_ = new double;
    pfSlice0P_pt_ = new double;
    pfSlice1P_pt_ = new double;
    pfSlice2P_pt_ = new double;
    pfSlice3P_pt_ = new double;
    pfSlice4P_pt_ = new double;
    pfSlice0N_pt_ = new double;
    pfSlice1N_pt_ = new double;
    pfSlice2N_pt_ = new double;
    pfSlice3N_pt_ = new double;
    pfSlice4N_pt_ = new double;
    npfSlice0P_ = new int;
    npfSlice1P_ = new int;
    npfSlice2P_ = new int;
    npfSlice3P_ = new int;
    npfSlice4P_ = new int;
    npfSlice0N_ = new int;
    npfSlice1N_ = new int;
    npfSlice2N_ = new int;
    npfSlice3N_ = new int;
    npfSlice4N_ = new int;
    nmuonSlice0P_ = new int;
    nmuonSlice1P_ = new int;
    nmuonSlice2P_ = new int;
    nmuonSlice3P_ = new int;
    nmuonSlice4P_ = new int;
    nmuonSlice0N_ = new int;
    nmuonSlice1N_ = new int;
    nmuonSlice2N_ = new int;
    nmuonSlice3N_ = new int;
    nmuonSlice4N_ = new int;
    nelectronSlice0P_ = new int;
    nelectronSlice1P_ = new int;
    nelectronSlice2P_ = new int;
    nelectronSlice3P_ = new int;
    nelectronSlice4P_ = new int;
    nelectronSlice0N_ = new int;
    nelectronSlice1N_ = new int;
    nelectronSlice2N_ = new int;
    nelectronSlice3N_ = new int;
    nelectronSlice4N_ = new int;
  }

  if(includePPSInfo_){
    protonsArm_ = new std::vector<int>;
    protonsStation_ = new std::vector<int>;
    protonsRP_ = new std::vector<int>;
    protonsX_ = new std::vector<double>;
    protonsXUnc_ = new std::vector<double>;
    protonsY_ = new std::vector<double>;
    protonsYUnc_ = new std::vector<double>;
    if(!includeProtonsReco_){
      protonsIsValid_ = new std::vector<bool>;
      protonsTx_ = new std::vector<double>;
      protonsTxUnc_ = new std::vector<double>;
      protonsTy_ = new std::vector<double>;
      protonsTyUnc_ = new std::vector<double>;
      protonsTime_ = new std::vector<double>;
      protonsTimeUnc_ = new std::vector<double>;
    }
  }

  if(includeProtonsReco_){
    singleProtonArm_ = new std::vector<int>;
    singleProtonStation_ = new std::vector<int>;
    singleProtonPot_ = new std::vector<int>;
    singleProtonXi_ = new std::vector<double>;
    singleProtonT_ = new std::vector<double>;
    singleProtonThetaX_ = new std::vector<double>;
    singleProtonThetaY_ = new std::vector<double>;
    multiProtonArm_ = new std::vector<int>;
    multiProtonXi_ = new std::vector<double>;
    multiProtonT_ = new std::vector<double>;
    multiProtonTime_ = new std::vector<double>;
    multiProtonTimeError_ = new std::vector<double>;
    multiProtonVz_ = new std::vector<double>;
    multiProtonThetaX_ = new std::vector<double>;
    multiProtonThetaY_ = new std::vector<double>;
  }

  tree_->Branch("run",run_,"run/I");
  tree_->Branch("event",ev_,"event/L");
  tree_->Branch("lumiblock",lumiblock_,"lumiblock/I");
  tree_->Branch("trigger",&trigger_);
  tree_->Branch("prescalesL1",&prescalesL1_);
  tree_->Branch("prescalesHLT",&prescalesHLT_);
  tree_->Branch("xangle",xangle_,"xangle/I");
  tree_->Branch("betastar",betastar_,"betastar/D");
  tree_->Branch("bunchesb1",bunchesb1_,"bunchesb1/I");
  tree_->Branch("bunchesb2",bunchesb2_,"bunchesb2/I");
  tree_->Branch("intensityb1",intensityb1_,"intensityb1/D");
  tree_->Branch("intensityb2",intensityb2_,"intensityb2/D");
  tree_->Branch("instlumi",instlumi_,"instlumi/D");

  if(includeVertices_){
    tree_->Branch("vertex_x",&vertex_x_);
    tree_->Branch("vertex_y",&vertex_y_);
    tree_->Branch("vertex_z",&vertex_z_);
    tree_->Branch("vertex_ntrack",&vertex_ntrack_);
    tree_->Branch("vertex_chi2",&vertex_chi2_);
    tree_->Branch("vertex_ndof",&vertex_ndof_);
  }

  if(enableMC_){ 
    tree_->Branch("PUInterac",PUinter_,"PUInterac/I");
    tree_->Branch("PUTrueInterac",PUtrueinter_,"PUTrueInterac/I");
    tree_->Branch("nGenParticles",nGenParticles_,"nGenParticles/I");
    if(includeElectrons_ || includeMuons_){
      tree_->Branch("nGenMuons",nGenMuons_,"nGenMuons/I");
      tree_->Branch("nGenElectrons",nGenElectrons_,"nGenElectrons/I");
      tree_->Branch("genleptons_pdgid",&genleptons_pdgid_);
      tree_->Branch("genleptons_status",&genleptons_status_);
      tree_->Branch("genleptons_energy",&genleptons_energy_);
      tree_->Branch("genleptons_pt",&genleptons_pt_);
      tree_->Branch("genleptons_eta",&genleptons_eta_);
      tree_->Branch("genleptons_phi",&genleptons_phi_);
      tree_->Branch("genleptons_px",&genleptons_px_);
      tree_->Branch("genleptons_py",&genleptons_py_);
      tree_->Branch("genleptons_pz",&genleptons_pz_);
      tree_->Branch("genleptons_charge",&genleptons_charge_);
      tree_->Branch("genleptons_vx",&genleptons_vx_);
      tree_->Branch("genleptons_vy",&genleptons_vy_);
      tree_->Branch("genleptons_vz",&genleptons_vz_);
    }
    if(includePhotons_){
      tree_->Branch("genphotons_pt",&genphotons_pt_);
      tree_->Branch("genphotons_eta",&genphotons_eta_);
      tree_->Branch("genphotons_phi",&genphotons_phi_);
      tree_->Branch("genphotons_energy",&genphotons_energy_);      
    }
    if(includeProtonsReco_){
      tree_->Branch("nGenProtons",nGenProtons_,"nGenProtons/I");
      tree_->Branch("genprotons_status",&genprotons_status_);
      tree_->Branch("genprotons_energy",&genprotons_energy_);
      tree_->Branch("genprotons_pt",&genprotons_pt_);
      tree_->Branch("genprotons_eta",&genprotons_eta_);
      tree_->Branch("genprotons_phi",&genprotons_phi_);
      tree_->Branch("genprotons_px",&genprotons_px_);
      tree_->Branch("genprotons_py",&genprotons_py_);
      tree_->Branch("genprotons_pz",&genprotons_pz_);
      tree_->Branch("genprotons_xi",&genprotons_xi_);      
    }
    if(includeMET_){
      tree_->Branch("genMissEt",genmisset_,"genMissEt/D");
      tree_->Branch("genMissEt_phi",genmisset_phi_,"genMissEt_phi/D");
    }
    if(includeJets_){
      tree_->Branch("genjetsak4_px",&genjetsak4_px_);
      tree_->Branch("genjetsak4_py",&genjetsak4_py_);
      tree_->Branch("genjetsak4_pz",&genjetsak4_pz_);
      tree_->Branch("genjetsak4_pt",&genjetsak4_pt_);
      tree_->Branch("genjetsak4_energy",&genjetsak4_energy_);
      tree_->Branch("genjetsak4_phi",&genjetsak4_phi_);
      tree_->Branch("genjetsak4_eta",&genjetsak4_eta_);
      tree_->Branch("genjetsak4_vz",&genjetsak4_vz_);
      tree_->Branch("genjetsak8_px",&genjetsak8_px_);
      tree_->Branch("genjetsak8_py",&genjetsak8_py_);
      tree_->Branch("genjetsak8_pz",&genjetsak8_pz_);
      tree_->Branch("genjetsak8_pt",&genjetsak8_pt_);
      tree_->Branch("genjetsak8_energy",&genjetsak8_energy_);
      tree_->Branch("genjetsak8_phi",&genjetsak8_phi_);
      tree_->Branch("genjetsak8_eta",&genjetsak8_eta_);
      tree_->Branch("genjetsak8_vz",&genjetsak8_vz_); 
    }
  }

  if(includeJets_){
    tree_->Branch("jetsak4_px",&jetsak4_px_);
    tree_->Branch("jetsak4_py",&jetsak4_py_);
    tree_->Branch("jetsak4_pz",&jetsak4_pz_);
    tree_->Branch("jetsak4_pt",&jetsak4_pt_);
    tree_->Branch("jetsak4_energy",&jetsak4_energy_);
    tree_->Branch("jetsak4_phi",&jetsak4_phi_);
    tree_->Branch("jetsak4_eta",&jetsak4_eta_);
    tree_->Branch("jetsak4_vx",&jetsak4_vx_);
    tree_->Branch("jetsak4_vy",&jetsak4_vy_);
    tree_->Branch("jetsak4_vz",&jetsak4_vz_);
    tree_->Branch("jetsak4_bdis",&jetsak4_bdis_);
    tree_->Branch("jetsak4_qgdis",&jetsak4_qgdis_);
    tree_->Branch("jetsak4_neutralEmFraction",&jetsak4_neutralEmFraction_);
    tree_->Branch("jetsak4_neutralHadFraction",&jetsak4_neutralHadFraction_);
    tree_->Branch("jetsak4_chargedEmFraction",&jetsak4_chargedEmFraction_);
    tree_->Branch("jetsak4_chargedHadFraction",&jetsak4_chargedHadFraction_);
    tree_->Branch("jetsak4_muonFraction",&jetsak4_muonFraction_);
    tree_->Branch("jetsak4_neutralMultiplicity",&jetsak4_neutralMultiplicity_);
    tree_->Branch("jetsak4_chargedMultiplicity",&jetsak4_chargedMultiplicity_);
    tree_->Branch("jetsak4_looseId",&jetsak4_looseJetId_);
    tree_->Branch("jetsak4_tightId",&jetsak4_tightJetId_);
    tree_->Branch("jetsak4_lepVeto",&jetsak4_lepVetoJetId_);
    tree_->Branch("jetsak4_puIdfdisc",&jetsak4_puIdfdisc_);
    tree_->Branch("jetsak4_puIdcbased",&jetsak4_puIdcbased_);
    tree_->Branch("jetsak4_puIdfid",&jetsak4_puIdfid_);
    tree_->Branch("jetsak8_px",&jetsak8_px_);
    tree_->Branch("jetsak8_py",&jetsak8_py_);
    tree_->Branch("jetsak8_pz",&jetsak8_pz_);
    tree_->Branch("jetsak8_pt",&jetsak8_pt_);
    tree_->Branch("jetsak8_energy",&jetsak8_energy_);
    tree_->Branch("jetsak8_phi",&jetsak8_phi_);
    tree_->Branch("jetsak8_eta",&jetsak8_eta_);
    tree_->Branch("jetsak8_vx",&jetsak8_vx_);
    tree_->Branch("jetsak8_vy",&jetsak8_vy_);
    tree_->Branch("jetsak8_vz",&jetsak8_vz_);
    tree_->Branch("jetsak8_bdis",&jetsak8_bdis_);
    tree_->Branch("jetsak8_neutralEmFraction",&jetsak8_neutralEmFraction_);
    tree_->Branch("jetsak8_neutralHadFraction",&jetsak8_neutralHadFraction_);
    tree_->Branch("jetsak8_chargedEmFraction",&jetsak8_chargedEmFraction_);
    tree_->Branch("jetsak8_chargedHadFraction",&jetsak8_chargedHadFraction_);
    tree_->Branch("jetsak8_muonFraction",&jetsak8_muonFraction_);
    tree_->Branch("jetsak8_neutralMultiplicity",&jetsak8_neutralMultiplicity_);
    tree_->Branch("jetsak8_chargedMultiplicity",&jetsak8_chargedMultiplicity_);
    tree_->Branch("jetsak8_looseId",&jetsak8_looseJetId_);
    tree_->Branch("jetsak8_tightId",&jetsak8_tightJetId_);
    tree_->Branch("jetsak8_lepVeto",&jetsak8_lepVetoJetId_);
  }

  if(includeElectrons_ || includeMuons_){
    tree_->Branch("nMuons",nMuons_,"nMuons/I");
    tree_->Branch("nMuons_looseId",nMuons_looseId_,"nMuons_looseId/I");
    tree_->Branch("nMuons_mediumId",nMuons_mediumId_,"nMuons_mediumId/I");
    tree_->Branch("nMuons_tightId",nMuons_tightId_,"nMuons_tightId/I");
    tree_->Branch("nMuons_pfIsoMedium",nMuons_pfIsoMedium_,"nMuons_pfIsoMedium/I");
    tree_->Branch("nMuons_miniIsoTight",nMuons_miniIsoTight_,"nMuons_miniIsoTight/I");
    tree_->Branch("nMuons_pfIsoVeryTight",nMuons_pfIsoVeryTight_,"nMuons_pfIsoVeryTight/I");
    tree_->Branch("nElectrons",nElectrons_,"nElectrons/I");
    tree_->Branch("nElectrons_looseId",nElectrons_looseId_,"nElectrons_looseId/I");
    tree_->Branch("nElectrons_mediumId",nElectrons_mediumId_,"nElectrons_mediumId/I");
    tree_->Branch("nElectrons_tightId",nElectrons_tightId_,"nElectrons_tightId/I");
    tree_->Branch("leptons_pdgid",&leptons_pdgid_);
    tree_->Branch("leptons_energy",&leptons_energy_);
    tree_->Branch("leptons_pt",&leptons_pt_);
    tree_->Branch("leptons_eta",&leptons_eta_);
    tree_->Branch("leptons_phi",&leptons_phi_);
    tree_->Branch("leptons_px",&leptons_px_);
    tree_->Branch("leptons_py",&leptons_py_);
    tree_->Branch("leptons_pz",&leptons_pz_);
    tree_->Branch("leptons_charge",&leptons_charge_);
    tree_->Branch("leptons_vx",&leptons_vx_);
    tree_->Branch("leptons_vy",&leptons_vy_);
    tree_->Branch("leptons_vz",&leptons_vz_);
    tree_->Branch("leptons_looseId",&leptons_looseId_);
    tree_->Branch("leptons_mediumId",&leptons_mediumId_);
    tree_->Branch("leptons_tightId",&leptons_tightId_);
    tree_->Branch("leptons_pfIsoMedium_",&leptons_pfIsoMedium_);
    tree_->Branch("leptons_miniIsoTight_",&leptons_miniIsoTight_);
    tree_->Branch("leptons_pfIsoVeryTight_",&leptons_pfIsoVeryTight_);
    tree_->Branch("leptons_pfIso_",&leptons_pfIso_);
    tree_->Branch("leptons_tkIso_",&leptons_tkIso_);
  }

  if(includePhotons_){
    tree_->Branch("photons_pt",&photons_pt_);
    tree_->Branch("photons_eta",&photons_eta_);
    tree_->Branch("photons_phi",&photons_phi_);
    tree_->Branch("photons_energy",&photons_energy_);
    tree_->Branch("photons_r9",&photons_r9_);
    tree_->Branch("photons_looseId",&photons_looseId_);
    tree_->Branch("photons_mediumId",&photons_mediumId_);
    tree_->Branch("photons_tightId",&photons_tightId_);
  }

  if(includeMET_){
    tree_->Branch("missEt",misset_,"missEt/D");
    tree_->Branch("missEt_phi",misset_phi_,"missEt_phi/D");
  }

  if(includePF_){
    tree_->Branch("nPF",nPF_,"nPF/I");
    tree_->Branch("sumEnergyPF",SumPF_energy_,"sumEnergyPF/D");
    tree_->Branch("pfHFHadEnergyPlus",pfHFHadPlus_energy_,"pfHFHadEnergyPlus/D");
    tree_->Branch("pfHFHadEnergyMinus",pfHFHadMinus_energy_,"pfHFHadEnergyMinus/D");
    tree_->Branch("pfHFHadPtPlus",pfHFHadPlus_pt_,"pfHFHadPtPlus/D");
    tree_->Branch("pfHFHadPtMinus",pfHFHadMinus_pt_,"pfHFHadPtMinus/D");
    tree_->Branch("pfHFEmEnergyPlus",pfHFEmPlus_energy_,"pfHFEmEnergyPlus/D");
    tree_->Branch("pfHFEmEnergyMinus",pfHFEmMinus_energy_,"pfHFEmEnergyMinus/D");
    tree_->Branch("pfHFEmPtPlus",pfHFEmPlus_pt_,"pfHFEmPtPlus/D");
    tree_->Branch("pfHFEmPtMinus",pfHFEmMinus_pt_,"pfHFEmPtMinus/D");
    tree_->Branch("pfEcalEnergyPlus",pfEcalPlus_energy_,"pfEcalEnergyPlus/D");
    tree_->Branch("pfEcalEnergyMinus",pfEcalMinus_energy_,"pfEcalEnergyMinus/D");
    tree_->Branch("pfHcalEnergyPlus",pfHcalPlus_energy_,"pfHcalEnergyPlus/D");
    tree_->Branch("pfHcalEnergyMinus",pfHcalMinus_energy_,"pfHcalEnergyMinus/D");
    tree_->Branch("pfHneutralEnergyPlus",pfHadNeutralPlus_energy_,"pfHneutralEnergyPlus/D");
    tree_->Branch("pfHneutralEnergyMinus",pfHadNeutralMinus_energy_,"pfHneutralEnergyMinus/D");
    tree_->Branch("nChargedPFMultiPV_Loose",nChargedPFMultiPV_Loose_,"nChargedPFMultiPV_Loose/I");
    tree_->Branch("nChargedPFMultiPV_Tight",nChargedPFMultiPV_Tight_,"nChargedPFMultiPV_Tight/I");
    tree_->Branch("nChargedPFMultiPV_UsedInFit",nChargedPFMultiPV_UsedInFit_,"nChargedPFMultiPV_UsedInFit/I");
    tree_->Branch("nChargedPFMultiPV_Tight_Fit",nChargedPFMultiPV_Tight_Fit_,"nChargedPFMultiPV_Tight_Fit/I");
    tree_->Branch("SumChargedPFMultiPV_pt_Loose",SumChargedPFMultiPV_pt_Loose_,"SumChargedPFMultiPV_pt_Loose/D");
    tree_->Branch("SumChargedPFMultiPV_pt_Tight",SumChargedPFMultiPV_pt_Tight_,"SumChargedPFMultiPV_pt_Tight/D");
    tree_->Branch("SumChargedPFMultiPV_pt_UsedInFit",SumChargedPFMultiPV_pt_UsedInFit_,"SumChargedPFMultiPV_pt_UsedInFit/D");
    tree_->Branch("SumChargedPFMultiPV_pt_Tight_Fit",SumChargedPFMultiPV_pt_Tight_Fit_,"SumChargedPFMultiPV_pt_Tight_Fit/D");
    tree_->Branch("nPFPlus",nPFPlus_,"nPFPlus/I");
    tree_->Branch("nPFMinus",nPFMinus_,"nPFMinus/I");
    tree_->Branch("nPFHFHadPlus",nPFHFHadPlus_,"nPFHFHadPlus/I");
    tree_->Branch("nPFHFHadMinus",nPFHFHadMinus_,"nPFHFHadMinus/I");
    tree_->Branch("nPFHFEmPlus",nPFHFEmPlus_,"nPFHFEmPlus/I");
    tree_->Branch("nPFHFEmMinus",nPFHFEmMinus_,"nPFHFEmMinus/I");
    tree_->Branch("pfSlice0P_energy",pfSlice0P_energy_,"pfSlice0P_energy/D");
    tree_->Branch("pfSlice1P_energy",pfSlice1P_energy_,"pfSlice1P_energy/D");
    tree_->Branch("pfSlice2P_energy",pfSlice2P_energy_,"pfSlice2P_energy/D");
    tree_->Branch("pfSlice3P_energy",pfSlice3P_energy_,"pfSlice3P_energy/D");
    tree_->Branch("pfSlice4P_energy",pfSlice4P_energy_,"pfSlice4P_energy/D");
    tree_->Branch("pfSlice0N_energy",pfSlice0N_energy_,"pfSlice0N_energy/D");
    tree_->Branch("pfSlice1N_energy",pfSlice1N_energy_,"pfSlice1N_energy/D");
    tree_->Branch("pfSlice2N_energy",pfSlice2N_energy_,"pfSlice2N_energy/D");
    tree_->Branch("pfSlice3N_energy",pfSlice3N_energy_,"pfSlice3N_energy/D");
    tree_->Branch("pfSlice4N_energy",pfSlice4N_energy_,"pfSlice4N_energy/D");
    tree_->Branch("pfSlice0P_pt",pfSlice0P_pt_,"pfSlice0P_pt/D");
    tree_->Branch("pfSlice1P_pt",pfSlice1P_pt_,"pfSlice1P_pt/D");
    tree_->Branch("pfSlice2P_pt",pfSlice2P_pt_,"pfSlice2P_pt/D");
    tree_->Branch("pfSlice3P_pt",pfSlice3P_pt_,"pfSlice3P_pt/D");
    tree_->Branch("pfSlice4P_pt",pfSlice4P_pt_,"pfSlice4P_pt/D");
    tree_->Branch("pfSlice0N_pt",pfSlice0N_pt_,"pfSlice0N_pt/D");
    tree_->Branch("pfSlice1N_pt",pfSlice1N_pt_,"pfSlice1N_pt/D");
    tree_->Branch("pfSlice2N_pt",pfSlice2N_pt_,"pfSlice2N_pt/D");
    tree_->Branch("pfSlice3N_pt",pfSlice3N_pt_,"pfSlice3N_pt/D");
    tree_->Branch("pfSlice4N_pt",pfSlice4N_pt_,"pfSlice4N_pt/D");
    tree_->Branch("npfSlice0P",npfSlice0P_,"npfSlice0P/I");
    tree_->Branch("npfSlice1P",npfSlice1P_,"npfSlice1P/I");
    tree_->Branch("npfSlice2P",npfSlice2P_,"npfSlice2P/I");
    tree_->Branch("npfSlice3P",npfSlice3P_,"npfSlice3P/I");
    tree_->Branch("npfSlice4P",npfSlice4P_,"npfSlice4P/I");
    tree_->Branch("npfSlice0N",npfSlice0N_,"npfSlice0N/I");
    tree_->Branch("npfSlice1N",npfSlice1N_,"npfSlice1N/I");
    tree_->Branch("npfSlice2N",npfSlice2N_,"npfSlice2N/I");
    tree_->Branch("npfSlice3N",npfSlice3N_,"npfSlice3N/I");
    tree_->Branch("npfSlice4N",npfSlice4N_,"npfSlice4N/I");
    tree_->Branch("nmuonSlice0P",nmuonSlice0P_,"nmuonSlice0P/I");
    tree_->Branch("nmuonSlice1P",nmuonSlice1P_,"nmuonSlice1P/I");
    tree_->Branch("nmuonSlice2P",nmuonSlice2P_,"nmuonSlice2P/I");
    tree_->Branch("nmuonSlice3P",nmuonSlice3P_,"nmuonSlice3P/I");
    tree_->Branch("nmuonSlice4P",nmuonSlice4P_,"nmuonSlice4P/I");
    tree_->Branch("nmuonSlice0N",nmuonSlice0N_,"nmuonSlice0N/I");
    tree_->Branch("nmuonSlice1N",nmuonSlice1N_,"nmuonSlice1N/I");
    tree_->Branch("nmuonSlice2N",nmuonSlice2N_,"nmuonSlice2N/I");
    tree_->Branch("nmuonSlice3N",nmuonSlice3N_,"nmuonSlice3N/I");
    tree_->Branch("nmuonSlice4N",nmuonSlice4N_,"nmuonSlice4N/I");
    tree_->Branch("nelectronSlice0P",nelectronSlice0P_,"nelectronSlice0P/I");
    tree_->Branch("nelectronSlice1P",nelectronSlice1P_,"nelectronSlice1P/I");
    tree_->Branch("nelectronSlice2P",nelectronSlice2P_,"nelectronSlice2P/I");
    tree_->Branch("nelectronSlice3P",nelectronSlice3P_,"nelectronSlice3P/I");
    tree_->Branch("nelectronSlice4P",nelectronSlice4P_,"nelectronSlice4P/I");
    tree_->Branch("nelectronSlice0N",nelectronSlice0N_,"nelectronSlice0N/I");
    tree_->Branch("nelectronSlice1N",nelectronSlice1N_,"nelectronSlice1N/I");
    tree_->Branch("nelectronSlice2N",nelectronSlice2N_,"nelectronSlice2N/I");
    tree_->Branch("nelectronSlice3N",nelectronSlice3N_,"nelectronSlice3N/I");
    tree_->Branch("nelectronSlice4N",nelectronSlice4N_,"nelectronSlice4N/I");
  }

  if(includePPSInfo_){
    tree_->Branch("protonsArm",&protonsArm_);
    tree_->Branch("protonsStation",&protonsStation_);
    tree_->Branch("protonsRP",&protonsRP_);
    tree_->Branch("protonsX",&protonsX_);
    tree_->Branch("protonsXUnc",&protonsXUnc_);
    tree_->Branch("protonsY",&protonsY_);
    tree_->Branch("protonsYUnc",&protonsYUnc_);
    if(!includeProtonsReco_){
      tree_->Branch("protonsIsValid",&protonsIsValid_);
      tree_->Branch("protonsTx",&protonsTx_);
      tree_->Branch("protonsTxUnc",&protonsTxUnc_);
      tree_->Branch("protonsTy",&protonsTy_);
      tree_->Branch("protonsTyUnc",&protonsTyUnc_);
      tree_->Branch("protonsTime",&protonsTime_);
      tree_->Branch("protonsTimeUnc",&protonsTimeUnc_);
    }
  }

  if(includeProtonsReco_){
    tree_->Branch("singleProtonArm",&singleProtonArm_);
    tree_->Branch("singleProtonStation",&singleProtonStation_);
    tree_->Branch("singleProtonPot",&singleProtonPot_);
    tree_->Branch("singleProtonXi",&singleProtonXi_);
    tree_->Branch("singleProtonT",&singleProtonT_);
    tree_->Branch("singleProtonThetaX",&singleProtonThetaX_);
    tree_->Branch("singleProtonThetaY",&singleProtonThetaY_);
    tree_->Branch("multiProtonArm",&multiProtonArm_);
    tree_->Branch("multiProtonXi",&multiProtonXi_);
    tree_->Branch("multiProtonT",&multiProtonT_);
    tree_->Branch("multiProtonTime",&multiProtonTime_);
    tree_->Branch("multiProtonTimeError",&multiProtonTimeError_);
    tree_->Branch("multiProtonVz",&multiProtonVz_);
    tree_->Branch("multiProtonThetaX",&multiProtonThetaX_);
    tree_->Branch("multiProtonThetaY",&multiProtonThetaY_);
  }

}

// ------------ method called once each job just after ending the event loop  ------------
  void 
MissingMassSearches::endJob() 
{

  delete run_;
  delete ev_;
  delete lumiblock_;

  delete trigger_;
  delete prescalesL1_;
  delete prescalesHLT_;

  delete xangle_;
  delete betastar_;
  delete bunchesb1_;
  delete bunchesb2_;
  delete intensityb1_;
  delete intensityb2_;
  delete instlumi_;

  delete PUinter_;
  delete PUtrueinter_;

  if(includeVertices_){
    delete vertex_x_;
    delete vertex_y_;
    delete vertex_z_;
    delete vertex_ntrack_;
    delete vertex_chi2_;
    delete vertex_ndof_;
  }

  if(enableMC_){
    delete nGenParticles_;
    if(includeElectrons_ || includeMuons_){
      delete nGenMuons_;
      delete nGenElectrons_;
      delete genleptons_status_;
      delete genleptons_pdgid_;
      delete genleptons_energy_;
      delete genleptons_pt_;
      delete genleptons_eta_;
      delete genleptons_phi_;
      delete genleptons_px_;
      delete genleptons_py_;
      delete genleptons_pz_;
      delete genleptons_charge_;
      delete genleptons_vx_;
      delete genleptons_vy_;
      delete genleptons_vz_;
    }
    if(includePhotons_){
      delete genphotons_pt_;
      delete genphotons_eta_;
      delete genphotons_phi_;
      delete genphotons_energy_;
    }
    if(includeProtonsReco_){
      delete nGenProtons_;
      delete genprotons_status_;
      delete genprotons_energy_;
      delete genprotons_pt_;
      delete genprotons_eta_;
      delete genprotons_phi_;
      delete genprotons_px_;
      delete genprotons_py_;
      delete genprotons_pz_;
      delete genprotons_xi_;
    }
    if(includeMET_){
      delete genmisset_;
      delete genmisset_phi_;
    }
    if(includeJets_){
      delete genjetsak4_px_;
      delete genjetsak4_py_;
      delete genjetsak4_pz_;
      delete genjetsak4_pt_;
      delete genjetsak4_energy_;
      delete genjetsak4_phi_;
      delete genjetsak4_eta_;
      delete genjetsak4_vz_;
      delete genjetsak8_px_;
      delete genjetsak8_py_;
      delete genjetsak8_pz_;
      delete genjetsak8_pt_;
      delete genjetsak8_energy_;
      delete genjetsak8_phi_;
      delete genjetsak8_eta_;
      delete genjetsak8_vz_;
    }
  }

  if(includeJets_){
    delete jetsak4_pt_;
    delete jetsak4_px_;
    delete jetsak4_py_;
    delete jetsak4_pz_;
    delete jetsak4_energy_;
    delete jetsak4_phi_;
    delete jetsak4_eta_;
    delete jetsak4_vx_;
    delete jetsak4_vy_;
    delete jetsak4_vz_;
    delete jetsak4_bdis_;
    delete jetsak4_qgdis_;
    delete jetsak4_neutralEmFraction_;
    delete jetsak4_neutralHadFraction_;
    delete jetsak4_chargedEmFraction_;
    delete jetsak4_chargedHadFraction_;
    delete jetsak4_muonFraction_;
    delete jetsak4_neutralMultiplicity_;
    delete jetsak4_chargedMultiplicity_;
    delete jetsak4_looseJetId_;
    delete jetsak4_tightJetId_;
    delete jetsak4_lepVetoJetId_;
    delete jetsak4_puIdfdisc_;
    delete jetsak4_puIdcbased_;
    delete jetsak4_puIdfid_;
    delete jetsak8_px_;
    delete jetsak8_py_;
    delete jetsak8_pz_;
    delete jetsak8_pt_;
    delete jetsak8_energy_;
    delete jetsak8_phi_;
    delete jetsak8_eta_;
    delete jetsak8_vx_;
    delete jetsak8_vy_;
    delete jetsak8_vz_;
    delete jetsak8_bdis_;
    delete jetsak8_neutralEmFraction_;
    delete jetsak8_neutralHadFraction_;
    delete jetsak8_chargedEmFraction_;
    delete jetsak8_chargedHadFraction_;
    delete jetsak8_muonFraction_;
    delete jetsak8_neutralMultiplicity_;
    delete jetsak8_chargedMultiplicity_;
    delete jetsak8_looseJetId_;
    delete jetsak8_tightJetId_;
    delete jetsak8_lepVetoJetId_;
  }

  if(includeElectrons_ || includeMuons_){
    delete nMuons_;
    delete nMuons_looseId_;
    delete nMuons_mediumId_;
    delete nMuons_tightId_;
    delete nMuons_pfIsoMedium_;
    delete nMuons_miniIsoTight_;
    delete nMuons_pfIsoVeryTight_;
    delete nElectrons_;
    delete nElectrons_looseId_;
    delete nElectrons_mediumId_;
    delete nElectrons_tightId_;
    delete leptons_pdgid_;
    delete leptons_energy_;
    delete leptons_pt_;
    delete leptons_eta_;
    delete leptons_phi_;
    delete leptons_px_;
    delete leptons_py_;
    delete leptons_pz_;
    delete leptons_charge_;
    delete leptons_vx_;
    delete leptons_vy_;
    delete leptons_vz_;
    delete leptons_looseId_;
    delete leptons_mediumId_;
    delete leptons_tightId_;
    delete leptons_pfIsoMedium_;
    delete leptons_miniIsoTight_;
    delete leptons_pfIsoVeryTight_;
    delete leptons_pfIso_;
    delete leptons_tkIso_;
  }

  if(includePhotons_){
    delete photons_pt_;
    delete photons_eta_;
    delete photons_phi_;
    delete photons_energy_;
    delete photons_r9_;
    delete photons_looseId_;
    delete photons_mediumId_;
    delete photons_tightId_;
  }

  if(includeMET_){
    delete misset_;
    delete misset_phi_;
  }

  if(includePF_){
    delete nPF_;
    delete SumPF_energy_;
    delete pfHFHadPlus_energy_;
    delete pfHFHadPlus_pt_;
    delete pfHFHadMinus_energy_;
    delete pfHFHadMinus_pt_;
    delete pfHFEmPlus_energy_;
    delete pfHFEmPlus_pt_;
    delete pfHFEmMinus_energy_;
    delete pfHFEmMinus_pt_;
    delete pfEcalMinus_energy_;
    delete pfEcalPlus_energy_;
    delete pfHcalMinus_energy_;
    delete pfHcalPlus_energy_;
    delete pfHadNeutralMinus_energy_;
    delete pfHadNeutralPlus_energy_;
    delete nChargedPFMultiPV_Loose_;
    delete nChargedPFMultiPV_Tight_;
    delete nChargedPFMultiPV_UsedInFit_;
    delete nChargedPFMultiPV_Tight_Fit_;
    delete SumChargedPFMultiPV_pt_Loose_;
    delete SumChargedPFMultiPV_pt_Tight_;
    delete SumChargedPFMultiPV_pt_UsedInFit_;  
    delete SumChargedPFMultiPV_pt_Tight_Fit_; 
    delete nPFPlus_;
    delete nPFMinus_;
    delete nPFHFHadPlus_;
    delete nPFHFHadMinus_;
    delete nPFHFEmPlus_;
    delete nPFHFEmMinus_;
    delete pfSlice0P_energy_;
    delete pfSlice1P_energy_;
    delete pfSlice2P_energy_;
    delete pfSlice3P_energy_;
    delete pfSlice4P_energy_;
    delete pfSlice0N_energy_;
    delete pfSlice1N_energy_;
    delete pfSlice2N_energy_;
    delete pfSlice3N_energy_;
    delete pfSlice4N_energy_;
    delete pfSlice0P_pt_;
    delete pfSlice1P_pt_;
    delete pfSlice2P_pt_;
    delete pfSlice3P_pt_;
    delete pfSlice4P_pt_;
    delete pfSlice0N_pt_;
    delete pfSlice1N_pt_;
    delete pfSlice2N_pt_;
    delete pfSlice3N_pt_;
    delete pfSlice4N_pt_;
    delete npfSlice0P_;
    delete npfSlice1P_;
    delete npfSlice2P_;
    delete npfSlice3P_;
    delete npfSlice4P_;
    delete npfSlice0N_;
    delete npfSlice1N_;
    delete npfSlice2N_;
    delete npfSlice3N_;
    delete npfSlice4N_;
    delete nmuonSlice0P_;
    delete nmuonSlice1P_;
    delete nmuonSlice2P_;
    delete nmuonSlice3P_;
    delete nmuonSlice4P_;
    delete nmuonSlice0N_;
    delete nmuonSlice1N_;
    delete nmuonSlice2N_;
    delete nmuonSlice3N_;
    delete nmuonSlice4N_;
    delete nelectronSlice0P_;
    delete nelectronSlice1P_;
    delete nelectronSlice2P_;
    delete nelectronSlice3P_;
    delete nelectronSlice4P_;
    delete nelectronSlice0N_;
    delete nelectronSlice1N_;
    delete nelectronSlice2N_;
    delete nelectronSlice3N_;
    delete nelectronSlice4N_;
  }

  if(includePPSInfo_){
    delete protonsArm_;
    delete protonsStation_;
    delete protonsRP_;
    delete protonsX_;
    delete protonsXUnc_;
    delete protonsY_;
    delete protonsYUnc_;
    if(!includeProtonsReco_){
      delete protonsIsValid_;
      delete protonsTx_;
      delete protonsTxUnc_;
      delete protonsTy_;
      delete protonsTyUnc_;
      delete protonsTime_;
      delete protonsTimeUnc_;
    }
  }

  if(includeProtonsReco_){
    delete singleProtonArm_;
    delete singleProtonStation_;
    delete singleProtonPot_;
    delete singleProtonXi_;
    delete singleProtonT_;
    delete singleProtonThetaX_;
    delete singleProtonThetaY_;
    delete multiProtonArm_;
    delete multiProtonXi_;
    delete multiProtonT_;
    delete multiProtonTime_;
    delete multiProtonTimeError_;
    delete multiProtonVz_;
    delete multiProtonThetaX_;
    delete multiProtonThetaY_;
  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MissingMassSearches::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ templates ------------

// ------------ DiSystemFilter
template <class T, class W>
bool MissingMassSearches::DiSystemFilter(T obj1_, W obj2_){
  bool selected = false;
  math::XYZTLorentzVector DiSys(0.,0.,0.,0.);
  DiSys = MissingMassSearches::DiSystem(obj1_, obj2_);
  if(DiSys.pt()>20.) selected = true;
  return selected;
}


// ------------ Invariant Mass
template <class T, class W>
math::XYZTLorentzVector MissingMassSearches::DiSystem(T obj1_, W obj2_){
  math::XYZTLorentzVector DiObj(0.,0.,0.,0.);
  DiObj += obj1_->p4();
  DiObj += obj2_->p4();
  return DiObj;
}

// ------------ Transverse Mass
template <class T, class W>
double MissingMassSearches::TransverseMass(T lepton_, W met_){

  double w_et = met_->et() + lepton_->pt();
  double w_px = met_->px() + lepton_->px();
  double w_py = met_->py() + lepton_->py();

  double tmass_ = w_et*w_et - w_px*w_px - w_py*w_py;
  tmass_ = (tmass_ > 0) ? sqrt(tmass_) : 0;

  return tmass_;

}

// ------------ Disentangle jets/leptons
template <class T, class W>
W MissingMassSearches::Unmatching(T lepton_, W jet_, double radius){

  W jetsCand_;
  jetsCand_.clear();

  // storing jets AK4 candidates which are not matching each lepton candidate (highest and second highest pt)
  try{
    for (std::size_t i = 0; i != jet_.size(); ++i) {
      if(!(lepton_.size()>1)) continue;
      if ((deltaR(lepton_[0]->eta(), lepton_[0]->phi(), jet_[i]->eta(), jet_[i]->phi()) > radius)&&(deltaR(lepton_[1]->eta(), lepton_[1]->phi(), jet_[i]->eta(), jet_[i]->phi())>radius)){
	jetsCand_.push_back(jet_[i]);
      }
    }
  }catch(...){}

  return jetsCand_;

}


template <class T, class W>
void MissingMassSearches::fillingJets(T jetak4_, W jetak8_){

  double NHFID, NEMFID, CHFID, MUFID, CEMFID, NumConstID, CHMID = -999;
  bool LooseJetID = false;
  bool TightJetID = false;
  bool tightLepVetoJetID = false;

  double vx_mean = -999.;
  double vy_mean = -999.;
  double vz_mean = -999.;
  double pt2_x = 0.;
  double pt2_y = 0.;
  double pt2_z = 0.;
  double pt2 = 0.;

  //#define sjets (typeid(jetak4_) == typeid(std::vector<const pat::Jet*>))
  //#define pfjets (!sjets)

  for ( const auto& jetsak4Evt: jetak4_) {
    if(debug_) std::cout << "<Candidate> Jets (AK4) pT: " << jetsak4Evt->pt() << " [GeV], eta: " << jetsak4Evt->eta() << ", phi: "<< jetsak4Evt->phi() << std::endl;
    (*jetsak4_px_).push_back(jetsak4Evt->px());
    (*jetsak4_py_).push_back(jetsak4Evt->py());
    (*jetsak4_pz_).push_back(jetsak4Evt->pz());
    (*jetsak4_pt_).push_back(jetsak4Evt->pt());
    (*jetsak4_energy_).push_back(jetsak4Evt->energy());
    (*jetsak4_phi_).push_back(jetsak4Evt->phi());
    (*jetsak4_eta_).push_back(jetsak4Evt->eta());
    if(tier_.find("mini") != std::string::npos || tier_.find("MINI") != std::string::npos){
      (*jetsak4_puIdfdisc_).push_back(jetsak4Evt->userFloat("pileupJetId:fullDiscriminant"));
      (*jetsak4_puIdcbased_).push_back(0);
      (*jetsak4_puIdfid_).push_back(jetsak4Evt->userInt("pileupJetId:fullId"));
      (*jetsak4_qgdis_).push_back(jetsak4Evt->userFloat("QGTagger:qgLikelihood"));
      (*jetsak4_bdis_).push_back(jetsak4Evt->bDiscriminator("pfDeepCSVJetTags:probb") + jetsak4Evt->bDiscriminator("pfDeepCSVJetTags:probbb"));
    }else{
      (*jetsak4_puIdfdisc_).push_back(jetsak4Evt->userFloat("AK4PFCHSpileupJetIdEvaluator:fullDiscriminant"));
      (*jetsak4_puIdcbased_).push_back(jetsak4Evt->userInt("AK4PFCHSpileupJetIdEvaluator:cutbasedId"));
      (*jetsak4_puIdfid_).push_back(jetsak4Evt->userInt("AK4PFCHSpileupJetIdEvaluator:fullId"));
      (*jetsak4_qgdis_).push_back(jetsak4Evt->userFloat("QGTaggerAK4PFCHS:qgLikelihood"));
      (*jetsak4_bdis_).push_back(jetsak4Evt->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    }
    if(jetsak4Evt->isJPTJet() || jetsak4Evt->isPFJet()){
      NHFID  = jetsak4Evt->neutralHadronEnergyFraction();
      NEMFID = jetsak4Evt->neutralEmEnergyFraction();
      CHFID  = jetsak4Evt->chargedHadronEnergyFraction();
      MUFID  = jetsak4Evt->muonEnergyFraction();
      CEMFID = jetsak4Evt->chargedEmEnergyFraction();
      NumConstID = jetsak4Evt->chargedMultiplicity()+jetsak4Evt->neutralMultiplicity();
      CHMID      = jetsak4Evt->chargedMultiplicity();

      (*jetsak4_neutralEmFraction_).push_back(jetsak4Evt->neutralEmEnergyFraction());
      (*jetsak4_neutralHadFraction_).push_back(jetsak4Evt->neutralHadronEnergyFraction());
      (*jetsak4_chargedEmFraction_).push_back(jetsak4Evt->chargedEmEnergyFraction());
      (*jetsak4_chargedHadFraction_).push_back(jetsak4Evt->chargedHadronEnergyFraction());
      (*jetsak4_muonFraction_).push_back(jetsak4Evt->muonEnergyFraction());
      (*jetsak4_neutralMultiplicity_).push_back(jetsak4Evt->neutralMultiplicity());
      (*jetsak4_chargedMultiplicity_).push_back(jetsak4Evt->chargedMultiplicity());

      if(debug_){
	std::cout << "\n --> ak4" << std::endl;
	std::cout << "Neutral EM f: " << jetsak4Evt->neutralEmEnergyFraction() << std::endl;
	std::cout << "Neutral Had f: " << jetsak4Evt->neutralHadronEnergyFraction() << std::endl;
	std::cout << "Charged EM f: " << jetsak4Evt->chargedEmEnergyFraction() << std::endl;
	std::cout << "Charged Had f: " << jetsak4Evt->chargedHadronEnergyFraction() << std::endl;
	std::cout << "Muon f: " << jetsak4Evt->muonEnergyFraction() << std::endl;
      }

      if(jetsak4Evt->eta()<=2.7&&jetsak4Evt->eta()>=-2.7){
	if(jetsak4Evt->eta()<=2.4&&jetsak4Evt->eta()>=-2.4){
	  LooseJetID = ((NHFID<0.99 && NEMFID<0.99 && NumConstID>1 && CHFID>0 && CHMID>0));
	  TightJetID = ((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 && CHFID>0 && CHMID>0));
	  tightLepVetoJetID=((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 && MUFID<0.8 && CHFID>0 && CHMID>0 && CEMFID<0.80));
	}else{
	  LooseJetID = ((NHFID<0.99 && NEMFID<0.99 && NumConstID>1));
	  TightJetID = ((NHFID<0.90 && NEMFID<0.90 && NumConstID>1));
	  tightLepVetoJetID=((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 &&MUFID<0.8));
	}
      }
      if((jetsak4Evt->eta()>2.7&&jetsak4Evt->eta()<=3)||(jetsak4Evt->eta()<-2.7&&jetsak4Evt->eta()>=-3)){
	LooseJetID = (NEMFID<0.99 && jetsak4Evt->neutralMultiplicity()>1);
	TightJetID = (NEMFID<0.99 && NEMFID>0.02 && jetsak4Evt->neutralMultiplicity()>2);
	tightLepVetoJetID=(NEMFID<0.90 && NEMFID>0.02 && jetsak4Evt->neutralMultiplicity()>2);
      }
      if(jetsak4Evt->eta()>3||jetsak4Evt->eta()<-3){
	LooseJetID = ((NHFID>0.02 && NEMFID<0.99 && jetsak4Evt->neutralMultiplicity()>5));
	TightJetID = ((NHFID>0.02 && NEMFID<0.90 && jetsak4Evt->neutralMultiplicity()>10));
	tightLepVetoJetID=((NHFID>0.02 && NEMFID<0.90 && jetsak4Evt->neutralMultiplicity()>15));
      }

      (*jetsak4_looseJetId_).push_back(LooseJetID);
      (*jetsak4_tightJetId_).push_back(TightJetID);
      (*jetsak4_lepVetoJetId_).push_back(tightLepVetoJetID);
    }else{
      (*jetsak4_looseJetId_).push_back(false);
      (*jetsak4_tightJetId_).push_back(false);
      (*jetsak4_lepVetoJetId_).push_back(false);
      (*jetsak4_neutralEmFraction_).push_back(-1);
      (*jetsak4_neutralHadFraction_).push_back(-1);
      (*jetsak4_chargedEmFraction_).push_back(-1);
      (*jetsak4_chargedHadFraction_).push_back(-1);
      (*jetsak4_muonFraction_).push_back(-1);
      (*jetsak4_neutralMultiplicity_).push_back(-1);
      (*jetsak4_chargedMultiplicity_).push_back(-1);

    }

    reco::CompositePtrCandidate::daughters pfconst = jetsak4Evt->daughterPtrVector();
    for (reco::CompositePtrCandidate::daughters::const_iterator itpf = pfconst.begin(); itpf != pfconst.end(); ++itpf) {
      pt2 += (*itpf)->pt()*(*itpf)->pt();
      pt2_x += (*itpf)->pt()*(*itpf)->pt()*(*itpf)->vx();
      pt2_y += (*itpf)->pt()*(*itpf)->pt()*(*itpf)->vy();
      pt2_z += (*itpf)->pt()*(*itpf)->pt()*(*itpf)->vz();
    }

    if (pt2 > 0.) {
      vx_mean = pt2_x/pt2;
      vy_mean = pt2_y/pt2;
      vz_mean = pt2_z/pt2;
    }

    (*jetsak4_vx_).push_back(vx_mean);
    (*jetsak4_vy_).push_back(vy_mean);
    (*jetsak4_vz_).push_back(vz_mean);

  }

  NHFID = -999;
  NEMFID = -999;
  CHFID = -999;
  MUFID = -999;
  CEMFID = -999;
  NumConstID = -999;
  CHMID = -999;
  LooseJetID = false;
  TightJetID = false;
  tightLepVetoJetID = false;

  vx_mean = -999.;
  vy_mean = -999.;
  vz_mean = -999.;
  pt2_x = 0.;
  pt2_y = 0.;
  pt2_z = 0.;
  pt2 = 0.;

  for ( const auto& jetsak8Evt: jetak8_) {
    if(debug_) std::cout << "<Candidate> Jets (AK8) pT: " << jetsak8Evt->pt() << " [GeV], eta: " << jetsak8Evt->eta() << ", phi: "<< jetsak8Evt->phi() << std::endl;
    (*jetsak8_px_).push_back(jetsak8Evt->px());
    (*jetsak8_py_).push_back(jetsak8Evt->py());
    (*jetsak8_pz_).push_back(jetsak8Evt->pz());
    (*jetsak8_pt_).push_back(jetsak8Evt->pt());
    (*jetsak8_energy_).push_back(jetsak8Evt->energy());
    (*jetsak8_phi_).push_back(jetsak8Evt->phi());
    (*jetsak8_eta_).push_back(jetsak8Evt->eta());

    //(*jetsak8_bdis_).push_back(jetsak8Evt->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    (*jetsak8_bdis_).push_back(jetsak8Evt->bDiscriminator("pfDeepCSVJetTags:probb") + jetsak8Evt->bDiscriminator("pfDeepCSVJetTags:probbb"));

    if(jetsak8Evt->isJPTJet() || jetsak8Evt->isPFJet()){
      NHFID  = jetsak8Evt->neutralHadronEnergyFraction();
      NEMFID = jetsak8Evt->neutralEmEnergyFraction();
      CHFID  = jetsak8Evt->chargedHadronEnergyFraction();
      MUFID  = jetsak8Evt->muonEnergyFraction();
      CEMFID = jetsak8Evt->chargedEmEnergyFraction();
      NumConstID = jetsak8Evt->chargedMultiplicity()+jetsak8Evt->neutralMultiplicity();
      CHMID      = jetsak8Evt->chargedMultiplicity();

      (*jetsak8_neutralEmFraction_).push_back(jetsak8Evt->neutralEmEnergyFraction());
      (*jetsak8_neutralHadFraction_).push_back(jetsak8Evt->neutralHadronEnergyFraction());
      (*jetsak8_chargedEmFraction_).push_back(jetsak8Evt->chargedEmEnergyFraction());
      (*jetsak8_chargedHadFraction_).push_back(jetsak8Evt->chargedHadronEnergyFraction());
      (*jetsak8_muonFraction_).push_back(jetsak8Evt->muonEnergyFraction());
      (*jetsak8_neutralMultiplicity_).push_back(jetsak8Evt->neutralMultiplicity());
      (*jetsak8_chargedMultiplicity_).push_back(jetsak8Evt->chargedMultiplicity());

      if(debug_){
	std::cout << "\n --> ak8" << std::endl;
	std::cout << "Neutral EM f: " << jetsak8Evt->neutralEmEnergyFraction() << std::endl;
	std::cout << "Neutral Had f: " << jetsak8Evt->neutralHadronEnergyFraction() << std::endl;
	std::cout << "Charged EM f: " << jetsak8Evt->chargedEmEnergyFraction() << std::endl;
	std::cout << "Charged Had f: " << jetsak8Evt->chargedHadronEnergyFraction() << std::endl;
	std::cout << "Muon f: " << jetsak8Evt->muonEnergyFraction() << std::endl;
      }

      if(jetsak8Evt->eta()<=2.7&&jetsak8Evt->eta()>=-2.7){
	if(jetsak8Evt->eta()<=2.4&&jetsak8Evt->eta()>=-2.4){
	  LooseJetID = ((NHFID<0.99 && NEMFID<0.99 && NumConstID>1 && CHFID>0 && CHMID>0));
	  TightJetID = ((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 && CHFID>0 && CHMID>0));
	  tightLepVetoJetID=((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 && MUFID<0.8 && CHFID>0 && CHMID>0 && CEMFID<0.80));
	}else{
	  LooseJetID = ((NHFID<0.99 && NEMFID<0.99 && NumConstID>1));
	  TightJetID = ((NHFID<0.90 && NEMFID<0.90 && NumConstID>1));
	  tightLepVetoJetID=((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 &&MUFID<0.8));
	}
      }
      if((jetsak8Evt->eta()>2.7&&jetsak8Evt->eta()<=3)||(jetsak8Evt->eta()<-2.7&&jetsak8Evt->eta()>=-3)){
	LooseJetID = (NEMFID<0.99 && jetsak8Evt->neutralMultiplicity()>1);
	TightJetID = (NEMFID<0.99 && NEMFID>0.02 && jetsak8Evt->neutralMultiplicity()>2);
	tightLepVetoJetID=(NEMFID<0.90 && NEMFID>0.02 && jetsak8Evt->neutralMultiplicity()>2);
      }
      if(jetsak8Evt->eta()>3||jetsak8Evt->eta()<-3){
	LooseJetID = ((NHFID>0.02 && NEMFID<0.99 && jetsak8Evt->neutralMultiplicity()>5));
	TightJetID = ((NHFID>0.02 && NEMFID<0.90 && jetsak8Evt->neutralMultiplicity()>10));
	tightLepVetoJetID=((NHFID>0.02 && NEMFID<0.90 && jetsak8Evt->neutralMultiplicity()>15));
      }

      (*jetsak8_looseJetId_).push_back(LooseJetID);
      (*jetsak8_tightJetId_).push_back(TightJetID);
      (*jetsak8_lepVetoJetId_).push_back(tightLepVetoJetID);
    }else{
      (*jetsak8_looseJetId_).push_back(false);
      (*jetsak8_tightJetId_).push_back(false);
      (*jetsak8_lepVetoJetId_).push_back(false);
      (*jetsak8_neutralEmFraction_).push_back(-1);
      (*jetsak8_neutralHadFraction_).push_back(-1);
      (*jetsak8_chargedEmFraction_).push_back(-1);
      (*jetsak8_chargedHadFraction_).push_back(-1);
      (*jetsak8_muonFraction_).push_back(-1);
      (*jetsak8_neutralMultiplicity_).push_back(-1);
      (*jetsak8_chargedMultiplicity_).push_back(-1);
    }

    reco::CompositePtrCandidate::daughters pfconst = jetsak8Evt->daughterPtrVector();
    for (reco::CompositePtrCandidate::daughters::const_iterator itpf = pfconst.begin(); itpf != pfconst.end(); ++itpf) {
      pt2 += (*itpf)->pt()*(*itpf)->pt();
      pt2_x += (*itpf)->pt()*(*itpf)->pt()*(*itpf)->vx();
      pt2_y += (*itpf)->pt()*(*itpf)->pt()*(*itpf)->vy();
      pt2_z += (*itpf)->pt()*(*itpf)->pt()*(*itpf)->vz();
    }

    if (pt2 > 0.) {
      vx_mean = pt2_x/pt2;
      vy_mean = pt2_y/pt2;
      vz_mean = pt2_z/pt2;
    }

    (*jetsak8_vx_).push_back(vx_mean);
    (*jetsak8_vy_).push_back(vy_mean);
    (*jetsak8_vz_).push_back(vz_mean);

  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(MissingMassSearches);
