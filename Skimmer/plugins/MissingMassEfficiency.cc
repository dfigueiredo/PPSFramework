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
// Original Authors:  Diego Figueiredo, Lorenzo Pagliai and Nicola Turini
//          Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

// user include files
#include "MissingMassEfficiency.h"

//
// constructors and destructor
//
MissingMassEfficiency::MissingMassEfficiency(const edm::ParameterSet& iConfig):
  debugging_               ( iConfig.getParameter<bool>( "debugging" ) ),
  verbosity_               ( iConfig.getParameter<int>( "verbosity" ) ),
  hltPrescaleProvider_     ( iConfig, consumesCollector(), *this),
  physics_                 ( iConfig.getParameter<std::string>( "physics" ) ),
  unmatching_              ( iConfig.getParameter<bool>( "unmatching" ) ),
  orderDiscriminator_      ( iConfig.getParameter<bool>( "orderDiscriminator" ) ),
  enableMC_                ( iConfig.getParameter<bool>( "enableMC" ) ),
  enableTrigger_           ( iConfig.getParameter<bool>( "enableTrigger" ) ),
  enablePrescales_         ( iConfig.getParameter<bool>( "enablePrescales" ) ),
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
  photonToken_         ( consumes<edm::View<pat::Photon>>                            ( iConfig.getParameter<edm::InputTag>( "photonTag" ) ) ),
  patmetToken_         ( consumes<edm::View<pat::MET>>                               ( iConfig.getParameter<edm::InputTag>( "patmetTag" ) ) ),
  vertexToken_         ( consumes<edm::View<reco::Vertex>>                           ( iConfig.getParameter<edm::InputTag>( "vertexTag" ) ) ),
  protonSingleToken_   ( consumes<reco::ForwardProtonCollection>                     ( iConfig.getParameter<edm::InputTag>( "protonSingleTag" ) ) ),
  protonMultiToken_    ( consumes<reco::ForwardProtonCollection>                     ( iConfig.getParameter<edm::InputTag>( "protonMultiTag" ) ) )
{

  //now do what ever initialization is needed
  usesResource("TFileService");

  std::cout<<"\n<CMSSW Plugin> Missing Mass Skimmer..." << std::endl;
  std::cout<<"\t-->Physics Option: " << physics_ <<  std::endl;
  std::cout<<"\t-->Unmatching: " << unmatching_ <<  std::endl;
  std::cout<<"\t-->Btag Discriminator (ordering the jets): " << orderDiscriminator_ <<  std::endl;
  std::cout<<"\t-->EnableMC: " << enableMC_ <<  std::endl;
  std::cout<<"\t-->EnableTrigger: " << enableTrigger_ <<  std::endl;
  std::cout<<"\t-->EnablePrescales: " << enablePrescales_ <<  std::endl;

}

MissingMassEfficiency::~MissingMassEfficiency()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ running over genparticles -----
void MissingMassEfficiency::fetchGenParticles(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  GenParticlesEvent *genpart;
  genpart = new GenParticlesEvent(iEvent, iSetup, GenPartToken_);
  genpartVec = genpart->GetGenParticles();

  GenProtonsEvent *genproton;
  genproton = new GenProtonsEvent(iEvent, iSetup, GenPartToken_);
  genprotonsVec = genproton->GetGenProtons();

}

// ------------- running over genmet ----------
void MissingMassEfficiency::fetchGenMET(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  GenMETEvent *genmissinget;
  genmissinget = new GenMETEvent(iEvent, iSetup, GenMETToken_);
  genmet = genmissinget->GetGenMET();

}

// ------------- running over genjets --------
void MissingMassEfficiency::fetchGenJets(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  GenJetsEvent *genjetAlgoA;
  genjetAlgoA = new GenJetsEvent(iEvent, iSetup, GenJetTokenA_);
  genjetsVecA = genjetAlgoA->GetGenJets();

  GenJetsEvent *genjetAlgoB;
  genjetAlgoB = new GenJetsEvent(iEvent, iSetup, GenJetTokenB_);
  genjetsVecB = genjetAlgoB->GetGenJets();

}

// ------------ running over muons  ------------
void MissingMassEfficiency::fetchMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  MuonsEvent *mu;
  mu = new MuonsEvent(iEvent, iSetup, muonToken_);
  muonsVec = mu->GetMuons();

}

// ------------ running over electrons  ------------
void MissingMassEfficiency::fetchElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  ElectronsEvent *ele;
  ele = new ElectronsEvent(iEvent, iSetup, eleToken_);
  electronsVec = ele->GetElectrons();

}

// ------------ running over jets  ------------
void MissingMassEfficiency::fetchJets(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  JetsEvent *jetAlgoA;
  jetAlgoA = new JetsEvent(iEvent, iSetup, patjetTokenA_);
  jetsVecA = jetAlgoA->GetJets();

  JetsEvent *jetAlgoB;
  jetAlgoB = new JetsEvent(iEvent, iSetup, patjetTokenB_);
  jetsVecB = jetAlgoB->GetJets();

}

// ------------ running over photons  ------------
void MissingMassEfficiency::fetchPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  PhotonsEvent *photon;
  photon = new PhotonsEvent(iEvent, iSetup, photonToken_);
  photonsVec = photon->GetPhotons();

}

// ------------ running over MET  ------------
void MissingMassEfficiency::fetchMET(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  METEvent *missinget;
  missinget = new METEvent(iEvent, iSetup, patmetToken_);
  met = missinget->GetMET();

}

// ------------ running over vertices  ------------
void MissingMassEfficiency::fetchVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  VerticesEvent *vtx;
  vtx = new VerticesEvent(iEvent, iSetup, vertexToken_);
  vtxVec = vtx->GetVertices();

}

// ------------ running over reco protons  ------------
void MissingMassEfficiency::fetchProtonsReco(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  ProtonsEvent *protonsSingle;
  protonsSingle = new ProtonsEvent(iEvent, iSetup, protonSingleToken_);
  protonSingleVec = protonsSingle->GetProtonsReco();

  ProtonsEvent *protonsMulti;
  protonsMulti = new ProtonsEvent(iEvent, iSetup, protonMultiToken_);
  protonMultiVec = protonsMulti->GetProtonsReco();

}

// ------------ saving trigger  ------------
bool MissingMassEfficiency::fetchTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup){

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

// ------------ method called for each event  ------------
  void
MissingMassEfficiency::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  cleaning();
  selectionHisto_->Fill("Total",1);

  // Trigger fired!
  if(fetchTrigger(iEvent,iSetup)){

    selectionHisto_->Fill("Triggered",1);

    // Filling Vectors
    if(enableMC_) {
      fetchGenParticles(iEvent, iSetup);
      fetchGenMET(iEvent, iSetup);
      fetchGenJets(iEvent, iSetup);
    }

    fetchMuons(iEvent, iSetup);
    fetchElectrons(iEvent, iSetup);
    fetchJets(iEvent, iSetup);
    fetchPhotons(iEvent, iSetup);
    fetchMET(iEvent, iSetup);
    fetchVertices(iEvent, iSetup);
    fetchProtonsReco(iEvent, iSetup); 
    selection(iEvent, iSetup); 

  } //end trigger

}

// ------------ method called for cleaning vectors memory ------------
void MissingMassEfficiency::cleaning(){

  genpartVec.clear();
  genprotonsVec.clear();
  genjetsVecA.clear();
  genjetsVecB.clear();
  jetsVecA.clear();
  jetsVecB.clear();
  electronsVec.clear();
  muonsVec.clear();
  photonsVec.clear();
  vtxVec.clear();
  protonSingleVec.clear();
  protonMultiVec.clear();
  jetsak4Cand.clear();
  jetsak8Cand.clear();


}

// ------------ method called for selecting events ------------
void MissingMassEfficiency::selection(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;

  std::sort(muonsVec.begin(), muonsVec.end(), orderPt());
  std::sort(electronsVec.begin(), electronsVec.end(), orderPt());
  std::sort(photonsVec.begin(), photonsVec.end(), orderPt());

  if(orderDiscriminator_){
    std::sort(jetsVecA.begin(), jetsVecA.end(), orderDiscriminator());
    std::sort(jetsVecB.begin(), jetsVecB.end(), orderDiscriminator());
  }else{
    std::sort(jetsVecA.begin(), jetsVecA.end(), orderPt());
    std::sort(jetsVecB.begin(), jetsVecB.end(), orderPt());
  }

  try{
    if(unmatching_){
      // storing jets AK4, AK8 candidates which are not matching each lepton candidate
      if(physics_ == "muons" || physics_ == "MUONS" || physics_ == "muon" || physics_=="MUON"){
	jetsak4Cand = Unmatching(muonsVec, jetsVecA, 0.5);
	jetsak8Cand = Unmatching(muonsVec, jetsVecB, 0.9);
      }else if(physics_ == "electrons" || physics_ == "ELECTRONS" || physics_ == "electron" || physics_=="ELECTRON"){
	jetsak4Cand = Unmatching(electronsVec, jetsVecA, 0.5);
	jetsak8Cand = Unmatching(electronsVec, jetsVecB, 0.9);
      }else{
	jetsak4Cand = jetsVecA;
	jetsak8Cand = jetsVecB;
      }
    }else{
      jetsak4Cand = jetsVecA;
      jetsak8Cand = jetsVecB;
    }
  }catch(...){}

  // Filling Candidate
  std::vector<double> cand_pt;
  std::vector<double> cand_eta;
  std::vector<double> cand_phi;
  std::vector<double> cand_energy;

  if(physics_ == "jet" || physics_ =="JET"){
    for(const auto& jetsak4Evt: jetsak4Cand){
      cand_pt.push_back(jetsak4Evt->pt());
      cand_eta.push_back(jetsak4Evt->eta());
      cand_phi.push_back(jetsak4Evt->phi());
      cand_energy.push_back(jetsak4Evt->energy());
    }
  }else if(physics_ == "electron" || physics_=="ELECTRON"){
    for(const auto& electronsEvt: electronsVec){
      cand_pt.push_back(electronsEvt->pt());
      cand_eta.push_back(electronsEvt->eta());
      cand_phi.push_back(electronsEvt->phi());
      cand_energy.push_back(electronsEvt->energy());
    }
  }else if(physics_=="muon" || physics_=="MUON"){
    for(const auto& muonsEvt: muonsVec){
      cand_pt.push_back(muonsEvt->pt());
      cand_eta.push_back(muonsEvt->eta());
      cand_phi.push_back(muonsEvt->phi());
      cand_energy.push_back(muonsEvt->energy());
    }
  }else if(physics_=="photon" || physics_=="PHOTON"){
    for(const auto& photonsEvt: photonsVec){
      cand_pt.push_back(photonsEvt->pt());
      cand_eta.push_back(photonsEvt->eta());
      cand_phi.push_back(photonsEvt->phi());
      cand_energy.push_back(photonsEvt->energy());
    }
  }else{
    throw cms::Exception("InputError") << "Please, insert a valid physics option (electrons, muons, jets or photons)" << std::endl;
    exit(EXIT_SUCCESS);
  }

  std::vector<double> protonArm0;
  std::vector<double> protonArm1;

  protonArm0.clear();
  protonArm1.clear();

  for(const auto& protonsEvt: protonMultiVec){
    const CTPPSDetId det_id_multi((*protonsEvt->contributingLocalTracks().begin())->getRPId());
    if(debugging_ && verbosity_>0) std::cout << "Proton Xi: " << protonsEvt->xi() << ", Arm: " << det_id_multi.arm() << std::endl;
    if(det_id_multi.arm()==0) protonArm0.push_back(protonsEvt->xi());
    if(det_id_multi.arm()==1) protonArm1.push_back(protonsEvt->xi());
  }

  std::sort(protonArm0.begin(), protonArm0.end());
  std::sort(protonArm1.begin(), protonArm1.end());

  if(protonArm0.size()>0&&protonArm1.size()>0){
    if(debugging_) std::cout << "Tagged [#arm0, #arm1]: (" << protonArm0.size() << ", " << protonArm1.size() << ")" << std::endl;

    for(const auto& protonsEvt: protonArm0){
      std::cout << "#xi, #arm0: " << protonArm0[protonsEvt] << std::endl;
    }

    for(const auto& protonsEvt: protonArm1){
      std::cout << "#xi, #arm1: " << protonArm1[protonsEvt] << std::endl;
    }

    p1Multi.SetPxPyPzE(0.,0., ECM*protonArm0[0]/2., ECM*protonArm0[0]/2.);
    p2Multi.SetPxPyPzE(0.,0., -ECM*protonArm1[0]/2., ECM*protonArm1[0]/2.);

    double DiffMass = 0.;
    double DiMass = 0.;
    double DiPt = 0.;
    double DiPhi = 0.;
    double DiEta = 0.;

    double DeltaPhi = 0.;
    double DeltaPt = 0.;

    double xDeltaPhi = 0.;
    double xDeltaPt = 0.;

    if(cand_pt.size()>1){

      cand1.SetPtEtaPhiE(cand_pt[0], cand_eta[0],cand_phi[0],cand_energy[0]);
      cand2.SetPtEtaPhiE(cand_pt[1], cand_eta[1],cand_phi[1],cand_energy[1]);
      disystem = cand1+cand2;

      DiMass = disystem.M();
      DiPt = disystem.Pt();
      DiPhi = disystem.Phi();
      DiEta = disystem.Eta();
      DeltaPhi = reduceRange(DiPhi - met->phi());
      DeltaPt = fabs(DiPt - met->et());

      diMassHisto_->Fill(DiMass);
      diPhiHisto_->Fill(DiPhi);
      diEtaHisto_->Fill(DiEta);
      diPtHisto_->Fill(DiPt);
      deltaPhiHisto_->Fill(DeltaPhi);
      deltaPtHisto_->Fill(DeltaPt);

      if(debugging_ && verbosity_>0){
	std::cout << "DiMass: " << disystem.M() << ", "<< DiMass << std::endl;
	std::cout << "DiPt: " << disystem.Pt() << " , "<< DiPt << std::endl;
	std::cout << "DiEta: " << disystem.Eta() << ", " << DiEta << std::endl;
	std::cout << "Delta Phi: " << DeltaPhi << std::endl;
	std::cout << "Delta Pt: " << DeltaPt << std::endl;
      }
      selectionHisto_->Fill("At least two jets",1);
    }

    for ( const auto& jetsak4Evt: jetsak4Cand) {
      if(debugging_ && verbosity_>0) std::cout << "Btagging discriminator: " << jetsak4Evt->bDiscriminator("pfDeepCSVJetTags:probb") + jetsak4Evt->bDiscriminator("pfDeepCSVJetTags:probbb") << ", pT[GeV] "<< jetsak4Evt->pt() << ", Parton Flavour: " << jetsak4Evt->partonFlavour() << std::endl;
    }

    xsystem = p1Multi + p2Multi - disystem;
    DiffMass = ECM*sqrt(protonArm0[0]*protonArm1[0]);

    xDeltaPhi = reduceRange(xsystem.Phi() - DiPhi);
    xDeltaPt = fabs(xsystem.Pt() - DiPt);

    std::cout << "X system (eta, phi, pt, M): (" << xsystem.Eta() << ", " << xsystem.Phi() << ", " << xsystem.Pt() << ", "<< xsystem.M() << ")" << std::endl;
    std::cout << "Diff Mass (xi1, xi2): " << DiffMass << std::endl;

    if(jetsak4Cand.size()>1){ 
      std::cout << "xDeltaPt: " << xDeltaPt << std::endl;
      std::cout << "DiPt: " << DiPt << std::endl;
      std::cout << "xPt: " << xsystem.Pt() << std::endl;
      std::cout << "xDeltaPhi: " << xDeltaPhi << std::endl;
      std::cout << "DiPhi: " << DiPhi << std::endl;
      std::cout << "xPhi: " << xsystem.Phi() << std::endl;
    }

    if(protonArm0.size()>0&&protonArm1.size()>0&&jetsak4Cand.size()>1){
      diffMassHisto_->Fill(DiffMass);
      xMassHisto_->Fill(xsystem.M());
      xPhiHisto_->Fill(xsystem.Phi());
      xEtaHisto_->Fill(xsystem.Eta());
      xPtHisto_->Fill(xsystem.Pt());
      xdeltaPhiHisto_->Fill(xDeltaPhi);
      xdeltaPtHisto_->Fill(xDeltaPt);
      xiArm0Histo_->Fill(protonArm0[0]);
      xiArm1Histo_->Fill(protonArm1[0]);
    }

    selectionHisto_->Fill("At least two protons",1);
  }

}

// ------------ method called when a run finishes ------------
void MissingMassEfficiency::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}

// ------------ method called when a new run number begins  ------------
void MissingMassEfficiency::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}

// ------------ method called once each job just before starting event loop  ------------
void MissingMassEfficiency::beginJob(){

  edm::Service<TFileService> fs;

  // Control Histograms
  selectionHisto_ = fs->make<TH1F>("Selection","Selection; Cuts; Events",1,0,1);
  selectionHisto_->SetCanExtend(TH1::kXaxis);

  diMassHisto_ = fs->make<TH1F>("DiSystemMass","Disystem Mass; Mass [GeV]; # Entries",500,0,3000);
  diPhiHisto_ = fs->make<TH1F>("DiSystemPhi","Disystem #phi; #phi; # Entries",100,-4,4);
  diEtaHisto_ = fs->make<TH1F>("DiSystemEta","Disystem #eta; #eta; # Entries",100,-6,6);
  diPtHisto_ = fs->make<TH1F>("DiSystemPt","Disystem p_{T}; p_{T}; # Entries",100,0,300);
  deltaPtHisto_ = fs->make<TH1F>("DeltaPt","#Delta p_{T} (#slash{E} - Disystem); #Delta p_{T}; # Entries",100,0,300);
  deltaPhiHisto_ = fs->make<TH1F>("DeltaPhi","#Delta#phi (#slash{E} - Disystem); #Delta#phi; # Entries",100,-4,4);

  xMassHisto_ = fs->make<TH1F>("XSystemMass","X system Mass; Mass [GeV]; # Entries",200,0,2500);
  xPhiHisto_ = fs->make<TH1F>("XSystemPhi","X system #phi; #phi; # Entries",100,-4,4);
  xEtaHisto_ = fs->make<TH1F>("XSystemEta","X system #eta; #eta; # Entries",100,-6,6);
  xPtHisto_ = fs->make<TH1F>("XSystemPt","X system p_{T}; p_{T}; # Entries",100,0,300);
  xdeltaPtHisto_ = fs->make<TH1F>("XDeltaPt","#Delta p_{T} (X - Disystem); #Delta p_{T}; # Entries",100,0,50);
  xdeltaPhiHisto_ = fs->make<TH1F>("XDeltaPhi","#Delta#phi (X - Disystem); #Delta#phi; # Entries",200,-4,4);

  xiArm0Histo_ = fs->make<TH1F>("XiArm0","#xi Arm 0; #xi; # Entries",100,0,0.5);
  xiArm1Histo_ = fs->make<TH1F>("XiArm1","#xi Arm 1; #xi; # Entries",100,0,0.5);
  diffMassHisto_ = fs->make<TH1F>("DiffractiveMass","Diffractive Mass; #sqrt{S#xi_{1}#xi_{2}} [GeV]; # Entries",200,0,2500);

}

// ------------ method called once each job just after ending the event loop  ------------
void MissingMassEfficiency::endJob(){}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MissingMassEfficiency::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ------------ templates ------------

// ------------ DiSystemFilter
template <class T, class W>
bool MissingMassEfficiency::DiSystemFilter(T obj1_, W obj2_){
  bool selected = false;
  math::XYZTLorentzVector DiSys(0.,0.,0.,0.);
  DiSys = MissingMassEfficiency::DiSystem(obj1_, obj2_);
  if(DiSys.pt()>20.) selected = true;
  return selected;
}


// ------------ Invariant Mass
template <class T, class W>
math::XYZTLorentzVector MissingMassEfficiency::DiSystem(T obj1_, W obj2_){
  math::XYZTLorentzVector DiObj(0.,0.,0.,0.);
  DiObj += obj1_->p4();
  DiObj += obj2_->p4();
  return DiObj;
}

// ------------ Transverse Mass
template <class T, class W>
double MissingMassEfficiency::TransverseMass(T lepton_, W met_){

  double w_et = met_->et() + lepton_->pt();
  double w_px = met_->px() + lepton_->px();
  double w_py = met_->py() + lepton_->py();

  double tmass_ = w_et*w_et - w_px*w_px - w_py*w_py;
  tmass_ = (tmass_ > 0) ? sqrt(tmass_) : 0;

  return tmass_;

}

// ------------ Disentangle jets/leptons
template <class T, class W>
W MissingMassEfficiency::Unmatching(T lepton_, W jet_, double radius){

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


// ------------ Angle
template<typename T>
constexpr T  MissingMassEfficiency::reduceRange(T x) {
  constexpr T o2pi = 1./(2.*M_PI);
  if (std::abs(x) <= T(M_PI)) return x;
  T n = std::round(x*o2pi);
  return x - n*T(2.*M_PI);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MissingMassEfficiency);
