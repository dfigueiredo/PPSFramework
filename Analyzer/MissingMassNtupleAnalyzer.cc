#define MissingMassNtupleAnalyzer_cxx

// ROOT header
#include <TSystem.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TRandom3.h"

// Header file for the classes stored in the TTree if any.
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <string>

// User Header
#include "statusbar.h"
#include "MissingMassNtupleAnalyzer.h"

// Code Constants
#define c_light 2.99792458e+8 // m/s
#define MASS_B 4.2 // GeV
#define MASS_MU 0.1057 // GeV
#define MASS_E 0.000511 // GeV
#define MASS_P 0.938272029 // GeV
#define pi 3.14159265359
#define ECM 13000.0 // GeV
#define ns_to_s_ 1e-9
#define m_to_cm_ 1e2

void MissingMassNtupleAnalyzer::Loop(char * era, char * mode, char * xa, char * jobid, char * outdir, char * datatype, bool createProtonFile, bool randomFlag, bool single, bool zerobias, bool protonsfilter, bool createEventFile)
{

  // debugging variables (printout values)
  bool debug = false;

  TString filenameout;
  TString filenameout_random_reco;
  TString filenameout_eventlist;
  TString filenameout_eventlist_kinematics;

  // saving TTree
  TFile* fileout; 
  TTree* tout;

  if (outdir!=NULL) {
    gSystem->MakeDirectory(outdir);
    if(jobid!=NULL){
      filenameout = std::string(outdir)+"/Ntuple_data_mode_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+"_"+std::string(jobid)+".root";
      filenameout_random_reco = std::string(outdir)+"/proton_reco_rphorizontal_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+"_"+std::string(jobid)+".txt";
      filenameout_eventlist = std::string(outdir)+"/event_list_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+"_"+std::string(jobid)+".txt";
      filenameout_eventlist_kinematics = std::string(outdir)+"/event_list_kinematics_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+"_"+std::string(jobid)+".txt";
    }else{
      filenameout = std::string(outdir)+"/Ntuple_data_mode_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+".root";
      filenameout_random_reco = std::string(outdir)+"/proton_reco_rphorizontal_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+".txt";
      filenameout_eventlist = std::string(outdir)+"/event_list_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+".txt";
      filenameout_eventlist_kinematics = std::string(outdir)+"/event_list_kinematics_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+".txt";
    }
  }else{
    if(jobid!=NULL){
      filenameout = "Ntuple_data_mode_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+"_"+std::string(jobid)+".root";
      filenameout_random_reco = "proton_reco_rphorizontal_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+"_"+std::string(jobid)+".txt";
      filenameout_eventlist = "event_list_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+"_"+std::string(jobid)+".txt";
      filenameout_eventlist_kinematics = "event_list_kinematics_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+"_"+std::string(jobid)+".txt";
    }else{
      filenameout = "Ntuple_data_mode_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+".root";
      filenameout_random_reco = "proton_reco_rphorizontal_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+".txt";
      filenameout_eventlist = "event_list_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+".txt";
      filenameout_eventlist_kinematics = "event_list_kinematics_"+std::string(mode)+"_era_"+std::string(era)+"_xa_"+std::string(xa)+".txt";
    }
  }

  fileout = new TFile(filenameout, "RECREATE"); 

  std::cout << "\t = Options =" << std::endl;
  std::cout << "\t\t Create Proton File: " << createProtonFile << std::endl;
  std::cout << "\t\t Create Event List File: " << createEventFile << std::endl;
  std::cout << "\t\t Random Protons: " << randomFlag << std::endl;
  std::cout << "\t\t Random Proton in a single arm: " << single << std::endl;
  std::cout << "\t\t Mode: " << mode << std::endl;
  std::cout << "\t\t Type: " << datatype << std::endl;
  std::cout << "\t\t ZeroBias: " << zerobias << std::endl;
  std::cout << "\t\t Protons Filter: " << protonsfilter << std::endl;
  std::cout << "\t\t Era: " << era << std::endl;
  std::cout << "\t\t X-Angle: " << xa << std::endl;
  if(jobid!=0) std::cout << "\t\t jobid: " << jobid << std::endl;
  std::cout << "\t\t output: " << filenameout << "\n" << std::endl;

  fileout->cd();

  int xangle_;
  int run_;
  int event_;
  int lumiblock_;
  int era_;
  int PUInterac_;
  int PUTrueInterac_;
  bool is45PUproton_ = false;
  bool is56PUproton_ = false;
  bool is2PUproton_ = false;

  int prescalesL1ZeroBias_;
  int prescalesL1ZeroBiasAfterTrain_;
  int prescalesL1ZeroBiasIsolatedBx_;
  int prescalesL1ZeroBiasAlignment_;
  int prescalesL1ZeroBiasBeamSpot_;
  int prescalesL1ZeroBiasUnpairedBptxMinus_;
  int prescalesL1ZeroBiasUnpairedBptxPlus_;
  int prescalesL1Physics_;
  int prescalesL1SingleMu_;

  bool triggerZeroBias = false;
  bool triggerZeroBiasAfterTrain = false;
  bool triggerZeroBiasIsolatedBx = false;
  bool triggerZeroBiasAlignment = false;
  bool triggerZeroBiasBeamSpot = false;
  bool triggerZeroBiasUnpairedBptxMinus = false;
  bool triggerZeroBiasUnpairedBptxPlus = false;
  bool triggerPhysics = false;
  bool triggerL1SingleMu = false;

  bool triggerIsoMu27 = false;
  bool triggerMu17TrkIsoMu8TrkIso = false;
  bool triggerMu17TrkIsoMu8TrkIsoMass8 = false;
  bool triggerMu17TrkIsoMu8TrkIsoMass3 = false;
  bool triggerDoubleMu43 = false;
  bool triggerDoubleMu48 = false;

  bool triggerEle27 = false;
  bool triggerEle23Ele12 = false;
  bool triggerEle23Ele12Dz = false;
  bool triggerDoubleEle33 = false;

  bool triggerMu23TrkIsoEle12 = false;
  bool triggerMu23TrkIsoEle12DZ = false;
  bool triggerMu8TrkIsoEle23DZ = false;

  bool triggerBTagMu5Ak4dijet20 = false;
  bool triggerBTagMu5Ak4dijet40 = false; 
  bool triggerBTagMu5Ak4dijet70 = false; 
  bool triggerBTagMu5Ak4dijet110 = false; 
  bool triggerBTagMu5Ak4dijet170 = false; 
  bool triggerBTagMu5Ak4dijet300 = false; 
  bool triggerBTagMu5Ak8dijet170 = false;
  bool triggerBTagMu5Ak8dijet300 = false;
  bool triggerPFHT380DoubleBTag = false;
  bool triggerPFHT430DoubleBTag = false;

  bool triggerPFMET100BTag = false;
  bool triggerPFMET110BTag = false;
  bool triggerPFMET120BTag = false;
  bool triggerPFMET130BTag = false;
  bool triggerPFMET140BTag = false;
  bool triggerEle15PFHT450BTag = false;
  bool triggerHT250BTagScouting = false;

  int nVertex = 0;
  int nTracksPerVertex = 0;
  double pvVertexX = 0.;
  double pvVertexY = 0.;
  double pvVertexZ = 0.;
  double pvVertexR = 0.;
  double pvVertexPhi = 0.;
  std::vector<double> allVertexZ;
  std::vector<std::pair<double, double> > MinimumDistance;

  int nGenLeptons = 0;
  int nGenParticles_ = 0;
  int nGenElectrons_ = 0;
  int nGenMuons_ = 0;
  int genleadingLeptonPDGId = 0;
  double genleadingLeptonEnergy = 0.;
  double genleadingLeptonPx = 0.;
  double genleadingLeptonPy = 0.;
  double genleadingLeptonPz = 0.;
  double genleadingLeptonPt = 0.;
  double genleadingLeptonEta = 0.;
  double genleadingLeptonPhi = 0.;
  int genleadingLeptonCharge = 0;
  double genleadingLeptonVx = 0.;
  double genleadingLeptonVy = 0.;
  double genleadingLeptonVz = 0.;
  double genleadingLeptonVr = 0.;
  double genleadingLeptonVphi = 0.;
  int gensecondLeptonPDGId = 0.;
  double gensecondLeptonEnergy = 0.;
  double gensecondLeptonPx = 0.;
  double gensecondLeptonPy = 0.;
  double gensecondLeptonPz = 0.;
  double gensecondLeptonPt = 0.;
  double gensecondLeptonEta = 0.;
  double gensecondLeptonPhi = 0.;
  int gensecondLeptonCharge = 0;
  double gensecondLeptonVx = 0.;
  double gensecondLeptonVy = 0.;
  double gensecondLeptonVz = 0.;
  double gensecondLeptonVr = 0.;
  double gensecondLeptonVphi = 0.;

  double genmissEt_ = 0.;
  double genmissEt_phi_ = 0.;

  int nGenJets = 0;
  double genleadingJetEnergy = 0.;
  double genleadingJetPx = 0.;
  double genleadingJetPy = 0.;
  double genleadingJetPz = 0.;
  double genleadingJetPt = 0.;
  double genleadingJetEta = 0.;
  double genleadingJetPhi = 0.;
  double genleadingJetVz = 0.;
  double gensecondJetEnergy = 0.;
  double gensecondJetPx = 0.;
  double gensecondJetPy = 0.;
  double gensecondJetPz = 0.;
  double gensecondJetPt = 0.;
  double gensecondJetEta = 0.;
  double gensecondJetPhi = 0.;
  double gensecondJetVz = 0.;

  int gendileptonCharge = 0;
  double gendileptonMass = 0.;
  double gendileptonEta = 0.;
  double gendileptonPhi = 0.;
  double gendileptonPt = 0.;
  double gendileptonRapidity = 0.;

  double gendijetMass = 0.;
  double gendijetEta = 0.;
  double gendijetPhi = 0.;
  double gendijetPt = 0.;
  double gendijetRapidity = 0.;

  int nLeptons = 0;
  int nElectrons_ = 0;
  int nMuons_ = 0;
  int nElectrons_looseId_ = 0;
  int nElectrons_mediumId_ = 0;
  int nElectrons_tightId_ = 0;
  int nMuons_looseId_ = 0;
  int nMuons_mediumId_ = 0;
  int nMuons_tightId_ = 0;
  int nChargedPFMultiPV_Loose_ = 0;
  int nChargedPFMultiPV_Tight_ = 0;
  int nChargedPFMultiPV_UsedInFit_ = 0;
  int nChargedPFMultiPV_Tight_Fit_ = 0;
  double SumChargedPFMultiPV_pt_Loose_ = 0;
  double SumChargedPFMultiPV_pt_Tight_ = 0;
  double SumChargedPFMultiPV_pt_UsedInFit_ = 0;
  double SumChargedPFMultiPV_pt_Tight_Fit_ = 0;
  int leadingLeptonPDGId = 0;
  double leadingLeptonEnergy = 0.;
  double leadingLeptonPx = 0.;
  double leadingLeptonPy = 0.;
  double leadingLeptonPz = 0.;
  double leadingLeptonPt = 0.;
  double leadingLeptonEta = 0.;
  double leadingLeptonPhi = 0.;
  int leadingLeptonCharge = 0;
  double leadingLeptonVx = 0.;
  double leadingLeptonVy = 0.;
  double leadingLeptonVz = 0.;
  double leadingLeptonVr = 0.;
  double leadingLeptonVphi = 0.;
  double leadingLeptonPFIso = 0.;
  double leadingLeptonTkIso = 0.;
  bool leadingLeptonLooseId = false;
  bool leadingLeptonMediumId = false;
  bool leadingLeptonTightId = false;
  bool leadingLeptonPfIsoMedium = false;
  bool leadingLeptonMiniIsoTight = false;
  bool leadingLeptonPfIsoVeryTight = false;
  int secondLeptonPDGId = 0;
  double secondLeptonEnergy = 0.;
  double secondLeptonPx = 0.;
  double secondLeptonPy = 0.;
  double secondLeptonPz = 0.;
  double secondLeptonPt = 0.;
  double secondLeptonEta = 0.;
  double secondLeptonPhi = 0.;
  int secondLeptonCharge = 0;
  double secondLeptonVx = 0.;
  double secondLeptonVy = 0.;
  double secondLeptonVz = 0.;
  double secondLeptonVr = 0.;
  double secondLeptonVphi = 0.;
  double secondLeptonPFIso = 0.;
  double secondLeptonTkIso = 0.;
  bool secondLeptonLooseId = false;
  bool secondLeptonMediumId = false;
  bool secondLeptonTightId = false;
  bool secondLeptonPfIsoMedium = false;
  bool secondLeptonMiniIsoTight = false;
  bool secondLeptonPfIsoVeryTight = false;
  double acoplanarity = 0.;

  int nJets = 0;
  int nJetsCandidatesLoose = 0;
  int nJetsCandidatesMedium = 0;
  int nJetsCandidatesTight = 0;
  double leadingJetEnergy = 0.;
  double leadingJetPx = 0.;
  double leadingJetPy = 0.;
  double leadingJetPz = 0.;
  double leadingJetPt = 0.;
  double leadingJetEta = 0.;
  double leadingJetPhi = 0.;
  double leadingJetVx = 0.;
  double leadingJetVy = 0.;
  double leadingJetVz = 0.;
  double leadingJetVr = 0.;
  double leadingJetVphi = 0.;
  double leadingJetBtag = 0.;
  double leadingJetQGdis = 0.;
  double leadingJetNeutralEmFrac = 0.;
  double leadingJetNeutralHadFrac = 0.;
  double leadingJetChargedEmFrac = 0.;
  double leadingJetChargedHadFrac = 0.;
  double leadingJetMuonFrac = 0.;
  int leadingJetNeutralMulti = 0;
  int leadingJetChargedMulti = 0;
  double leadingJetpuIdfdisc = 0.;
  int leadingJetpuIdcbased = 0;
  int leadingJetpuIdfid = 0;
  bool leadingJetLooseId = false;
  bool leadingJetTightId = false;
  bool leadingJetLepVeto = false;
  double Rmpf = 0;
  double JZBalance = 0;
  double JZBalance2 = 0;

  double secondJetEnergy = 0.;
  double secondJetPx = 0.;
  double secondJetPy = 0.;
  double secondJetPz = 0.;
  double secondJetPt = 0.;
  double secondJetEta = 0.;
  double secondJetPhi = 0.;
  double secondJetVx = 0.;
  double secondJetVy = 0.;
  double secondJetVz = 0.;
  double secondJetVr = 0.;
  double secondJetVphi = 0.;
  double secondJetBtag = 0.;
  double secondJetQGdis = 0.;
  double secondJetNeutralEmFrac = 0.;
  double secondJetNeutralHadFrac = 0.;
  double secondJetChargedEmFrac = 0.;
  double secondJetChargedHadFrac = 0.;
  double secondJetMuonFrac = 0.;
  int secondJetNeutralMulti = 0;
  int secondJetChargedMulti = 0;
  double secondJetpuIdfdisc = 0.;
  int secondJetpuIdcbased = 0;
  int secondJetpuIdfid = 0;
  bool secondJetLooseId = false;
  bool secondJetTightId = false;
  bool secondJetLepVeto = false;

  int leptonmetCharge = 0;
  double leptonmetMassT = 0.;
  double leptonmetEta = 0.;
  double leptonmetPhi = 0.;
  double leptonmetPt = 0.;
  double leptonmetRapidity = 0.;

  int dileptonCharge = 0;
  double dileptonMass = 0.;
  double dileptonEta = 0.;
  double dileptonPhi = 0.;
  double dileptonPt = 0.;
  double dileptonRapidity = 0.;

  double dijetMass = 0.;
  double dijetEta = 0.;
  double dijetPhi = 0.;
  double dijetPt = 0.;
  double dijetRapidity = 0.;

  int leptonsystemCharge = 0;
  double leptonsystemMass = 0.;
  double leptonsystemEta = 0.;
  double leptonsystemPhi = 0.;
  double leptonsystemPt = 0.;
  double leptonsystemRapidity = 0.;

  double jetsystemMass = 0.;
  double jetsystemEta = 0.;
  double jetsystemPhi = 0.;
  double jetsystemPt = 0.;
  double jetsystemRapidity = 0.;

  double jetcandidatesystemMass = 0.;
  double jetcandidatesystemEta = 0.;
  double jetcandidatesystemPhi = 0.;
  double jetcandidatesystemPt = 0.;
  double jetcandidatesystemRapidity = 0.;

  double missingMassDijetRP210 = 0.;
  double missingEtaDijetRP210 = 0.;
  double missingPhiDijetRP210 = 0.;
  double missingPtDijetRP210 = 0.;
  double missingRapidityDijetRP210 = 0.;

  double missingMassDijetRP220 = 0.;
  double missingEtaDijetRP220 = 0.;
  double missingPhiDijetRP220 = 0.;
  double missingPtDijetRP220 = 0.;
  double missingRapidityDijetRP220 = 0.;

  double missingMassDijetMulti = 0.;
  double missingEtaDijetMulti = 0.;
  double missingPhiDijetMulti = 0.;
  double missingPtDijetMulti = 0.;
  double missingRapidityDijetMulti = 0.;

  double missingMassDileptonRP210 = 0.;
  double missingEtaDileptonRP210 = 0.;
  double missingPhiDileptonRP210 = 0.;
  double missingPtDileptonRP210 = 0.;
  double missingRapidityDileptonRP210 = 0.;

  double missingMassDileptonRP220 = 0.;
  double missingEtaDileptonRP220 = 0.;
  double missingPhiDileptonRP220 = 0.;
  double missingPtDileptonRP220 = 0.;
  double missingRapidityDileptonRP220 = 0.;

  double missingMassDileptonMulti = 0.;
  double missingEtaDileptonMulti = 0.;
  double missingPhiDileptonMulti = 0.;
  double missingPtDileptonMulti = 0.;
  double missingRapidityDileptonMulti = 0.;

  double missingMassDileptonJet1RP220 = 0.;
  double missingEtaDileptonJet1RP220 = 0.;
  double missingPhiDileptonJet1RP220 = 0.;
  double missingPtDileptonJet1RP220 = 0.;
  double missingRapidityDileptonJet1RP220 = 0.;

  double missingMassDileptonJet1Multi = 0.;
  double missingEtaDileptonJet1Multi = 0.;
  double missingPhiDileptonJet1Multi = 0.;
  double missingPtDileptonJet1Multi = 0.;
  double missingRapidityDileptonJet1Multi = 0.;

  double missingMassDileptonDijetRP210 = 0.;
  double missingEtaDileptonDijetRP210 = 0.;
  double missingPhiDileptonDijetRP210 = 0.;
  double missingPtDileptonDijetRP210 = 0.;
  double missingRapidityDileptonDijetRP210 = 0.;

  double missingMassDileptonDijetRP220 = 0.;
  double missingEtaDileptonDijetRP220 = 0.;
  double missingPhiDileptonDijetRP220 = 0.;
  double missingPtDileptonDijetRP220 = 0.;
  double missingRapidityDileptonDijetRP220 = 0.;

  double missingMassDileptonDijetMulti = 0.;
  double missingEtaDileptonDijetMulti = 0.;
  double missingPhiDileptonDijetMulti = 0.;
  double missingPtDileptonDijetMulti = 0.;
  double missingRapidityDileptonDijetMulti = 0.;

  double missEt_ = 0.;
  double missEt_phi_ = 0.;
  double dphi_ = 0.;

  double diffMassRP210 = -1.;
  double diffMassRP220 = -1.;
  double diffMassMulti = -1.;

  double proton_pz_rp210Arm45 = 0.;
  double proton_pz_rp210Arm56 = 0.;
  double proton_pz_rp220Arm45 = 0.;
  double proton_pz_rp220Arm56 = 0.;
  double proton_pz_multiArm45 = 0.;
  double proton_pz_multiArm56 = 0.;

  int nprotonRP210_sec45 = 0;
  int nprotonRP210_sec56 = 0;
  int nprotonRP220_sec45 = 0;
  int nprotonRP220_sec56 = 0;
  int nprotonMulti_sec45 = 0;
  int nprotonMulti_sec56 = 0;
  int nprotonRP210_sec45_Reduced = 0;
  int nprotonRP210_sec56_Reduced = 0;
  int nprotonRP220_sec45_Reduced = 0;
  int nprotonRP220_sec56_Reduced = 0;

  double xi_rp210_Arm45 = 0.;
  double t_rp210_Arm45 = 0.;
  double xi_rp210_Arm56 = 0.;
  double t_rp210_Arm56 = 0.;
  double xi_rp220_Arm45 = 0.;
  double t_rp220_Arm45 = 0.;
  double xi_rp220_Arm56 = 0.;
  double t_rp220_Arm56 = 0.;

  double xi_multiArm45 = 0.;
  double time_multiArm45 = 0.;
  double timeunc_multiArm45 = 0.;
  double thx_multiArm45 = 0.;
  double thy_multiArm45 = 0.;
  double xi_multiArm56 = 0.;
  double time_multiArm56 = 0.;
  double timeunc_multiArm56 = 0.;
  double thx_multiArm56 = 0.;
  double thy_multiArm56 = 0.;

  double t45andt56 = 0.;
  double vzpps = 0.;

  double distanceVertexZPPSandCMS = 0.;
  double selectedVertexZ = 0.;

  bool isprotonRP210 = false;
  bool isprotonRP220 = false;
  bool isprotonRP210_Reduced = false;
  bool isprotonRP220_Reduced = false;
  bool isprotonMulti = false;

  tout = new TTree("Events", "Events");
  tout->Branch("era",&era_,"era/I");
  tout->Branch("run",&run_,"run/I");
  tout->Branch("event",&event_,"event/I");
  tout->Branch("lumiblock",&lumiblock,"lumiblock/I");
  tout->Branch("xangle",&xangle_,"xangle/I");
  if(zerobias){
    tout->Branch("triggerZeroBias",&triggerZeroBias,"triggerZeroBias/B");
    tout->Branch("triggerZeroBiasAfterTrain",&triggerZeroBiasAfterTrain,"triggerZeroBiasAfterTrain/B");
    tout->Branch("triggerZeroBiasIsolatedBx",&triggerZeroBiasIsolatedBx,"triggerZeroBiasIsolatedBx/B");
    tout->Branch("triggerZeroBiasAlignment",&triggerZeroBiasAlignment,"triggerZeroBiasAlignment/B");
    tout->Branch("triggerZeroBiasBeamSpot",&triggerZeroBiasBeamSpot,"triggerZeroBiasBeamSpot/B");
    tout->Branch("triggerZeroBiasUnpairedBptxMinus",&triggerZeroBiasUnpairedBptxMinus,"triggerZeroBiasUnpairedBptxMinus/B");
    tout->Branch("triggerZeroBiasUnpairedBptxPlus",&triggerZeroBiasUnpairedBptxPlus,"triggerZeroBiasUnpairedBptxPlus/B");
    tout->Branch("triggerPhysics",&triggerPhysics,"triggerPhysics/B");
    tout->Branch("triggerL1SingleMu",&triggerL1SingleMu,"triggerL1SingleMu/B");
    tout->Branch("prescalesL1ZeroBias",&prescalesL1ZeroBias_,"prescalesL1ZeroBias/I");
    tout->Branch("prescalesL1ZeroBiasAfterTrain",&prescalesL1ZeroBiasAfterTrain_,"prescalesL1ZeroBiasAfterTrain/I");
    tout->Branch("prescalesL1ZeroBiasIsolatedBx",&prescalesL1ZeroBiasIsolatedBx_,"prescalesL1ZeroBiasIsolatedBx/I");
    tout->Branch("prescalesL1ZeroBiasAlignment",&prescalesL1ZeroBiasAlignment_,"prescalesL1ZeroBiasAlignment/I");
    tout->Branch("prescalesL1ZeroBiasBeamSpot",&prescalesL1ZeroBiasBeamSpot_,"prescalesL1ZeroBeamSpot/I");
    tout->Branch("prescalesL1ZeroBiasUnpairedBptxMinus",&prescalesL1ZeroBiasUnpairedBptxMinus_,"prescalesL1ZeroBeamUnpairedBptxMinus/I");
    tout->Branch("prescalesL1ZeroBiasUnpairedBptxPlus",&prescalesL1ZeroBiasUnpairedBptxPlus_,"prescalesL1ZeroBeamUnpairedBptxPlus/I");
    tout->Branch("prescalesL1Physics",&prescalesL1Physics_,"prescalesL1Physics/I");
    tout->Branch("prescalesL1SingleMu",&prescalesL1SingleMu_,"prescalesL1SingleMu/I");
  }else{
    //if((strcmp(mode, "Muon")==0 || strcmp(mode, "MC_Muon")==0) && strcmp(datatype, "SD")!=0 && strcmp(datatype, "Signal")!=0){
    tout->Branch("triggerIsoMu27",&triggerIsoMu27,"triggerIsoMu27/B");
    tout->Branch("triggerMu17TrkIsoMu8TrkIso",&triggerMu17TrkIsoMu8TrkIso,"triggerMu17TrkIsoMu8TrkIso/B");
    tout->Branch("triggerMu17TrkIsoMu8TrkIsoMass3",&triggerMu17TrkIsoMu8TrkIsoMass3,"triggerMu17TrkIsoMu8TrkIsoMass3/B");
    tout->Branch("triggerMu17TrkIsoMu8TrkIsoMass8",&triggerMu17TrkIsoMu8TrkIsoMass8,"triggerMu17TrkIsoMu8TrkIsoMass8/B");
    tout->Branch("triggerDoubleMu43",&triggerDoubleMu43,"triggerDoubleMu43/B");
    tout->Branch("triggerDoubleMu48",&triggerDoubleMu48,"triggerDoubleMu48/B");
    //}
    //else if(strcmp(mode, "Electron")==0 || strcmp(mode, "MC_Electron")==0 && strcmp(datatype, "SD")!=0 && strcmp(datatype, "Signal")!=0){
    tout->Branch("triggerEle27",&triggerEle27,"triggerEle27/B");
    tout->Branch("triggerEle23Ele12",&triggerEle23Ele12,"triggerEle23Ele12/B");
    tout->Branch("triggerEle23Ele12Dz",&triggerEle23Ele12Dz,"triggerEle23Ele12Dz/B");
    tout->Branch("triggerDoubleEle33",&triggerDoubleEle33,"triggerDoubleEle33/B");
    tout->Branch("triggerBTagMu5Ak4dijet20",&triggerBTagMu5Ak4dijet20,"triggerBTagMu5Ak4dijet20/B");
    tout->Branch("triggerBTagMu5Ak4dijet40",&triggerBTagMu5Ak4dijet40,"triggerBTagMu5Ak4dijet40/B");
    tout->Branch("triggerBTagMu5Ak4dijet70",&triggerBTagMu5Ak4dijet70,"triggerBTagMu5Ak4dijet70/B");
    tout->Branch("triggerBTagMu5Ak4dijet110",&triggerBTagMu5Ak4dijet110,"triggerBTagMu5Ak4dijet110/B");
    tout->Branch("triggerBTagMu5Ak4dijet170",&triggerBTagMu5Ak4dijet170,"triggerBTagMu5Ak4dijet170/B");
    tout->Branch("triggerBTagMu5Ak4dijet300",&triggerBTagMu5Ak4dijet300,"triggerBTagMu5Ak4dijet300/B");
    tout->Branch("triggerBTagMu5Ak8dijet170",&triggerBTagMu5Ak8dijet170,"triggerBTagMu5Ak8dijet170/B");
    tout->Branch("triggerBTagMu5Ak8dijet300",&triggerBTagMu5Ak8dijet300,"triggerBTagMu5Ak8dijet300/B");
    tout->Branch("triggerPFHT380DoubleBTag",&triggerPFHT380DoubleBTag,"triggerPFHT380DoubleBTag/B");
    tout->Branch("triggerPFHT430DoubleBTag",&triggerPFHT430DoubleBTag,"triggerPFHT430DoubleBTag/B");
    tout->Branch("triggerPFMET100BTag",&triggerPFMET100BTag,"triggerPFMET100BTag/B");
    tout->Branch("triggerPFMET110BTag",&triggerPFMET110BTag,"triggerPFMET110BTag/B");
    tout->Branch("triggerPFMET120BTag",&triggerPFMET120BTag,"triggerPFMET120BTag/B");
    tout->Branch("triggerPFMET130BTag",&triggerPFMET130BTag,"triggerPFMET130BTag/B");
    tout->Branch("triggerPFMET140BTag",&triggerPFMET140BTag,"triggerPFMET140BTag/B");
    tout->Branch("triggerEle15PFHT450BTag",&triggerEle15PFHT450BTag,"triggerEle15PFHT450BTag/B");
    tout->Branch("triggerHT250BTagScouting",&triggerHT250BTagScouting,"triggerHT250BTagScouting/B");
    //}
    //else if(strcmp(mode, "Bjets")==0 || strcmp(mode, "MC_Bjets")==0 && strcmp(datatype, "SD")!=0 && strcmp(datatype, "Signal")!=0){
    tout->Branch("triggerIsoMu27",&triggerIsoMu27,"triggerIsoMu27/B");
    tout->Branch("triggerMu17TrkIsoMu8TrkIso",&triggerMu17TrkIsoMu8TrkIso,"triggerMu17TrkIsoMu8TrkIso/B");
    tout->Branch("triggerMu17TrkIsoMu8TrkIsoMass3",&triggerMu17TrkIsoMu8TrkIsoMass3,"triggerMu17TrkIsoMu8TrkIsoMass3/B");
    tout->Branch("triggerMu17TrkIsoMu8TrkIsoMass8",&triggerMu17TrkIsoMu8TrkIsoMass8,"triggerMu17TrkIsoMu8TrkIsoMass8/B");
    tout->Branch("triggerDoubleMu43",&triggerDoubleMu43,"triggerDoubleMu43/B");
    tout->Branch("triggerDoubleMu48",&triggerDoubleMu48,"triggerDoubleMu48/B");
    tout->Branch("triggerBTagMu5Ak4dijet20",&triggerBTagMu5Ak4dijet20,"triggerBTagMu5Ak4dijet20/B");
    tout->Branch("triggerBTagMu5Ak4dijet40",&triggerBTagMu5Ak4dijet40,"triggerBTagMu5Ak4dijet40/B");
    tout->Branch("triggerBTagMu5Ak4dijet70",&triggerBTagMu5Ak4dijet70,"triggerBTagMu5Ak4dijet70/B");
    tout->Branch("triggerBTagMu5Ak4dijet110",&triggerBTagMu5Ak4dijet110,"triggerBTagMu5Ak4dijet110/B");
    tout->Branch("triggerBTagMu5Ak4dijet170",&triggerBTagMu5Ak4dijet170,"triggerBTagMu5Ak4dijet170/B");
    tout->Branch("triggerBTagMu5Ak4dijet300",&triggerBTagMu5Ak4dijet300,"triggerBTagMu5Ak4dijet300/B");
    tout->Branch("triggerBTagMu5Ak8dijet170",&triggerBTagMu5Ak8dijet170,"triggerBTagMu5Ak8dijet170/B");
    tout->Branch("triggerBTagMu5Ak8dijet300",&triggerBTagMu5Ak8dijet300,"triggerBTagMu5Ak8dijet300/B");
    tout->Branch("triggerPFHT380DoubleBTag",&triggerPFHT380DoubleBTag,"triggerPFHT380DoubleBTag/B");
    tout->Branch("triggerPFHT430DoubleBTag",&triggerPFHT430DoubleBTag,"triggerPFHT430DoubleBTag/B");
    //}
    //else if(strcmp(mode, "EMu")==0){
    //  tout->Branch("triggerMu23TrkIsoEle12",&triggerMu23TrkIsoEle12,"triggerMu23TrkIsoEle12/B");
    //  tout->Branch("triggerMu23TrkIsoEle12DZ",&triggerMu23TrkIsoEle12DZ,"triggerMu23TrkIsoEle12DZ/B");
    //  tout->Branch("triggerMu8TrkIsoEle23DZ",&triggerMu8TrkIsoEle23DZ,"triggerMu8TrkIsoEle23DZ/B");
    //}
  }

  // In miniAOD this is empty (slimmed Vertex do not have the track ref.)
  tout->Branch("nTracksPerVertex",&nTracksPerVertex,"nTracksPerVertex/I");
  tout->Branch("nVertex",&nVertex,"nVertex/I");
  tout->Branch("pvVertexX",&pvVertexX,"pvVertexX/D");
  tout->Branch("pvVertexY",&pvVertexY,"pvVertexY/D");
  tout->Branch("pvVertexZ",&pvVertexZ,"pvVertexZ/D");
  tout->Branch("pvVertexR",&pvVertexR,"pvVertexR/D");
  tout->Branch("pvVertexPhi",&pvVertexPhi,"pvVertexPhi/D");
  tout->Branch("allVertexZ",&allVertexZ);
  if (strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0){
    tout->Branch("PUInterac",&PUInterac,"PUInterac/I");
    tout->Branch("PUTrueInterac",&PUTrueInterac,"PUTrueInterac/I");
    tout->Branch("is45PUproton",&is45PUproton_,"is45PUproton/B");
    tout->Branch("is56PUproton",&is56PUproton_,"is56PUproton/B");
    tout->Branch("is2PUproton",&is2PUproton_,"is2PUproton/B");
    tout->Branch("nGenLeptons",&nGenLeptons,"nGenLeptons/I");
    tout->Branch("nGenParticles",&nGenParticles,"nGenParticles/I");
    tout->Branch("nGenElectrons",&nGenElectrons,"nGenElectrons/I");
    tout->Branch("nGenMuons",&nGenMuons,"nGenMuons/I");
    tout->Branch("genleadingLeptonPDGId",&genleadingLeptonPDGId,"genleadingLeptonPDGId/I");
    tout->Branch("genleadingLeptonEnergy",&genleadingLeptonEnergy,"genleadingLeptonEnergy/D");
    tout->Branch("genleadingLeptonPx",&genleadingLeptonPx,"genleadingLeptonPx/D");
    tout->Branch("genleadingLeptonPy",&genleadingLeptonPy,"genleadingLeptonPy/D");
    tout->Branch("genleadingLeptonPz",&genleadingLeptonPz,"genleadingLeptonPz/D");
    tout->Branch("genleadingLeptonPt",&genleadingLeptonPt,"genleadingLeptonPt/D");
    tout->Branch("genleadingLeptonEta",&genleadingLeptonEta,"genleadingLeptonEta/D");
    tout->Branch("genleadingLeptonPhi",&genleadingLeptonPhi,"genleadingLeptonPhi/D");
    tout->Branch("genleadingLeptonCharge",&genleadingLeptonCharge,"genleadingLeptonCharge/I");
    tout->Branch("genleadingLeptonVx",&genleadingLeptonVx,"genleadingLeptonVx/D");
    tout->Branch("genleadingLeptonVy",&genleadingLeptonVy,"genleadingLeptonVy/D");
    tout->Branch("genleadingLeptonVz",&genleadingLeptonVz,"genleadingLeptonVz/D");
    tout->Branch("genleadingLeptonVr",&genleadingLeptonVr,"genleadingLeptonVr/D");
    tout->Branch("genleadingLeptonVphi",&genleadingLeptonVphi,"genleadingLeptonVphi/D");
    tout->Branch("gensecondLeptonPDGId",&gensecondLeptonPDGId,"gensecondLeptonPDGId/I");
    tout->Branch("gensecondLeptonEnergy",&gensecondLeptonEnergy,"gensecondLeptonEnergy/D");
    tout->Branch("gensecondLeptonPx",&gensecondLeptonPx,"gensecondLeptonPx/D");
    tout->Branch("gensecondLeptonPy",&gensecondLeptonPy,"gensecondLeptonPy/D");
    tout->Branch("gensecondLeptonPz",&gensecondLeptonPz,"gensecondLeptonPz/D");
    tout->Branch("gensecondLeptonPt",&gensecondLeptonPt,"gensecondLeptonPt/D");
    tout->Branch("gensecondLeptonEta",&gensecondLeptonEta,"gensecondLeptonEta/D");
    tout->Branch("gensecondLeptonPhi",&gensecondLeptonPhi,"gensecondLeptonPhi/D");
    tout->Branch("gensecondLeptonCharge",&gensecondLeptonCharge,"gensecondLeptonCharge/I");
    tout->Branch("gensecondLeptonVx",&gensecondLeptonVx,"gensecondLeptonVx/D");
    tout->Branch("gensecondLeptonVy",&gensecondLeptonVy,"gensecondLeptonVy/D");
    tout->Branch("gensecondLeptonVz",&gensecondLeptonVz,"gensecondLeptonVz/D");
    tout->Branch("gensecondLeptonVr",&gensecondLeptonVr,"gensecondLeptonVr/D");
    tout->Branch("gensecondLeptonVphi",&gensecondLeptonVphi,"gensecondLeptonVphi/D");
    tout->Branch("genmissEt",&genmissEt_,"genmissEt/D");
    tout->Branch("genmissEt_phi",&genmissEt_phi_,"genmissEt_phi/D");
    tout->Branch("nGenJets",&nGenJets,"nGenJets/I");
    tout->Branch("genleadingJetEnergy",&genleadingJetEnergy,"genleadingJetEnergy/D");
    tout->Branch("genleadingJetPx",&genleadingJetPx,"genleadingJetPx/D");
    tout->Branch("genleadingJetPy",&genleadingJetPy,"genleadingJetPy/D");
    tout->Branch("genleadingJetPz",&genleadingJetPz,"genleadingJetPz/D");
    tout->Branch("genleadingJetPt",&genleadingJetPt,"genleadingJetPt/D");
    tout->Branch("genleadingJetEta",&genleadingJetEta,"genleadingJetEta/D");
    tout->Branch("genleadingJetPhi",&genleadingJetPhi,"genleadingJetPhi/D");
    tout->Branch("genleadingJetVz",&genleadingJetVz,"genleadingJetVz/D");
    tout->Branch("gensecondJetEnergy",&gensecondJetEnergy,"gensecondJetEnergy/D");
    tout->Branch("gensecondJetPx",&gensecondJetPx,"gensecondJetPx/D");
    tout->Branch("gensecondJetPy",&gensecondJetPy,"gensecondJetPy/D");
    tout->Branch("gensecondJetPz",&gensecondJetPz,"gensecondJetPz/D");
    tout->Branch("gensecondJetPt",&gensecondJetPt,"gensecondJetPt/D");
    tout->Branch("gensecondJetEta",&gensecondJetEta,"gensecondJetEta/D");
    tout->Branch("gensecondJetPhi",&gensecondJetPhi,"gensecondJetPhi/D");
    tout->Branch("gensecondJetVz",&gensecondJetVz,"gensecondJetVz/D");
    tout->Branch("gendileptonCharge",&gendileptonCharge,"gendileptonCharge/I");
    tout->Branch("gendileptonMass",&gendileptonMass,"gendileptonMass/D");
    tout->Branch("gendileptonEta",&gendileptonEta,"gendileptonEta/D");
    tout->Branch("gendileptonPhi",&gendileptonPhi,"gendileptonPhi/D");
    tout->Branch("gendileptonPt",&gendileptonPt,"gendileptonPt/D");
    tout->Branch("gendileptonRapidity",&gendileptonRapidity,"gendileptonRapidity/D");
    tout->Branch("gendijetMass",&gendijetMass,"gendijetMass/D");
    tout->Branch("gendijetEta",&gendijetEta,"gendijetEta/D");
    tout->Branch("gendijetPhi",&gendijetPhi,"gendijetPhi/D");
    tout->Branch("gendijetPt",&gendijetPt,"gendijetPt/D");
    tout->Branch("gendijetRapidity",&gendijetRapidity,"gendijetRapidity/D");
  }
  tout->Branch("nLeptons",&nLeptons,"nLeptons/I");
  tout->Branch("nElectrons",&nElectrons_,"nElectrons/I");
  tout->Branch("nElectronsLooseId",&nElectrons_looseId_,"nElectronsLooseId/I");
  tout->Branch("nElectronsMediumId",&nElectrons_mediumId_,"nElectronsMediumId/I");
  tout->Branch("nElectronsTightId",&nElectrons_tightId_,"nElectronsTightId/I");
  tout->Branch("nMuons",&nMuons_,"nMuons/I");
  tout->Branch("nMuonsLooseId",&nMuons_looseId_,"nMuonsLooseId/I");
  tout->Branch("nMuonsMediumId",&nMuons_mediumId_,"nMuonsMediumId/I");
  tout->Branch("nMuonsTightId",&nMuons_tightId_,"nMuonsTightId/I");
  tout->Branch("nChargedPFMultiPV_Loose",&nChargedPFMultiPV_Loose_,"nChargedPFMultiPV_Loose/I");
  tout->Branch("nChargedPFMultiPV_Tight",&nChargedPFMultiPV_Tight_,"nChargedPFMultiPV_Tight/I");
  tout->Branch("nChargedPFMultiPV_UsedInFit",&nChargedPFMultiPV_UsedInFit_,"nChargedPFMultiPV_UsedInFit/I");
  tout->Branch("nChargedPFMultiPV_Tight_Fit",&nChargedPFMultiPV_Tight_Fit_,"nChargedPFMultiPV_Tight_Fit/I");
  tout->Branch("SumChargedPFMultiPV_pt_Loose",&SumChargedPFMultiPV_pt_Loose_,"SumChargedPFMultiPV_pt_Loose/D");
  tout->Branch("SumChargedPFMultiPV_pt_Tight",&SumChargedPFMultiPV_pt_Tight_,"SumChargedPFMultiPV_pt_Tight/D");
  tout->Branch("SumChargedPFMultiPV_pt_UsedInFit",&SumChargedPFMultiPV_pt_UsedInFit_,"SumChargedPFMultiPV_pt_UsedInFit/D");
  tout->Branch("SumChargedPFMultiPV_pt_Tight_Fit",&SumChargedPFMultiPV_pt_Tight_Fit_,"SumChargedPFMultiPV_pt_Tight_Fit/D");
  tout->Branch("leadingLeptonPDGId",&leadingLeptonPDGId,"leadingLeptonPDGId/I");
  tout->Branch("leadingLeptonEnergy",&leadingLeptonEnergy,"leadingLeptonEnergy/D");
  tout->Branch("leadingLeptonPx",&leadingLeptonPx,"leadingLeptonPx/D");
  tout->Branch("leadingLeptonPy",&leadingLeptonPy,"leadingLeptonPy/D");
  tout->Branch("leadingLeptonPz",&leadingLeptonPz,"leadingLeptonPz/D");
  tout->Branch("leadingLeptonPt",&leadingLeptonPt,"leadingLeptonPt/D");
  tout->Branch("leadingLeptonEta",&leadingLeptonEta,"leadingLeptonEta/D");
  tout->Branch("leadingLeptonPhi",&leadingLeptonPhi,"leadingLeptonPhi/D");
  tout->Branch("leadingLeptonCharge",&leadingLeptonCharge,"leadingLeptonCharge/I");
  tout->Branch("leadingLeptonVx",&leadingLeptonVx,"leadingLeptonVx/D");
  tout->Branch("leadingLeptonVy",&leadingLeptonVy,"leadingLeptonVy/D");
  tout->Branch("leadingLeptonVz",&leadingLeptonVz,"leadingLeptonVz/D");
  tout->Branch("leadingLeptonVr",&leadingLeptonVr,"leadingLeptonVr/D");
  tout->Branch("leadingLeptonVphi",&leadingLeptonVphi,"leadingLeptonVphi/D");
  tout->Branch("leadingLeptonPFIso",&leadingLeptonPFIso,"leadingLeptonPFIso/D");
  tout->Branch("leadingLeptonTkIso",&leadingLeptonTkIso,"leadingLeptonTkIso/D");
  tout->Branch("leadingLeptonLooseId",&leadingLeptonLooseId,"leadingLeptonLooseId/B");
  tout->Branch("leadingLeptonMediumId",&leadingLeptonMediumId,"leadingLeptonMediumId/B");
  tout->Branch("leadingLeptonTightId",&leadingLeptonTightId,"leadingLeptonTightId/B");
  tout->Branch("leadingLeptonPfIsoMedium",&leadingLeptonPfIsoMedium,"leadingLeptonPfIsoMedium/B");
  tout->Branch("leadingLeptonMiniIsoTight",&leadingLeptonMiniIsoTight,"leadingLeptonMiniIsoTight/B");
  tout->Branch("leadingLeptonPfIsoVeryTight",&leadingLeptonPfIsoVeryTight,"leadingLeptonPfIsoVeryTight/B");
  tout->Branch("secondLeptonPDGId",&secondLeptonPDGId,"secondLeptonPDGId/I");
  tout->Branch("secondLeptonEnergy",&secondLeptonEnergy,"secondLeptonEnergy/D");
  tout->Branch("secondLeptonPx",&secondLeptonPx,"secondLeptonPx/D");
  tout->Branch("secondLeptonPy",&secondLeptonPy,"secondLeptonPy/D");
  tout->Branch("secondLeptonPz",&secondLeptonPz,"secondLeptonPz/D");
  tout->Branch("secondLeptonPt",&secondLeptonPt,"secondLeptonPt/D");
  tout->Branch("secondLeptonEta",&secondLeptonEta,"secondLeptonEta/D");
  tout->Branch("secondLeptonPhi",&secondLeptonPhi,"secondLeptonPhi/D");
  tout->Branch("secondLeptonCharge",&secondLeptonCharge,"secondLeptonCharge/I");
  tout->Branch("secondLeptonVx",&secondLeptonVx,"secondLeptonVx/D");
  tout->Branch("secondLeptonVy",&secondLeptonVy,"secondLeptonVy/D");
  tout->Branch("secondLeptonVz",&secondLeptonVz,"secondLeptonVz/D");
  tout->Branch("secondLeptonVr",&secondLeptonVr,"secondLeptonVr/D");
  tout->Branch("secondLeptonVphi",&secondLeptonVphi,"secondLeptonVphi/D");
  tout->Branch("secondLeptonPFIso",&secondLeptonPFIso,"secondLeptonPFIso/D");
  tout->Branch("secondLeptonTkIso",&secondLeptonTkIso,"secondLeptonTkIso/D");
  tout->Branch("secondLeptonLooseId",&secondLeptonLooseId,"secondLeptonLooseId/B");
  tout->Branch("secondLeptonMediumId",&secondLeptonMediumId,"secondLeptonMediumId/B");
  tout->Branch("secondLeptonTightId",&secondLeptonTightId,"secondLeptonTightId/B");
  tout->Branch("secondLeptonPfIsoMedium",&secondLeptonPfIsoMedium,"secondLeptonPfIsoMedium/B");
  tout->Branch("secondLeptonMiniIsoTight",&secondLeptonMiniIsoTight,"secondLeptonMiniIsoTight/B");
  tout->Branch("secondLeptonPfIsoVeryTight",&secondLeptonPfIsoVeryTight,"secondLeptonPfIsoVeryTight/B");
  tout->Branch("acoplanarity",&acoplanarity,"acoplanarity/D");
  tout->Branch("nJets",&nJets,"nJets/I");
  tout->Branch("leadingJetEnergy",&leadingJetEnergy,"leadingJetEnergy/D");
  tout->Branch("leadingJetPx",&leadingJetPx,"leadingJetPx/D");
  tout->Branch("leadingJetPy",&leadingJetPy,"leadingJetPy/D");
  tout->Branch("leadingJetPz",&leadingJetPz,"leadingJetPz/D");
  tout->Branch("leadingJetPt",&leadingJetPt,"leadingJetPt/D");
  tout->Branch("leadingJetEta",&leadingJetEta,"leadingJetEta/D");
  tout->Branch("leadingJetPhi",&leadingJetPhi,"leadingJetPhi/D");
  tout->Branch("leadingJetVx",&leadingJetVx,"leadingJetVx/D");
  tout->Branch("leadingJetVy",&leadingJetVy,"leadingJetVy/D");
  tout->Branch("leadingJetVz",&leadingJetVz,"leadingJetVz/D");
  tout->Branch("leadingJetVr",&leadingJetVr,"leadingJetVr/D");
  tout->Branch("leadingJetVphi",&leadingJetVphi,"leadingJetVphi/D");
  tout->Branch("leadingJetBtag",&leadingJetBtag,"leadingJetBtag/D");
  tout->Branch("leadingJetQGdis",&leadingJetQGdis,"leadingJetQGdis/D");
  tout->Branch("leadingJetNeutralEmFrac",&leadingJetNeutralEmFrac,"leadingJetNeutralEmFrac/D");
  tout->Branch("leadingJetNeutralHadFrac",&leadingJetNeutralHadFrac,"leadingJetNeutralHadFrac/D");
  tout->Branch("leadingJetChargedEmFrac",&leadingJetChargedEmFrac,"leadingJetChargedEmFrac/D");
  tout->Branch("leadingJetChargedHadFrac",&leadingJetChargedHadFrac,"leadingJetChargedHadFrac/D");
  tout->Branch("leadingJetMuonFrac",&leadingJetMuonFrac,"leadingJetMuonFrac/D");
  tout->Branch("leadingJetNeutralMulti",&leadingJetNeutralMulti,"leadingJetNeutralMulti/I");
  tout->Branch("leadingJetChargedMulti",&leadingJetChargedMulti,"leadingJetChargedMulti/I");
  tout->Branch("leadingJetpuIdfdisc",&leadingJetpuIdfdisc,"leadingJetpuIdfdisc/D");
  tout->Branch("leadingJetpuIdcbased",&leadingJetpuIdcbased,"leadingJetpuIdcbased/I");
  tout->Branch("leadingJetpuIdfid",&leadingJetpuIdfid,"leadingJetpuIdfid/I");
  tout->Branch("leadingJetLooseId",&leadingJetLooseId,"leadingJetLooseId/B");
  tout->Branch("leadingJetTightId",&leadingJetTightId,"leadingJetTightId/B");
  tout->Branch("leadingJetLepVeto",&leadingJetLepVeto,"leadingJetLepVeto/B");
  tout->Branch("secondJetEnergy",&secondJetEnergy,"secondJetEnergy/D");
  tout->Branch("secondJetPx",&secondJetPx,"secondJetPx/D");
  tout->Branch("secondJetPy",&secondJetPy,"secondJetPy/D");
  tout->Branch("secondJetPz",&secondJetPz,"secondJetPz/D");
  tout->Branch("secondJetPt",&secondJetPt,"secondJetPt/D");
  tout->Branch("secondJetEta",&secondJetEta,"secondJetEta/D");
  tout->Branch("secondJetPhi",&secondJetPhi,"secondJetPhi/D");
  tout->Branch("secondJetVx",&secondJetVx,"secondJetVx/D");
  tout->Branch("secondJetVy",&secondJetVy,"secondJetVy/D");
  tout->Branch("secondJetVz",&secondJetVz,"secondJetVz/D");
  tout->Branch("secondJetVr",&secondJetVr,"secondJetVr/D");
  tout->Branch("secondJetVphi",&secondJetVphi,"secondJetVphi/D");
  tout->Branch("secondJetBtag",&secondJetBtag,"secondJetBtag/D");
  tout->Branch("secondJetQGdis",&secondJetQGdis,"secondJetQGdis/D");
  tout->Branch("secondJetNeutralEmFrac",&secondJetNeutralEmFrac,"secondJetNeutralEmFrac/D");
  tout->Branch("secondJetNeutralHadFrac",&secondJetNeutralHadFrac,"secondJetNeutralHadFrac/D");
  tout->Branch("secondJetChargedEmFrac",&secondJetChargedEmFrac,"secondJetChargedEmFrac/D");
  tout->Branch("secondJetChargedHadFrac",&secondJetChargedHadFrac,"secondJetChargedHadFrac/D");
  tout->Branch("secondJetMuonFrac",&secondJetMuonFrac,"secondJetMuonFrac/D");
  tout->Branch("secondJetNeutralMulti",&secondJetNeutralMulti,"secondJetNeutralMulti/I");
  tout->Branch("secondJetChargedMulti",&secondJetChargedMulti,"secondJetChargedMulti/I");
  tout->Branch("secondJetpuIdfdisc",&secondJetpuIdfdisc,"secondJetpuIdfdisc/D");
  tout->Branch("secondJetpuIdcbased",&secondJetpuIdcbased,"secondJetpuIdcbased/I");
  tout->Branch("secondJetpuIdfid",&secondJetpuIdfid,"secondJetpuIdfid/I");
  tout->Branch("secondJetLooseId",&secondJetLooseId,"secondJetLooseId/B");
  tout->Branch("secondJetTightId",&secondJetTightId,"secondJetTightId/B");
  tout->Branch("secondJetLepVeto",&secondJetLepVeto,"secondJetLepVeto/B");
  tout->Branch("leptonmetCharge",&leptonmetCharge,"leptonmetCharge/I");
  tout->Branch("leptonmetMassT",&leptonmetMassT,"leptonmetMassT/D");
  tout->Branch("leptonmetEta",&leptonmetEta,"leptonmetEta/D");
  tout->Branch("leptonmetPhi",&leptonmetPhi,"leptonmetPhi/D");
  tout->Branch("leptonmetPt",&leptonmetPt,"leptonmetPt/D");
  tout->Branch("leptonmetRapidity",&leptonmetRapidity,"leptonmetRapidity/D");
  tout->Branch("dileptonCharge",&dileptonCharge,"dileptonCharge/I");
  tout->Branch("dileptonMass",&dileptonMass,"dileptonMass/D");
  tout->Branch("dileptonEta",&dileptonEta,"dileptonEta/D");
  tout->Branch("dileptonPhi",&dileptonPhi,"dileptonPhi/D");
  tout->Branch("dileptonPt",&dileptonPt,"dileptonPt/D");
  tout->Branch("dileptonRapidity",&dileptonRapidity,"dileptonRapidity/D");
  tout->Branch("dijetMass",&dijetMass,"dijetMass/D");
  tout->Branch("dijetEta",&dijetEta,"dijetEta/D");
  tout->Branch("dijetPhi",&dijetPhi,"dijetPhi/D");
  tout->Branch("dijetPt",&dijetPt,"dijetPt/D");
  tout->Branch("dijetRapidity",&dijetRapidity,"dijetRapidity/D");
  tout->Branch("leptonsystemCharge",&leptonsystemCharge,"leptonsystemCharge/I");
  tout->Branch("leptonsystemMass",&leptonsystemMass,"leptonsystemMass/D");
  tout->Branch("leptonsystemEta",&leptonsystemEta,"leptonsystemEta/D");
  tout->Branch("leptonsystemPhi",&leptonsystemPhi,"leptonsystemPhi/D");
  tout->Branch("leptonsystemPt",&leptonsystemPt,"leptonsystemPt/D");
  tout->Branch("leptonsystemRapidity",&leptonsystemRapidity,"leptonsystemRapidity/D");
  tout->Branch("jetsystemMass",&jetsystemMass,"jetsystemMass/D");
  tout->Branch("jetsystemEta",&jetsystemEta,"jetsystemEta/D");
  tout->Branch("jetsystemPhi",&jetsystemPhi,"jetsystemPhi/D");
  tout->Branch("jetsystemPt",&jetsystemPt,"jetsystemPt/D");
  tout->Branch("jetsystemRapidity",&jetsystemRapidity,"jetsystemRapidity/D");
  tout->Branch("nJetsCandidatesLoose",&nJetsCandidatesLoose,"nJetsCandidatesLoose/I");
  tout->Branch("nJetsCandidatesMedium",&nJetsCandidatesMedium,"nJetsCandidatesMedium/I");
  tout->Branch("nJetsCandidatesTight",&nJetsCandidatesTight,"nJetsCandidatesTight/I");
  tout->Branch("Rmpf",&Rmpf,"Rmpf/D");
  tout->Branch("JZBalance",&JZBalance,"JZBalance/D");
  tout->Branch("JZBalance2",&JZBalance2,"JZBalance2/D");
  tout->Branch("jetcandidatesystemMass",&jetcandidatesystemMass,"jetcandidatesystemMass/D");
  tout->Branch("jetcandidatesystemEta",&jetcandidatesystemEta,"jetcandidatesystemEta/D");
  tout->Branch("jetcandidatesystemPhi",&jetcandidatesystemPhi,"jetcandidatesystemPhi/D");
  tout->Branch("jetcandidatesystemPt",&jetcandidatesystemPt,"jetcandidatesystemPt/D");
  tout->Branch("jetcandidatesystemRapidity",&jetcandidatesystemRapidity,"jetcandidatesystemRapidity/D");
  tout->Branch("missingMassDijetRP210",&missingMassDijetRP210,"missingMassDijetRP210/D");
  tout->Branch("missingEtaDijetRP210",&missingEtaDijetRP210,"missingEtaDijetRP210/D");
  tout->Branch("missingPhiDijetRP210",&missingPhiDijetRP210,"missingPhiDijetRP210/D");
  tout->Branch("missingPtDijetRP210",&missingPtDijetRP210,"missingPtDijetRP210/D");
  tout->Branch("missingRapidityDijetRP210",&missingRapidityDijetRP210,"missingRapidityDijetRP210/D");
  tout->Branch("missingMassDijetRP220",&missingMassDijetRP220,"missingMassDijetRP220/D");
  tout->Branch("missingEtaDijetRP220",&missingEtaDijetRP220,"missingEtaDijetRP220/D");
  tout->Branch("missingPhiDijetRP220",&missingPhiDijetRP220,"missingPhiDijetRP220/D");
  tout->Branch("missingPtDijetRP220",&missingPtDijetRP220,"missingPtDijetRP220/D");
  tout->Branch("missingRapidityDijetRP220",&missingRapidityDijetRP220,"missingRapidityDijetRP220/D");
  tout->Branch("missingMassDijetMulti",&missingMassDijetMulti,"missingMassDijetMulti/D");
  tout->Branch("missingEtaDijetMulti",&missingEtaDijetMulti,"missingEtaDijetMulti/D");
  tout->Branch("missingPhiDijetMulti",&missingPhiDijetMulti,"missingPhiDijetMulti/D");
  tout->Branch("missingPtDijetMulti",&missingPtDijetMulti,"missingPtDijetMulti/D");
  tout->Branch("missingRapidityDijetMulti",&missingRapidityDijetMulti,"missingRapidityDijetMulti/D");
  tout->Branch("missingMassDileptonRP210",&missingMassDileptonRP210,"missingMassDileptonRP210/D");
  tout->Branch("missingEtaDileptonRP210",&missingEtaDileptonRP210,"missingEtaDileptonRP210/D");
  tout->Branch("missingPhiDileptonRP210",&missingPhiDileptonRP210,"missingPhiDileptonRP210/D");
  tout->Branch("missingPtDileptonRP210",&missingPtDileptonRP210,"missingPtDileptonRP210/D");
  tout->Branch("missingRapidityDileptonRP210",&missingRapidityDileptonRP210,"missingRapidityDileptonRP210/D");
  tout->Branch("missingMassDileptonRP220",&missingMassDileptonRP220,"missingMassDileptonRP220/D");
  tout->Branch("missingEtaDileptonRP220",&missingEtaDileptonRP220,"missingEtaDileptonRP220/D");
  tout->Branch("missingPhiDileptonRP220",&missingPhiDileptonRP220,"missingPhiDileptonRP220/D");
  tout->Branch("missingPtDileptonRP220",&missingPtDileptonRP220,"missingPtDileptonRP220/D");
  tout->Branch("missingRapidityDileptonRP220",&missingRapidityDileptonRP220,"missingRapidityDileptonRP220/D");
  tout->Branch("missingMassDileptonMulti",&missingMassDileptonMulti,"missingMassDileptonMulti/D");
  tout->Branch("missingEtaDileptonMulti",&missingEtaDileptonMulti,"missingEtaDileptonMulti/D");
  tout->Branch("missingPhiDileptonMulti",&missingPhiDileptonMulti,"missingPhiDileptonMulti/D");
  tout->Branch("missingPtDileptonMulti",&missingPtDileptonMulti,"missingPtDileptonMulti/D");
  tout->Branch("missingRapidityDileptonMulti",&missingRapidityDileptonMulti,"missingRapidityDileptonMulti/D");
  tout->Branch("missingMassDileptonJet1RP220",&missingMassDileptonJet1RP220,"missingMassDileptonJet1RP220/D");
  tout->Branch("missingEtaDileptonJet1RP220",&missingEtaDileptonJet1RP220,"missingEtaDileptonJet1RP220/D");
  tout->Branch("missingPhiDileptonJet1RP220",&missingPhiDileptonJet1RP220,"missingPhiDileptonJet1RP220/D");
  tout->Branch("missingPtDileptonJet1RP220",&missingPtDileptonJet1RP220,"missingPtDileptonJet1RP220/D");
  tout->Branch("missingRapidityDileptonJet1RP220",&missingRapidityDileptonJet1RP220,"missingRapidityDileptonJet1RP220/D");
  tout->Branch("missingMassDileptonJet1Multi",&missingMassDileptonJet1Multi,"missingMassDileptonJet1Multi/D");
  tout->Branch("missingEtaDileptonJet1Multi",&missingEtaDileptonJet1Multi,"missingEtaDileptonJet1Multi/D");
  tout->Branch("missingPhiDileptonJet1Multi",&missingPhiDileptonJet1Multi,"missingPhiDileptonJet1Multi/D");
  tout->Branch("missingPtDileptonJet1Multi",&missingPtDileptonJet1Multi,"missingPtDileptonJet1Multi/D");
  tout->Branch("missingRapidityDileptonJet1Multi",&missingRapidityDileptonJet1Multi,"missingRapidityDileptonJet1Multi/D");
  tout->Branch("missingMassDileptonDijetRP210",&missingMassDileptonDijetRP210,"missingMassDileptonDijetRP210/D");
  tout->Branch("missingEtaDileptonDijetRP210",&missingEtaDileptonDijetRP210,"missingEtaDileptonDijetRP210/D");
  tout->Branch("missingPhiDileptonDijetRP210",&missingPhiDileptonDijetRP210,"missingPhiDileptonDijetRP210/D");
  tout->Branch("missingPtDileptonDijetRP210",&missingPtDileptonDijetRP210,"missingPtDileptonDijetRP210/D");
  tout->Branch("missingRapidityDileptonDijetRP210",&missingRapidityDileptonDijetRP210,"missingRapidityDileptonDijetRP210/D");
  tout->Branch("missingMassDileptonDijetRP220",&missingMassDileptonDijetRP220,"missingMassDileptonDijetRP220/D");
  tout->Branch("missingEtaDileptonDijetRP220",&missingEtaDileptonDijetRP220,"missingEtaDileptonDijetRP220/D");
  tout->Branch("missingPhiDileptonDijetRP220",&missingPhiDileptonDijetRP220,"missingPhiDileptonDijetRP220/D");
  tout->Branch("missingPtDileptonDijetRP220",&missingPtDileptonDijetRP220,"missingPtDileptonDijetRP220/D");
  tout->Branch("missingRapidityDileptonDijetRP220",&missingRapidityDileptonDijetRP220,"missingRapidityDileptonDijetRP220/D");
  tout->Branch("missingMassDileptonDijetMulti",&missingMassDileptonDijetMulti,"missingMassDileptonDijetMulti/D");
  tout->Branch("missingEtaDileptonDijetMulti",&missingEtaDileptonDijetMulti,"missingEtaDileptonDijetMulti/D");
  tout->Branch("missingPhiDileptonDijetMulti",&missingPhiDileptonDijetMulti,"missingPhiDileptonDijetMulti/D");
  tout->Branch("missingPtDileptonDijetMulti",&missingPtDileptonDijetMulti,"missingPtDileptonDijetMulti/D");
  tout->Branch("missingRapidityDileptonDijetMulti",&missingRapidityDileptonDijetMulti,"missingRapidityDileptonDijetMulti/D");
  tout->Branch("missEt",&missEt_,"missEt/D");
  tout->Branch("missEt_phi",&missEt_phi_,"missEt_phi/D");
  tout->Branch("dphi",&dphi_,"dphi/D");
  tout->Branch("diffMassRP210",&diffMassRP210,"diffMassRP210/D");
  tout->Branch("diffMassRP220",&diffMassRP220,"diffMassRP220/D");
  tout->Branch("diffMassMulti",&diffMassMulti,"diffMassMulti/D");
  tout->Branch("protonsrp210Arm45Pz",&proton_pz_rp210Arm45,"protonsrp210Arm45Pz/D");
  tout->Branch("protonsrp210Arm56Pz",&proton_pz_rp210Arm56,"protonsrp210Arm56Pz/D");
  tout->Branch("protonsrp220Arm45Pz",&proton_pz_rp220Arm45,"protonsrp220Arm45Pz/D");
  tout->Branch("protonsrp220Arm56Pz",&proton_pz_rp220Arm56,"protonsrp220Arm56Pz/D");
  tout->Branch("protonsmultiArm45Pz",&proton_pz_multiArm45,"protonsmultiArm45Pz/D");
  tout->Branch("protonsmultiArm56Pz",&proton_pz_multiArm56,"protonsmultiArm56Pz/D");
  tout->Branch("nprotonRP210_Arm45",&nprotonRP210_sec45,"nprotonRP210_Arm45/I");
  tout->Branch("nprotonRP210_Arm56",&nprotonRP210_sec56,"nprotonRP210_Arm56/I");
  tout->Branch("nprotonRP220_Arm45",&nprotonRP220_sec45,"nprotonRP220_Arm45/I");
  tout->Branch("nprotonRP220_Arm56",&nprotonRP220_sec56,"nprotonRP220_Arm56/I");
  tout->Branch("nprotonMultiArm45",&nprotonMulti_sec45,"nprotonMultiArm45/I");
  tout->Branch("nprotonMultiArm56",&nprotonMulti_sec56,"nprotonMultiArm56/I");
  tout->Branch("isprotonRP210",&isprotonRP210,"isprotonRP210/B");
  tout->Branch("isprotonRP220",&isprotonRP220,"isprotonRP220/B");
  tout->Branch("isprotonMulti",&isprotonMulti,"isprotonMulti/B");
  tout->Branch("isprotonRP210_Reduced ",&isprotonRP210_Reduced,"isprotonRP210_Reduced/B");
  tout->Branch("isprotonRP220_Reduced ",&isprotonRP220_Reduced,"isprotonRP220_Reduced/B");
  tout->Branch("xi_rp210_Arm45",&xi_rp210_Arm45,"xi_rp210_Arm45/D");
  tout->Branch("xi_rp210_Arm56",&xi_rp210_Arm56,"xi_rp210_Arm56/D");
  tout->Branch("xi_rp220_Arm45",&xi_rp220_Arm45,"xi_rp220_Arm45/D");
  tout->Branch("xi_rp220_Arm56",&xi_rp220_Arm56,"xi_rp220_Arm56/D");
  tout->Branch("xi_multiArm45",&xi_multiArm45,"xi_multiArm45/D");
  tout->Branch("time_multiArm45",&time_multiArm45,"time_multiArm45/D");
  tout->Branch("timeunc_multiArm45",&timeunc_multiArm45,"timeunc_multiArm45/D");
  tout->Branch("thx_multiArm45",&thx_multiArm45,"thx_multiArm45/D");
  tout->Branch("thy_multiArm45",&thy_multiArm45,"thy_multiArm45/D");
  tout->Branch("xi_multiArm56",&xi_multiArm56,"xi_multiArm56/D");
  tout->Branch("time_multiArm56",&time_multiArm56,"time_multiArm56/D");
  tout->Branch("timeunc_multiArm56",&timeunc_multiArm56,"timeunc_multiArm56/D");
  tout->Branch("thx_multiArm56",&thx_multiArm56,"thx_multiArm56/D");
  tout->Branch("thy_multiArm56",&thy_multiArm56,"thy_multiArm56/D");
  tout->Branch("t45andt56",&t45andt56,"t45andt56/D");
  tout->Branch("vzpps",&vzpps,"vzpps/D");
  tout->Branch("distanceVertexZPPSandCMS",&distanceVertexZPPSandCMS,"distanceVertexZPPSandCMS/D");
  tout->Branch("selectedVertexZ",&selectedVertexZ,"selectedVertexZ/D");

  std::ofstream file_protons_reco;
  if(createProtonFile){
    file_protons_reco.open (filenameout_random_reco);
  }

  std::ofstream file_eventlist;
  std::ofstream file_eventlist_kinematics;
  if(createEventFile){
    file_eventlist.open (filenameout_eventlist);
    file_eventlist_kinematics.open (filenameout_eventlist_kinematics);
  }

  if (fChain == 0) return;

  TRandom3 rand(0); 
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  Long64_t rangemin = 0;
  Long64_t rangemax = nentries;

  if(strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0){

    Double_t frac[20];
    if(strcmp(mode, "MC_Muon")==0){ //Fractions for 1 Multi-RP and 1 pixels
      frac[0]  = 0.009941;   //0.0135092;       
      frac[1]  = 0.011229;   //0.0150113;
      frac[2]  = 0.016891;   //0.0222738;
      frac[3]  = 0.018763;   //0.0244954;
      frac[4]  = 0.065784;   //0.0847967;
      frac[5]  = 0.097749;   //0.119802;
      frac[6]  = 0.064967;   //0.0747702;
      frac[7]  = 0.093342;   //0.103186;
      frac[8]  = 0.019649;   //0.021513;
      frac[9]  = 0.041326;   //0.0417492;
      frac[10] = 0.041221;   //0.0389873;
      frac[11] = 0.042031;   //0.0375865;
      frac[12] = 0.023791;   //0.0297421;
      frac[13] = 0.019781;   //0.023547;
      frac[14] = 0.000317;   //0.000372487;
      frac[15] = 0.008085;   //0.00945544;
      frac[16] = 0.185997;   //0.174744;
      frac[17] = 0.098481;   //0.0764684;
      frac[18] = 0.013582;   //0.0125681;
      frac[19] = 0.104090;   //0.0754216;   
    }
    else if (strcmp(mode, "MC_Electron")==0){
      frac[0]  = 0.015610;   //0.0204228;     
      frac[1]  = 0.022636;   //0.0292415;
      frac[2]  = 0.031372;   //0.0399177;
      frac[3]  = 0.039083;   //0.0491031;
      frac[4]  = 0.072814;   //0.090098;
      frac[5]  = 0.119773;   //0.140086;
      frac[6]  = 0.090694;   //0.0986203;
      frac[7]  = 0.144755;   //0.151347;
      frac[8]  = 0.013952;   //0.014786;
      frac[9]  = 0.031176;   //0.0300992;
      frac[10] = 0.032689;   //0.0298066;
      frac[11] = 0.038521;   //0.0299309;
      frac[12] = 0.017496;   //0.0211467;
      frac[13] = 0.015451;   //0.0176158;
      frac[14] = 0.000225;   //0.000250771;
      frac[15] = 0.005792;   //0.00657242;
      frac[16] = 0.131171;   //0.117324;
      frac[17] = 0.072942;   //0.0537374;
      frac[18] = 0.009866;   //0.00885834;
      frac[19] = 0.074057;   //0.0510357;   
    }
    if(randomFlag && strcmp(datatype, "DY")==0){
      if(strcmp(era, "B")==0 && strcmp(xa, "120")==0){
	rangemin=0;
	rangemax=floor(nentries*frac[0]);  
      }
      else if(strcmp(era, "B")==0 && strcmp(xa, "130")==0){
	rangemin=floor(nentries*frac[0]);
	rangemax=rangemin + floor(nentries*frac[1]);
      }
      else if(strcmp(era, "B")==0 && strcmp(xa, "140")==0){
	rangemin=floor(nentries*(frac[0] + frac[1]));
	rangemax=rangemin + floor(nentries*frac[2]);
      }
      else if(strcmp(era, "B")==0 && strcmp(xa, "150")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2]));
	rangemax=rangemin + floor(nentries*frac[3]);
      }
      else if(strcmp(era, "C")==0 && strcmp(xa, "120")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3]));
	rangemax=rangemin + floor(nentries*frac[4]);
      }
      else if(strcmp(era, "C")==0 && strcmp(xa, "130")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4]));
	rangemax=rangemin + floor(nentries*frac[5]);
      }
      else if(strcmp(era, "C")==0 && strcmp(xa, "140")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5]));
	rangemax=rangemin + floor(nentries*frac[6]);
      }
      else if(strcmp(era, "C")==0 && strcmp(xa, "150")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6]));
	rangemax=rangemin + floor(nentries*frac[7]);
      }
      else if(strcmp(era, "D")==0 && strcmp(xa, "120")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7]));
	rangemax=rangemin + floor(nentries*frac[8]);    
      }
      else if(strcmp(era, "D")==0 && strcmp(xa, "130")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8]));
	rangemax=rangemin + floor(nentries*frac[9]);
      }
      else if(strcmp(era, "D")==0 && strcmp(xa, "140")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9]));
	rangemax=rangemin + floor(nentries*frac[10]);
      }
      else if(strcmp(era, "D")==0 && strcmp(xa, "150")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9] + frac[10]));
	rangemax=rangemin + floor(nentries*frac[11]);
      }
      else if(strcmp(era, "E")==0 && strcmp(xa, "120")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9] + frac[10] + frac[11]));
	rangemax=rangemin + floor(nentries*frac[12]);
      }
      else if(strcmp(era, "E")==0 && strcmp(xa, "130")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9] + frac[10] + frac[11] + frac[12]));
	rangemax=rangemin + floor(nentries*frac[13]);
      }
      else if(strcmp(era, "E")==0 && strcmp(xa, "140")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9] + frac[10] + frac[11] + frac[12] + frac[13]));
	rangemax=rangemin + floor(nentries*frac[14]);
	std::cout << rangemin << std::endl;
	std::cout << rangemax << std::endl;
      }
      else if(strcmp(era, "E")==0 && strcmp(xa, "150")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9] + frac[10] + frac[11] + frac[12] + frac[13] + frac[14]));
	rangemax=rangemin + floor(nentries*frac[15]);
      }
      else if(strcmp(era, "F")==0 && strcmp(xa, "120")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9] + frac[10] + frac[11] + frac[12] + frac[13] + frac[14] + frac[15]));
	rangemax=rangemin + floor(nentries*frac[16]);
      }
      else if(strcmp(era, "F")==0 && strcmp(xa, "130")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9] + frac[10] + frac[11] + frac[12] + frac[13] + frac[14] + frac[15] + frac[16]));
	rangemax=rangemin + floor(nentries*frac[17]);
      }
      else if(strcmp(era, "F")==0 && strcmp(xa, "140")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9] + frac[10] + frac[11] + frac[12] + frac[13] + frac[14] + frac[15] + frac[16] + frac[17]));
	rangemax=rangemin + floor(nentries*frac[18]);
      }
      else if(strcmp(era, "F")==0 && strcmp(xa, "150")==0){
	rangemin=floor(nentries*(frac[0] + frac[1] + frac[2] + frac[3] + frac[4] + frac[5] + frac[6] + frac[7] + frac[8] + frac[9] + frac[10] + frac[11] + frac[12] + frac[13] + frac[14] + frac[15] + frac[16] + frac[17] + frac[18]));
	rangemax=rangemin + floor(nentries*frac[19]);
      }
      else{
	std::cout << "Era or X-angle not present in 2017 dataset! " << std::endl;
      }
    }
  }

  for (Long64_t jentry=rangemin; jentry<rangemax;jentry++) {

    xangle_ = 0;
    run_ = 0;
    event_ = 0;
    lumiblock_ = 0;
    era_ = -1;

    prescalesL1ZeroBias_ = 1;
    prescalesL1ZeroBiasAfterTrain_ = 1;
    prescalesL1ZeroBiasIsolatedBx_ = 1;
    prescalesL1ZeroBiasAlignment_ = 1;
    prescalesL1ZeroBiasBeamSpot_ = 1;
    prescalesL1ZeroBiasUnpairedBptxMinus_ = 1;
    prescalesL1ZeroBiasUnpairedBptxPlus_ = 1;
    prescalesL1Physics_ = 1;
    prescalesL1SingleMu_ = 1;

    triggerZeroBias = false;
    triggerZeroBiasAfterTrain = false;
    triggerZeroBiasIsolatedBx = false;
    triggerZeroBiasAlignment = false;
    triggerZeroBiasBeamSpot = false;
    triggerZeroBiasUnpairedBptxMinus = false;
    triggerZeroBiasUnpairedBptxPlus = false;
    triggerPhysics = false;
    triggerL1SingleMu = false;

    triggerIsoMu27 = false;
    triggerMu17TrkIsoMu8TrkIso = false;
    triggerMu17TrkIsoMu8TrkIsoMass8 = false;
    triggerMu17TrkIsoMu8TrkIsoMass3 = false;
    triggerDoubleMu43 = false;
    triggerDoubleMu48 = false;

    triggerEle27 = false;
    triggerEle23Ele12 = false;
    triggerEle23Ele12Dz = false;
    triggerDoubleEle33 = false;

    triggerBTagMu5Ak4dijet20 = false;
    triggerBTagMu5Ak4dijet40 = false;
    triggerBTagMu5Ak4dijet70 = false;
    triggerBTagMu5Ak4dijet110 = false;
    triggerBTagMu5Ak4dijet170 = false;
    triggerBTagMu5Ak4dijet300 = false;
    triggerBTagMu5Ak8dijet170 = false;
    triggerBTagMu5Ak8dijet300 = false;
    triggerPFHT380DoubleBTag = false;
    triggerPFHT430DoubleBTag = false;

    triggerPFMET100BTag = false;
    triggerPFMET110BTag = false;
    triggerPFMET120BTag = false;
    triggerPFMET130BTag = false;
    triggerPFMET140BTag = false;
    triggerEle15PFHT450BTag = false;
    triggerHT250BTagScouting = false;

    triggerMu23TrkIsoEle12 = false;
    triggerMu23TrkIsoEle12DZ = false;
    triggerMu8TrkIsoEle23DZ = false;

    nVertex = 0;
    nTracksPerVertex = 0;
    pvVertexX = 0.;
    pvVertexY = 0.;
    pvVertexZ = 0.;
    pvVertexR = 0.;
    pvVertexPhi = 0.;

    allVertexZ.clear();
    MinimumDistance.clear();

    if(strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0){   
      PUInterac_ = 0;
      PUTrueInterac_ = 0;
      is45PUproton_ = false;
      is56PUproton_ = false;
      is2PUproton_ = false;
      nGenLeptons = 0;
      nGenParticles = 0;
      nGenElectrons = 0;
      nGenMuons = 0;
      genleadingLeptonPDGId = 0;
      genleadingLeptonEnergy = 0.;
      genleadingLeptonPx = 0.;
      genleadingLeptonPy = 0.;
      genleadingLeptonPz = 0.;
      genleadingLeptonPt = 0.;
      genleadingLeptonEta = 0.;
      genleadingLeptonPhi = 0.;
      genleadingLeptonCharge = 0.;
      genleadingLeptonVx = 0.;
      genleadingLeptonVy = 0.;
      genleadingLeptonVz = 0.;
      genleadingLeptonVr = 0.;
      genleadingLeptonVphi = 0.;

      gensecondLeptonPDGId = 0;
      gensecondLeptonEnergy = 0.;
      gensecondLeptonPx = 0.;
      gensecondLeptonPy = 0.;
      gensecondLeptonPz = 0.;
      gensecondLeptonPt = 0.;
      gensecondLeptonEta = 0.;
      gensecondLeptonPhi = 0.;
      gensecondLeptonCharge = 0.;
      gensecondLeptonVx = 0.;
      gensecondLeptonVy = 0.;
      gensecondLeptonVz = 0.;
      gensecondLeptonVr = 0.;
      gensecondLeptonVphi = 0.;

      genmissEt_ = 0.;
      genmissEt_phi_ = 0.;

      nGenJets = 0;

      genleadingJetEnergy = 0.;
      genleadingJetPx = 0.;
      genleadingJetPy = 0.;
      genleadingJetPz = 0.;
      genleadingJetPt = 0.;
      genleadingJetEta = 0.;
      genleadingJetPhi = 0.;
      genleadingJetVz = 0.;

      gensecondJetEnergy = 0.;
      gensecondJetPx = 0.;
      gensecondJetPy = 0.;
      gensecondJetPz = 0.;
      gensecondJetPt = 0.;
      gensecondJetEta = 0.;
      gensecondJetPhi = 0.;
      gensecondJetVz = 0.;

      gendileptonCharge = 0.;
      gendileptonMass = 0.;
      gendileptonEta = 0.;
      gendileptonPhi = 0.;
      gendileptonPt = 0.;
      gendileptonRapidity = 0.;

      gendijetMass = 0.;
      gendijetEta = 0.;
      gendijetPhi = 0.;
      gendijetPt = 0.;
      gendijetRapidity = 0.;
    }

    nLeptons = 0;
    nElectrons = 0;
    nMuons = 0;
    nElectrons_looseId = 0;
    nElectrons_mediumId = 0;
    nElectrons_tightId = 0;
    nMuons_looseId = 0;
    nMuons_mediumId = 0;
    nMuons_tightId = 0;
    nChargedPFMultiPV_Loose = 0;
    nChargedPFMultiPV_Tight = 0;
    nChargedPFMultiPV_UsedInFit = 0;
    nChargedPFMultiPV_Tight_Fit = 0;
    SumChargedPFMultiPV_pt_Loose = 0;
    SumChargedPFMultiPV_pt_Tight = 0;
    SumChargedPFMultiPV_pt_UsedInFit = 0;
    SumChargedPFMultiPV_pt_Tight_Fit = 0;
    leadingLeptonPDGId = 0;
    leadingLeptonEnergy = 0.;
    leadingLeptonPx = 0.;
    leadingLeptonPy = 0.;
    leadingLeptonPz = 0.;
    leadingLeptonPt = 0.;
    leadingLeptonEta = 0.;
    leadingLeptonPhi = 0.;
    leadingLeptonCharge = 0.;
    leadingLeptonVx = 0.;
    leadingLeptonVy = 0.;
    leadingLeptonVz = 0.;
    leadingLeptonVr = 0.;
    leadingLeptonVphi = 0.;
    leadingLeptonPFIso = 0.;
    leadingLeptonTkIso = 0.;
    leadingLeptonLooseId = false;
    leadingLeptonMediumId = false;
    leadingLeptonTightId = false;
    leadingLeptonPfIsoMedium = false;
    leadingLeptonMiniIsoTight = false;
    leadingLeptonPfIsoVeryTight = false;

    secondLeptonPDGId = 0;
    secondLeptonEnergy = 0.;
    secondLeptonPx = 0.;
    secondLeptonPy = 0.;
    secondLeptonPz = 0.;
    secondLeptonPt = 0.;
    secondLeptonEta = 0.;
    secondLeptonPhi = 0.;
    secondLeptonCharge = 0.;
    secondLeptonVx = 0.;
    secondLeptonVy = 0.;
    secondLeptonVz = 0.;
    secondLeptonVr = 0.;
    secondLeptonVphi = 0.;
    secondLeptonPFIso = 0.;
    secondLeptonTkIso = 0.;
    secondLeptonLooseId = false;
    secondLeptonMediumId = false;
    secondLeptonTightId = false;
    secondLeptonPfIsoMedium = false;
    secondLeptonMiniIsoTight = false;
    secondLeptonPfIsoVeryTight = false;

    acoplanarity = 0;

    nJets = 0;
    nJetsCandidatesLoose = 0;
    nJetsCandidatesMedium = 0;
    nJetsCandidatesTight = 0;

    leadingJetEnergy = 0.;
    leadingJetPx = 0.;
    leadingJetPy = 0.;
    leadingJetPz = 0.;
    leadingJetPt = 0.;
    leadingJetEta = 0.;
    leadingJetPhi = 0.;
    leadingJetVx = 0.;
    leadingJetVy = 0.;
    leadingJetVz = 0.;
    leadingJetVr = 0.;
    leadingJetVphi = 0.;
    leadingJetBtag = 0.;
    leadingJetQGdis = 0.;
    leadingJetNeutralEmFrac = 0.;
    leadingJetNeutralHadFrac = 0.;
    leadingJetChargedEmFrac = 0.;
    leadingJetChargedHadFrac = 0.;
    leadingJetMuonFrac = 0.;
    leadingJetNeutralMulti = 0;
    leadingJetChargedMulti = 0;
    leadingJetpuIdfdisc = 0.;
    leadingJetpuIdcbased = 0;
    leadingJetpuIdfid = 0;
    leadingJetLooseId = false;
    leadingJetTightId = false;
    leadingJetLepVeto = false;

    secondJetEnergy = 0.;
    secondJetPx = 0.;
    secondJetPy = 0.;
    secondJetPz = 0.;
    secondJetPt = 0.;
    secondJetEta = 0.;
    secondJetPhi = 0.;
    secondJetVx = 0.;
    secondJetVy = 0.;
    secondJetVz = 0.;
    secondJetVr = 0.;
    secondJetVphi = 0.;
    secondJetBtag = 0.;
    secondJetQGdis = 0.; 
    secondJetNeutralEmFrac = 0.; 
    secondJetNeutralHadFrac = 0.; 
    secondJetChargedEmFrac = 0.; 
    secondJetChargedHadFrac = 0.; 
    secondJetMuonFrac = 0.; 
    secondJetNeutralMulti = 0; 
    secondJetChargedMulti = 0; 
    secondJetpuIdfdisc = 0.; 
    secondJetpuIdcbased = 0; 
    secondJetpuIdfid = 0; 
    secondJetLooseId = false;
    secondJetTightId = false;
    secondJetLepVeto = false;

    Rmpf = 0.;
    JZBalance = 0.;
    JZBalance2 = 0.;

    leptonmetCharge = 0.;
    leptonmetMassT = 0.;
    leptonmetEta = 0.;
    leptonmetPhi = 0.;
    leptonmetPt = 0.;
    leptonmetRapidity = 0.;

    dileptonCharge = 0.;
    dileptonMass = 0.;
    dileptonEta = 0.;
    dileptonPhi = 0.;
    dileptonPt = 0.;
    dileptonRapidity = 0.;

    dijetMass = 0.;
    dijetEta = 0.;
    dijetPhi = 0.;
    dijetPt = 0.;
    dijetRapidity = 0.;

    jetsystemMass = 0.;
    jetsystemEta = 0.;
    jetsystemPhi = 0.;
    jetsystemPt = 0.;
    jetsystemRapidity = 0.;

    jetcandidatesystemMass = 0.;
    jetcandidatesystemEta = 0.;
    jetcandidatesystemPhi = 0.;
    jetcandidatesystemPt = 0.;
    jetcandidatesystemRapidity = 0.;

    missingMassDijetRP210 = 0.;
    missingMassDijetRP210 = 0.;
    missingEtaDijetRP210 = 0.;
    missingPhiDijetRP210 = 0.;
    missingPtDijetRP210 = 0.;
    missingRapidityDijetRP210 = 0.;

    missingMassDijetRP220 = 0.;
    missingMassDijetRP220 = 0.;
    missingEtaDijetRP220 = 0.;
    missingPhiDijetRP220 = 0.;
    missingPtDijetRP220 = 0.;
    missingRapidityDijetRP220 = 0.;

    missingMassDijetMulti = 0.;
    missingMassDijetMulti = 0.;
    missingEtaDijetMulti = 0.;
    missingPhiDijetMulti = 0.;
    missingPtDijetMulti = 0.;
    missingRapidityDijetMulti = 0.;

    missingMassDileptonRP210 = 0.;
    missingMassDileptonRP210 = 0.;
    missingEtaDileptonRP210 = 0.;
    missingPhiDileptonRP210 = 0.;
    missingPtDileptonRP210 = 0.;
    missingRapidityDileptonRP210 = 0.;

    missingMassDileptonRP220 = 0.;
    missingMassDileptonRP220 = 0.;
    missingEtaDileptonRP220 = 0.;
    missingPhiDileptonRP220 = 0.;
    missingPtDileptonRP220 = 0.;
    missingRapidityDileptonRP220 = 0.;

    missingMassDileptonMulti = 0.;
    missingMassDileptonMulti = 0.;
    missingEtaDileptonMulti = 0.;
    missingPhiDileptonMulti = 0.;
    missingPtDileptonMulti = 0.;
    missingRapidityDileptonMulti = 0.;

    missingMassDileptonJet1RP220 = 0.;
    missingMassDileptonJet1RP220 = 0.;
    missingEtaDileptonJet1RP220 = 0.;
    missingPhiDileptonJet1RP220 = 0.;
    missingPtDileptonJet1RP220 = 0.;
    missingRapidityDileptonJet1RP220 = 0.;

    missingMassDileptonJet1Multi = 0.;
    missingMassDileptonJet1Multi = 0.;
    missingEtaDileptonJet1Multi = 0.;
    missingPhiDileptonJet1Multi = 0.;
    missingPtDileptonJet1Multi = 0.;
    missingRapidityDileptonJet1Multi = 0.;

    missingMassDileptonDijetRP210 = 0.;
    missingMassDileptonDijetRP210 = 0.;
    missingEtaDileptonDijetRP210 = 0.;
    missingPhiDileptonDijetRP210 = 0.;
    missingPtDileptonDijetRP210 = 0.;
    missingRapidityDileptonDijetRP210 = 0.;

    missingMassDileptonDijetRP220 = 0.;
    missingMassDileptonDijetRP220 = 0.;
    missingEtaDileptonDijetRP220 = 0.;
    missingPhiDileptonDijetRP220 = 0.;
    missingPtDileptonDijetRP220 = 0.;
    missingRapidityDileptonDijetRP220 = 0.;

    missingMassDileptonDijetMulti = 0.;
    missingMassDileptonDijetMulti = 0.;
    missingEtaDileptonDijetMulti = 0.;
    missingPhiDileptonDijetMulti = 0.;
    missingPtDileptonDijetMulti = 0.;
    missingRapidityDileptonDijetMulti = 0.;

    missEt_ = 0.;
    missEt_phi_ = 0.;
    dphi_ = 0.;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    if(!debug) loadBar(jentry, nentries);
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //std::cout << "nbytes: \t" << nbytes << std::endl;

    int random_idx;
    if(randomFlag){
      random_idx = rand.Integer(rndXi_RP210_sec45.size()); // rndTrack45, rndTrack56 or rndMulti have the same size.
      //std::cout << "random idx: \t" << random_idx << std::endl;
    }

    xi_rp210_Arm45 = 0;
    xi_rp210_Arm56 = 0;
    xi_rp220_Arm45 = 0;
    xi_rp220_Arm56 = 0;

    xi_multiArm45 = 0;
    time_multiArm45 = 0;
    timeunc_multiArm45 = 0;
    thx_multiArm45 = 0;
    thy_multiArm45 = 0;
    xi_multiArm56 = 0;
    time_multiArm56 = 0;
    timeunc_multiArm56 = 0;
    thx_multiArm56 = 0;
    thy_multiArm56 = 0;

    isprotonRP210 = false;
    isprotonRP220 = false;
    isprotonMulti = false;
    isprotonRP210_Reduced = false;
    isprotonRP220_Reduced = false;

    nprotonRP210_sec45 = 0;
    nprotonRP210_sec56 = 0;
    nprotonRP220_sec45 = 0;
    nprotonRP220_sec56 = 0;
    nprotonRP210_sec45_Reduced = 0;
    nprotonRP220_sec45_Reduced = 0;
    nprotonRP210_sec56_Reduced = 0;
    nprotonRP220_sec56_Reduced = 0;
    nprotonMulti_sec45 = 0;
    nprotonMulti_sec56 = 0;

    if((strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0)){
      // Era flag
      int random_era;
      if(std::string(era)=="A"||std::string(era)=="a") era_ = 0;
      else if(std::string(era)=="B"||std::string(era)=="b") era_=1;
      else if(std::string(era)=="C"||std::string(era)=="c") era_=2;
      else if(std::string(era)=="D"||std::string(era)=="d") era_=3;
      else if(std::string(era)=="E"||std::string(era)=="e") era_=4;
      else if(std::string(era)=="F"||std::string(era)=="f") era_=5;
      else if(std::string(era)=="preTS2"||std::string(era)=="pre TS2"){
	random_era = rand.Integer(3) + 1; //preTS2 refers to era B, C and era D
	era_ = random_era;
      }
      else if(std::string(era)=="postTS2"||std::string(era)=="post TS2"){
	random_era = rand.Integer(2) + 4; //post TS2 refers to era E and F
	era_ = random_era;
      }
      else era_=-1;

      if(std::string(xa) == "120") xangle = 120;
      else if(std::string(xa) == "130") xangle = 130;
      else if(std::string(xa) == "140") xangle = 140;
      else if(std::string(xa) == "150") xangle = 150;
    }

    // Loop over protons with single RP
    for(std::vector<int>::size_type i = 0; i != singleProtonArm->size(); i++){
      if(singleProtonStation->at(i)==0&&singleProtonPot->at(i)==3){
	if(singleProtonArm->at(i)==0){
	  xi_rp210_Arm45 = singleProtonXi->at(i);
	  nprotonRP210_sec45++;
	  if(xi_rp210_Arm45 > 0.04) ++nprotonRP210_sec45_Reduced;
	}
	if(singleProtonArm->at(i)==1){
	  xi_rp210_Arm56 = singleProtonXi->at(i);
	  nprotonRP210_sec56++;
	  if(xi_rp210_Arm56 > 0.04) ++nprotonRP210_sec56_Reduced;
	}
      }
      if(singleProtonStation->at(i)==2&&singleProtonPot->at(i)==3){
	if(singleProtonArm->at(i)==0){
	  xi_rp220_Arm45 = singleProtonXi->at(i);
	  nprotonRP220_sec45++;
	  if(xi_rp220_Arm45 > 0.04) ++nprotonRP220_sec45_Reduced;
	}
	if(singleProtonArm->at(i)==1){
	  xi_rp220_Arm56 = singleProtonXi->at(i);
	  nprotonRP220_sec56++;
	  if(xi_rp220_Arm56 > 0.04) ++nprotonRP220_sec56_Reduced;
	}
      }
    }

    // Loop over protons with multiple RP
    for(std::vector<int>::size_type i = 0; i != multiProtonArm->size(); i++){
      if(multiProtonArm->at(i)==0){
	xi_multiArm45 = multiProtonXi->at(i);
	time_multiArm45 = multiProtonTime->at(i);
	timeunc_multiArm45 = multiProtonTimeError->at(i);
	thx_multiArm45 = multiProtonThetaX->at(i);
	thy_multiArm45 = multiProtonThetaY->at(i);
	nprotonMulti_sec45++;
      }
      if(multiProtonArm->at(i)==1){
	xi_multiArm56 = multiProtonXi->at(i);
	time_multiArm56 = multiProtonTime->at(i);
	timeunc_multiArm56 = multiProtonTimeError->at(i);
	thx_multiArm56 = multiProtonThetaX->at(i);
	thy_multiArm56 = multiProtonThetaY->at(i);
	nprotonMulti_sec56++;
      }
    }

    if(nprotonRP210_sec45==1 && nprotonRP210_sec56==1) isprotonRP210 = true;
    if(nprotonRP220_sec45==1 && nprotonRP220_sec56==1) isprotonRP220 = true;
    if(nprotonRP210_sec45_Reduced==1 && nprotonRP210_sec56_Reduced==1) isprotonRP210_Reduced = true;
    if(nprotonRP220_sec45_Reduced==1 && nprotonRP220_sec56_Reduced==1) isprotonRP220_Reduced = true;
    if(nprotonMulti_sec45==1 && nprotonMulti_sec56==1) isprotonMulti = true;
    if(!protonsfilter){
      isprotonRP210 = true;
      isprotonRP220 = true;
      isprotonMulti = true;
      isprotonRP210_Reduced = true;
      isprotonRP220_Reduced = true;
    }

    if((strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0) || (strcmp(mode, "Muon")==0 || strcmp(mode, "Electron")==0 || strcmp(mode, "Bjets")==0)){
      if(randomFlag){
	if(strcmp(datatype, "DY")==0){
	  xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx];
	  xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	  xi_multiArm45 = rndXi_Multi_sec45[random_idx];
	  xi_rp210_Arm56 = rndXi_RP210_sec56[random_idx];
	  xi_rp220_Arm56 = rndXi_RP220_sec56[random_idx];
	  xi_multiArm56 = rndXi_Multi_sec56[random_idx];
	  is2PUproton_ = true;
	  isprotonRP210 = true;
	  isprotonMulti = true;	  
	}
	else if(strcmp(datatype, "Signal")==0 || strcmp(datatype, "SD")==0){
	  //Correcting for the efficiency
	  //Reading Efficiency File
	  TFile *f_eff = TFile::Open("PreliminaryEfficiencies_October92019_1D2DMultiTrack.root");

	  Double_t rad_eff_45[60];
	  Double_t rad_eff_56[60];
	  TH1D *h_rad_45;
	  TH1D *h_rad_56;

	  if (era_ == 1){
	    if (xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017B/h45_2017B_120_1D");     
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017B/h56_2017B_120_1D");    
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017B/h45_2017B_130_1D");  
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017B/h56_2017B_130_1D"); 
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }   
	    }
	    else if (xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017B/h45_2017B_140_1D");  
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017B/h56_2017B_140_1D"); 
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 150){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017B/h45_2017B_150_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017B/h56_2017B_150_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	  }
	  else if (era_ == 2){
	    if (xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017C/h45_2017C_120_1D");  
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017C/h56_2017C_120_1D"); 
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017C/h45_2017C_130_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017C/h56_2017C_130_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017C/h45_2017C_140_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017C/h56_2017C_140_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 150){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017C/h45_2017C_150_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017C/h56_2017C_150_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	  }
	  else if (era_ == 3){
	    if (xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017D/h45_2017D_120_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017D/h56_2017D_120_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017D/h45_2017D_130_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017D/h56_2017D_130_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017D/h45_2017D_140_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017D/h56_2017D_140_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 150){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017D/h45_2017D_150_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017D/h56_2017D_150_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	  }
	  else if (era_ == 4){
	    if (xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017E/h45_2017E_120_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017E/h56_2017E_120_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017E/h45_2017E_130_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017E/h56_2017E_130_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017E/h45_2017E_140_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017E/h56_2017E_140_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 150){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017E/h45_2017E_150_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017E/h56_2017E_150_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	  }
	  else if (era_ == 5){
	    if (xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017F/h45_2017F_120_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017F/h56_2017F_120_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017F/h45_2017F_130_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017F/h56_2017F_130_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017F/h45_2017F_140_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017F/h56_2017F_140_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (xangle == 150){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017F/h45_2017F_150_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017F/h56_2017F_150_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	  }

	  f_eff->Close();
	  delete f_eff;

	  double random_eff_45;
	  double random_eff_56;
	  if (nprotonRP210_sec45>1 || nprotonRP220_sec56>1){
	    std::cout << "More than 1 strip hit" << std::endl;
	  }
	  if (nprotonRP210_sec45 > 0){
	    random_eff_45 = rand.Rndm();	  
	    for (int j = 0; j < 60; j++){         
	      if(xi_rp210_Arm45 > (j * 0.005) && xi_rp210_Arm45 < ((j * 0.005)+0.005)){
		if(debug){
		  std::cout << "####### 45 #######"<< std::endl;
		  std::cout << "era\t" << era_ << std::endl;
		  std::cout << "xi_45\t" << xi_rp210_Arm45 << std::endl;
		  std::cout << "random_45\t" << random_eff_45 << std::endl;
		  std::cout << "rad_eff_45\t" << rad_eff_45[j] << std::endl;
		  std::cout << "##################"<< std::endl;
		}
		if (random_eff_45 > rad_eff_45[j]){   
		  nprotonRP210_sec45 = 0;
		  nprotonMulti_sec45 = 0;
		  xi_rp210_Arm45 = 0;
		  xi_multiArm45 = 0;
		  isprotonRP210 = false;
		  isprotonMulti = false;
		  goto end45;
		}  
		else {goto end45;}
	      }
	    }
	  }
end45:;

      if (nprotonRP210_sec56 > 0){
	random_eff_56 = rand.Rndm();
	for (int j = 0; j < 60; j++){
	  if(xi_rp210_Arm56 > (j * 0.005) && xi_rp210_Arm56 < ((j * 0.005)+0.005)){
	    if(debug){
	      std::cout << "####### 56 #######"<< std::endl;
	      std::cout << "era\t" << era_ << std::endl;
	      std::cout << "xi_56\t" << xi_rp210_Arm56 << std::endl;
	      std::cout << "random_56\t" << random_eff_56 << std::endl;
	      std::cout << "rad_eff_56\t" << rad_eff_56[j] << std::endl;
	      std::cout << "##################"<< std::endl;
	    }
	    if (random_eff_56 > rad_eff_56[j]){
	      nprotonRP210_sec56 = 0;
	      nprotonMulti_sec56 = 0;
	      isprotonRP210 = false;
	      isprotonMulti = false;
	      xi_rp210_Arm56 = 0;
	      xi_multiArm56 = 0;
	      goto end56;
	    }
	    else {goto end56;}
	  }
	}
      }
end56:;

      //Out of acceptance background
      if(nprotonRP210_sec45 == 0 && nprotonRP210_sec56 > 0){
	xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx];
	xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	xi_multiArm45  = rndXi_Multi_sec45[random_idx];
	nprotonRP210_sec45 = 1;
	nprotonRP220_sec45++;
	nprotonMulti_sec45 = 1;
	is45PUproton_ = true;
	isprotonRP210 = true;
      }
      else if(nprotonRP210_sec56 == 0 && nprotonRP210_sec45 > 0){
	xi_rp210_Arm56 = rndXi_RP210_sec56[random_idx];
	xi_rp220_Arm56 = rndXi_RP220_sec56[random_idx];
	xi_multiArm56  = rndXi_Multi_sec56[random_idx];
	nprotonRP210_sec56 = 1;
	nprotonRP220_sec56++;
	nprotonMulti_sec56 = 1;
	is56PUproton_ = true;
	isprotonRP210 = true;
      }
      else if(nprotonRP210_sec56 == 0 && nprotonRP210_sec45 == 0){
	xi_rp210_Arm56 = rndXi_RP210_sec56[random_idx];
	xi_rp220_Arm56 = rndXi_RP220_sec56[random_idx];
	xi_multiArm56 = rndXi_Multi_sec56[random_idx];
	xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx];
	xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	xi_multiArm45 = rndXi_Multi_sec45[random_idx];
	nprotonRP210_sec45 = 1;
	nprotonRP220_sec45++;
	nprotonMulti_sec45 = 1;
	nprotonRP210_sec56 = 1;
	nprotonRP220_sec56++;
	nprotonMulti_sec56 = 1;
	is2PUproton_ = true;
	isprotonRP210 = true;
      }
      else if(nprotonRP210_sec56 == 1 && nprotonRP210_sec45 == 1){
	is2PUproton_ = false;
	is45PUproton_ = false;
	is56PUproton_ = false;
	isprotonRP210 = true;	    
      }
	}
	if(nprotonRP210_sec45==1 && nprotonRP210_sec56==1) isprotonRP210 = true;
	if(nprotonRP220_sec45==1 && nprotonRP220_sec56==1) isprotonRP220 = true;
	if(nprotonMulti_sec45==1 && nprotonMulti_sec56==1) isprotonMulti = true;
	if(single){
	  xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx]; 
	  xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	  xi_multiArm45 = rndXi_Multi_sec45[random_idx];
	  isprotonRP210 = true;
	  isprotonRP220 = true;
	  isprotonMulti = true;
	}
	else{
	  xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx];
	  xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	  xi_multiArm45 = rndXi_Multi_sec45[random_idx];
	  xi_rp210_Arm56 = rndXi_RP210_sec56[random_idx];
	  xi_rp220_Arm56 = rndXi_RP220_sec56[random_idx];
	  xi_multiArm56 = rndXi_Multi_sec56[random_idx];
	  isprotonRP210 = true;
	  isprotonRP220 = true;
	  isprotonMulti = true;
	}
	if(debug){
	  std::cout << "proton pixels Xi, 210, arm 45: " << xi_rp210_Arm45 << std::endl;
	  std::cout << "proton pixels Xi, 210, arm 56: " << xi_rp210_Arm56 << std::endl;
	  std::cout << "proton pixels Xi, 220, arm 45: " << xi_rp220_Arm45 << std::endl;
	  std::cout << "proton pixels Xi, 220, arm 56: " << xi_rp220_Arm56 << std::endl;
	  std::cout << "proton multi Xi, arm 45: " << xi_multiArm45 << std::endl;
	  std::cout << "proton multi Xi, arm 56: " << xi_multiArm56 << std::endl;
	}
      }

      TLorentzVector p1RP210;
      TLorentzVector p2RP210;
      TLorentzVector p1RP220;
      TLorentzVector p2RP220;
      TLorentzVector p1Multi;
      TLorentzVector p2Multi;

      TLorentzVector genlepton1;
      TLorentzVector genlepton2;
      TLorentzVector gendilepton;
      TLorentzVector genmet;

      TLorentzVector gendijet;
      TLorentzVector genjet1;
      TLorentzVector genjet2;

      TLorentzVector lepton1;
      TLorentzVector lepton2;
      TLorentzVector leptonmet;
      TLorentzVector met;
      TLorentzVector dilepton;

      TLorentzVector lepton;
      TLorentzVector leptonsystem;
      TLorentzVector jet;
      TLorentzVector jetsystem;

      TLorentzVector jetcandidate;
      TLorentzVector jetcandidatesystem;

      TLorentzVector dijet;
      TLorentzVector jet1;
      TLorentzVector jet2;

      TLorentzVector X_dijetRP210;
      TLorentzVector X_dijetRP220;
      TLorentzVector X_dijetMulti;

      TLorentzVector X_dileptonRP210;
      TLorentzVector X_dileptonRP220;
      TLorentzVector X_dileptonMulti;

      TLorentzVector X_dileptonjet1RP210;
      TLorentzVector X_dileptonjet1RP220;
      TLorentzVector X_dileptonjet1Multi;

      TLorentzVector X_dileptondijetRP210;
      TLorentzVector X_dileptondijetRP220;
      TLorentzVector X_dileptondijetMulti;

      p1RP210.SetPxPyPzE(0.,0., ECM*xi_rp210_Arm45/2., ECM*xi_rp210_Arm45/2.);
      p2RP210.SetPxPyPzE(0.,0., -ECM*xi_rp210_Arm56/2., ECM*xi_rp210_Arm56/2.);

      p1RP220.SetPxPyPzE(0.,0., ECM*xi_rp220_Arm45/2., ECM*xi_rp220_Arm45/2.);
      p2RP220.SetPxPyPzE(0.,0., -ECM*xi_rp220_Arm56/2., ECM*xi_rp220_Arm56/2.);

      p1Multi.SetPxPyPzE(0.,0., ECM*xi_multiArm45/2., ECM*xi_multiArm45/2.);
      p2Multi.SetPxPyPzE(0.,0., -ECM*xi_multiArm56/2., ECM*xi_multiArm56/2.);

      diffMassRP210 = ECM*sqrt(xi_rp210_Arm45*xi_rp210_Arm56); 
      diffMassRP220 = ECM*sqrt(xi_rp220_Arm45*xi_rp220_Arm56); 
      diffMassMulti = ECM*sqrt(xi_multiArm45*xi_multiArm56); 

      proton_pz_rp210Arm45 = (1.-(xi_rp210_Arm45))*(ECM/2.);
      proton_pz_rp210Arm56 = -(1.-(xi_rp210_Arm56))*(ECM/2.);

      proton_pz_rp220Arm45 = (1.-(xi_rp220_Arm45))*(ECM/2.);
      proton_pz_rp220Arm56 = -(1.-(xi_rp220_Arm56))*(ECM/2.);

      proton_pz_multiArm45 = (1.-(xi_multiArm45))*(ECM/2.);
      proton_pz_multiArm56 = -(1.-(xi_multiArm56))*(ECM/2.);

      run_ = run;
      event_ = event;
      lumiblock_ = lumiblock;
      xangle_ = xangle;

      std::vector<int> index_leptons;
      index_leptons.clear();
      for (auto i: sort_indexes(*leptons_pt)) {
	index_leptons.push_back(i);
      }

      // Filling TTree                                                           
      if(strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0){
	PUInterac_ = PUInterac;
	PUTrueInterac_ = PUTrueInterac;
	genmet.SetPtEtaPhiM(genmissEt, 0, genmissEt_phi, 0);
	genmet.SetPz(0);
	genmet.SetE(genmissEt);

	if(genleptons_pt->size()>0){
	  genlepton1.SetPtEtaPhiE(genleptons_pt->at(0), genleptons_eta->at(0), genleptons_phi->at(0), genleptons_energy->at(0));
	}
	if(genleptons_pt->size()>1){
	  genlepton2.SetPtEtaPhiE(genleptons_pt->at(1), genleptons_eta->at(1), genleptons_phi->at(1), genleptons_energy->at(1));
	}
	gendilepton = genlepton1+genlepton2;

	if(genjetsak4_pt->size()>1 && jetsak4_pt->size()>1){
	  if(MCMatch(genjetsak4_phi->at(0), genjetsak4_eta->at(0), genjetsak4_pt->at(0), jetsak4_phi->at(0), jetsak4_eta->at(0), jetsak4_pt->at(0)) && MCMatch(genjetsak4_phi->at(1), genjetsak4_eta->at(1), genjetsak4_pt->at(1), jetsak4_phi->at(1), jetsak4_eta->at(1), jetsak4_pt->at(1))){
	    genjet1.SetPtEtaPhiE(genjetsak4_pt->at(0), genjetsak4_eta->at(0), genjetsak4_phi->at(0), genjetsak4_energy->at(0));
	    genjet2.SetPtEtaPhiE(genjetsak4_pt->at(1), genjetsak4_eta->at(1), genjetsak4_phi->at(1), genjetsak4_energy->at(1));
	    gendijet = genjet1 + genjet2;
	  }
	}
      }

      if(createProtonFile) {
	if(isprotonRP210 && isprotonRP220 && isprotonMulti && xangle == atoi(xa)){
	  if((strcmp(mode, "Muon")==0 || strcmp(mode, "Electron")==0) && dileptonPt>10){
	    file_protons_reco << xi_rp210_Arm45 << "\t" << xi_rp210_Arm56 << "\t" << xi_rp220_Arm45 << "\t"<< xi_rp220_Arm56 << "\t" << xi_multiArm45 << "\t" << xi_multiArm56 << "\n"; //write to file
	  }
	}
      }

      // HLT Muons: 'HLT_IsoMu27_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*'
      // HLT Electrons: 'HLT_Ele27_WPTight_Gsf_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_DoubleEle33_CaloIdL_MW_v*'
      // HLT_ZeroBias: 'HLT_ZeroBias_v*', 'HLT_ZeroBias_FirstBXAfterTrain_v*','HLT_ZeroBias_IsolatedBunches_v*','HLT_ZeroBias_Alignment_v*','HLT_ZeroBias_Beamspot_v*','HLT_L1UnpairedBunchBptxMinus_v*','HLT_L1UnpairedBunchBptxPlus_v*','HLT_Physics_v*','HLT_L1SingleMu_v*'
      // If -1: not present in the path, 0: not triggered and 1 triggered.

      if(zerobias){
	if(trigger->at(0)==1){
	  triggerZeroBias = true;
	  prescalesL1ZeroBias_ = prescalesL1->at(0);
	}
	if(trigger->at(1)==1){
	  triggerZeroBiasAfterTrain = true;
	  prescalesL1ZeroBiasAfterTrain_ = prescalesL1->at(1);
	}
	if(trigger->at(2)==1){
	  triggerZeroBiasIsolatedBx = true;
	  prescalesL1ZeroBiasIsolatedBx_ = prescalesL1->at(2);
	}
	if(trigger->at(3)==1){
	  triggerZeroBiasAlignment = true;
	  prescalesL1ZeroBiasAlignment_ = prescalesL1->at(3);
	}
	if(trigger->at(4)==1){
	  triggerZeroBiasBeamSpot = true;
	  prescalesL1ZeroBiasBeamSpot_ = prescalesL1->at(4);
	}
	if(trigger->at(5)==1){
	  triggerZeroBiasUnpairedBptxMinus = true;
	  prescalesL1ZeroBiasUnpairedBptxMinus_ = prescalesL1->at(5);
	}
	if(trigger->at(6)==1){
	  triggerZeroBiasUnpairedBptxPlus = true;
	  prescalesL1ZeroBiasUnpairedBptxPlus_ = prescalesL1->at(6);
	}
	if(trigger->at(7)==1){
	  triggerPhysics = true;
	  prescalesL1Physics_ = prescalesL1->at(7);
	}
	if(trigger->at(8)==1){
	  triggerL1SingleMu = true;
	  prescalesL1SingleMu_ = prescalesL1->at(8);
	}
      }else{
	if((strcmp(mode, "Muon")==0 || strcmp(mode, "MC_Muon")==0) && strcmp(datatype, "SD")!=0 && strcmp(datatype, "Signal")!=0){
	  if(trigger->at(0)==1) triggerIsoMu27 = true;
	  if(trigger->at(1)==1) triggerMu17TrkIsoMu8TrkIso = true;
	  if(trigger->at(2)==1) triggerMu17TrkIsoMu8TrkIsoMass8 = true;
	  if(trigger->at(3)==1) triggerMu17TrkIsoMu8TrkIsoMass3 = true;
	  if(trigger->at(4)==1) triggerDoubleMu43 = true;
	  if(trigger->at(5)==1) triggerDoubleMu48 = true;
	}else if(strcmp(mode, "Electron")==0 || strcmp(mode, "MC_Electron")==0 && strcmp(datatype, "SD")!=0 && strcmp(datatype, "Signal")!=0){
	  if(trigger->at(0)==1) triggerEle27 = true;
	  if(trigger->at(1)==1) triggerEle23Ele12 = true;
	  if(trigger->at(2)==1) triggerEle23Ele12Dz = true;
	  if(trigger->at(3)==1) triggerDoubleEle33 = true;
	  if(trigger->at(4)==1) triggerBTagMu5Ak4dijet20 = true;
	  if(trigger->at(5)==1) triggerBTagMu5Ak4dijet40 = true;
	  if(trigger->at(6)==1) triggerBTagMu5Ak4dijet70 = true;
	  if(trigger->at(7)==1) triggerBTagMu5Ak4dijet110 = true;
	  if(trigger->at(8)==1) triggerBTagMu5Ak4dijet170 = true;
	  if(trigger->at(9)==1) triggerBTagMu5Ak4dijet300 = true;
	  if(trigger->at(10)==1) triggerBTagMu5Ak8dijet170 = true;
	  if(trigger->at(11)==1) triggerBTagMu5Ak8dijet300 = true;
	  if(trigger->at(12)==1) triggerPFHT380DoubleBTag = true;
	  if(trigger->at(13)==1) triggerPFHT430DoubleBTag = true;
	  if(trigger->at(14)==1) triggerPFMET100BTag = true;
	  if(trigger->at(15)==1) triggerPFMET110BTag = true;
	  if(trigger->at(16)==1) triggerPFMET120BTag = true;
	  if(trigger->at(17)==1) triggerPFMET130BTag = true;
	  if(trigger->at(18)==1) triggerPFMET140BTag = true;
	  if(trigger->at(19)==1) triggerEle15PFHT450BTag = true;
	  if(trigger->at(20)==1) triggerHT250BTagScouting = true;
	}else if(strcmp(mode, "Bjets")==0 || (strcmp(datatype, "SD")!=0 && strcmp(datatype, "Signal")!=0)){
	  if(trigger->at(0)==1) triggerIsoMu27 = true;
	  if(trigger->at(1)==1) triggerMu17TrkIsoMu8TrkIso = true;
	  if(trigger->at(2)==1) triggerMu17TrkIsoMu8TrkIsoMass8 = true;
	  if(trigger->at(3)==1) triggerMu17TrkIsoMu8TrkIsoMass3 = true;
	  if(trigger->at(4)==1) triggerDoubleMu43 = true;
	  if(trigger->at(5)==1) triggerDoubleMu48 = true;
	  if(trigger->at(6)==1) triggerBTagMu5Ak4dijet20 = true;
	  if(trigger->at(7)==1) triggerBTagMu5Ak4dijet40 = true;
	  if(trigger->at(8)==1) triggerBTagMu5Ak4dijet70 = true;
	  if(trigger->at(9)==1) triggerBTagMu5Ak4dijet110 = true;
	  if(trigger->at(10)==1) triggerBTagMu5Ak4dijet170 = true;
	  if(trigger->at(11)==1) triggerBTagMu5Ak4dijet300 = true;
	  if(trigger->at(12)==1) triggerBTagMu5Ak8dijet170 = true;
	  if(trigger->at(13)==1) triggerBTagMu5Ak8dijet300 = true;
	  if(trigger->at(14)==1) triggerPFHT380DoubleBTag = true;
	  if(trigger->at(15)==1) triggerPFHT430DoubleBTag = true;
	}else if(strcmp(mode, "EMu")==0){
	  if(trigger->at(0)==1) triggerMu23TrkIsoEle12 = true;        
	  if(trigger->at(1)==1) triggerMu23TrkIsoEle12DZ = true;
	  if(trigger->at(2)==1) triggerMu8TrkIsoEle23DZ = true;
	}else{
	  std::cout << "\nNo Mode option!\n" << std::endl;
	  exit(EXIT_FAILURE);
	}
      }

      /*
	 if(zerobias){
	 if(!(vertex_z->size()>0)) continue;
	 }else{
	 if(!(leptons_pt->size()>1&&jetsak4_pt->size()>1)) continue;
	 }
	 */

      if(!(isprotonRP210 || isprotonRP220 || isprotonMulti)) continue;
      if(!(vertex_z->size()>0)) continue;
      if(!(missEt>30)) continue;

      met.SetPtEtaPhiM(missEt, 0, missEt_phi, 0);
      met.SetPz(0);
      met.SetE(missEt);

      if(leptons_pt->size()>0){
	lepton1.SetPtEtaPhiE(leptons_pt->at(index_leptons.at(0)), leptons_eta->at(index_leptons.at(0)), leptons_phi->at(index_leptons.at(0)), leptons_energy->at(index_leptons.at(0)));
      }
      if(leptons_pt->size()>1){
	lepton2.SetPtEtaPhiE(leptons_pt->at(index_leptons.at(1)), leptons_eta->at(index_leptons.at(1)), leptons_phi->at(index_leptons.at(1)), leptons_energy->at(index_leptons.at(1)));
      }
      leptonmet = lepton1+met;
      dilepton = lepton1+lepton2;

      acoplanarity = 1. - fabs(lepton1.Phi()-lepton2.Phi())/pi;

      // Leptons and Jets are already sorted by pT
      for(std::vector<int>::size_type i = 0; i != leptons_pt->size(); i++){
	lepton.SetPtEtaPhiE(leptons_pt->at(i), leptons_eta->at(i), leptons_phi->at(i), leptons_energy->at(i));
	leptonsystem+=lepton;
	leptonsystemCharge+=leptons_charge->at(i);
      }

      nJetsCandidatesLoose = 0;
      nJetsCandidatesMedium = 0;
      nJetsCandidatesTight = 0;
      if(leptons_pt->size()>1){
	for(std::vector<int>::size_type i = 0; i != jetsak4_pt->size(); i++){
	  jet.SetPtEtaPhiE(jetsak4_pt->at(i), jetsak4_eta->at(i), jetsak4_phi->at(i), jetsak4_energy->at(i));
	  jetsystem+=jet;
	  // Counting Jets Close to the Leading and SubLeading Lepton (200 microns)
	  if(fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(0)))<0.02&&fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(1)))<0.02&&jetsak4_pt->at(i)>15.){
	    nJetsCandidatesTight++;
	    jetcandidate.SetPtEtaPhiE(jetsak4_pt->at(i), jetsak4_eta->at(i), jetsak4_phi->at(i), jetsak4_energy->at(i));
	    jetcandidatesystem+=jetcandidate;
	  }
	  if(fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(0)))<0.05&&fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(1)))<0.05&&jetsak4_pt->at(i)>15.){
	    nJetsCandidatesMedium++;
	  }
	  if(fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(0)))<0.1&&fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(1)))<0.1&&jetsak4_pt->at(i)>15.){
	    nJetsCandidatesLoose++;
	  }
	}
      }

      dilepton.SetPtEtaPhiM(dilepton.Pt(), dilepton.Eta(), dilepton.Phi(), dilepton.M());
      met.SetPtEtaPhiE(missEt, 0, missEt_phi, 0);
      Rmpf = 1.+(met.Px() * dilepton.Px() + met.Py() * dilepton.Py()) / (dilepton.Pt()*dilepton.Pt());

      if(jetsak4_pt->size()>0){
	jet1.SetPtEtaPhiE(jetsak4_pt->at(0), jetsak4_eta->at(0), jetsak4_phi->at(0), jetsak4_energy->at(0));
	JZBalance = jet1.Pt() - dilepton.Pt();
      }

      if(jetsak4_pt->size()>1){
	jet2.SetPtEtaPhiE(jetsak4_pt->at(1), jetsak4_eta->at(1), jetsak4_phi->at(1), jetsak4_energy->at(1));
	dijet = jet1 + jet2;
	JZBalance2 = (dijet).Pt() - dilepton.Pt();
      }

      // Invisible particle associated with leptonsystem

      X_dijetRP210 = p1RP210 + p2RP210 - dijet;
      X_dijetRP220 = p1RP220 + p2RP220 - dijet;
      X_dijetMulti = p1Multi + p2Multi - dijet;

      X_dileptonRP210 = p1RP210 + p2RP210 - dilepton;
      X_dileptonRP220 = p1RP220 + p2RP220 - dilepton;
      X_dileptonMulti = p1Multi + p2Multi - dilepton;

      X_dileptonjet1RP210 = p1RP210 + p2RP210 - dilepton - jet1;
      X_dileptonjet1RP220 = p1RP220 + p2RP220 - dilepton - jet1;
      X_dileptonjet1Multi =  p1Multi + p2Multi - dilepton - jet1;

      X_dileptondijetRP210 =  p1RP210 + p2RP210 - dilepton - dijet;
      X_dileptondijetRP220 =  p1RP220 + p2RP220 - dilepton - dijet;
      X_dileptondijetMulti =  p1Multi + p2Multi - dilepton - dijet;

      if(debug){
	std::cout << "\n\t == Event " << jentry << " ==" << std::endl;
	for(std::vector<int>::size_type i = 0; i != trigger->size(); i++){
	  std::cout << "\t\t --> Trigger (" << i << "): " << trigger->at(i) << std::endl;
	}
	std::cout << "\t\t --> # Jets: " << jetsak4_pt->size() << std::endl;
	std::cout << "\t\t --> # Fat Jets: " << jetsak8_phi->size() << std::endl;
	std::cout << "\t\t --> # Leptons: " << leptons_pt->size() << std::endl;
	std::cout << "\t\t --> lepton system Mass [GeV]: " << leptonsystem.M() << std::endl;
	std::cout << "\t\t --> jet system Mass [GeV]: " << jetsystem.M() << std::endl;
	std::cout << "\t\t --> candidate jet system Mass [GeV]: " << jetcandidatesystem.M() << std::endl;
      }

      nVertex = vertex_z->size();
      nTracksPerVertex = vertex_ntrack->at(0);
      pvVertexX = vertex_x->at(0);
      pvVertexY = vertex_y->at(0);
      pvVertexZ = vertex_z->at(0);
      pvVertexR = sqrt(vertex_x->at(0) * vertex_x->at(0) + vertex_y->at(0) * vertex_y->at(0));
      pvVertexPhi = atan2(vertex_y->at(0), vertex_x->at(0));

      for (auto evtVz : *vertex_z){
	allVertexZ.push_back(evtVz);
      }

      if(strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0){
	nGenLeptons = genleptons_pt->size();
	nGenParticles_ = nGenParticles;
	nGenElectrons_ = nGenElectrons;
	nGenMuons_ = nGenMuons;
	if(genleptons_pt->size()>0){
	  genleadingLeptonPDGId = genleptons_pdgid->at(0);
	  genleadingLeptonEnergy = genleptons_energy->at(0);
	  genleadingLeptonPx = genleptons_px->at(0);
	  genleadingLeptonPy = genleptons_py->at(0);
	  genleadingLeptonPz = genleptons_pz->at(0);
	  genleadingLeptonPt = genleptons_pt->at(0);
	  genleadingLeptonEta = genleptons_eta->at(0);
	  genleadingLeptonPhi = genleptons_phi->at(0);
	  genleadingLeptonCharge = genleptons_charge->at(0);
	  genleadingLeptonVx = genleptons_vx->at(0);
	  genleadingLeptonVy = genleptons_vy->at(0);
	  genleadingLeptonVz = genleptons_vz->at(0);
	  genleadingLeptonVr = sqrt(genleptons_vx->at(0) * genleptons_vx->at(0) + genleptons_vy->at(0) * genleptons_vy->at(0));
	  genleadingLeptonVphi = atan2(genleptons_vy->at(0), genleptons_vx->at(0));
	}

	if(genleptons_pt->size()>1){
	  gensecondLeptonPDGId = genleptons_pdgid->at(1);
	  gensecondLeptonEnergy = genleptons_energy->at(1);
	  gensecondLeptonPx = genleptons_px->at(1);
	  gensecondLeptonPy = genleptons_py->at(1);
	  gensecondLeptonPz = genleptons_pz->at(1);
	  gensecondLeptonPt = genleptons_pt->at(1);
	  gensecondLeptonEta = genleptons_eta->at(1);
	  gensecondLeptonPhi = genleptons_phi->at(1);
	  gensecondLeptonCharge = genleptons_charge->at(1);
	  gensecondLeptonVx = genleptons_vx->at(1);
	  gensecondLeptonVy = genleptons_vy->at(1);
	  gensecondLeptonVz = genleptons_vz->at(1);
	  gensecondLeptonVr = sqrt(genleptons_vx->at(1) * genleptons_vx->at(1) + genleptons_vy->at(1) * genleptons_vy->at(1));
	  gensecondLeptonVphi = atan2(genleptons_vy->at(1), genleptons_vx->at(1));
	}

	if(genleptons_pt->size()>1){
	  gendileptonCharge = genleptons_charge->at(0) + genleptons_charge->at(1);
	  gendileptonMass = gendilepton.M();
	  gendileptonEta = gendilepton.Eta();
	  gendileptonPhi = gendilepton.Phi();
	  gendileptonPt = gendilepton.Pt();
	  gendileptonRapidity = gendilepton.Rapidity();
	}

	genmissEt_ = genmissEt;
	genmissEt_phi_ = genmissEt_phi;

	nGenJets = jetsak4_pt->size();
	if(genjetsak4_pt->size()>0 && jetsak4_pt->size()>0){
	  if(MCMatch(genjetsak4_phi->at(0), genjetsak4_eta->at(0), genjetsak4_pt->at(0), jetsak4_phi->at(0), jetsak4_eta->at(0), jetsak4_pt->at(0))){
	    genleadingJetEnergy = genjetsak4_energy->at(0);
	    genleadingJetPx = genjetsak4_px->at(0);
	    genleadingJetPy = genjetsak4_py->at(0);
	    genleadingJetPz = genjetsak4_pz->at(0);
	    genleadingJetPt = genjetsak4_pt->at(0);
	    genleadingJetEta = genjetsak4_eta->at(0);
	    genleadingJetPhi = genjetsak4_phi->at(0);
	    genleadingJetVz = genjetsak4_vz->at(0);
	  }
	}
	if(genjetsak4_pt->size()>1 && jetsak4_pt->size()>1){
	  if(MCMatch(genjetsak4_phi->at(1), genjetsak4_eta->at(1), genjetsak4_pt->at(1), jetsak4_phi->at(1), jetsak4_eta->at(1), jetsak4_pt->at(1))){
	    gensecondJetEnergy = genjetsak4_energy->at(1);
	    gensecondJetPx = genjetsak4_px->at(1);
	    gensecondJetPy = genjetsak4_py->at(1);
	    gensecondJetPz = genjetsak4_pz->at(1);
	    gensecondJetPt = genjetsak4_pt->at(1);
	    gensecondJetEta = genjetsak4_eta->at(1);
	    gensecondJetPhi = genjetsak4_phi->at(1);
	    gensecondJetVz = genjetsak4_vz->at(1);
	  }
	}

	gendijetMass = gendijet.M();
	gendijetEta = gendijet.Eta();
	gendijetPhi = gendijet.Phi();
	gendijetPt = gendijet.Pt();
	gendijetRapidity = gendijet.Rapidity();
      }

      nLeptons = leptons_pt->size();

      if(leptons_pt->size()>0){
	leadingLeptonPDGId = leptons_pdgid->at(index_leptons.at(0));
	leadingLeptonEnergy = leptons_energy->at(index_leptons.at(0));
	leadingLeptonPx = leptons_px->at(index_leptons.at(0));
	leadingLeptonPy = leptons_py->at(index_leptons.at(0));
	leadingLeptonPz = leptons_pz->at(index_leptons.at(0));
	leadingLeptonPt = leptons_pt->at(index_leptons.at(0));
	leadingLeptonEta = leptons_eta->at(index_leptons.at(0));
	leadingLeptonPhi = leptons_phi->at(index_leptons.at(0));
	leadingLeptonCharge = leptons_charge->at(index_leptons.at(0));
	leadingLeptonVx = leptons_vx->at(index_leptons.at(0));
	leadingLeptonVy = leptons_vy->at(index_leptons.at(0));	
	leadingLeptonVz = leptons_vz->at(index_leptons.at(0));
	leadingLeptonVr = sqrt(leptons_vx->at(index_leptons.at(0)) * leptons_vx->at(index_leptons.at(0)) + leptons_vy->at(index_leptons.at(0)) * leptons_vy->at(index_leptons.at(0)));
	leadingLeptonVphi = atan2(leptons_vy->at(index_leptons.at(0)), leptons_vx->at(index_leptons.at(0)));
	leadingLeptonPFIso = leptons_pfIso_->at(index_leptons.at(0)); 
	leadingLeptonTkIso = leptons_tkIso_->at(index_leptons.at(0));
	leadingLeptonLooseId = leptons_looseId->at(index_leptons.at(0));
	leadingLeptonMediumId = leptons_mediumId->at(index_leptons.at(0));
	leadingLeptonTightId = leptons_tightId->at(index_leptons.at(0));
	leadingLeptonPfIsoMedium = leptons_pfIsoMedium_->at(index_leptons.at(0));
	leadingLeptonMiniIsoTight = leptons_miniIsoTight_->at(index_leptons.at(0));
	leadingLeptonPfIsoVeryTight = leptons_pfIsoVeryTight_->at(index_leptons.at(0));
      }

      if(leptons_pt->size()>1){
	secondLeptonPDGId = leptons_pdgid->at(index_leptons.at(1));
	secondLeptonEnergy = leptons_energy->at(index_leptons.at(1));
	secondLeptonPx = leptons_px->at(index_leptons.at(1));
	secondLeptonPy = leptons_py->at(index_leptons.at(1));
	secondLeptonPz = leptons_pz->at(index_leptons.at(1));
	secondLeptonPt = leptons_pt->at(index_leptons.at(1));
	secondLeptonEta = leptons_eta->at(index_leptons.at(1));
	secondLeptonPhi = leptons_phi->at(index_leptons.at(1));
	secondLeptonCharge = leptons_charge->at(index_leptons.at(1));
	secondLeptonVx = leptons_vx->at(index_leptons.at(1));
	secondLeptonVy = leptons_vy->at(index_leptons.at(1)); 
	secondLeptonVz = leptons_vz->at(index_leptons.at(1));
	secondLeptonVr = sqrt(leptons_vx->at(index_leptons.at(1)) * leptons_vx->at(index_leptons.at(1)) + leptons_vy->at(index_leptons.at(1)) * leptons_vy->at(index_leptons.at(1)));
	secondLeptonVphi = atan2(leptons_vy->at(index_leptons.at(1)), leptons_vx->at(index_leptons.at(1)));
	secondLeptonPFIso = leptons_pfIso_->at(index_leptons.at(1));
	secondLeptonTkIso = leptons_tkIso_->at(index_leptons.at(1));
	secondLeptonLooseId = leptons_looseId->at(index_leptons.at(1));
	secondLeptonMediumId = leptons_mediumId->at(index_leptons.at(1));
	secondLeptonTightId = leptons_tightId->at(index_leptons.at(1));
	secondLeptonPfIsoMedium = leptons_pfIsoMedium_->at(index_leptons.at(1));
	secondLeptonMiniIsoTight = leptons_miniIsoTight_->at(index_leptons.at(1));
	secondLeptonPfIsoVeryTight = leptons_pfIsoVeryTight_->at(index_leptons.at(1));
      }     

      nJets = jetsak4_pt->size();
      if(jetsak4_pt->size()>0){
	leadingJetEnergy = jetsak4_energy->at(0);
	leadingJetPx = jetsak4_px->at(0);
	leadingJetPy = jetsak4_py->at(0);
	leadingJetPz = jetsak4_pz->at(0);
	leadingJetPt = jetsak4_pt->at(0);
	leadingJetEta = jetsak4_eta->at(0);
	leadingJetPhi = jetsak4_phi->at(0);
	leadingJetVx = jetsak4_vx->at(0);
	leadingJetVy = jetsak4_vy->at(0);
	leadingJetVz = jetsak4_vz->at(0);
	leadingJetVr = sqrt(jetsak4_vx->at(0) * jetsak4_vx->at(0) + jetsak4_vy->at(0) * jetsak4_vy->at(0));
	leadingJetVphi = atan2(jetsak4_vy->at(0), jetsak4_vx->at(0));
	leadingJetBtag = jetsak4_bdis->at(0);
	leadingJetQGdis = jetsak4_qgdis->at(0);
	leadingJetNeutralEmFrac = jetsak4_neutralemfrac->at(0);
	leadingJetNeutralHadFrac = jetsak4_neutralhadfrac->at(0);
	leadingJetChargedEmFrac = jetsak4_chargedemfrac->at(0);
	leadingJetChargedHadFrac = jetsak4_chargedhadfrac->at(0);
	leadingJetMuonFrac = jetsak4_muonfrac->at(0);
	leadingJetNeutralMulti = jetsak4_neutralmulti->at(0);
	leadingJetChargedMulti = jetsak4_chargedmulti->at(0);
	leadingJetpuIdfdisc = jetsak4_puIdfdisc->at(0);
	leadingJetpuIdcbased = jetsak4_puIdcbased->at(0);
	leadingJetpuIdfid = jetsak4_puIdfid->at(0);
	leadingJetLooseId = jetsak4_looseId->at(0);
	leadingJetTightId = jetsak4_tightId->at(0);
	leadingJetLepVeto = jetsak4_lepVeto->at(0);
      }
      if(jetsak4_pt->size()>1){
	secondJetEnergy = jetsak4_energy->at(1);
	secondJetPx = jetsak4_px->at(1);
	secondJetPy = jetsak4_py->at(1);
	secondJetPz = jetsak4_pz->at(1);
	secondJetPt = jetsak4_pt->at(1);
	secondJetEta = jetsak4_eta->at(1);
	secondJetPhi = jetsak4_phi->at(1);
	secondJetVx = jetsak4_vx->at(1);
	secondJetVy = jetsak4_vy->at(1);
	secondJetVz = jetsak4_vz->at(1);
	secondJetVr = sqrt(jetsak4_vx->at(1) * jetsak4_vx->at(1) + jetsak4_vy->at(1) * jetsak4_vy->at(1));
	secondJetVphi = atan2(jetsak4_vy->at(1), jetsak4_vx->at(1));
	secondJetBtag = jetsak4_bdis->at(1);
	secondJetQGdis = jetsak4_qgdis->at(1);
	secondJetNeutralEmFrac = jetsak4_neutralemfrac->at(1);
	secondJetNeutralHadFrac = jetsak4_neutralhadfrac->at(1);
	secondJetChargedEmFrac = jetsak4_chargedemfrac->at(1);
	secondJetChargedHadFrac = jetsak4_chargedhadfrac->at(1);
	secondJetMuonFrac = jetsak4_muonfrac->at(1);
	secondJetNeutralMulti = jetsak4_neutralmulti->at(1);
	secondJetChargedMulti = jetsak4_chargedmulti->at(1);
	secondJetpuIdfdisc = jetsak4_puIdfdisc->at(1);
	secondJetpuIdcbased = jetsak4_puIdcbased->at(1);
	secondJetpuIdfid = jetsak4_puIdfid->at(1);
	secondJetLooseId = jetsak4_looseId->at(1);
	secondJetTightId = jetsak4_tightId->at(1);
	secondJetLepVeto = jetsak4_lepVeto->at(1);
      }

      dijetMass = dijet.M();
      dijetEta = dijet.Eta();
      dijetPhi = dijet.Phi();
      dijetPt = dijet.Pt();
      dijetRapidity = dijet.Rapidity();

      jetsystemMass = jetsystem.M();
      jetsystemEta = jetsystem.Eta();
      jetsystemPhi = jetsystem.Phi();
      jetsystemPt = jetsystem.Pt();
      jetsystemRapidity = jetsystem.Rapidity();

      jetcandidatesystemMass = jetcandidatesystem.M();
      jetcandidatesystemEta = jetcandidatesystem.Eta();
      jetcandidatesystemPhi = jetcandidatesystem.Phi();
      jetcandidatesystemPt = jetcandidatesystem.Pt();
      jetcandidatesystemRapidity = jetcandidatesystem.Rapidity();

      if(leptons_pt->size()>0){
	leptonmetCharge = leptons_charge->at(index_leptons.at(0));
	leptonmetMassT = leptonmet.M();
	leptonmetEta = leptonmet.Eta();
	leptonmetPhi = leptonmet.Phi();
	leptonmetPt = leptonmet.Pt();
	leptonmetRapidity = leptonmet.Rapidity();
      }

      if(leptons_pt->size()>1){
	dileptonCharge = leptons_charge->at(index_leptons.at(0)) + leptons_charge->at(index_leptons.at(1));
	dileptonMass = dilepton.M();
	dileptonEta = dilepton.Eta();
	dileptonPhi = dilepton.Phi();
	dileptonPt = dilepton.Pt();
	dileptonRapidity = dilepton.Rapidity();
      }

      leptonsystemMass = leptonsystem.M();
      leptonsystemEta = leptonsystem.Eta();
      leptonsystemPhi = leptonsystem.Phi();
      leptonsystemPt = leptonsystem.Pt();
      leptonsystemRapidity = leptonsystem.Rapidity();

      missingMassDijetRP210 = X_dijetRP210.M();
      missingEtaDijetRP210 = X_dijetRP210.Eta();
      missingPhiDijetRP210 = X_dijetRP210.Phi();
      missingPtDijetRP210 = X_dijetRP210.Pt();
      missingRapidityDijetRP210 = X_dijetRP210.Rapidity();

      missingMassDijetRP220 = X_dijetRP220.M();
      missingEtaDijetRP220 = X_dijetRP220.Eta();
      missingPhiDijetRP220 = X_dijetRP220.Phi();
      missingPtDijetRP220 = X_dijetRP220.Pt();
      missingRapidityDijetRP220 = X_dijetRP220.Rapidity();

      missingMassDijetMulti = X_dijetMulti.M();
      missingEtaDijetMulti = X_dijetMulti.Eta();
      missingPhiDijetMulti = X_dijetMulti.Phi();
      missingPtDijetMulti = X_dijetMulti.Pt();
      missingRapidityDijetMulti = X_dijetMulti.Rapidity();

      missingMassDileptonRP210 = X_dileptonRP210.M();
      missingEtaDileptonRP210 = X_dileptonRP210.Eta();
      missingPhiDileptonRP210 = X_dileptonRP210.Phi();
      missingPtDileptonRP210 = X_dileptonRP210.Pt();
      missingRapidityDileptonRP210 = X_dileptonRP210.Rapidity();

      missingMassDileptonRP220 = X_dileptonRP220.M();
      missingEtaDileptonRP220 = X_dileptonRP220.Eta();
      missingPhiDileptonRP220 = X_dileptonRP220.Phi();
      missingPtDileptonRP220 = X_dileptonRP220.Pt();
      missingRapidityDileptonRP220 = X_dileptonRP220.Rapidity();

      missingMassDileptonMulti = X_dileptonMulti.M();
      missingEtaDileptonMulti = X_dileptonMulti.Eta();
      missingPhiDileptonMulti = X_dileptonMulti.Phi();
      missingPtDileptonMulti = X_dileptonMulti.Pt();
      missingRapidityDileptonMulti = X_dileptonMulti.Rapidity();

      missingMassDileptonJet1RP220 = X_dileptonjet1RP220.M();
      missingEtaDileptonJet1RP220 = X_dileptonjet1RP220.Eta();
      missingPhiDileptonJet1RP220 = X_dileptonjet1RP220.Phi();
      missingPtDileptonJet1RP220 = X_dileptonjet1RP220.Pt();
      missingRapidityDileptonJet1RP220 = X_dileptonjet1RP220.Rapidity();

      missingMassDileptonJet1Multi = X_dileptonjet1Multi.M();
      missingEtaDileptonJet1Multi = X_dileptonjet1Multi.Eta();
      missingPhiDileptonJet1Multi = X_dileptonjet1Multi.Phi();
      missingPtDileptonJet1Multi = X_dileptonjet1Multi.Pt();
      missingRapidityDileptonJet1Multi = X_dileptonjet1Multi.Rapidity();

      missingMassDileptonDijetRP210 = X_dileptondijetRP210.M();
      missingEtaDileptonDijetRP210 = X_dileptondijetRP210.Eta();
      missingPhiDileptonDijetRP210 = X_dileptondijetRP210.Phi();
      missingPtDileptonDijetRP210 = X_dileptondijetRP210.Pt();
      missingRapidityDileptonDijetRP210 = X_dileptondijetRP210.Rapidity();

      missingMassDileptonDijetRP220 = X_dileptondijetRP220.M();
      missingEtaDileptonDijetRP220 = X_dileptondijetRP220.Eta();
      missingPhiDileptonDijetRP220 = X_dileptondijetRP220.Phi();
      missingPtDileptonDijetRP220 = X_dileptondijetRP220.Pt();
      missingRapidityDileptonDijetRP220 = X_dileptondijetRP220.Rapidity();

      missingMassDileptonDijetMulti = X_dileptondijetMulti.M();
      missingEtaDileptonDijetMulti = X_dileptondijetMulti.Eta();
      missingPhiDileptonDijetMulti = X_dileptondijetMulti.Phi();
      missingPtDileptonDijetMulti = X_dileptondijetMulti.Pt();
      missingRapidityDileptonDijetMulti = X_dileptondijetMulti.Rapidity();

      missEt_ = missEt;
      missEt_phi_ = missEt_phi;

      nElectrons_ = nElectrons;
      nMuons_ = nMuons;
      nElectrons_looseId_ = nElectrons_looseId;
      nElectrons_mediumId_ = nElectrons_mediumId;
      nElectrons_tightId_ = nElectrons_tightId;
      nMuons_looseId_ = nMuons_looseId;
      nMuons_mediumId_ = nMuons_mediumId;
      nMuons_tightId_ = nMuons_tightId; 
      nChargedPFMultiPV_Loose_ = nChargedPFMultiPV_Loose;
      nChargedPFMultiPV_Tight_ = nChargedPFMultiPV_Tight;
      nChargedPFMultiPV_UsedInFit_ = nChargedPFMultiPV_UsedInFit;
      nChargedPFMultiPV_Tight_Fit_ = nChargedPFMultiPV_Tight_Fit;
      SumChargedPFMultiPV_pt_Loose_ = SumChargedPFMultiPV_pt_Loose; 
      SumChargedPFMultiPV_pt_Tight_ = SumChargedPFMultiPV_pt_Tight;
      SumChargedPFMultiPV_pt_UsedInFit_ = SumChargedPFMultiPV_pt_UsedInFit;
      SumChargedPFMultiPV_pt_Tight_Fit_ = SumChargedPFMultiPV_pt_Tight_Fit;

      double dphi = fabs(dilepton.Phi()-jetsystem.Phi());
      dphi = (dphi<pi) ? dphi : 2.*pi-dphi;
      dphi_ = dphi;

      if(isprotonMulti &&  timeunc_multiArm45>0 && timeunc_multiArm56>0){
	t45andt56 = c_light*ns_to_s_*m_to_cm_*(time_multiArm56 + time_multiArm45);
	vzpps = (c_light/2)*ns_to_s_*m_to_cm_*(time_multiArm56-time_multiArm45);
	for (auto evtVz : *vertex_z){
	  allVertexZ.push_back(evtVz);
	  MinimumDistance.push_back(std::make_pair(fabs(vzpps-evtVz), evtVz));
	}
      }else{
	t45andt56 = -999;
	vzpps = -999;
      }

      // Finding the minimum distance between CMS vtx and PPS vtx
      std::sort(MinimumDistance.begin(), MinimumDistance.end());

      if(MinimumDistance.size()>0){
	distanceVertexZPPSandCMS = MinimumDistance[0].first;
	selectedVertexZ = MinimumDistance[0].second;
      }else{
	distanceVertexZPPSandCMS = -999.;
	selectedVertexZ = -999.;
      }

      tout->Fill();

      if(createProtonFile) {
	if(isprotonMulti && xangle == atoi(xa)){
	  if((strcmp(mode, "Muon")==0 || strcmp(mode, "Electron")==0  || strcmp(mode, "Bjets")==0) && dileptonPt<10){ 
	    file_protons_reco << xi_rp210_Arm45 << "\t" << xi_rp210_Arm56 << "\t" << xi_rp220_Arm45 << "\t"<< xi_rp220_Arm56 << "\t" << xi_multiArm45 << "\t" << xi_multiArm56 << "\n"; //write to file
	  }
	}
      }

      // Format file suggested by https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPickEvents#Run_edm_PickEvents_py_with_CRAB
      if(createEventFile){
	file_eventlist << run_ <<":"<< lumiblock_ << ":" << event_ << "\n"; //write to file
	if(isprotonRP210 && isprotonRP220 && dileptonMass>88 && dileptonMass<94 &&run_==297101) file_eventlist_kinematics << run_ <<"\t"<< lumiblock_ << "\t" << event_ << "\t" << leadingLeptonPt << "\t"<< leadingLeptonEta << "\t" << leadingLeptonPhi << "\t" << leadingLeptonVz <<"\t" << secondLeptonPt << "\t" << secondLeptonEta << "\t" << secondLeptonPhi << "\t" << secondLeptonVz << "\t" << dileptonPt << "\t" << dileptonEta << "\t"<< dileptonPhi << "\t" << dileptonMass << "\t" << nMuons_looseId << "\t" << xi_rp210_Arm45 << "\t" << xi_rp210_Arm56 << "\t" << xi_rp220_Arm45 << "\t"<< xi_rp220_Arm56 << "\t" << xi_multiArm45 << "\t" << xi_multiArm56 << "\n"; // write to file
      }

    } // proton bracket

  } // end of the event loop

  std::cout << "\n\n<END> good luck!\n" << std::endl;
  fileout->cd();
  tout->Write();
  fileout->Close();
  if(createProtonFile) {
    file_protons_reco.close();
  }
  if(createEventFile){
    file_eventlist.close();
    file_eventlist_kinematics.close();
  }

}

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
  {
    return *itr;
  }
  return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

void randomProtons(TString file_random){

  rndXi_RP210_sec45.clear();
  rndXi_RP210_sec56.clear();
  rndXi_RP220_sec45.clear();
  rndXi_RP220_sec56.clear();
  rndXi_Multi_sec45.clear();
  rndXi_Multi_sec56.clear();

  std::ifstream in;
  in.open(file_random);

  if (!in.good()){
    std::cout << "\n\t --> Random Reco files are not present! Please, produce them first running the option --protonfile\n" << std::endl;
    exit(EXIT_FAILURE);
  }

  Float_t xi_rp210_sec45, xi_rp210_sec56;
  Float_t xi_rp220_sec45, xi_rp220_sec56;
  Float_t xi_multi_sec45, xi_multi_sec56;
  Int_t nlines = 0;

  while (1) {
    in >> xi_rp210_sec45 >> xi_rp210_sec56 >> xi_rp220_sec45 >> xi_rp220_sec56 >> xi_multi_sec45 >> xi_multi_sec56;
    if (!in.good()) break;

    rndXi_RP210_sec45.push_back(xi_rp210_sec45);
    rndXi_RP210_sec56.push_back(xi_rp210_sec56);
    rndXi_RP220_sec45.push_back(xi_rp220_sec45);
    rndXi_RP220_sec56.push_back(xi_rp220_sec56);
    rndXi_Multi_sec45.push_back(xi_multi_sec45);
    rndXi_Multi_sec56.push_back(xi_multi_sec56);

    nlines++;
  }
  printf(" found %d points\n",nlines);

  in.close();

}

bool MissingMassNtupleAnalyzer::MCMatch(double phi1, double eta1, double pt1, double phi2, double eta2, double pt2){

  bool selected = false;
  double dp = std::fabs(phi1 - phi2);
  double deta = std::fabs(eta1 - eta2);
  double dpt = std::fabs(pt1 - pt2);

  std::cout << "dp: " << dp << ", deta: " << deta << ", dpt: " << dpt << std::endl;
  //std::cout << "phi1: " << phi1 << ", eta1: " << eta1 << ", pt1: " << pt1 << ", phi2: " << phi2 << ", eta2: " << eta2 << ", pt2: " << pt2 << std::endl;

  if(dpt<0.2*pt1) selected = true;

  if(selected){
    std::cout << "phi1: " << phi1 << ", eta1: " << eta1 << ", pt1: " << pt1 << ", phi2: " << phi2 << ", eta2: " << eta2 << ", pt2: " << pt2 << std::endl;
    std::cout << "\n\nSelected? " << selected << std::endl;
  }

  return selected;

}

int main(int argc, char * argv[])
{

  std::cout << "\n===================" << std::endl;
  std::cout << "Missing Mass Search" << std::endl;
  std::cout << "===================" << std::endl;

  TTree* tree = NULL;
  if(cmdOptionExists(argv, argv+argc, "--h")||cmdOptionExists(argv, argv+argc, "--help"))
  {
    std::cout << "\n== Help ==\n" << std::endl;
    std::cout << "\t --f filename.root (input file)" << std::endl;
    std::cout << "\t --mode Muon, Electron or Bjets" << std::endl;
    std::cout << "\t --datatype Signal (or DY or SD or Data)" << std::endl;
    std::cout << "\t --xa 120, 130, 140 or 150 (2017 conditions only), needed when running on MC" << std::endl;
    std::cout << "\t --jobid job1 (tag to be added in the outputfile, _option_)" << std::endl;
    std::cout << "\t --outdir Output (Output dir for jobs, _option_)" << std::endl;
    std::cout << "\t --protonfile (it generates a text file with info from protons)" << std::endl;
    std::cout << "\t --eventfile (it generates a text file with run:ls:event_number format)" << std::endl;
    std::cout << "\t --random (performing analysis using random protons)" << std::endl;
    std::cout << "\t --single (fixing arm 45 with random protons.)" << std::endl;
    std::cout << "\t --zerobias (option to run with zerobias triggers. It includes the option --noprotonsfilter)" << std::endl;
    std::cout << "\t --noprotonsfilter (option to run without proton selection)\n" << std::endl;
    std::cout << ">> Important: options \"--protonfile\" and \"--random\" _can not_ be used together! In case you do not have the text file with the proton x[mm] hit, first run --protonfile to create the files and after, run --random.\n\n" << std::endl;
    return 0;
  }

  char * filename = getCmdOption(argv, argv + argc, "--f");
  char * era = getCmdOption(argv, argv + argc, "--era");
  char * xa = getCmdOption(argv, argv + argc, "--xa");
  char * mode = getCmdOption(argv, argv + argc, "--mode");
  char * datatype = getCmdOption(argv, argv + argc, "--datatype");
  char * jobid = getCmdOption(argv, argv + argc, "--jobid");
  char * outdir = getCmdOption(argv, argv + argc, "--outdir");

  if((cmdOptionExists(argv, argv+argc, "--protonfile") || cmdOptionExists(argv, argv+argc, "--random") || cmdOptionExists(argv, argv+argc, "--single"))&&(cmdOptionExists(argv, argv+argc, "--zerobias")||cmdOptionExists(argv, argv+argc, "--noprotonsfilter"))){
    std::cout << "\n\t ---> Options (--protonfile or --random or --single) can not be used together with (--noprotonsfilter and --zerobias) parameter! Please try again!\n" << std::endl;
    return 0;
  }

  if (outdir && strstr(outdir,"--")!=0) {
    std::cout << "\n\t ---> Missing --outdir parameter! Please try again!\n" << std::endl;
    return 0;
  }

  if (jobid && strstr(jobid,"--")!=0) {
    std::cout << "\n\t ---> Missing --jobid parameter! Please try again!\n" << std::endl;
    return 0;
  }

  if (filename && strstr(filename,"--")!=0) {
    std::cout << "\n\t ---> Missing --filename parameter! Please try again!\n" << std::endl;
    return 0;
  }

  if (era && strstr(era,"--")!=0) {
    std::cout << "\n\t ---> Missing --era parameter! Please try again!\n" << std::endl;
    return 0;
  } 

  if (mode && strstr(mode,"--")!=0) {
    std::cout << "\n\t ---> Missing --mode parameter! Please try again!\n" << std::endl;
    return 0;
  }

  if (datatype && strstr(datatype,"--")!=0) {
    std::cout << "\n\t ---> Missing --datatype parameter! Please try again!\n" << std::endl;
    return 0;
  }

  if (xa && strstr(xa,"--")!=0) {
    std::cout << "\n\t ---> Missing --xa parameter! Please try again!\n" << std::endl;
    return 0;
  }

  if(cmdOptionExists(argv, argv+argc, "--protonfile") && cmdOptionExists(argv, argv+argc, "--random")){
    std::cout << "\n\t --> Please, use only --protonfile or --random option.\n" << std::endl;
    return 0;
  }

  if (filename && era && mode)
  {
    TTree* tree = NULL;
    TFile filecheck(filename);
    if(filecheck.IsZombie()){
      std::cout << "\n\t ---> Corrupted file! Please try another file!\n" << std::endl;
      return 0;
    }

    if (!filecheck.GetDirectory("missing_mass")){
      std::cout << "\n---------------------------------------------------" << std::endl;
      std::cout << " There is no directory/path "<< std::endl;
      std::cout << " in the file." << std::endl;
      std::cout << "---------------------------------------------------\n" << std::endl;
      return 0;
    }

    std::cout << "Reading file " << filename << std::endl;

    TFile* file = TFile::Open(filename,"READ");
    tree = (TTree*) file->Get( "missing_mass/analyzer" );

    bool createProtonFile = false;
    bool createEventFile = false;
    bool randomFlag = false;
    bool single = false;
    bool zerobias = false;
    bool protonsfilter = true;

    if(strcmp(era, "A")==0 || strcmp(era, "B")==0 || strcmp(era, "C")==0 || strcmp(era, "D")==0 || strcmp(era, "E")==0 || strcmp(era,"F")==0 || strcmp(era, "preTS2")==0 || strcmp(era, "postTS2")==0 || strcmp(era, "a")==0 || strcmp(era, "b")==0 || strcmp(era, "c")==0 || strcmp(era, "d")==0 || strcmp(era, "e")==0 || strcmp(era,"f")==0 || strcmp(era, "pre TS2")==0 || strcmp(era, "post TS2")==0){}
    else{
      std::cout << "\nNo --era option! It should be A, B, C, D, E, F, preTS2 or postTS2.\n" << std::endl;
      return 0;
    }

    if(cmdOptionExists(argv, argv+argc, "--protonfile")){
      createProtonFile = true;
    }

    if(cmdOptionExists(argv, argv+argc, "--eventfile")){
      createEventFile = true;
    }

    if(cmdOptionExists(argv, argv+argc, "--single")) single = true;
    if(cmdOptionExists(argv, argv+argc, "--zerobias")){
      zerobias = true;
    }
    if(cmdOptionExists(argv, argv+argc, "--noprotonsfilter")) protonsfilter = false;

    if(cmdOptionExists(argv, argv+argc, "--random")){
      randomFlag = true;
      if(strcmp(mode, "MC_Muon")==0){
	if(strcmp(era, "B")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_B_xa_120.txt");
	}
	else if(strcmp(era, "B")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_B_xa_130.txt");
	}
	else if(strcmp(era, "B")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_B_xa_140.txt"); 
	}
	else if(strcmp(era, "B")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_B_xa_150.txt"); 
	}
	else if(strcmp(era, "C")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_C_xa_120.txt");
	}
	else if(strcmp(era, "C")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_C_xa_130.txt");
	}
	else if(strcmp(era, "C")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_C_xa_140.txt");
	}
	else if(strcmp(era, "C")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_C_xa_150.txt");
	}
	else if(strcmp(era, "D")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_D_xa_120.txt");
	}
	else if(strcmp(era, "D")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_D_xa_130.txt");
	}
	else if(strcmp(era, "D")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_D_xa_140.txt");
	}
	else if(strcmp(era, "D")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_D_xa_150.txt");
	}
	else if(strcmp(era, "preTS2")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_preTS2_xa_120.txt");
	}
	else if(strcmp(era, "preTS2")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_preTS2_xa_130.txt");
	}
	else if(strcmp(era, "preTS2")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_preTS2_xa_140.txt");
	}
	else if(strcmp(era, "preTS2")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_preTS2_xa_150.txt");
	}
	else if(strcmp(era, "E")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_E_xa_120.txt");
	}
	else if(strcmp(era, "E")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_E_xa_130.txt");      
	}
	else if(strcmp(era, "E")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_E_xa_140.txt");   
	}
	else if(strcmp(era, "E")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_E_xa_150.txt"); 
	}
	else if(strcmp(era, "F")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_F_xa_120.txt"); 
	}
	else if(strcmp(era, "F")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_F_xa_130.txt");
	}
	else if(strcmp(era, "F")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_F_xa_140.txt");      
	}
	else if(strcmp(era, "F")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_F_xa_150.txt"); 
	}
	else if(strcmp(era, "postTS2")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_postTS2_xa_120.txt");
	}
	else if(strcmp(era, "postTS2")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_postTS2_xa_130.txt");
	}
	else if(strcmp(era, "postTS2")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_postTS2_xa_140.txt");
	}
	else if(strcmp(era, "postTS2")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Muon_era_postTS2_xa_150.txt");
	}
	else{
	  std::cout << "Era or X-angle not present in 2017 dataset!" << std::endl;
	}
      }

      else if(strcmp(mode, "MC_Electron")==0){ 
	if(strcmp(era, "B")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_B_xa_120.txt"); 
	}
	else if(strcmp(era, "B")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_B_xa_130.txt");
	}
	else if(strcmp(era, "B")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_B_xa_140.txt"); 
	}
	else if(strcmp(era, "B")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_B_xa_150.txt"); 
	}
	else if(strcmp(era, "C")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_C_xa_120.txt");
	}
	else if(strcmp(era, "C")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_C_xa_130.txt");
	}
	else if(strcmp(era, "C")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_C_xa_140.txt");
	}
	else if(strcmp(era, "C")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_C_xa_150.txt");
	}
	else if(strcmp(era, "D")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_D_xa_120.txt");
	}
	else if(strcmp(era, "D")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_D_xa_130.txt");
	}
	else if(strcmp(era, "D")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_D_xa_140.txt");
	}
	else if(strcmp(era, "D")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_D_xa_150.txt");
	}
	else if(strcmp(era, "preTS2")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_preTS2_xa_120.txt");
	}
	else if(strcmp(era, "preTS2")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_preTS2_xa_130.txt");
	}
	else if(strcmp(era, "preTS2")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_preTS2_xa_140.txt");
	}
	else if(strcmp(era, "preTS2")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_preTS2_xa_150.txt");
	}
	else if(strcmp(era, "E")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_E_xa_120.txt");
	}
	else if(strcmp(era, "E")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_E_xa_130.txt");      
	}
	else if(strcmp(era, "E")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_E_xa_140.txt");   
	}
	else if(strcmp(era, "E")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_E_xa_150.txt"); 
	}
	else if(strcmp(era, "F")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_F_xa_120.txt"); 
	}
	else if(strcmp(era, "F")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_F_xa_130.txt");
	}
	else if(strcmp(era, "F")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_F_xa_140.txt");      
	}
	else if(strcmp(era, "F")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_F_xa_150.txt"); 
	}
	else if(strcmp(era, "postTS2")==0 && strcmp(xa, "120")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_postTS2_xa_120.txt");
	}
	else if(strcmp(era, "postTS2")==0 && strcmp(xa, "130")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_postTS2_xa_130.txt");
	}
	else if(strcmp(era, "postTS2")==0 && strcmp(xa, "140")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_postTS2_xa_140.txt");
	}
	else if(strcmp(era, "postTS2")==0 && strcmp(xa, "150")==0){
	  randomProtons("proton_reco_rphorizontal_Electron_era_postTS2_xa_150.txt");
	}
	else{
	  std::cout << "Era or X-angle not present in 2017 dataset!" << std::endl;
	}
      }

      //Not correlated events, for muons the list of protons in electrons samples have been taken.
      else if(strcmp(mode, "Muon")==0||strcmp(mode, "muon")==0){
	if(strcmp(era, "A")==0 || strcmp(era, "a")==0){
	  randomProtons("proton_reco_muona.txt");
	}
	else if(strcmp(era, "B")==0 || strcmp(era, "b")==0){
	  randomProtons("proton_reco_muonb.txt");
	}
	else if(strcmp(era, "C")==0 || strcmp(era, "c")==0){
	  randomProtons("proton_reco_muonc.txt");
	}
	else if(strcmp(era, "D")==0 || strcmp(era, "d")==0){
	  randomProtons("proton_reco_muond.txt");
	}
	else if(strcmp(era, "E")==0 || strcmp(era, "e")==0){
	  randomProtons("proton_reco_muone.txt");
	}
	else if(strcmp(era, "F")==0 || strcmp(era, "f")==0){
	  randomProtons("proton_reco_muonf.txt");
	}
	else{
	  std::cout << "\n\t --> Please, insert --era B, C, D, E or F\n" << std::endl;
	  return 0;
	} 
      }

      //Not correlated events, for muons the list of protons in electrons samples have been taken.
      else if(strcmp(mode, "Bjets")==0||strcmp(mode, "bjets")==0){
	if(strcmp(era, "A")==0 || strcmp(era, "a")==0){
	  randomProtons("proton_reco_bjetsa.txt");
	}
	else if(strcmp(era, "B")==0 || strcmp(era, "b")==0){
	  randomProtons("proton_reco_bjetsb.txt");
	}
	else if(strcmp(era, "C")==0 || strcmp(era, "c")==0){
	  randomProtons("proton_reco_bjetsc.txt");
	}
	else if(strcmp(era, "D")==0 || strcmp(era, "d")==0){
	  randomProtons("proton_reco_bjetsd.txt");
	}
	else if(strcmp(era, "E")==0 || strcmp(era, "e")==0){
	  randomProtons("proton_reco_bjetse.txt");
	}
	else if(strcmp(era, "F")==0 || strcmp(era, "f")==0){
	  randomProtons("proton_reco_bjetsf.txt");
	}
	else{
	  std::cout << "\n\t --> Please, insert --era B, C, D, E or F\n" << std::endl;
	  return 0;
	}
      }

      //Not correlated events, for electrons the list of protons in the muons samples have been taken.
      else if(strcmp(mode, "Electron")==0||strcmp(mode, "electron")==0){
	if(strcmp(era, "A")==0 || strcmp(era, "a")==0){
	  randomProtons("proton_reco_electrona.txt");
	}
	else if(strcmp(era, "B")==0 || strcmp(era, "b")==0){
	  randomProtons("proton_reco_electronb.txt");
	}
	else if(strcmp(era, "C")==0 || strcmp(era, "c")==0){
	  randomProtons("proton_reco_electronc.txt");
	}
	else if(strcmp(era, "D")==0 || strcmp(era, "d")==0){
	  randomProtons("proton_reco_electrond.txt");
	}
	else if(strcmp(era, "E")==0 || strcmp(era, "e")==0){
	  randomProtons("proton_reco_electrone.txt");
	}
	else if(strcmp(era, "F")==0 || strcmp(era, "f")==0){
	  randomProtons("proton_reco_electronf.txt");
	}
	else{
	  std::cout << "\n\t --> Please, insert --era C or D\n" << std::endl;
	  return 0;
	} 
      }
      else{
	std::cout << "\n\t --> Please, insert --mode Muon, Electron, BJets, MC_Muon, MC_Electron\n" << std::endl;
	return 0;
      }
    }

    // Accessing Missing Mass Object
    MissingMassNtupleAnalyzer m(tree); 
    m.Loop(era, mode, xa, jobid, outdir, datatype, createProtonFile, randomFlag, single, zerobias, protonsfilter, createEventFile);
  }else{
    std::cout << "\n\t --> Please, insert --f filename.root and --era A, B, C, D, E and F --mode Muon (or Electron)\n" << std::endl;
    return 0;
  }

  return 0;
}

