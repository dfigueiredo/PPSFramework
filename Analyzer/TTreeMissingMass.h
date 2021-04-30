#include <iostream>
#include <list>
#include <TH1.h>
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>

class TTreeMissingMass {

  private:
    TTree* out;
    TFile* fileout;

  public:
    virtual ~TTreeMissingMass() {}
    void CreateTTree();
    void Fill();
    void Storing();
    void Clearing();

    bool switchZeroBias = false;
    bool switchMC = false;
    bool switchMuon = false;
    bool switchElectron = false;
    bool switchBjets = false;
    bool switchElectronMuon = false;

    int xangle = 0;
    int run = 0;
    int event = 0;
    int lumiblock = 0;
    int era = 0;
    int PUInterac = 0;
    int PUTrueInterac = 0;

    bool is45PUproton = false;
    bool is56PUproton = false;
    bool is2PUproton = false;

    int prescalesL1ZeroBias = 0;
    int prescalesL1ZeroBiasAfterTrain = 0;
    int prescalesL1ZeroBiasIsolatedBx = 0;
    int prescalesL1ZeroBiasAlignment = 0;
    int prescalesL1ZeroBiasBeamSpot = 0;
    int prescalesL1ZeroBiasUnpairedBptxMinus = 0;
    int prescalesL1ZeroBiasUnpairedBptxPlus = 0;
    int prescalesL1Physics = 0;
    int prescalesL1SingleMu = 0;

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
    int nGenParticles = 0;
    int nGenElectrons = 0;
    int nGenMuons = 0;
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

    double genmissEt = 0.;
    double genmissEt_phi = 0.;

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
    int nElectrons = 0;
    int nMuons = 0;
    int nElectrons_looseId = 0;
    int nElectrons_mediumId = 0;
    int nElectrons_tightId = 0;
    int nMuons_looseId = 0;
    int nMuons_mediumId = 0;
    int nMuons_tightId = 0;
    int nChargedPFMultiPV_Loose = 0;
    int nChargedPFMultiPV_Tight = 0;
    int nChargedPFMultiPV_UsedInFit = 0;
    int nChargedPFMultiPV_Tight_Fit = 0;
    double SumChargedPFMultiPV_pt_Loose = 0;
    double SumChargedPFMultiPV_pt_Tight = 0;
    double SumChargedPFMultiPV_pt_UsedInFit = 0;
    double SumChargedPFMultiPV_pt_Tight_Fit = 0;

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

    double missEt = 0.;
    double missEt_phi = 0.;
    double dphi = 0.;

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

};

void TTreeMissingMass::Clearing(){

  xangle = 0;
  run = 0;
  event = 0;
  lumiblock = 0;
  era = -1;

  prescalesL1ZeroBias = 1;
  prescalesL1ZeroBiasAfterTrain = 1;
  prescalesL1ZeroBiasIsolatedBx = 1;
  prescalesL1ZeroBiasAlignment = 1;
  prescalesL1ZeroBiasBeamSpot = 1;
  prescalesL1ZeroBiasUnpairedBptxMinus = 1;
  prescalesL1ZeroBiasUnpairedBptxPlus = 1;
  prescalesL1Physics = 1;
  prescalesL1SingleMu = 1;

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

  PUInterac = 0;
  PUTrueInterac = 0;
  is45PUproton = false;
  is56PUproton = false;
  is2PUproton = false;
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

  genmissEt = 0.;
  genmissEt_phi = 0.;

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

  missEt = 0.;
  missEt_phi = 0.;
  dphi = 0.;

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

}

void TTreeMissingMass::CreateTTree(){

  fileout = new TFile("testfile.root", "RECREATE");
  out = new TTree("Events", "Events");
  out->Branch("era",&era,"era/I");
  out->Branch("run",&run,"run/I");
  out->Branch("event",&event,"event/I");
  out->Branch("lumiblock",&lumiblock,"lumiblock/I");
  out->Branch("xangle",&xangle,"xangle/I");

  if(switchZeroBias){
    out->Branch("triggerZeroBias",&triggerZeroBias,"triggerZeroBias/B");
    out->Branch("triggerZeroBiasAfterTrain",&triggerZeroBiasAfterTrain,"triggerZeroBiasAfterTrain/B");
    out->Branch("triggerZeroBiasIsolatedBx",&triggerZeroBiasIsolatedBx,"triggerZeroBiasIsolatedBx/B");
    out->Branch("triggerZeroBiasAlignment",&triggerZeroBiasAlignment,"triggerZeroBiasAlignment/B");
    out->Branch("triggerZeroBiasBeamSpot",&triggerZeroBiasBeamSpot,"triggerZeroBiasBeamSpot/B");
    out->Branch("triggerZeroBiasUnpairedBptxMinus",&triggerZeroBiasUnpairedBptxMinus,"triggerZeroBiasUnpairedBptxMinus/B");
    out->Branch("triggerZeroBiasUnpairedBptxPlus",&triggerZeroBiasUnpairedBptxPlus,"triggerZeroBiasUnpairedBptxPlus/B");
    out->Branch("triggerPhysics",&triggerPhysics,"triggerPhysics/B");
    out->Branch("triggerL1SingleMu",&triggerL1SingleMu,"triggerL1SingleMu/B");
    out->Branch("prescalesL1ZeroBias",&prescalesL1ZeroBias,"prescalesL1ZeroBias/I");
    out->Branch("prescalesL1ZeroBiasAfterTrain",&prescalesL1ZeroBiasAfterTrain,"prescalesL1ZeroBiasAfterTrain/I");
    out->Branch("prescalesL1ZeroBiasIsolatedBx",&prescalesL1ZeroBiasIsolatedBx,"prescalesL1ZeroBiasIsolatedBx/I");
    out->Branch("prescalesL1ZeroBiasAlignment",&prescalesL1ZeroBiasAlignment,"prescalesL1ZeroBiasAlignment/I");
    out->Branch("prescalesL1ZeroBiasBeamSpot",&prescalesL1ZeroBiasBeamSpot,"prescalesL1ZeroBeamSpot/I");
    out->Branch("prescalesL1ZeroBiasUnpairedBptxMinus",&prescalesL1ZeroBiasUnpairedBptxMinus,"prescalesL1ZeroBeamUnpairedBptxMinus/I");
    out->Branch("prescalesL1ZeroBiasUnpairedBptxPlus",&prescalesL1ZeroBiasUnpairedBptxPlus,"prescalesL1ZeroBeamUnpairedBptxPlus/I");
    out->Branch("prescalesL1Physics",&prescalesL1Physics,"prescalesL1Physics/I");
    out->Branch("prescalesL1SingleMu",&prescalesL1SingleMu,"prescalesL1SingleMu/I");
  }
  if(switchMuon){
    out->Branch("triggerIsoMu27",&triggerIsoMu27,"triggerIsoMu27/B");
    out->Branch("triggerMu17TrkIsoMu8TrkIso",&triggerMu17TrkIsoMu8TrkIso,"triggerMu17TrkIsoMu8TrkIso/B");
    out->Branch("triggerMu17TrkIsoMu8TrkIsoMass3",&triggerMu17TrkIsoMu8TrkIsoMass3,"triggerMu17TrkIsoMu8TrkIsoMass3/B");
    out->Branch("triggerMu17TrkIsoMu8TrkIsoMass8",&triggerMu17TrkIsoMu8TrkIsoMass8,"triggerMu17TrkIsoMu8TrkIsoMass8/B");
    out->Branch("triggerDoubleMu43",&triggerDoubleMu43,"triggerDoubleMu43/B");
    out->Branch("triggerDoubleMu48",&triggerDoubleMu48,"triggerDoubleMu48/B");
  }
  if(switchElectron){
    out->Branch("triggerEle27",&triggerEle27,"triggerEle27/B");
    out->Branch("triggerEle23Ele12",&triggerEle23Ele12,"triggerEle23Ele12/B");
    out->Branch("triggerEle23Ele12Dz",&triggerEle23Ele12Dz,"triggerEle23Ele12Dz/B");
    out->Branch("triggerDoubleEle33",&triggerDoubleEle33,"triggerDoubleEle33/B");
    out->Branch("triggerBTagMu5Ak4dijet20",&triggerBTagMu5Ak4dijet20,"triggerBTagMu5Ak4dijet20/B");
    out->Branch("triggerBTagMu5Ak4dijet40",&triggerBTagMu5Ak4dijet40,"triggerBTagMu5Ak4dijet40/B");
    out->Branch("triggerBTagMu5Ak4dijet70",&triggerBTagMu5Ak4dijet70,"triggerBTagMu5Ak4dijet70/B");
    out->Branch("triggerBTagMu5Ak4dijet110",&triggerBTagMu5Ak4dijet110,"triggerBTagMu5Ak4dijet110/B");
    out->Branch("triggerBTagMu5Ak4dijet170",&triggerBTagMu5Ak4dijet170,"triggerBTagMu5Ak4dijet170/B");
    out->Branch("triggerBTagMu5Ak4dijet300",&triggerBTagMu5Ak4dijet300,"triggerBTagMu5Ak4dijet300/B");
    out->Branch("triggerBTagMu5Ak8dijet170",&triggerBTagMu5Ak8dijet170,"triggerBTagMu5Ak8dijet170/B");
    out->Branch("triggerBTagMu5Ak8dijet300",&triggerBTagMu5Ak8dijet300,"triggerBTagMu5Ak8dijet300/B");
    out->Branch("triggerPFHT380DoubleBTag",&triggerPFHT380DoubleBTag,"triggerPFHT380DoubleBTag/B");
    out->Branch("triggerPFHT430DoubleBTag",&triggerPFHT430DoubleBTag,"triggerPFHT430DoubleBTag/B");
    out->Branch("triggerPFMET100BTag",&triggerPFMET100BTag,"triggerPFMET100BTag/B");
    out->Branch("triggerPFMET110BTag",&triggerPFMET110BTag,"triggerPFMET110BTag/B");
    out->Branch("triggerPFMET120BTag",&triggerPFMET120BTag,"triggerPFMET120BTag/B");
    out->Branch("triggerPFMET130BTag",&triggerPFMET130BTag,"triggerPFMET130BTag/B");
    out->Branch("triggerPFMET140BTag",&triggerPFMET140BTag,"triggerPFMET140BTag/B");
    out->Branch("triggerEle15PFHT450BTag",&triggerEle15PFHT450BTag,"triggerEle15PFHT450BTag/B");
    out->Branch("triggerHT250BTagScouting",&triggerHT250BTagScouting,"triggerHT250BTagScouting/B");
  }
  if(switchBjets){
    out->Branch("triggerIsoMu27",&triggerIsoMu27,"triggerIsoMu27/B");
    out->Branch("triggerMu17TrkIsoMu8TrkIso",&triggerMu17TrkIsoMu8TrkIso,"triggerMu17TrkIsoMu8TrkIso/B");
    out->Branch("triggerMu17TrkIsoMu8TrkIsoMass3",&triggerMu17TrkIsoMu8TrkIsoMass3,"triggerMu17TrkIsoMu8TrkIsoMass3/B");
    out->Branch("triggerMu17TrkIsoMu8TrkIsoMass8",&triggerMu17TrkIsoMu8TrkIsoMass8,"triggerMu17TrkIsoMu8TrkIsoMass8/B");
    out->Branch("triggerDoubleMu43",&triggerDoubleMu43,"triggerDoubleMu43/B");
    out->Branch("triggerDoubleMu48",&triggerDoubleMu48,"triggerDoubleMu48/B");
    out->Branch("triggerBTagMu5Ak4dijet20",&triggerBTagMu5Ak4dijet20,"triggerBTagMu5Ak4dijet20/B");
    out->Branch("triggerBTagMu5Ak4dijet40",&triggerBTagMu5Ak4dijet40,"triggerBTagMu5Ak4dijet40/B");
    out->Branch("triggerBTagMu5Ak4dijet70",&triggerBTagMu5Ak4dijet70,"triggerBTagMu5Ak4dijet70/B");
    out->Branch("triggerBTagMu5Ak4dijet110",&triggerBTagMu5Ak4dijet110,"triggerBTagMu5Ak4dijet110/B");
    out->Branch("triggerBTagMu5Ak4dijet170",&triggerBTagMu5Ak4dijet170,"triggerBTagMu5Ak4dijet170/B");
    out->Branch("triggerBTagMu5Ak4dijet300",&triggerBTagMu5Ak4dijet300,"triggerBTagMu5Ak4dijet300/B");
    out->Branch("triggerBTagMu5Ak8dijet170",&triggerBTagMu5Ak8dijet170,"triggerBTagMu5Ak8dijet170/B");
    out->Branch("triggerBTagMu5Ak8dijet300",&triggerBTagMu5Ak8dijet300,"triggerBTagMu5Ak8dijet300/B");
    out->Branch("triggerPFHT380DoubleBTag",&triggerPFHT380DoubleBTag,"triggerPFHT380DoubleBTag/B");
    out->Branch("triggerPFHT430DoubleBTag",&triggerPFHT430DoubleBTag,"triggerPFHT430DoubleBTag/B");
  }
  if(switchElectronMuon){
    out->Branch("triggerMu23TrkIsoEle12",&triggerMu23TrkIsoEle12,"triggerMu23TrkIsoEle12/B");
    out->Branch("triggerMu23TrkIsoEle12DZ",&triggerMu23TrkIsoEle12DZ,"triggerMu23TrkIsoEle12DZ/B");
    out->Branch("triggerMu8TrkIsoEle23DZ",&triggerMu8TrkIsoEle23DZ,"triggerMu8TrkIsoEle23DZ/B");
  }

  // In miniAOD this is empty (slimmed Vertex do not have the track ref.)
  out->Branch("nTracksPerVertex",&nTracksPerVertex,"nTracksPerVertex/I");
  out->Branch("nVertex",&nVertex,"nVertex/I");
  out->Branch("pvVertexX",&pvVertexX,"pvVertexX/D");
  out->Branch("pvVertexY",&pvVertexY,"pvVertexY/D");
  out->Branch("pvVertexZ",&pvVertexZ,"pvVertexZ/D");
  out->Branch("pvVertexR",&pvVertexR,"pvVertexR/D");
  out->Branch("pvVertexPhi",&pvVertexPhi,"pvVertexPhi/D");
  out->Branch("allVertexZ",&allVertexZ);
  if(switchMC){
    out->Branch("PUInterac",&PUInterac,"PUInterac/I");
    out->Branch("PUTrueInterac",&PUTrueInterac,"PUTrueInterac/I");
    out->Branch("is45PUproton",&is45PUproton,"is45PUproton/B");
    out->Branch("is56PUproton",&is56PUproton,"is56PUproton/B");
    out->Branch("is2PUproton",&is2PUproton,"is2PUproton/B");
    out->Branch("nGenLeptons",&nGenLeptons,"nGenLeptons/I");
    out->Branch("nGenParticles",&nGenParticles,"nGenParticles/I");
    out->Branch("nGenElectrons",&nGenElectrons,"nGenElectrons/I");
    out->Branch("nGenMuons",&nGenMuons,"nGenMuons/I");
    out->Branch("genleadingLeptonPDGId",&genleadingLeptonPDGId,"genleadingLeptonPDGId/I");
    out->Branch("genleadingLeptonEnergy",&genleadingLeptonEnergy,"genleadingLeptonEnergy/D");
    out->Branch("genleadingLeptonPx",&genleadingLeptonPx,"genleadingLeptonPx/D");
    out->Branch("genleadingLeptonPy",&genleadingLeptonPy,"genleadingLeptonPy/D");
    out->Branch("genleadingLeptonPz",&genleadingLeptonPz,"genleadingLeptonPz/D");
    out->Branch("genleadingLeptonPt",&genleadingLeptonPt,"genleadingLeptonPt/D");
    out->Branch("genleadingLeptonEta",&genleadingLeptonEta,"genleadingLeptonEta/D");
    out->Branch("genleadingLeptonPhi",&genleadingLeptonPhi,"genleadingLeptonPhi/D");
    out->Branch("genleadingLeptonCharge",&genleadingLeptonCharge,"genleadingLeptonCharge/I");
    out->Branch("genleadingLeptonVx",&genleadingLeptonVx,"genleadingLeptonVx/D");
    out->Branch("genleadingLeptonVy",&genleadingLeptonVy,"genleadingLeptonVy/D");
    out->Branch("genleadingLeptonVz",&genleadingLeptonVz,"genleadingLeptonVz/D");
    out->Branch("genleadingLeptonVr",&genleadingLeptonVr,"genleadingLeptonVr/D");
    out->Branch("genleadingLeptonVphi",&genleadingLeptonVphi,"genleadingLeptonVphi/D");
    out->Branch("gensecondLeptonPDGId",&gensecondLeptonPDGId,"gensecondLeptonPDGId/I");
    out->Branch("gensecondLeptonEnergy",&gensecondLeptonEnergy,"gensecondLeptonEnergy/D");
    out->Branch("gensecondLeptonPx",&gensecondLeptonPx,"gensecondLeptonPx/D");
    out->Branch("gensecondLeptonPy",&gensecondLeptonPy,"gensecondLeptonPy/D");
    out->Branch("gensecondLeptonPz",&gensecondLeptonPz,"gensecondLeptonPz/D");
    out->Branch("gensecondLeptonPt",&gensecondLeptonPt,"gensecondLeptonPt/D");
    out->Branch("gensecondLeptonEta",&gensecondLeptonEta,"gensecondLeptonEta/D");
    out->Branch("gensecondLeptonPhi",&gensecondLeptonPhi,"gensecondLeptonPhi/D");
    out->Branch("gensecondLeptonCharge",&gensecondLeptonCharge,"gensecondLeptonCharge/I");
    out->Branch("gensecondLeptonVx",&gensecondLeptonVx,"gensecondLeptonVx/D");
    out->Branch("gensecondLeptonVy",&gensecondLeptonVy,"gensecondLeptonVy/D");
    out->Branch("gensecondLeptonVz",&gensecondLeptonVz,"gensecondLeptonVz/D");
    out->Branch("gensecondLeptonVr",&gensecondLeptonVr,"gensecondLeptonVr/D");
    out->Branch("gensecondLeptonVphi",&gensecondLeptonVphi,"gensecondLeptonVphi/D");
    out->Branch("genmissEt",&genmissEt,"genmissEt/D");
    out->Branch("genmissEt_phi",&genmissEt_phi,"genmissEt_phi/D");
    out->Branch("nGenJets",&nGenJets,"nGenJets/I");
    out->Branch("genleadingJetEnergy",&genleadingJetEnergy,"genleadingJetEnergy/D");
    out->Branch("genleadingJetPx",&genleadingJetPx,"genleadingJetPx/D");
    out->Branch("genleadingJetPy",&genleadingJetPy,"genleadingJetPy/D");
    out->Branch("genleadingJetPz",&genleadingJetPz,"genleadingJetPz/D");
    out->Branch("genleadingJetPt",&genleadingJetPt,"genleadingJetPt/D");
    out->Branch("genleadingJetEta",&genleadingJetEta,"genleadingJetEta/D");
    out->Branch("genleadingJetPhi",&genleadingJetPhi,"genleadingJetPhi/D");
    out->Branch("genleadingJetVz",&genleadingJetVz,"genleadingJetVz/D");
    out->Branch("gensecondJetEnergy",&gensecondJetEnergy,"gensecondJetEnergy/D");
    out->Branch("gensecondJetPx",&gensecondJetPx,"gensecondJetPx/D");
    out->Branch("gensecondJetPy",&gensecondJetPy,"gensecondJetPy/D");
    out->Branch("gensecondJetPz",&gensecondJetPz,"gensecondJetPz/D");
    out->Branch("gensecondJetPt",&gensecondJetPt,"gensecondJetPt/D");
    out->Branch("gensecondJetEta",&gensecondJetEta,"gensecondJetEta/D");
    out->Branch("gensecondJetPhi",&gensecondJetPhi,"gensecondJetPhi/D");
    out->Branch("gensecondJetVz",&gensecondJetVz,"gensecondJetVz/D");
    out->Branch("gendileptonCharge",&gendileptonCharge,"gendileptonCharge/I");
    out->Branch("gendileptonMass",&gendileptonMass,"gendileptonMass/D");
    out->Branch("gendileptonEta",&gendileptonEta,"gendileptonEta/D");
    out->Branch("gendileptonPhi",&gendileptonPhi,"gendileptonPhi/D");
    out->Branch("gendileptonPt",&gendileptonPt,"gendileptonPt/D");
    out->Branch("gendileptonRapidity",&gendileptonRapidity,"gendileptonRapidity/D");
    out->Branch("gendijetMass",&gendijetMass,"gendijetMass/D");
    out->Branch("gendijetEta",&gendijetEta,"gendijetEta/D");
    out->Branch("gendijetPhi",&gendijetPhi,"gendijetPhi/D");
    out->Branch("gendijetPt",&gendijetPt,"gendijetPt/D");
    out->Branch("gendijetRapidity",&gendijetRapidity,"gendijetRapidity/D");
  }
  out->Branch("nLeptons",&nLeptons,"nLeptons/I");
  out->Branch("nElectrons",&nElectrons,"nElectrons/I");
  out->Branch("nElectronsLooseId",&nElectrons_looseId,"nElectronsLooseId/I");
  out->Branch("nElectronsMediumId",&nElectrons_mediumId,"nElectronsMediumId/I");
  out->Branch("nElectronsTightId",&nElectrons_tightId,"nElectronsTightId/I");
  out->Branch("nMuons",&nMuons,"nMuons/I");
  out->Branch("nMuonsLooseId",&nMuons_looseId,"nMuonsLooseId/I");
  out->Branch("nMuonsMediumId",&nMuons_mediumId,"nMuonsMediumId/I");
  out->Branch("nMuonsTightId",&nMuons_tightId,"nMuonsTightId/I");
  out->Branch("nChargedPFMultiPV_Loose",&nChargedPFMultiPV_Loose,"nChargedPFMultiPV_Loose/I");
  out->Branch("nChargedPFMultiPV_Tight",&nChargedPFMultiPV_Tight,"nChargedPFMultiPV_Tight/I");
  out->Branch("nChargedPFMultiPV_UsedInFit",&nChargedPFMultiPV_UsedInFit,"nChargedPFMultiPV_UsedInFit/I");
  out->Branch("nChargedPFMultiPV_Tight_Fit",&nChargedPFMultiPV_Tight_Fit,"nChargedPFMultiPV_Tight_Fit/I");
  out->Branch("SumChargedPFMultiPV_pt_Loose",&SumChargedPFMultiPV_pt_Loose,"SumChargedPFMultiPV_pt_Loose/D");
  out->Branch("SumChargedPFMultiPV_pt_Tight",&SumChargedPFMultiPV_pt_Tight,"SumChargedPFMultiPV_pt_Tight/D");
  out->Branch("SumChargedPFMultiPV_pt_UsedInFit",&SumChargedPFMultiPV_pt_UsedInFit,"SumChargedPFMultiPV_pt_UsedInFit/D");
  out->Branch("SumChargedPFMultiPV_pt_Tight_Fit",&SumChargedPFMultiPV_pt_Tight_Fit,"SumChargedPFMultiPV_pt_Tight_Fit/D");
  out->Branch("leadingLeptonPDGId",&leadingLeptonPDGId,"leadingLeptonPDGId/I");
  out->Branch("leadingLeptonEnergy",&leadingLeptonEnergy,"leadingLeptonEnergy/D");
  out->Branch("leadingLeptonPx",&leadingLeptonPx,"leadingLeptonPx/D");
  out->Branch("leadingLeptonPy",&leadingLeptonPy,"leadingLeptonPy/D");
  out->Branch("leadingLeptonPz",&leadingLeptonPz,"leadingLeptonPz/D");
  out->Branch("leadingLeptonPt",&leadingLeptonPt,"leadingLeptonPt/D");
  out->Branch("leadingLeptonEta",&leadingLeptonEta,"leadingLeptonEta/D");
  out->Branch("leadingLeptonPhi",&leadingLeptonPhi,"leadingLeptonPhi/D");
  out->Branch("leadingLeptonCharge",&leadingLeptonCharge,"leadingLeptonCharge/I");
  out->Branch("leadingLeptonVx",&leadingLeptonVx,"leadingLeptonVx/D");
  out->Branch("leadingLeptonVy",&leadingLeptonVy,"leadingLeptonVy/D");
  out->Branch("leadingLeptonVz",&leadingLeptonVz,"leadingLeptonVz/D");
  out->Branch("leadingLeptonVr",&leadingLeptonVr,"leadingLeptonVr/D");
  out->Branch("leadingLeptonVphi",&leadingLeptonVphi,"leadingLeptonVphi/D");
  out->Branch("leadingLeptonPFIso",&leadingLeptonPFIso,"leadingLeptonPFIso/D");
  out->Branch("leadingLeptonTkIso",&leadingLeptonTkIso,"leadingLeptonTkIso/D");
  out->Branch("leadingLeptonLooseId",&leadingLeptonLooseId,"leadingLeptonLooseId/B");
  out->Branch("leadingLeptonMediumId",&leadingLeptonMediumId,"leadingLeptonMediumId/B");
  out->Branch("leadingLeptonTightId",&leadingLeptonTightId,"leadingLeptonTightId/B");
  out->Branch("leadingLeptonPfIsoMedium",&leadingLeptonPfIsoMedium,"leadingLeptonPfIsoMedium/B");
  out->Branch("leadingLeptonMiniIsoTight",&leadingLeptonMiniIsoTight,"leadingLeptonMiniIsoTight/B");
  out->Branch("leadingLeptonPfIsoVeryTight",&leadingLeptonPfIsoVeryTight,"leadingLeptonPfIsoVeryTight/B");
  out->Branch("secondLeptonPDGId",&secondLeptonPDGId,"secondLeptonPDGId/I");
  out->Branch("secondLeptonEnergy",&secondLeptonEnergy,"secondLeptonEnergy/D");
  out->Branch("secondLeptonPx",&secondLeptonPx,"secondLeptonPx/D");
  out->Branch("secondLeptonPy",&secondLeptonPy,"secondLeptonPy/D");
  out->Branch("secondLeptonPz",&secondLeptonPz,"secondLeptonPz/D");
  out->Branch("secondLeptonPt",&secondLeptonPt,"secondLeptonPt/D");
  out->Branch("secondLeptonEta",&secondLeptonEta,"secondLeptonEta/D");
  out->Branch("secondLeptonPhi",&secondLeptonPhi,"secondLeptonPhi/D");
  out->Branch("secondLeptonCharge",&secondLeptonCharge,"secondLeptonCharge/I");
  out->Branch("secondLeptonVx",&secondLeptonVx,"secondLeptonVx/D");
  out->Branch("secondLeptonVy",&secondLeptonVy,"secondLeptonVy/D");
  out->Branch("secondLeptonVz",&secondLeptonVz,"secondLeptonVz/D");
  out->Branch("secondLeptonVr",&secondLeptonVr,"secondLeptonVr/D");
  out->Branch("secondLeptonVphi",&secondLeptonVphi,"secondLeptonVphi/D");
  out->Branch("secondLeptonPFIso",&secondLeptonPFIso,"secondLeptonPFIso/D");
  out->Branch("secondLeptonTkIso",&secondLeptonTkIso,"secondLeptonTkIso/D");
  out->Branch("secondLeptonLooseId",&secondLeptonLooseId,"secondLeptonLooseId/B");
  out->Branch("secondLeptonMediumId",&secondLeptonMediumId,"secondLeptonMediumId/B");
  out->Branch("secondLeptonTightId",&secondLeptonTightId,"secondLeptonTightId/B");
  out->Branch("secondLeptonPfIsoMedium",&secondLeptonPfIsoMedium,"secondLeptonPfIsoMedium/B");
  out->Branch("secondLeptonMiniIsoTight",&secondLeptonMiniIsoTight,"secondLeptonMiniIsoTight/B");
  out->Branch("secondLeptonPfIsoVeryTight",&secondLeptonPfIsoVeryTight,"secondLeptonPfIsoVeryTight/B");
  out->Branch("acoplanarity",&acoplanarity,"acoplanarity/D");
  out->Branch("nJets",&nJets,"nJets/I");
  out->Branch("leadingJetEnergy",&leadingJetEnergy,"leadingJetEnergy/D");
  out->Branch("leadingJetPx",&leadingJetPx,"leadingJetPx/D");
  out->Branch("leadingJetPy",&leadingJetPy,"leadingJetPy/D");
  out->Branch("leadingJetPz",&leadingJetPz,"leadingJetPz/D");
  out->Branch("leadingJetPt",&leadingJetPt,"leadingJetPt/D");
  out->Branch("leadingJetEta",&leadingJetEta,"leadingJetEta/D");
  out->Branch("leadingJetPhi",&leadingJetPhi,"leadingJetPhi/D");
  out->Branch("leadingJetVx",&leadingJetVx,"leadingJetVx/D");
  out->Branch("leadingJetVy",&leadingJetVy,"leadingJetVy/D");
  out->Branch("leadingJetVz",&leadingJetVz,"leadingJetVz/D");
  out->Branch("leadingJetVr",&leadingJetVr,"leadingJetVr/D");
  out->Branch("leadingJetVphi",&leadingJetVphi,"leadingJetVphi/D");
  out->Branch("leadingJetBtag",&leadingJetBtag,"leadingJetBtag/D");
  out->Branch("leadingJetQGdis",&leadingJetQGdis,"leadingJetQGdis/D");
  out->Branch("leadingJetNeutralEmFrac",&leadingJetNeutralEmFrac,"leadingJetNeutralEmFrac/D");
  out->Branch("leadingJetNeutralHadFrac",&leadingJetNeutralHadFrac,"leadingJetNeutralHadFrac/D");
  out->Branch("leadingJetChargedEmFrac",&leadingJetChargedEmFrac,"leadingJetChargedEmFrac/D");
  out->Branch("leadingJetChargedHadFrac",&leadingJetChargedHadFrac,"leadingJetChargedHadFrac/D");
  out->Branch("leadingJetMuonFrac",&leadingJetMuonFrac,"leadingJetMuonFrac/D");
  out->Branch("leadingJetNeutralMulti",&leadingJetNeutralMulti,"leadingJetNeutralMulti/I");
  out->Branch("leadingJetChargedMulti",&leadingJetChargedMulti,"leadingJetChargedMulti/I");
  out->Branch("leadingJetpuIdfdisc",&leadingJetpuIdfdisc,"leadingJetpuIdfdisc/D");
  out->Branch("leadingJetpuIdcbased",&leadingJetpuIdcbased,"leadingJetpuIdcbased/I");
  out->Branch("leadingJetpuIdfid",&leadingJetpuIdfid,"leadingJetpuIdfid/I");
  out->Branch("leadingJetLooseId",&leadingJetLooseId,"leadingJetLooseId/B");
  out->Branch("leadingJetTightId",&leadingJetTightId,"leadingJetTightId/B");
  out->Branch("leadingJetLepVeto",&leadingJetLepVeto,"leadingJetLepVeto/B");
  out->Branch("secondJetEnergy",&secondJetEnergy,"secondJetEnergy/D");
  out->Branch("secondJetPx",&secondJetPx,"secondJetPx/D");
  out->Branch("secondJetPy",&secondJetPy,"secondJetPy/D");
  out->Branch("secondJetPz",&secondJetPz,"secondJetPz/D");
  out->Branch("secondJetPt",&secondJetPt,"secondJetPt/D");
  out->Branch("secondJetEta",&secondJetEta,"secondJetEta/D");
  out->Branch("secondJetPhi",&secondJetPhi,"secondJetPhi/D");
  out->Branch("secondJetVx",&secondJetVx,"secondJetVx/D");
  out->Branch("secondJetVy",&secondJetVy,"secondJetVy/D");
  out->Branch("secondJetVz",&secondJetVz,"secondJetVz/D");
  out->Branch("secondJetVr",&secondJetVr,"secondJetVr/D");
  out->Branch("secondJetVphi",&secondJetVphi,"secondJetVphi/D");
  out->Branch("secondJetBtag",&secondJetBtag,"secondJetBtag/D");
  out->Branch("secondJetQGdis",&secondJetQGdis,"secondJetQGdis/D");
  out->Branch("secondJetNeutralEmFrac",&secondJetNeutralEmFrac,"secondJetNeutralEmFrac/D");
  out->Branch("secondJetNeutralHadFrac",&secondJetNeutralHadFrac,"secondJetNeutralHadFrac/D");
  out->Branch("secondJetChargedEmFrac",&secondJetChargedEmFrac,"secondJetChargedEmFrac/D");
  out->Branch("secondJetChargedHadFrac",&secondJetChargedHadFrac,"secondJetChargedHadFrac/D");
  out->Branch("secondJetMuonFrac",&secondJetMuonFrac,"secondJetMuonFrac/D");
  out->Branch("secondJetNeutralMulti",&secondJetNeutralMulti,"secondJetNeutralMulti/I");
  out->Branch("secondJetChargedMulti",&secondJetChargedMulti,"secondJetChargedMulti/I");
  out->Branch("secondJetpuIdfdisc",&secondJetpuIdfdisc,"secondJetpuIdfdisc/D");
  out->Branch("secondJetpuIdcbased",&secondJetpuIdcbased,"secondJetpuIdcbased/I");
  out->Branch("secondJetpuIdfid",&secondJetpuIdfid,"secondJetpuIdfid/I");
  out->Branch("secondJetLooseId",&secondJetLooseId,"secondJetLooseId/B");
  out->Branch("secondJetTightId",&secondJetTightId,"secondJetTightId/B");
  out->Branch("secondJetLepVeto",&secondJetLepVeto,"secondJetLepVeto/B");
  out->Branch("leptonmetCharge",&leptonmetCharge,"leptonmetCharge/I");
  out->Branch("leptonmetMassT",&leptonmetMassT,"leptonmetMassT/D");
  out->Branch("leptonmetEta",&leptonmetEta,"leptonmetEta/D");
  out->Branch("leptonmetPhi",&leptonmetPhi,"leptonmetPhi/D");
  out->Branch("leptonmetPt",&leptonmetPt,"leptonmetPt/D");
  out->Branch("leptonmetRapidity",&leptonmetRapidity,"leptonmetRapidity/D");
  out->Branch("dileptonCharge",&dileptonCharge,"dileptonCharge/I");
  out->Branch("dileptonMass",&dileptonMass,"dileptonMass/D");
  out->Branch("dileptonEta",&dileptonEta,"dileptonEta/D");
  out->Branch("dileptonPhi",&dileptonPhi,"dileptonPhi/D");
  out->Branch("dileptonPt",&dileptonPt,"dileptonPt/D");
  out->Branch("dileptonRapidity",&dileptonRapidity,"dileptonRapidity/D");
  out->Branch("dijetMass",&dijetMass,"dijetMass/D");
  out->Branch("dijetEta",&dijetEta,"dijetEta/D");
  out->Branch("dijetPhi",&dijetPhi,"dijetPhi/D");
  out->Branch("dijetPt",&dijetPt,"dijetPt/D");
  out->Branch("dijetRapidity",&dijetRapidity,"dijetRapidity/D");
  out->Branch("leptonsystemCharge",&leptonsystemCharge,"leptonsystemCharge/I");
  out->Branch("leptonsystemMass",&leptonsystemMass,"leptonsystemMass/D");
  out->Branch("leptonsystemEta",&leptonsystemEta,"leptonsystemEta/D");
  out->Branch("leptonsystemPhi",&leptonsystemPhi,"leptonsystemPhi/D");
  out->Branch("leptonsystemPt",&leptonsystemPt,"leptonsystemPt/D");
  out->Branch("leptonsystemRapidity",&leptonsystemRapidity,"leptonsystemRapidity/D");
  out->Branch("jetsystemMass",&jetsystemMass,"jetsystemMass/D");
  out->Branch("jetsystemEta",&jetsystemEta,"jetsystemEta/D");
  out->Branch("jetsystemPhi",&jetsystemPhi,"jetsystemPhi/D");
  out->Branch("jetsystemPt",&jetsystemPt,"jetsystemPt/D");
  out->Branch("jetsystemRapidity",&jetsystemRapidity,"jetsystemRapidity/D");
  out->Branch("nJetsCandidatesLoose",&nJetsCandidatesLoose,"nJetsCandidatesLoose/I");
  out->Branch("nJetsCandidatesMedium",&nJetsCandidatesMedium,"nJetsCandidatesMedium/I");
  out->Branch("nJetsCandidatesTight",&nJetsCandidatesTight,"nJetsCandidatesTight/I");
  out->Branch("Rmpf",&Rmpf,"Rmpf/D");
  out->Branch("JZBalance",&JZBalance,"JZBalance/D");
  out->Branch("JZBalance2",&JZBalance2,"JZBalance2/D");
  out->Branch("jetcandidatesystemMass",&jetcandidatesystemMass,"jetcandidatesystemMass/D");
  out->Branch("jetcandidatesystemEta",&jetcandidatesystemEta,"jetcandidatesystemEta/D");
  out->Branch("jetcandidatesystemPhi",&jetcandidatesystemPhi,"jetcandidatesystemPhi/D");
  out->Branch("jetcandidatesystemPt",&jetcandidatesystemPt,"jetcandidatesystemPt/D");
  out->Branch("jetcandidatesystemRapidity",&jetcandidatesystemRapidity,"jetcandidatesystemRapidity/D");
  out->Branch("missingMassDijetRP210",&missingMassDijetRP210,"missingMassDijetRP210/D");
  out->Branch("missingEtaDijetRP210",&missingEtaDijetRP210,"missingEtaDijetRP210/D");
  out->Branch("missingPhiDijetRP210",&missingPhiDijetRP210,"missingPhiDijetRP210/D");
  out->Branch("missingPtDijetRP210",&missingPtDijetRP210,"missingPtDijetRP210/D");
  out->Branch("missingRapidityDijetRP210",&missingRapidityDijetRP210,"missingRapidityDijetRP210/D");
  out->Branch("missingMassDijetRP220",&missingMassDijetRP220,"missingMassDijetRP220/D");
  out->Branch("missingEtaDijetRP220",&missingEtaDijetRP220,"missingEtaDijetRP220/D");
  out->Branch("missingPhiDijetRP220",&missingPhiDijetRP220,"missingPhiDijetRP220/D");
  out->Branch("missingPtDijetRP220",&missingPtDijetRP220,"missingPtDijetRP220/D");
  out->Branch("missingRapidityDijetRP220",&missingRapidityDijetRP220,"missingRapidityDijetRP220/D");
  out->Branch("missingMassDijetMulti",&missingMassDijetMulti,"missingMassDijetMulti/D");
  out->Branch("missingEtaDijetMulti",&missingEtaDijetMulti,"missingEtaDijetMulti/D");
  out->Branch("missingPhiDijetMulti",&missingPhiDijetMulti,"missingPhiDijetMulti/D");
  out->Branch("missingPtDijetMulti",&missingPtDijetMulti,"missingPtDijetMulti/D");
  out->Branch("missingRapidityDijetMulti",&missingRapidityDijetMulti,"missingRapidityDijetMulti/D");
  out->Branch("missingMassDileptonRP210",&missingMassDileptonRP210,"missingMassDileptonRP210/D");
  out->Branch("missingEtaDileptonRP210",&missingEtaDileptonRP210,"missingEtaDileptonRP210/D");
  out->Branch("missingPhiDileptonRP210",&missingPhiDileptonRP210,"missingPhiDileptonRP210/D");
  out->Branch("missingPtDileptonRP210",&missingPtDileptonRP210,"missingPtDileptonRP210/D");
  out->Branch("missingRapidityDileptonRP210",&missingRapidityDileptonRP210,"missingRapidityDileptonRP210/D");
  out->Branch("missingMassDileptonRP220",&missingMassDileptonRP220,"missingMassDileptonRP220/D");
  out->Branch("missingEtaDileptonRP220",&missingEtaDileptonRP220,"missingEtaDileptonRP220/D");
  out->Branch("missingPhiDileptonRP220",&missingPhiDileptonRP220,"missingPhiDileptonRP220/D");
  out->Branch("missingPtDileptonRP220",&missingPtDileptonRP220,"missingPtDileptonRP220/D");
  out->Branch("missingRapidityDileptonRP220",&missingRapidityDileptonRP220,"missingRapidityDileptonRP220/D");
  out->Branch("missingMassDileptonMulti",&missingMassDileptonMulti,"missingMassDileptonMulti/D");
  out->Branch("missingEtaDileptonMulti",&missingEtaDileptonMulti,"missingEtaDileptonMulti/D");
  out->Branch("missingPhiDileptonMulti",&missingPhiDileptonMulti,"missingPhiDileptonMulti/D");
  out->Branch("missingPtDileptonMulti",&missingPtDileptonMulti,"missingPtDileptonMulti/D");
  out->Branch("missingRapidityDileptonMulti",&missingRapidityDileptonMulti,"missingRapidityDileptonMulti/D");
  out->Branch("missingMassDileptonJet1RP220",&missingMassDileptonJet1RP220,"missingMassDileptonJet1RP220/D");
  out->Branch("missingEtaDileptonJet1RP220",&missingEtaDileptonJet1RP220,"missingEtaDileptonJet1RP220/D");
  out->Branch("missingPhiDileptonJet1RP220",&missingPhiDileptonJet1RP220,"missingPhiDileptonJet1RP220/D");
  out->Branch("missingPtDileptonJet1RP220",&missingPtDileptonJet1RP220,"missingPtDileptonJet1RP220/D");
  out->Branch("missingRapidityDileptonJet1RP220",&missingRapidityDileptonJet1RP220,"missingRapidityDileptonJet1RP220/D");
  out->Branch("missingMassDileptonJet1Multi",&missingMassDileptonJet1Multi,"missingMassDileptonJet1Multi/D");
  out->Branch("missingEtaDileptonJet1Multi",&missingEtaDileptonJet1Multi,"missingEtaDileptonJet1Multi/D");
  out->Branch("missingPhiDileptonJet1Multi",&missingPhiDileptonJet1Multi,"missingPhiDileptonJet1Multi/D");
  out->Branch("missingPtDileptonJet1Multi",&missingPtDileptonJet1Multi,"missingPtDileptonJet1Multi/D");
  out->Branch("missingRapidityDileptonJet1Multi",&missingRapidityDileptonJet1Multi,"missingRapidityDileptonJet1Multi/D");
  out->Branch("missingMassDileptonDijetRP210",&missingMassDileptonDijetRP210,"missingMassDileptonDijetRP210/D");
  out->Branch("missingEtaDileptonDijetRP210",&missingEtaDileptonDijetRP210,"missingEtaDileptonDijetRP210/D");
  out->Branch("missingPhiDileptonDijetRP210",&missingPhiDileptonDijetRP210,"missingPhiDileptonDijetRP210/D");
  out->Branch("missingPtDileptonDijetRP210",&missingPtDileptonDijetRP210,"missingPtDileptonDijetRP210/D");
  out->Branch("missingRapidityDileptonDijetRP210",&missingRapidityDileptonDijetRP210,"missingRapidityDileptonDijetRP210/D");
  out->Branch("missingMassDileptonDijetRP220",&missingMassDileptonDijetRP220,"missingMassDileptonDijetRP220/D");
  out->Branch("missingEtaDileptonDijetRP220",&missingEtaDileptonDijetRP220,"missingEtaDileptonDijetRP220/D");
  out->Branch("missingPhiDileptonDijetRP220",&missingPhiDileptonDijetRP220,"missingPhiDileptonDijetRP220/D");
  out->Branch("missingPtDileptonDijetRP220",&missingPtDileptonDijetRP220,"missingPtDileptonDijetRP220/D");
  out->Branch("missingRapidityDileptonDijetRP220",&missingRapidityDileptonDijetRP220,"missingRapidityDileptonDijetRP220/D");
  out->Branch("missingMassDileptonDijetMulti",&missingMassDileptonDijetMulti,"missingMassDileptonDijetMulti/D");
  out->Branch("missingEtaDileptonDijetMulti",&missingEtaDileptonDijetMulti,"missingEtaDileptonDijetMulti/D");
  out->Branch("missingPhiDileptonDijetMulti",&missingPhiDileptonDijetMulti,"missingPhiDileptonDijetMulti/D");
  out->Branch("missingPtDileptonDijetMulti",&missingPtDileptonDijetMulti,"missingPtDileptonDijetMulti/D");
  out->Branch("missingRapidityDileptonDijetMulti",&missingRapidityDileptonDijetMulti,"missingRapidityDileptonDijetMulti/D");
  out->Branch("missEt",&missEt,"missEt/D");
  out->Branch("missEt_phi",&missEt_phi,"missEt_phi/D");
  out->Branch("dphi",&dphi,"dphi/D");
  out->Branch("diffMassRP210",&diffMassRP210,"diffMassRP210/D");
  out->Branch("diffMassRP220",&diffMassRP220,"diffMassRP220/D");
  out->Branch("diffMassMulti",&diffMassMulti,"diffMassMulti/D");
  out->Branch("protonsrp210Arm45Pz",&proton_pz_rp210Arm45,"protonsrp210Arm45Pz/D");
  out->Branch("protonsrp210Arm56Pz",&proton_pz_rp210Arm56,"protonsrp210Arm56Pz/D");
  out->Branch("protonsrp220Arm45Pz",&proton_pz_rp220Arm45,"protonsrp220Arm45Pz/D");
  out->Branch("protonsrp220Arm56Pz",&proton_pz_rp220Arm56,"protonsrp220Arm56Pz/D");
  out->Branch("protonsmultiArm45Pz",&proton_pz_multiArm45,"protonsmultiArm45Pz/D");
  out->Branch("protonsmultiArm56Pz",&proton_pz_multiArm56,"protonsmultiArm56Pz/D");
  out->Branch("nprotonRP210_Arm45",&nprotonRP210_sec45,"nprotonRP210_Arm45/I");
  out->Branch("nprotonRP210_Arm56",&nprotonRP210_sec56,"nprotonRP210_Arm56/I");
  out->Branch("nprotonRP220_Arm45",&nprotonRP220_sec45,"nprotonRP220_Arm45/I");
  out->Branch("nprotonRP220_Arm56",&nprotonRP220_sec56,"nprotonRP220_Arm56/I");
  out->Branch("nprotonMultiArm45",&nprotonMulti_sec45,"nprotonMultiArm45/I");
  out->Branch("nprotonMultiArm56",&nprotonMulti_sec56,"nprotonMultiArm56/I");
  out->Branch("isprotonRP210",&isprotonRP210,"isprotonRP210/B");
  out->Branch("isprotonRP220",&isprotonRP220,"isprotonRP220/B");
  out->Branch("isprotonMulti",&isprotonMulti,"isprotonMulti/B");
  out->Branch("isprotonRP210_Reduced ",&isprotonRP210_Reduced,"isprotonRP210_Reduced/B");
  out->Branch("isprotonRP220_Reduced ",&isprotonRP220_Reduced,"isprotonRP220_Reduced/B");
  out->Branch("xi_rp210_Arm45",&xi_rp210_Arm45,"xi_rp210_Arm45/D");
  out->Branch("xi_rp210_Arm56",&xi_rp210_Arm56,"xi_rp210_Arm56/D");
  out->Branch("xi_rp220_Arm45",&xi_rp220_Arm45,"xi_rp220_Arm45/D");
  out->Branch("xi_rp220_Arm56",&xi_rp220_Arm56,"xi_rp220_Arm56/D");
  out->Branch("xi_multiArm45",&xi_multiArm45,"xi_multiArm45/D");
  out->Branch("time_multiArm45",&time_multiArm45,"time_multiArm45/D");
  out->Branch("timeunc_multiArm45",&timeunc_multiArm45,"timeunc_multiArm45/D");
  out->Branch("thx_multiArm45",&thx_multiArm45,"thx_multiArm45/D");
  out->Branch("thy_multiArm45",&thy_multiArm45,"thy_multiArm45/D");
  out->Branch("xi_multiArm56",&xi_multiArm56,"xi_multiArm56/D");
  out->Branch("time_multiArm56",&time_multiArm56,"time_multiArm56/D");
  out->Branch("timeunc_multiArm56",&timeunc_multiArm56,"timeunc_multiArm56/D");
  out->Branch("thx_multiArm56",&thx_multiArm56,"thx_multiArm56/D");
  out->Branch("thy_multiArm56",&thy_multiArm56,"thy_multiArm56/D");
  out->Branch("t45andt56",&t45andt56,"t45andt56/D");
  out->Branch("vzpps",&vzpps,"vzpps/D");
  out->Branch("distanceVertexZPPSandCMS",&distanceVertexZPPSandCMS,"distanceVertexZPPSandCMS/D");
  out->Branch("selectedVertexZ",&selectedVertexZ,"selectedVertexZ/D");

}


void TTreeMissingMass::Fill(){
  out->Fill();
}

void TTreeMissingMass::Storing(){
  fileout->cd();
  out->Write();
  fileout->Close();
}


