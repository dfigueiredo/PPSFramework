//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct 14 14:42:49 2018 by ROOT version 6.10/09
// from TTree analyzer/analyzer
// found on file: /eos/cms/store/group/phys_pps/MissingMassSearch/MuonB-MiniAOD-v2/DoubleMuon/MuonB-MiniAOD-v2/181013_221704/0000/output_30.root
//////////////////////////////////////////////////////////

#ifndef MissingMassNtupleAnalyzer_h
#define MissingMassNtupleAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

// Header file for the classes stored in the TTree if any.
#include <vector>

std::vector<double> rndXi_RP210_sec45;
std::vector<double> rndXi_RP210_sec56;
std::vector<double> rndXi_RP220_sec45;
std::vector<double> rndXi_RP220_sec56;
std::vector<double> rndXi_Multi_sec45;
std::vector<double> rndXi_Multi_sec56;

class MissingMassNtupleAnalyzer {
	public :
		TTree          *fChain;   //!pointer to the analyzed TTree or TChain
		Int_t           fCurrent; //!current Tree number in a TChain

		// Fixed size dimensions of array or collections stored in the TTree if any.
		static constexpr Int_t kMaxleptons_pfIsoMedium = 1;
		static constexpr Int_t kMaxleptons_miniIsoTight = 1;
		static constexpr Int_t kMaxleptons_pfIsoVeryTight = 1;
		static constexpr Int_t kMaxleptons_pfIso = 1;
		static constexpr Int_t kMaxleptons_tkIso = 1;

		// Declaration of leaf types
		Int_t           xangle;
		Int_t           run;
		Long64_t        event;
		Int_t           lumiblock;
		Int_t		PUInterac;
		Int_t		PUTrueInterac;
		std::vector<int>     *trigger;
		std::vector<int>     *prescalesL1;
		std::vector<int>     *prescalesHLT;
		std::vector<double>   *vertex_x;
		std::vector<double>   *vertex_y;
		std::vector<double>   *vertex_z;
		std::vector<int>      *vertex_ntrack;
		std::vector<double>   *vertex_chi2;
		std::vector<double>   *vertex_ndof;
		
	        Int_t 	        nGenParticles;
	        Int_t           nGenElectrons;
	        Int_t           nGenMuons;

	        std::vector<int>    *genleptons_pdgid;
	        std::vector<double> *genleptons_energy;
	        std::vector<double> *genleptons_pt;
	        std::vector<double> *genleptons_eta;
	        std::vector<double> *genleptons_phi;
	        std::vector<double> *genleptons_px;
	        std::vector<double> *genleptons_py;
	        std::vector<double> *genleptons_pz;
	        std::vector<double> *genleptons_charge;
	        std::vector<double> *genleptons_vx;
	        std::vector<double> *genleptons_vy;
	        std::vector<double> *genleptons_vz;

	        std::vector<double> *genphotons_pt;
	        std::vector<double> *genphotons_eta;
	        std::vector<double> *genphotons_phi;
	        std::vector<double> *genphotons_energy;

	        double_t genmissEt;
	        double_t genmissEt_phi;

	        std::vector<double> *genjetsak4_px;
	        std::vector<double> *genjetsak4_py;
	        std::vector<double> *genjetsak4_pz;
	        std::vector<double> *genjetsak4_pt;
	        std::vector<double> *genjetsak4_energy;
	        std::vector<double> *genjetsak4_phi;
	        std::vector<double> *genjetsak4_eta;
	        std::vector<double> *genjetsak4_vz;

	        std::vector<double> *genjetsak8_px;
	        std::vector<double> *genjetsak8_py;
	        std::vector<double> *genjetsak8_pz;
	        std::vector<double> *genjetsak8_pt;
	        std::vector<double> *genjetsak8_energy;
	        std::vector<double> *genjetsak8_phi;
	        std::vector<double> *genjetsak8_eta;
	        std::vector<double> *genjetsak8_vz;	

		std::vector<double>   *jetsak4_px;
		std::vector<double>   *jetsak4_py;
		std::vector<double>   *jetsak4_pz;
		std::vector<double>   *jetsak4_pt;
		std::vector<double>   *jetsak4_energy;
		std::vector<double>   *jetsak4_phi;
		std::vector<double>   *jetsak4_eta;
		std::vector<double>   *jetsak4_vx;
		std::vector<double>   *jetsak4_vy;
		std::vector<double>   *jetsak4_vz;
		std::vector<double>   *jetsak4_bdis;
		std::vector<double>   *jetsak4_qgdis;
		std::vector<double>   *jetsak4_neutralemfrac;
		std::vector<double>   *jetsak4_neutralhadfrac;
		std::vector<double>   *jetsak4_chargedemfrac;
		std::vector<double>   *jetsak4_chargedhadfrac;
		std::vector<double>   *jetsak4_muonfrac;
		std::vector<int>      *jetsak4_neutralmulti;
		std::vector<int>      *jetsak4_chargedmulti;
		std::vector<double>   *jetsak4_puIdfdisc;
		std::vector<int>      *jetsak4_puIdcbased;
		std::vector<int>      *jetsak4_puIdfid;
		std::vector<bool>     *jetsak4_looseId;
		std::vector<bool>     *jetsak4_tightId;
		std::vector<bool>     *jetsak4_lepVeto;
		std::vector<double>   *jetsak8_px;
		std::vector<double>   *jetsak8_py;
		std::vector<double>   *jetsak8_pz;
		std::vector<double>   *jetsak8_pt;
		std::vector<double>   *jetsak8_energy;
		std::vector<double>   *jetsak8_phi;
		std::vector<double>   *jetsak8_eta;
		std::vector<double>   *jetsak8_vx;
		std::vector<double>   *jetsak8_vy;
		std::vector<double>   *jetsak8_vz;
		std::vector<double>   *jetsak8_bdis;
		std::vector<bool>    *jetsak8_looseId;
		std::vector<bool>    *jetsak8_tightId;
		std::vector<bool>    *jetsak8_lepVeto;

		Int_t           nElectrons;
		Int_t           nElectrons_looseId;
		Int_t           nElectrons_mediumId;
		Int_t           nElectrons_tightId;
		Int_t           nMuons;
		Int_t           nMuons_looseId;
		Int_t           nMuons_mediumId;
		Int_t           nMuons_tightId;
		Int_t 	        nChargedPFMultiPV_Loose;
		Int_t		nChargedPFMultiPV_Tight;
		Int_t		nChargedPFMultiPV_UsedInFit;
		Int_t		nChargedPFMultiPV_Tight_Fit;
		Double_t	SumChargedPFMultiPV_pt_Loose;
		Double_t        SumChargedPFMultiPV_pt_Tight;
		Double_t        SumChargedPFMultiPV_pt_UsedInFit;
		Double_t	SumChargedPFMultiPV_pt_Tight_Fit;
		std::vector<int>      *leptons_pdgid;
		std::vector<double>   *leptons_energy;
		std::vector<double>   *leptons_pt;
		std::vector<double>   *leptons_eta;
		std::vector<double>   *leptons_phi;
		std::vector<double>   *leptons_px;
		std::vector<double>   *leptons_py;
		std::vector<double>   *leptons_pz;
		std::vector<double>   *leptons_charge;
		std::vector<double>   *leptons_vx;
		std::vector<double>   *leptons_vy;
		std::vector<double>   *leptons_vz;
		std::vector<bool>    *leptons_looseId;
		std::vector<bool>    *leptons_mediumId;
		std::vector<bool>    *leptons_tightId;
		std::vector<bool>    *leptons_pfIsoMedium_;
		std::vector<bool>    *leptons_miniIsoTight_;
		std::vector<bool>    *leptons_pfIsoVeryTight_;
		std::vector<double>   *leptons_pfIso_;
		std::vector<double>   *leptons_tkIso_;
		double_t         missEt;
		double_t         missEt_phi;
		std::vector<bool>    *protonsIsValid;
		std::vector<int>     *protonsArm;
		std::vector<int>     *protonsStation;
		std::vector<int>     *protonsRP;
		std::vector<double>   *protonsX;
		std::vector<double>   *protonsXUnc;
		std::vector<double>   *protonsY;
		std::vector<double>   *protonsYUnc;

		std::vector<int> *singleProtonArm;
		std::vector<int> *singleProtonStation;
		std::vector<int> *singleProtonPot;
		std::vector<double> *singleProtonXi;
		std::vector<double> *singleProtonThetaX;
		std::vector<double> *singleProtonThetaY;

		std::vector<int> *multiProtonArm;
		std::vector<double> *multiProtonXi;
		std::vector<double> *multiProtonTime;
		std::vector<double> *multiProtonTimeError;
		std::vector<double> *multiProtonThetaX;
		std::vector<double> *multiProtonThetaY;

		// List of branches
		TBranch        *b_xangle;   //!
		TBranch        *b_run;   //!
		TBranch        *b_event;   //!
		TBranch        *b_lumiblock;   //!
		TBranch        *b_trigger;   //!
		TBranch        *b_prescalesL1;   //!
		TBranch        *b_prescalesHLT;   //!
		TBranch        *b_vertex_x;   //!
		TBranch        *b_vertex_y;   //!
		TBranch        *b_vertex_z;   //!
		TBranch        *b_vertex_ntrack;   //!
		TBranch        *b_vertex_chi2;   //!
		TBranch        *b_vertex_ndof;   //!

	        TBranch        *b_PUInterac;
	        TBranch        *b_PUTrueInterac;
	        TBranch        *b_nGenParticles;
	        TBranch        *b_nGenElectrons;
	        TBranch        *b_nGenMuons;

	        TBranch        *b_genleptons_pdgid;
	        TBranch        *b_genleptons_energy;
	        TBranch        *b_genleptons_pt;
	        TBranch        *b_genleptons_eta;
	        TBranch        *b_genleptons_phi;
	        TBranch        *b_genleptons_px;
	        TBranch        *b_genleptons_py;
	        TBranch        *b_genleptons_pz;
	        TBranch        *b_genleptons_charge;
	        TBranch        *b_genleptons_vx;
	        TBranch        *b_genleptons_vy;
	        TBranch        *b_genleptons_vz;

	        TBranch        *b_genphotons_pt;
	        TBranch        *b_genphotons_eta;
	        TBranch        *b_genphotons_phi;
	        TBranch        *b_genphotons_energy;

	        TBranch        *b_genmissEt;
	        TBranch        *b_genmissEt_phi;

	        TBranch        *b_genjetsak4_px;
	        TBranch        *b_genjetsak4_py;
	        TBranch        *b_genjetsak4_pz;
	        TBranch        *b_genjetsak4_pt;
	        TBranch        *b_genjetsak4_energy;
	        TBranch        *b_genjetsak4_phi;
	        TBranch        *b_genjetsak4_eta;
	        TBranch        *b_genjetsak4_vz;
	        TBranch        *b_genjetsak8_px;
	        TBranch        *b_genjetsak8_py;
	        TBranch        *b_genjetsak8_pz;
	        TBranch        *b_genjetsak8_pt;
	        TBranch        *b_genjetsak8_energy;
	        TBranch        *b_genjetsak8_phi;
	        TBranch        *b_genjetsak8_eta;
	        TBranch        *b_genjetsak8_vz;
		
		TBranch        *b_jetsak4_px;   //!
		TBranch        *b_jetsak4_py;   //!
		TBranch        *b_jetsak4_pz;   //!
		TBranch        *b_jetsak4_pt;   //!
		TBranch        *b_jetsak4_energy;   //!
		TBranch        *b_jetsak4_phi;   //!
		TBranch        *b_jetsak4_eta;   //!
		TBranch        *b_jetsak4_vx; 
		TBranch        *b_jetsak4_vy; 
		TBranch        *b_jetsak4_vz;   //!
		TBranch        *b_jetsak4_bdis;   //
		TBranch        *b_jetsak4_qgdis;
		TBranch        *b_jetsak4_neutralemfrac;
		TBranch        *b_jetsak4_neutralhadfrac;
		TBranch        *b_jetsak4_chargedemfrac;
		TBranch        *b_jetsak4_chargedhadfrac;
		TBranch        *b_jetsak4_muonfrac;
		TBranch        *b_jetsak4_neutralmulti;
		TBranch        *b_jetsak4_chargedmulti;
		TBranch        *b_jetsak4_puIdfdisc;
		TBranch        *b_jetsak4_puIdcbased;
		TBranch        *b_jetsak4_puIdfid;
		TBranch        *b_jetsak4_looseId;   //!
		TBranch        *b_jetsak4_tightId;   //!
		TBranch        *b_jetsak4_lepVeto;   //!
		TBranch        *b_jetsak8_px;   //!
		TBranch        *b_jetsak8_py;   //!
		TBranch        *b_jetsak8_pz;   //!
		TBranch        *b_jetsak8_pt;   //!
		TBranch        *b_jetsak8_energy;   //!
		TBranch        *b_jetsak8_phi;   //!
		TBranch        *b_jetsak8_eta;   //!
		TBranch        *b_jetsak8_vx;
		TBranch        *b_jetsak8_vy;
		TBranch        *b_jetsak8_vz;   //!*/
		TBranch        *b_jetsak8_bdis;   //!
		TBranch        *b_jetsak8_looseId;   //!
		TBranch        *b_jetsak8_tightId;   //!
		TBranch        *b_jetsak8_lepVeto;   //!

		TBranch        *b_nElectrons;   //!
		TBranch        *b_nElectrons_looseId;   //!
		TBranch        *b_nElectrons_mediumId;   //!
		TBranch        *b_nElectrons_tightId;   //!
		TBranch        *b_nMuons;   //!
		TBranch        *b_nMuons_looseId;   //!
		TBranch        *b_nMuons_mediumId;   //!
		TBranch        *b_nMuons_tightId;   //!
 		TBranch	       *b_nChargedPFMultiPV_Loose;
 		TBranch        *b_nChargedPFMultiPV_Tight;
 		TBranch	       *b_nChargedPFMultiPV_UsedInFit;
 		TBranch	       *b_nChargedPFMultiPV_Tight_Fit;
                TBranch        *b_SumChargedPFMultiPV_pt_Loose;
                TBranch        *b_SumChargedPFMultiPV_pt_Tight;
                TBranch        *b_SumChargedPFMultiPV_pt_UsedInFit;
                TBranch        *b_SumChargedPFMultiPV_pt_Tight_Fit;
		TBranch        *b_leptons_pdgid;   //!
		TBranch        *b_leptons_energy;   //!
		TBranch        *b_leptons_pt;   //!
		TBranch        *b_leptons_eta;   //!
		TBranch        *b_leptons_phi;   //!
		TBranch        *b_leptons_px;   //!
		TBranch        *b_leptons_py;   //!
		TBranch        *b_leptons_pz;   //!
		TBranch        *b_leptons_charge;   //!
		TBranch        *b_leptons_vx;   //!
		TBranch        *b_leptons_vy;   //!
		TBranch        *b_leptons_vz;   //!
		TBranch        *b_leptons_looseId;   //!
		TBranch        *b_leptons_mediumId;   //!
		TBranch        *b_leptons_tightId;   //!
		TBranch        *b_leptons_pfIsoMedium_;   //!
		TBranch        *b_leptons_miniIsoTight_;   //!
		TBranch        *b_leptons_pfIsoVeryTight_;   //!
		TBranch        *b_leptons_pfIso_;   //!
		TBranch        *b_leptons_tkIso_;   //!
		TBranch        *b_missEt;   //!
		TBranch        *b_missEt_phi;   //!
		TBranch        *b_protonsArm;   //!
		TBranch        *b_protonsStation;   //!
		TBranch        *b_protonsRP;   //!
		TBranch        *b_protonsX;   //!
		TBranch        *b_protonsXUnc;   //!
		TBranch        *b_protonsY;   //!
		TBranch        *b_protonsYUnc;   //!

		TBranch        *b_singleProtonArm;
		TBranch        *b_singleProtonStation;
		TBranch        *b_singleProtonPot;
		TBranch        *b_singleProtonXi;
		TBranch        *b_singleProtonThetaX;
		TBranch        *b_singleProtonThetaY;

		TBranch        *b_multiProtonArm;
		TBranch        *b_multiProtonXi;
		TBranch        *b_multiProtonTime;
		TBranch        *b_multiProtonTimeError;
		TBranch        *b_multiProtonThetaX;
		TBranch        *b_multiProtonThetaY;

		MissingMassNtupleAnalyzer(TTree *tree=0);
		virtual ~MissingMassNtupleAnalyzer();
		virtual Int_t    Cut(Long64_t entry);
		virtual Int_t    GetEntry(Long64_t entry);
		virtual Long64_t LoadTree(Long64_t entry);
		virtual void     Init(TTree *tree);
		virtual void     Loop(char*, char*, char*, char*, char*, char*, bool, bool, bool, bool, bool, bool);
		virtual Bool_t   Notify();
		virtual void     Show(Long64_t entry = -1);

		char* getCmdOption(char**, char**, const std::string&);
		bool cmdOptionExists(char**, char**, const std::string&);
		void randomProtons(TString);
		bool MCMatch(double, double, double, double, double, double);
};

#endif

#ifdef MissingMassNtupleAnalyzer_cxx
MissingMassNtupleAnalyzer::MissingMassNtupleAnalyzer(TTree *tree) : fChain(0) 
{
	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	if (tree == 0) {
		std::cout << "\n\t --> Please, insert a valid TTree file\n" << std::endl;
		exit(EXIT_FAILURE);
	}
	Init(tree);
}

MissingMassNtupleAnalyzer::~MissingMassNtupleAnalyzer()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t MissingMassNtupleAnalyzer::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t MissingMassNtupleAnalyzer::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (fChain->GetTreeNumber() != fCurrent) {
		fCurrent = fChain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void MissingMassNtupleAnalyzer::Init(TTree *tree)
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set object pointer
	trigger = 0;
	prescalesL1 = 0;
	prescalesHLT = 0;
	vertex_x = 0;
	vertex_y = 0;
	vertex_z = 0;
	vertex_ntrack = 0;
	vertex_chi2 = 0;
	vertex_ndof = 0;

	PUInterac = 0;
	PUTrueInterac = 0;

	nGenParticles = 0;
	nGenElectrons = 0;
	nGenMuons = 0;

        genleptons_pdgid = 0;
        genleptons_energy = 0;
        genleptons_pt = 0;
        genleptons_eta = 0;
        genleptons_phi = 0;
        genleptons_px = 0;
        genleptons_py = 0;
        genleptons_pz = 0;
        genleptons_charge = 0;
        genleptons_vx = 0;
        genleptons_vy = 0;
        genleptons_vz = 0;
        
        genphotons_pt = 0;
        genphotons_eta = 0;
        genphotons_phi = 0;
        genphotons_energy = 0;
        
        genmissEt = 0;
        genmissEt_phi = 0;
        
        genjetsak4_px = 0;
        genjetsak4_py = 0;
        genjetsak4_pz = 0;
        genjetsak4_pt = 0;
        genjetsak4_energy = 0;
        genjetsak4_phi = 0;
        genjetsak4_eta = 0;
        genjetsak4_vz = 0;
        genjetsak8_px = 0;
        genjetsak8_py = 0;
        genjetsak8_pz = 0;
        genjetsak8_pt = 0;
        genjetsak8_energy = 0;
        genjetsak8_phi = 0;
        genjetsak8_eta = 0;
        genjetsak8_vz = 0;
	
	jetsak4_px = 0;
	jetsak4_py = 0;
	jetsak4_pz = 0;
	jetsak4_pt = 0;
	jetsak4_energy = 0;
	jetsak4_phi = 0;
	jetsak4_eta = 0;
	jetsak4_vx = 0;
	jetsak4_vy = 0;
	jetsak4_vz = 0;
	jetsak4_bdis = 0;
	jetsak4_qgdis = 0;
	jetsak4_neutralemfrac = 0;
	jetsak4_neutralhadfrac = 0;
	jetsak4_chargedemfrac = 0;
	jetsak4_chargedhadfrac = 0;
	jetsak4_muonfrac = 0;
	jetsak4_neutralmulti = 0;
	jetsak4_chargedmulti = 0;
	jetsak4_puIdfdisc = 0;
	jetsak4_puIdcbased = 0;
	jetsak4_puIdfid = 0;
	jetsak4_looseId = 0;
	jetsak4_tightId = 0;
	jetsak4_lepVeto = 0;
	jetsak8_px = 0;
	jetsak8_py = 0;
	jetsak8_pz = 0;
	jetsak8_pt = 0;
	jetsak8_energy = 0;
	jetsak8_phi = 0;
	jetsak8_eta = 0;
	jetsak8_vx = 0;
	jetsak8_vy = 0;
	jetsak8_vz = 0;
	jetsak8_bdis = 0;
	jetsak8_looseId = 0;
	jetsak8_tightId = 0;
	jetsak8_lepVeto = 0;

	nElectrons = 0;
	nElectrons_looseId = 0;
	nElectrons_mediumId = 0;
	nElectrons_tightId = 0;
	nMuons = 0;
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
	leptons_pdgid = 0;
	leptons_energy = 0;
	leptons_pt = 0;
	leptons_eta = 0;
	leptons_phi = 0;
	leptons_px = 0;
	leptons_py = 0;
	leptons_pz = 0;
	leptons_charge = 0;
	leptons_vx = 0;
	leptons_vy = 0;
	leptons_vz = 0;
	leptons_looseId = 0;
	leptons_mediumId = 0;
	leptons_tightId = 0;
	leptons_pfIsoMedium_ = 0;
	leptons_miniIsoTight_ = 0;
	leptons_pfIsoVeryTight_ = 0;
	leptons_pfIso_ = 0;
	leptons_tkIso_ = 0;
	protonsArm = 0;
	protonsStation = 0;
	protonsRP = 0;
	protonsX = 0;
	protonsXUnc = 0;
	protonsY = 0;
	protonsYUnc = 0;
	singleProtonArm = 0;
	singleProtonStation = 0;
	singleProtonPot = 0;
	singleProtonXi = 0;
	singleProtonThetaX = 0;
	singleProtonThetaY = 0;
	multiProtonArm = 0;
	multiProtonXi = 0;
	multiProtonTime = 0;
	multiProtonTimeError = 0;
	multiProtonThetaX = 0;
	multiProtonThetaY = 0;

	// Set branch addresses and branch pointers
	if (!tree) return;
	fChain = tree;
	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("xangle", &xangle, &b_xangle);
	fChain->SetBranchAddress("run", &run, &b_run);
	fChain->SetBranchAddress("event", &event, &b_event);
	fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
	fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
	fChain->SetBranchAddress("prescalesL1", &prescalesL1, &b_prescalesL1);
	fChain->SetBranchAddress("prescalesHLT", &prescalesHLT, &b_prescalesHLT);
	fChain->SetBranchAddress("vertex_x", &vertex_x, &b_vertex_x);
	fChain->SetBranchAddress("vertex_y", &vertex_y, &b_vertex_y);
	fChain->SetBranchAddress("vertex_z", &vertex_z, &b_vertex_z);
	fChain->SetBranchAddress("vertex_ntrack", &vertex_ntrack, &b_vertex_ntrack);
	fChain->SetBranchAddress("vertex_chi2", &vertex_chi2, &b_vertex_chi2);
	fChain->SetBranchAddress("vertex_ndof", &vertex_ndof, &b_vertex_ndof);

/*
        fChain->SetBranchAddress("PUInterac", &PUInterac, &b_PUInterac);
        fChain->SetBranchAddress("PUTrueInterac", &PUTrueInterac, &b_PUTrueInterac);
        fChain->SetBranchAddress("nGenParticles", &nGenParticles, &b_nGenParticles);
        fChain->SetBranchAddress("nGenMuons", &nGenMuons, &b_nGenMuons);
        fChain->SetBranchAddress("nGenElectrons", &nGenElectrons, &b_nGenElectrons);
        fChain->SetBranchAddress("genleptons_pdgid", &genleptons_pdgid, &b_genleptons_pdgid);
        fChain->SetBranchAddress("genleptons_energy", &genleptons_energy, &b_genleptons_energy);
        fChain->SetBranchAddress("genleptons_pt", &genleptons_pt, &b_genleptons_pt);
        fChain->SetBranchAddress("genleptons_eta", &genleptons_eta, &b_genleptons_eta);
        fChain->SetBranchAddress("genleptons_phi", &genleptons_phi, &b_genleptons_phi);
        fChain->SetBranchAddress("genleptons_px", &genleptons_px, &b_genleptons_px);
        fChain->SetBranchAddress("genleptons_py", &genleptons_py, &b_genleptons_py);
        fChain->SetBranchAddress("genleptons_pz", &genleptons_pz, &b_genleptons_pz);
        fChain->SetBranchAddress("genleptons_charge", &genleptons_charge, &b_genleptons_charge);
        fChain->SetBranchAddress("genleptons_vx", &genleptons_vx, &b_genleptons_vx);
        fChain->SetBranchAddress("genleptons_vy", &genleptons_vy, &b_genleptons_vy);
        fChain->SetBranchAddress("genleptons_vz", &genleptons_vz, &b_genleptons_vz);
        fChain->SetBranchAddress("genphotons_pt", &genphotons_pt, &b_genphotons_pt);
        fChain->SetBranchAddress("genphotons_eta", &genphotons_eta, &b_genphotons_eta);
        fChain->SetBranchAddress("genphotons_phi", &genphotons_phi, &b_genphotons_phi);
        fChain->SetBranchAddress("genphotons_energy", &genphotons_energy, &b_genphotons_energy);
        fChain->SetBranchAddress("genMissEt", &genmissEt, &b_genmissEt);
        fChain->SetBranchAddress("genMissEt_phi", &genmissEt_phi, &b_genmissEt_phi);
        fChain->SetBranchAddress("genjetsak4_px", &genjetsak4_px, &b_genjetsak4_px);
        fChain->SetBranchAddress("genjetsak4_py", &genjetsak4_py, &b_genjetsak4_py);
        fChain->SetBranchAddress("genjetsak4_pz", &genjetsak4_pz, &b_genjetsak4_pz);
        fChain->SetBranchAddress("genjetsak4_pt", &genjetsak4_pt, &b_genjetsak4_pt);
        fChain->SetBranchAddress("genjetsak4_energy", &genjetsak4_energy, &b_genjetsak4_energy);
        fChain->SetBranchAddress("genjetsak4_phi", &genjetsak4_phi, &b_genjetsak4_phi);
        fChain->SetBranchAddress("genjetsak4_eta", &genjetsak4_eta, &b_genjetsak4_eta);
        fChain->SetBranchAddress("genjetsak4_vz", &genjetsak4_vz, &b_genjetsak4_vz);
        fChain->SetBranchAddress("genjetsak8_px", &genjetsak8_px, &b_genjetsak8_px);
        fChain->SetBranchAddress("genjetsak8_py", &genjetsak8_py, &b_genjetsak8_py);
        fChain->SetBranchAddress("genjetsak8_pz", &genjetsak8_pz, &b_genjetsak8_pz);
        fChain->SetBranchAddress("genjetsak8_pt", &genjetsak8_pt, &b_genjetsak8_pt);
        fChain->SetBranchAddress("genjetsak8_energy", &genjetsak8_energy, &b_genjetsak8_energy);
        fChain->SetBranchAddress("genjetsak8_phi", &genjetsak8_phi, &b_genjetsak8_phi);
        fChain->SetBranchAddress("genjetsak8_eta", &genjetsak8_eta, &b_genjetsak8_eta);
        fChain->SetBranchAddress("genjetsak8_vz", &genjetsak8_vz, &b_genjetsak8_vz);
*/

	fChain->SetBranchAddress("jetsak4_px", &jetsak4_px, &b_jetsak4_px);
	fChain->SetBranchAddress("jetsak4_py", &jetsak4_py, &b_jetsak4_py);
	fChain->SetBranchAddress("jetsak4_pz", &jetsak4_pz, &b_jetsak4_pz);
	fChain->SetBranchAddress("jetsak4_pt", &jetsak4_pt, &b_jetsak4_pt);
	fChain->SetBranchAddress("jetsak4_energy", &jetsak4_energy, &b_jetsak4_energy);
	fChain->SetBranchAddress("jetsak4_phi", &jetsak4_phi, &b_jetsak4_phi);
	fChain->SetBranchAddress("jetsak4_eta", &jetsak4_eta, &b_jetsak4_eta);
	fChain->SetBranchAddress("jetsak4_vx", &jetsak4_vx, &b_jetsak4_vx);
	fChain->SetBranchAddress("jetsak4_vy", &jetsak4_vy, &b_jetsak4_vy);
	fChain->SetBranchAddress("jetsak4_vz", &jetsak4_vz, &b_jetsak4_vz);
	fChain->SetBranchAddress("jetsak4_bdis", &jetsak4_bdis, &b_jetsak4_bdis);
	fChain->SetBranchAddress("jetsak4_qgdis",&jetsak4_qgdis);
	fChain->SetBranchAddress("jetsak4_neutralEmFraction",&jetsak4_neutralemfrac);
	fChain->SetBranchAddress("jetsak4_neutralHadFraction",&jetsak4_neutralhadfrac);
	fChain->SetBranchAddress("jetsak4_chargedEmFraction",&jetsak4_chargedemfrac);
	fChain->SetBranchAddress("jetsak4_chargedHadFraction",&jetsak4_chargedhadfrac);
	fChain->SetBranchAddress("jetsak4_muonFraction",&jetsak4_muonfrac);
	fChain->SetBranchAddress("jetsak4_neutralMultiplicity",&jetsak4_neutralmulti);
	fChain->SetBranchAddress("jetsak4_chargedMultiplicity",&jetsak4_chargedmulti);
	fChain->SetBranchAddress("jetsak4_puIdfdisc",&jetsak4_puIdfdisc);
	fChain->SetBranchAddress("jetsak4_puIdcbased",&jetsak4_puIdcbased);
	fChain->SetBranchAddress("jetsak4_puIdfid",&jetsak4_puIdfid);
	fChain->SetBranchAddress("jetsak4_looseId", &jetsak4_looseId, &b_jetsak4_looseId);
	fChain->SetBranchAddress("jetsak4_tightId", &jetsak4_tightId, &b_jetsak4_tightId);
	fChain->SetBranchAddress("jetsak4_lepVeto", &jetsak4_lepVeto, &b_jetsak4_lepVeto);
	fChain->SetBranchAddress("jetsak8_px", &jetsak8_px, &b_jetsak8_px);
	fChain->SetBranchAddress("jetsak8_py", &jetsak8_py, &b_jetsak8_py);
	fChain->SetBranchAddress("jetsak8_pz", &jetsak8_pz, &b_jetsak8_pz);
	fChain->SetBranchAddress("jetsak8_pt", &jetsak8_pt, &b_jetsak8_pt);
	fChain->SetBranchAddress("jetsak8_energy", &jetsak8_energy, &b_jetsak8_energy);
	fChain->SetBranchAddress("jetsak8_phi", &jetsak8_phi, &b_jetsak8_phi);
	fChain->SetBranchAddress("jetsak8_eta", &jetsak8_eta, &b_jetsak8_eta);
	fChain->SetBranchAddress("jetsak8_vx", &jetsak8_vx, &b_jetsak8_vx);
	fChain->SetBranchAddress("jetsak8_vy", &jetsak8_vy, &b_jetsak8_vy);
	fChain->SetBranchAddress("jetsak8_vz", &jetsak8_vz, &b_jetsak8_vz);
	fChain->SetBranchAddress("jetsak8_bdis", &jetsak8_bdis, &b_jetsak8_bdis);
	fChain->SetBranchAddress("jetsak8_looseId", &jetsak8_looseId, &b_jetsak8_looseId);
	fChain->SetBranchAddress("jetsak8_tightId", &jetsak8_tightId, &b_jetsak8_tightId);
	fChain->SetBranchAddress("jetsak8_lepVeto", &jetsak8_lepVeto, &b_jetsak8_lepVeto);
	fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
	fChain->SetBranchAddress("nElectrons_looseId", &nElectrons_looseId, &b_nElectrons_looseId);
	fChain->SetBranchAddress("nElectrons_mediumId", &nElectrons_mediumId, &b_nElectrons_mediumId);
	fChain->SetBranchAddress("nElectrons_tightId", &nElectrons_tightId, &b_nElectrons_tightId);
	fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
	fChain->SetBranchAddress("nMuons_looseId", &nMuons_looseId, &b_nMuons_looseId);
	fChain->SetBranchAddress("nMuons_mediumId", &nMuons_mediumId, &b_nMuons_mediumId);
	fChain->SetBranchAddress("nMuons_tightId", &nMuons_tightId, &b_nMuons_tightId);
 	fChain->SetBranchAddress("nChargedPFMultiPV_Loose", &nChargedPFMultiPV_Loose, &b_nChargedPFMultiPV_Loose);
 	fChain->SetBranchAddress("nChargedPFMultiPV_Tight", &nChargedPFMultiPV_Tight, &b_nChargedPFMultiPV_Tight);
 	fChain->SetBranchAddress("nChargedPFMultiPV_UsedInFit", &nChargedPFMultiPV_UsedInFit, &b_nChargedPFMultiPV_UsedInFit);
 	fChain->SetBranchAddress("nChargedPFMultiPV_Tight_Fit", &nChargedPFMultiPV_Tight_Fit, &b_nChargedPFMultiPV_Tight_Fit);
        fChain->SetBranchAddress("SumChargedPFMultiPV_pt_Loose", &SumChargedPFMultiPV_pt_Loose, &b_SumChargedPFMultiPV_pt_Loose);
        fChain->SetBranchAddress("SumChargedPFMultiPV_pt_Tight", &SumChargedPFMultiPV_pt_Tight, &b_SumChargedPFMultiPV_pt_Tight);
        fChain->SetBranchAddress("SumChargedPFMultiPV_pt_UsedInFit", &SumChargedPFMultiPV_pt_UsedInFit, &b_SumChargedPFMultiPV_pt_UsedInFit);
        fChain->SetBranchAddress("SumChargedPFMultiPV_pt_Tight_Fit", &SumChargedPFMultiPV_pt_Tight_Fit, &b_SumChargedPFMultiPV_pt_Tight_Fit);
 	fChain->SetBranchAddress("leptons_pdgid", &leptons_pdgid, &b_leptons_pdgid);
 	fChain->SetBranchAddress("leptons_energy", &leptons_energy, &b_leptons_energy);
	fChain->SetBranchAddress("leptons_pt", &leptons_pt, &b_leptons_pt);
	fChain->SetBranchAddress("leptons_eta", &leptons_eta, &b_leptons_eta);
	fChain->SetBranchAddress("leptons_phi", &leptons_phi, &b_leptons_phi);
	fChain->SetBranchAddress("leptons_px", &leptons_px, &b_leptons_px);
	fChain->SetBranchAddress("leptons_py", &leptons_py, &b_leptons_py);
	fChain->SetBranchAddress("leptons_pz", &leptons_pz, &b_leptons_pz);
	fChain->SetBranchAddress("leptons_charge", &leptons_charge, &b_leptons_charge);
	fChain->SetBranchAddress("leptons_vx", &leptons_vx, &b_leptons_vx);
	fChain->SetBranchAddress("leptons_vy", &leptons_vy, &b_leptons_vy);
	fChain->SetBranchAddress("leptons_vz", &leptons_vz, &b_leptons_vz);
	fChain->SetBranchAddress("leptons_looseId", &leptons_looseId, &b_leptons_looseId);
	fChain->SetBranchAddress("leptons_mediumId", &leptons_mediumId, &b_leptons_mediumId);
	fChain->SetBranchAddress("leptons_tightId", &leptons_tightId, &b_leptons_tightId);
	fChain->SetBranchAddress("leptons_pfIsoMedium_", &leptons_pfIsoMedium_, &b_leptons_pfIsoMedium_);
	fChain->SetBranchAddress("leptons_miniIsoTight_", &leptons_miniIsoTight_, &b_leptons_miniIsoTight_);
	fChain->SetBranchAddress("leptons_pfIsoVeryTight_", &leptons_pfIsoVeryTight_, &b_leptons_pfIsoVeryTight_);
	fChain->SetBranchAddress("leptons_pfIso_", &leptons_pfIso_, &b_leptons_pfIso_);
	fChain->SetBranchAddress("leptons_tkIso_", &leptons_tkIso_, &b_leptons_tkIso_);
	fChain->SetBranchAddress("missEt", &missEt, &b_missEt);
	fChain->SetBranchAddress("missEt_phi", &missEt_phi, &b_missEt_phi);
	fChain->SetBranchAddress("protonsArm", &protonsArm, &b_protonsArm);
	fChain->SetBranchAddress("protonsStation", &protonsStation, &b_protonsStation);
	fChain->SetBranchAddress("protonsRP", &protonsRP, &b_protonsRP);
	fChain->SetBranchAddress("protonsX", &protonsX, &b_protonsX);
	fChain->SetBranchAddress("protonsXUnc", &protonsXUnc, &b_protonsXUnc);
	fChain->SetBranchAddress("protonsY", &protonsY, &b_protonsY);
	fChain->SetBranchAddress("protonsYUnc", &protonsYUnc, &b_protonsYUnc);
	fChain->SetBranchAddress("singleProtonArm", &singleProtonArm, &b_singleProtonArm);
	fChain->SetBranchAddress("singleProtonStation", &singleProtonStation, &b_singleProtonStation);
	fChain->SetBranchAddress("singleProtonPot", &singleProtonPot, &b_singleProtonPot);
	fChain->SetBranchAddress("singleProtonXi", &singleProtonXi, &b_singleProtonXi);
	fChain->SetBranchAddress("singleProtonThetaX", &singleProtonThetaX, &b_singleProtonThetaX);
	fChain->SetBranchAddress("singleProtonThetaY", &singleProtonThetaY, &b_singleProtonThetaY);
	fChain->SetBranchAddress("multiProtonArm", &multiProtonArm, &b_multiProtonArm);
	fChain->SetBranchAddress("multiProtonXi", &multiProtonXi, &b_multiProtonXi);
	fChain->SetBranchAddress("multiProtonTime", &multiProtonTime, &b_multiProtonTime);
	fChain->SetBranchAddress("multiProtonTimeError", &multiProtonTimeError, &b_multiProtonTimeError);
	fChain->SetBranchAddress("multiProtonThetaX", &multiProtonThetaX, &b_multiProtonThetaX);
	fChain->SetBranchAddress("multiProtonThetaY", &multiProtonThetaY, &b_multiProtonThetaY);
	Notify();
}

Bool_t MissingMassNtupleAnalyzer::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}

void MissingMassNtupleAnalyzer::Show(Long64_t entry)
{
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t MissingMassNtupleAnalyzer::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

/*
void sortArr(std::vector<double> &v) 
{ 
  
    vector<pair<int, int> > vp; 
  
    // Inserting element in pair vector 
    // to keep track of previous indexes 
    for (int i = 0; i < n; ++i) { 
        vp.push_back(make_pair(arr[i], i)); 
    } 
  
    // Sorting pair vector 
    sort(vp.begin(), vp.end()); 
  
    // Displaying sorted element 
    // with previous indexes 
    // corresponding to each element 
    cout << "Element\t"
         << "index" << endl; 
    for (int i = 0; i < vp.size(); i++) { 
        cout << vp[i].first << "\t"
             << vp[i].second << endl; 
    } 
} 
*/

#endif // #ifdef MissingMassNtupleAnalyzer_cxx

