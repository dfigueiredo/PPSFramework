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
#include "ProtonsRandom.h"
#include "FilesInput.h"
#include "TTreeMissingMass.h"

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

  //ProtonsRandom protonsR;
  //protonsR.createHistogram();

  FilesInput FileInput;
  std::list<TString> filestr1, filestr2, filestr3;

  filestr1 = {TString(outdir), TString(mode), TString(era)};
  filestr2 = {TString(outdir), TString(mode), TString(jobid)};
  filestr3 = {TString(outdir), TString(mode), TString(xa)};

  TString filenameout1 = FileInput.GetFileName(filestr1); 
  TString filenameout2 = FileInput.GetFileName(filestr2); 
  TString filenameout3 = FileInput.GetFileName(filestr3); 

  std::cout << "filenameout1: " << filenameout1 << std::endl;
  std::cout << "filenameout2: " << filenameout2 << std::endl;
  std::cout << "filenameout3: " << filenameout3 << std::endl;

  // debugging variables (printout values)
  bool debug = false;

  TString filenameout;
  TString filenameout_random_reco;
  TString filenameout_eventlist;
  TString filenameout_eventlist_kinematics;

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

  TTreeMissingMass ttreeAnalysis;
  ttreeAnalysis.switchZeroBias = false;
  ttreeAnalysis.CreateTTree();

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

    ttreeAnalysis.Clearing();

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

    if((strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0)){
      // Era flag
      int random_era;
      if(std::string(era)=="A"||std::string(era)=="a") ttreeAnalysis.era = 0;
      else if(std::string(era)=="B"||std::string(era)=="b") ttreeAnalysis.era=1;
      else if(std::string(era)=="C"||std::string(era)=="c") ttreeAnalysis.era=2;
      else if(std::string(era)=="D"||std::string(era)=="d") ttreeAnalysis.era=3;
      else if(std::string(era)=="E"||std::string(era)=="e") ttreeAnalysis.era=4;
      else if(std::string(era)=="F"||std::string(era)=="f") ttreeAnalysis.era=5;
      else if(std::string(era)=="preTS2"||std::string(era)=="pre TS2"){
	random_era = rand.Integer(3) + 1; //preTS2 refers to era B, C and era D
	ttreeAnalysis.era = random_era;
      }
      else if(std::string(era)=="postTS2"||std::string(era)=="post TS2"){
	random_era = rand.Integer(2) + 4; //post TS2 refers to era E and F
	ttreeAnalysis.era = random_era;
      }
      else ttreeAnalysis.era=-1;

      if(std::string(xa) == "120") ttreeAnalysis.xangle = 120;
      else if(std::string(xa) == "130") ttreeAnalysis.xangle = 130;
      else if(std::string(xa) == "140") ttreeAnalysis.xangle = 140;
      else if(std::string(xa) == "150") ttreeAnalysis.xangle = 150;
    }

    // Loop over protons with single RP
    for(std::vector<int>::size_type i = 0; i != singleProtonArm->size(); i++){
      if(singleProtonStation->at(i)==0&&singleProtonPot->at(i)==3){
	if(singleProtonArm->at(i)==0){
	  ttreeAnalysis.xi_rp210_Arm45 = singleProtonXi->at(i);
	  ttreeAnalysis.nprotonRP210_sec45++;
	  if(ttreeAnalysis.xi_rp210_Arm45 > 0.04) ++ttreeAnalysis.nprotonRP210_sec45_Reduced;
	}
	if(singleProtonArm->at(i)==1){
	  ttreeAnalysis.xi_rp210_Arm56 = singleProtonXi->at(i);
	  ttreeAnalysis.nprotonRP210_sec56++;
	  if(ttreeAnalysis.xi_rp210_Arm56 > 0.04) ++ttreeAnalysis.nprotonRP210_sec56_Reduced;
	}
      }
      if(singleProtonStation->at(i)==2&&singleProtonPot->at(i)==3){
	if(singleProtonArm->at(i)==0){
	  ttreeAnalysis.xi_rp220_Arm45 = singleProtonXi->at(i);
	  ttreeAnalysis.nprotonRP220_sec45++;
	  if(ttreeAnalysis.xi_rp220_Arm45 > 0.04) ++ttreeAnalysis.nprotonRP220_sec45_Reduced;
	}
	if(singleProtonArm->at(i)==1){
	  ttreeAnalysis.xi_rp220_Arm56 = singleProtonXi->at(i);
	  ttreeAnalysis.nprotonRP220_sec56++;
	  if(ttreeAnalysis.xi_rp220_Arm56 > 0.04) ++ttreeAnalysis.nprotonRP220_sec56_Reduced;
	}
      }
    }

    // Loop over protons with multiple RP
    for(std::vector<int>::size_type i = 0; i != multiProtonArm->size(); i++){
      if(multiProtonArm->at(i)==0){
	ttreeAnalysis.xi_multiArm45 = multiProtonXi->at(i);
	ttreeAnalysis.time_multiArm45 = multiProtonTime->at(i);
	ttreeAnalysis.timeunc_multiArm45 = multiProtonTimeError->at(i);
	ttreeAnalysis.thx_multiArm45 = multiProtonThetaX->at(i);
	ttreeAnalysis.thy_multiArm45 = multiProtonThetaY->at(i);
	ttreeAnalysis.nprotonMulti_sec45++;
      }
      if(multiProtonArm->at(i)==1){
	ttreeAnalysis.xi_multiArm56 = multiProtonXi->at(i);
	ttreeAnalysis.time_multiArm56 = multiProtonTime->at(i);
	ttreeAnalysis.timeunc_multiArm56 = multiProtonTimeError->at(i);
	ttreeAnalysis.thx_multiArm56 = multiProtonThetaX->at(i);
	ttreeAnalysis.thy_multiArm56 = multiProtonThetaY->at(i);
	ttreeAnalysis.nprotonMulti_sec56++;
      }
    }

    if(ttreeAnalysis.nprotonRP210_sec45==1 && ttreeAnalysis.nprotonRP210_sec56==1) ttreeAnalysis.isprotonRP210 = true;
    if(ttreeAnalysis.nprotonRP220_sec45==1 && ttreeAnalysis.nprotonRP220_sec56==1) ttreeAnalysis.isprotonRP220 = true;
    if(ttreeAnalysis.nprotonRP210_sec45_Reduced==1 && ttreeAnalysis.nprotonRP210_sec56_Reduced==1) ttreeAnalysis.isprotonRP210_Reduced = true;
    if(ttreeAnalysis.nprotonRP220_sec45_Reduced==1 && ttreeAnalysis.nprotonRP220_sec56_Reduced==1) ttreeAnalysis.isprotonRP220_Reduced = true;
    if(ttreeAnalysis.nprotonMulti_sec45==1 && ttreeAnalysis.nprotonMulti_sec56==1) ttreeAnalysis.isprotonMulti = true;
    if(!protonsfilter){
      ttreeAnalysis.isprotonRP210 = true;
      ttreeAnalysis.isprotonRP220 = true;
      ttreeAnalysis.isprotonMulti = true;
      ttreeAnalysis.isprotonRP210_Reduced = true;
      ttreeAnalysis.isprotonRP220_Reduced = true;
    }

    if((strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0) || (strcmp(mode, "Muon")==0 || strcmp(mode, "Electron")==0 || strcmp(mode, "Bjets")==0)){
      if(randomFlag){
	if(strcmp(datatype, "DY")==0){
	  ttreeAnalysis.xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx];
	  ttreeAnalysis.xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	  ttreeAnalysis.xi_multiArm45 = rndXi_Multi_sec45[random_idx];
	  ttreeAnalysis.xi_rp210_Arm56 = rndXi_RP210_sec56[random_idx];
	  ttreeAnalysis.xi_rp220_Arm56 = rndXi_RP220_sec56[random_idx];
	  ttreeAnalysis.xi_multiArm56 = rndXi_Multi_sec56[random_idx];
	  ttreeAnalysis.is2PUproton = true;
	  ttreeAnalysis.isprotonRP210 = true;
	  ttreeAnalysis.isprotonMulti = true;	  
	}
	else if(strcmp(datatype, "Signal")==0 || strcmp(datatype, "SD")==0){
	  //Correcting for the efficiency
	  //Reading Efficiency File
	  TFile *f_eff = TFile::Open("PreliminaryEfficiencies_October92019_1D2DMultiTrack.root");

	  Double_t rad_eff_45[60];
	  Double_t rad_eff_56[60];
	  TH1D *h_rad_45;
	  TH1D *h_rad_56;

	  if (ttreeAnalysis.era == 1){
	    if (ttreeAnalysis.xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017B/h45_2017B_120_1D");     
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017B/h56_2017B_120_1D");    
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017B/h45_2017B_130_1D");  
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017B/h56_2017B_130_1D"); 
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }   
	    }
	    else if (ttreeAnalysis.xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017B/h45_2017B_140_1D");  
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017B/h56_2017B_140_1D"); 
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 150){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017B/h45_2017B_150_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017B/h56_2017B_150_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	  }
	  else if (ttreeAnalysis.era == 2){
	    if (ttreeAnalysis.xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017C/h45_2017C_120_1D");  
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017C/h56_2017C_120_1D"); 
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017C/h45_2017C_130_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017C/h56_2017C_130_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017C/h45_2017C_140_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017C/h56_2017C_140_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 150){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017C/h45_2017C_150_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017C/h56_2017C_150_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	  }
	  else if (ttreeAnalysis.era == 3){
	    if (ttreeAnalysis.xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017D/h45_2017D_120_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017D/h56_2017D_120_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017D/h45_2017D_130_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017D/h56_2017D_130_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017D/h45_2017D_140_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017D/h56_2017D_140_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 150){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017D/h45_2017D_150_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017D/h56_2017D_150_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	  }
	  else if (ttreeAnalysis.era == 4){
	    if (ttreeAnalysis.xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017E/h45_2017E_120_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017E/h56_2017E_120_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017E/h45_2017E_130_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017E/h56_2017E_130_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017E/h45_2017E_140_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017E/h56_2017E_140_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 150){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017E/h45_2017E_150_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017E/h56_2017E_150_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	  }
	  else if (ttreeAnalysis.era == 5){
	    if (ttreeAnalysis.xangle == 120){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017F/h45_2017F_120_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017F/h56_2017F_120_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 130){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017F/h45_2017F_130_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017F/h56_2017F_130_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 140){
	      h_rad_45 = (TH1D*)f_eff->Get("Strips/2017/2017F/h45_2017F_140_1D");
	      h_rad_56 = (TH1D*)f_eff->Get("Strips/2017/2017F/h56_2017F_140_1D");
	      for(int i=0; i < 60; i++){
		rad_eff_45[i] = h_rad_45->GetBinContent(i+1);
		rad_eff_56[i] = h_rad_56->GetBinContent(i+1);
	      }
	    }
	    else if (ttreeAnalysis.xangle == 150){
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
	  if (ttreeAnalysis.nprotonRP210_sec45>1 || ttreeAnalysis.nprotonRP220_sec56>1){
	    std::cout << "More than 1 strip hit" << std::endl;
	  }
	  if (ttreeAnalysis.nprotonRP210_sec45 > 0){
	    random_eff_45 = rand.Rndm();	  
	    for (int j = 0; j < 60; j++){         
	      if(ttreeAnalysis.xi_rp210_Arm45 > (j * 0.005) && ttreeAnalysis.xi_rp210_Arm45 < ((j * 0.005)+0.005)){
		if(debug){
		  std::cout << "####### 45 #######"<< std::endl;
		  std::cout << "era\t" << ttreeAnalysis.era << std::endl;
		  std::cout << "xi_45\t" << ttreeAnalysis.xi_rp210_Arm45 << std::endl;
		  std::cout << "random_45\t" << random_eff_45 << std::endl;
		  std::cout << "rad_eff_45\t" << rad_eff_45[j] << std::endl;
		  std::cout << "##################"<< std::endl;
		}
		if (random_eff_45 > rad_eff_45[j]){   
		  ttreeAnalysis.nprotonRP210_sec45 = 0;
		  ttreeAnalysis.nprotonMulti_sec45 = 0;
		  ttreeAnalysis.xi_rp210_Arm45 = 0;
		  ttreeAnalysis.xi_multiArm45 = 0;
		  ttreeAnalysis.isprotonRP210 = false;
		  ttreeAnalysis.isprotonMulti = false;
		  goto end45;
		}  
		else {goto end45;}
	      }
	    }
	  }
end45:;

      if (ttreeAnalysis.nprotonRP210_sec56 > 0){
	random_eff_56 = rand.Rndm();
	for (int j = 0; j < 60; j++){
	  if(ttreeAnalysis.xi_rp210_Arm56 > (j * 0.005) && ttreeAnalysis.xi_rp210_Arm56 < ((j * 0.005)+0.005)){
	    if(debug){
	      std::cout << "####### 56 #######"<< std::endl;
	      std::cout << "era\t" << ttreeAnalysis.era << std::endl;
	      std::cout << "xi_56\t" << ttreeAnalysis.xi_rp210_Arm56 << std::endl;
	      std::cout << "random_56\t" << random_eff_56 << std::endl;
	      std::cout << "rad_eff_56\t" << rad_eff_56[j] << std::endl;
	      std::cout << "##################"<< std::endl;
	    }
	    if (random_eff_56 > rad_eff_56[j]){
	      ttreeAnalysis.nprotonRP210_sec56 = 0;
	      ttreeAnalysis.nprotonMulti_sec56 = 0;
	      ttreeAnalysis.isprotonRP210 = false;
	      ttreeAnalysis.isprotonMulti = false;
	      ttreeAnalysis.xi_rp210_Arm56 = 0;
	      ttreeAnalysis.xi_multiArm56 = 0;
	      goto end56;
	    }
	    else {goto end56;}
	  }
	}
      }
end56:;

      //Out of acceptance background
      if(ttreeAnalysis.nprotonRP210_sec45 == 0 && ttreeAnalysis.nprotonRP210_sec56 > 0){
	ttreeAnalysis.xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx];
	ttreeAnalysis.xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	ttreeAnalysis.xi_multiArm45  = rndXi_Multi_sec45[random_idx];
	ttreeAnalysis.nprotonRP210_sec45 = 1;
	ttreeAnalysis.nprotonRP220_sec45++;
	ttreeAnalysis.nprotonMulti_sec45 = 1;
	ttreeAnalysis.is45PUproton = true;
	ttreeAnalysis.isprotonRP210 = true;
      }
      else if(ttreeAnalysis.nprotonRP210_sec56 == 0 && ttreeAnalysis.nprotonRP210_sec45 > 0){
	ttreeAnalysis.xi_rp210_Arm56 = rndXi_RP210_sec56[random_idx];
	ttreeAnalysis.xi_rp220_Arm56 = rndXi_RP220_sec56[random_idx];
	ttreeAnalysis.xi_multiArm56  = rndXi_Multi_sec56[random_idx];
	ttreeAnalysis.nprotonRP210_sec56 = 1;
	ttreeAnalysis.nprotonRP220_sec56++;
	ttreeAnalysis.nprotonMulti_sec56 = 1;
	ttreeAnalysis.is56PUproton = true;
	ttreeAnalysis.isprotonRP210 = true;
      }
      else if(ttreeAnalysis.nprotonRP210_sec56 == 0 && ttreeAnalysis.nprotonRP210_sec45 == 0){
	ttreeAnalysis.xi_rp210_Arm56 = rndXi_RP210_sec56[random_idx];
	ttreeAnalysis.xi_rp220_Arm56 = rndXi_RP220_sec56[random_idx];
	ttreeAnalysis.xi_multiArm56 = rndXi_Multi_sec56[random_idx];
	ttreeAnalysis.xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx];
	ttreeAnalysis.xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	ttreeAnalysis.xi_multiArm45 = rndXi_Multi_sec45[random_idx];
	ttreeAnalysis.nprotonRP210_sec45 = 1;
	ttreeAnalysis.nprotonRP220_sec45++;
	ttreeAnalysis.nprotonMulti_sec45 = 1;
	ttreeAnalysis.nprotonRP210_sec56 = 1;
	ttreeAnalysis.nprotonRP220_sec56++;
	ttreeAnalysis.nprotonMulti_sec56 = 1;
	ttreeAnalysis.is2PUproton = true;
	ttreeAnalysis.isprotonRP210 = true;
      }
      else if(ttreeAnalysis.nprotonRP210_sec56 == 1 && ttreeAnalysis.nprotonRP210_sec45 == 1){
	ttreeAnalysis.is2PUproton = false;
	ttreeAnalysis.is45PUproton = false;
	ttreeAnalysis.is56PUproton = false;
	ttreeAnalysis.isprotonRP210 = true;	    
      }
	}
	if(ttreeAnalysis.nprotonRP210_sec45==1 && ttreeAnalysis.nprotonRP210_sec56==1) ttreeAnalysis.isprotonRP210 = true;
	if(ttreeAnalysis.nprotonRP220_sec45==1 && ttreeAnalysis.nprotonRP220_sec56==1) ttreeAnalysis.isprotonRP220 = true;
	if(ttreeAnalysis.nprotonMulti_sec45==1 && ttreeAnalysis.nprotonMulti_sec56==1) ttreeAnalysis.isprotonMulti = true;
	if(single){
	  ttreeAnalysis.xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx]; 
	  ttreeAnalysis.xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	  ttreeAnalysis.xi_multiArm45 = rndXi_Multi_sec45[random_idx];
	  ttreeAnalysis.isprotonRP210 = true;
	  ttreeAnalysis.isprotonRP220 = true;
	  ttreeAnalysis.isprotonMulti = true;
	}
	else{
	  ttreeAnalysis.xi_rp210_Arm45 = rndXi_RP210_sec45[random_idx];
	  ttreeAnalysis.xi_rp220_Arm45 = rndXi_RP220_sec45[random_idx];
	  ttreeAnalysis.xi_multiArm45 = rndXi_Multi_sec45[random_idx];
	  ttreeAnalysis.xi_rp210_Arm56 = rndXi_RP210_sec56[random_idx];
	  ttreeAnalysis.xi_rp220_Arm56 = rndXi_RP220_sec56[random_idx];
	  ttreeAnalysis.xi_multiArm56 = rndXi_Multi_sec56[random_idx];
	  ttreeAnalysis.isprotonRP210 = true;
	  ttreeAnalysis.isprotonRP220 = true;
	  ttreeAnalysis.isprotonMulti = true;
	}
	if(debug){
	  std::cout << "proton pixels Xi, 210, arm 45: " << ttreeAnalysis.xi_rp210_Arm45 << std::endl;
	  std::cout << "proton pixels Xi, 210, arm 56: " << ttreeAnalysis.xi_rp210_Arm56 << std::endl;
	  std::cout << "proton pixels Xi, 220, arm 45: " << ttreeAnalysis.xi_rp220_Arm45 << std::endl;
	  std::cout << "proton pixels Xi, 220, arm 56: " << ttreeAnalysis.xi_rp220_Arm56 << std::endl;
	  std::cout << "proton multi Xi, arm 45: " << ttreeAnalysis.xi_multiArm45 << std::endl;
	  std::cout << "proton multi Xi, arm 56: " << ttreeAnalysis.xi_multiArm56 << std::endl;
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

      p1RP210.SetPxPyPzE(0.,0., ECM*ttreeAnalysis.xi_rp210_Arm45/2., ECM*ttreeAnalysis.xi_rp210_Arm45/2.);
      p2RP210.SetPxPyPzE(0.,0., -ECM*ttreeAnalysis.xi_rp210_Arm56/2., ECM*ttreeAnalysis.xi_rp210_Arm56/2.);

      p1RP220.SetPxPyPzE(0.,0., ECM*ttreeAnalysis.xi_rp220_Arm45/2., ECM*ttreeAnalysis.xi_rp220_Arm45/2.);
      p2RP220.SetPxPyPzE(0.,0., -ECM*ttreeAnalysis.xi_rp220_Arm56/2., ECM*ttreeAnalysis.xi_rp220_Arm56/2.);

      p1Multi.SetPxPyPzE(0.,0., ECM*ttreeAnalysis.xi_multiArm45/2., ECM*ttreeAnalysis.xi_multiArm45/2.);
      p2Multi.SetPxPyPzE(0.,0., -ECM*ttreeAnalysis.xi_multiArm56/2., ECM*ttreeAnalysis.xi_multiArm56/2.);

      ttreeAnalysis.diffMassRP210 = ECM*sqrt(ttreeAnalysis.xi_rp210_Arm45*ttreeAnalysis.xi_rp210_Arm56); 
      ttreeAnalysis.diffMassRP220 = ECM*sqrt(ttreeAnalysis.xi_rp220_Arm45*ttreeAnalysis.xi_rp220_Arm56); 
      ttreeAnalysis.diffMassMulti = ECM*sqrt(ttreeAnalysis.xi_multiArm45*ttreeAnalysis.xi_multiArm56); 

      ttreeAnalysis.proton_pz_rp210Arm45 = (1.-(ttreeAnalysis.xi_rp210_Arm45))*(ECM/2.);
      ttreeAnalysis.proton_pz_rp210Arm56 = -(1.-(ttreeAnalysis.xi_rp210_Arm56))*(ECM/2.);

      ttreeAnalysis.proton_pz_rp220Arm45 = (1.-(ttreeAnalysis.xi_rp220_Arm45))*(ECM/2.);
      ttreeAnalysis.proton_pz_rp220Arm56 = -(1.-(ttreeAnalysis.xi_rp220_Arm56))*(ECM/2.);

      ttreeAnalysis.proton_pz_multiArm45 = (1.-(ttreeAnalysis.xi_multiArm45))*(ECM/2.);
      ttreeAnalysis.proton_pz_multiArm56 = -(1.-(ttreeAnalysis.xi_multiArm56))*(ECM/2.);

      ttreeAnalysis.run = run;
      ttreeAnalysis.event = event;
      ttreeAnalysis.lumiblock = lumiblock;
      ttreeAnalysis.xangle = xangle;

      //TTreeTest.test1_ = 10;
      //TTreeTest.test2_ = 20;
      //TTreeTest.test2_ = 20;

      std::vector<int> index_leptons;
      index_leptons.clear();
      for (auto i: sort_indexes(*leptons_pt)) {
	index_leptons.push_back(i);
      }

      // Filling TTree                                                           
      if(strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0){
	ttreeAnalysis.PUInterac = PUInterac;
	ttreeAnalysis.PUTrueInterac = PUTrueInterac;
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

      // HLT Muons: 'HLT_IsoMu27_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*'
      // HLT Electrons: 'HLT_Ele27_WPTight_Gsf_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*', 'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*', 'HLT_DoubleEle33_CaloIdL_MW_v*'
      // HLT_ZeroBias: 'HLT_ZeroBias_v*', 'HLT_ZeroBias_FirstBXAfterTrain_v*','HLT_ZeroBias_IsolatedBunches_v*','HLT_ZeroBias_Alignment_v*','HLT_ZeroBias_Beamspot_v*','HLT_L1UnpairedBunchBptxMinus_v*','HLT_L1UnpairedBunchBptxPlus_v*','HLT_Physics_v*','HLT_L1SingleMu_v*'
      // If -1: not present in the path, 0: not triggered and 1 triggered.

      if(zerobias){
	if(trigger->at(0)==1){
	  ttreeAnalysis.triggerZeroBias = true;
	  ttreeAnalysis.prescalesL1ZeroBias = prescalesL1->at(0);
	}
	if(trigger->at(1)==1){
	  ttreeAnalysis.triggerZeroBiasAfterTrain = true;
	  ttreeAnalysis.prescalesL1ZeroBiasAfterTrain = prescalesL1->at(1);
	}
	if(trigger->at(2)==1){
	  ttreeAnalysis.triggerZeroBiasIsolatedBx = true;
	  ttreeAnalysis.prescalesL1ZeroBiasIsolatedBx = prescalesL1->at(2);
	}
	if(trigger->at(3)==1){
	  ttreeAnalysis.triggerZeroBiasAlignment = true;
	  ttreeAnalysis.prescalesL1ZeroBiasAlignment = prescalesL1->at(3);
	}
	if(trigger->at(4)==1){
	  ttreeAnalysis.triggerZeroBiasBeamSpot = true;
	  ttreeAnalysis.prescalesL1ZeroBiasBeamSpot = prescalesL1->at(4);
	}
	if(trigger->at(5)==1){
	  ttreeAnalysis.triggerZeroBiasUnpairedBptxMinus = true;
	  ttreeAnalysis.prescalesL1ZeroBiasUnpairedBptxMinus = prescalesL1->at(5);
	}
	if(trigger->at(6)==1){
	  ttreeAnalysis.triggerZeroBiasUnpairedBptxPlus = true;
	  ttreeAnalysis.prescalesL1ZeroBiasUnpairedBptxPlus = prescalesL1->at(6);
	}
	if(trigger->at(7)==1){
	  ttreeAnalysis.triggerPhysics = true;
	  ttreeAnalysis.prescalesL1Physics = prescalesL1->at(7);
	}
	if(trigger->at(8)==1){
	  ttreeAnalysis.triggerL1SingleMu = true;
	  ttreeAnalysis.prescalesL1SingleMu = prescalesL1->at(8);
	}
      }else{
	if((strcmp(mode, "Muon")==0 || strcmp(mode, "MC_Muon")==0) && strcmp(datatype, "SD")!=0 && strcmp(datatype, "Signal")!=0){
	  if(trigger->at(0)==1) ttreeAnalysis.triggerIsoMu27 = true;
	  if(trigger->at(1)==1) ttreeAnalysis.triggerMu17TrkIsoMu8TrkIso = true;
	  if(trigger->at(2)==1) ttreeAnalysis.triggerMu17TrkIsoMu8TrkIsoMass8 = true;
	  if(trigger->at(3)==1) ttreeAnalysis.triggerMu17TrkIsoMu8TrkIsoMass3 = true;
	  if(trigger->at(4)==1) ttreeAnalysis.triggerDoubleMu43 = true;
	  if(trigger->at(5)==1) ttreeAnalysis.triggerDoubleMu48 = true;
	}else if(strcmp(mode, "Electron")==0 || strcmp(mode, "MC_Electron")==0 && strcmp(datatype, "SD")!=0 && strcmp(datatype, "Signal")!=0){
	  if(trigger->at(0)==1) ttreeAnalysis.triggerEle27 = true;
	  if(trigger->at(1)==1) ttreeAnalysis.triggerEle23Ele12 = true;
	  if(trigger->at(2)==1) ttreeAnalysis.triggerEle23Ele12Dz = true;
	  if(trigger->at(3)==1) ttreeAnalysis.triggerDoubleEle33 = true;
	  if(trigger->at(4)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet20 = true;
	  if(trigger->at(5)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet40 = true;
	  if(trigger->at(6)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet70 = true;
	  if(trigger->at(7)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet110 = true;
	  if(trigger->at(8)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet170 = true;
	  if(trigger->at(9)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet300 = true;
	  if(trigger->at(10)==1) ttreeAnalysis.triggerBTagMu5Ak8dijet170 = true;
	  if(trigger->at(11)==1) ttreeAnalysis.triggerBTagMu5Ak8dijet300 = true;
	  if(trigger->at(12)==1) ttreeAnalysis.triggerPFHT380DoubleBTag = true;
	  if(trigger->at(13)==1) ttreeAnalysis.triggerPFHT430DoubleBTag = true;
	  if(trigger->at(14)==1) ttreeAnalysis.triggerPFMET100BTag = true;
	  if(trigger->at(15)==1) ttreeAnalysis.triggerPFMET110BTag = true;
	  if(trigger->at(16)==1) ttreeAnalysis.triggerPFMET120BTag = true;
	  if(trigger->at(17)==1) ttreeAnalysis.triggerPFMET130BTag = true;
	  if(trigger->at(18)==1) ttreeAnalysis.triggerPFMET140BTag = true;
	  if(trigger->at(19)==1) ttreeAnalysis.triggerEle15PFHT450BTag = true;
	  if(trigger->at(20)==1) ttreeAnalysis.triggerHT250BTagScouting = true;
	}else if(strcmp(mode, "Bjets")==0 || (strcmp(datatype, "SD")!=0 && strcmp(datatype, "Signal")!=0)){
	  if(trigger->at(0)==1) ttreeAnalysis.triggerIsoMu27 = true;
	  if(trigger->at(1)==1) ttreeAnalysis.triggerMu17TrkIsoMu8TrkIso = true;
	  if(trigger->at(2)==1) ttreeAnalysis.triggerMu17TrkIsoMu8TrkIsoMass8 = true;
	  if(trigger->at(3)==1) ttreeAnalysis.triggerMu17TrkIsoMu8TrkIsoMass3 = true;
	  if(trigger->at(4)==1) ttreeAnalysis.triggerDoubleMu43 = true;
	  if(trigger->at(5)==1) ttreeAnalysis.triggerDoubleMu48 = true;
	  if(trigger->at(6)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet20 = true;
	  if(trigger->at(7)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet40 = true;
	  if(trigger->at(8)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet70 = true;
	  if(trigger->at(9)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet110 = true;
	  if(trigger->at(10)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet170 = true;
	  if(trigger->at(11)==1) ttreeAnalysis.triggerBTagMu5Ak4dijet300 = true;
	  if(trigger->at(12)==1) ttreeAnalysis.triggerBTagMu5Ak8dijet170 = true;
	  if(trigger->at(13)==1) ttreeAnalysis.triggerBTagMu5Ak8dijet300 = true;
	  if(trigger->at(14)==1) ttreeAnalysis.triggerPFHT380DoubleBTag = true;
	  if(trigger->at(15)==1) ttreeAnalysis.triggerPFHT430DoubleBTag = true;
	}else if(strcmp(mode, "EMu")==0){
	  if(trigger->at(0)==1) ttreeAnalysis.triggerMu23TrkIsoEle12 = true;        
	  if(trigger->at(1)==1) ttreeAnalysis.triggerMu23TrkIsoEle12DZ = true;
	  if(trigger->at(2)==1) ttreeAnalysis.triggerMu8TrkIsoEle23DZ = true;
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

      if(!(ttreeAnalysis.isprotonRP210 || ttreeAnalysis.isprotonRP220 || ttreeAnalysis.isprotonMulti)) continue;
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

      ttreeAnalysis.acoplanarity = 1. - fabs(lepton1.Phi()-lepton2.Phi())/pi;

      // Leptons and Jets are already sorted by pT
      for(std::vector<int>::size_type i = 0; i != leptons_pt->size(); i++){
	lepton.SetPtEtaPhiE(leptons_pt->at(i), leptons_eta->at(i), leptons_phi->at(i), leptons_energy->at(i));
	leptonsystem+=lepton;
	ttreeAnalysis.leptonsystemCharge+=leptons_charge->at(i);
      }

      ttreeAnalysis.nJetsCandidatesLoose = 0;
      ttreeAnalysis.nJetsCandidatesMedium = 0;
      ttreeAnalysis.nJetsCandidatesTight = 0;
      if(leptons_pt->size()>1){
	for(std::vector<int>::size_type i = 0; i != jetsak4_pt->size(); i++){
	  jet.SetPtEtaPhiE(jetsak4_pt->at(i), jetsak4_eta->at(i), jetsak4_phi->at(i), jetsak4_energy->at(i));
	  jetsystem+=jet;
	  // Counting Jets Close to the Leading and SubLeading Lepton (200 microns)
	  if(fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(0)))<0.02&&fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(1)))<0.02&&jetsak4_pt->at(i)>15.){
	    ttreeAnalysis.nJetsCandidatesTight++;
	    jetcandidate.SetPtEtaPhiE(jetsak4_pt->at(i), jetsak4_eta->at(i), jetsak4_phi->at(i), jetsak4_energy->at(i));
	    jetcandidatesystem+=jetcandidate;
	  }
	  if(fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(0)))<0.05&&fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(1)))<0.05&&jetsak4_pt->at(i)>15.){
	    ttreeAnalysis.nJetsCandidatesMedium++;
	  }
	  if(fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(0)))<0.1&&fabs(jetsak4_vz->at(i)-leptons_vz->at(index_leptons.at(1)))<0.1&&jetsak4_pt->at(i)>15.){
	    ttreeAnalysis.nJetsCandidatesLoose++;
	  }
	}
      }

      dilepton.SetPtEtaPhiM(dilepton.Pt(), dilepton.Eta(), dilepton.Phi(), dilepton.M());
      met.SetPtEtaPhiE(missEt, 0, missEt_phi, 0);
      ttreeAnalysis.Rmpf = 1.+(met.Px() * dilepton.Px() + met.Py() * dilepton.Py()) / (dilepton.Pt()*dilepton.Pt());

      if(jetsak4_pt->size()>0){
	jet1.SetPtEtaPhiE(jetsak4_pt->at(0), jetsak4_eta->at(0), jetsak4_phi->at(0), jetsak4_energy->at(0));
	ttreeAnalysis.JZBalance = jet1.Pt() - dilepton.Pt();
      }

      if(jetsak4_pt->size()>1){
	jet2.SetPtEtaPhiE(jetsak4_pt->at(1), jetsak4_eta->at(1), jetsak4_phi->at(1), jetsak4_energy->at(1));
	dijet = jet1 + jet2;
	ttreeAnalysis.JZBalance2 = (dijet).Pt() - dilepton.Pt();
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

      ttreeAnalysis.nVertex = vertex_z->size();
      ttreeAnalysis.nTracksPerVertex = vertex_ntrack->at(0);
      ttreeAnalysis.pvVertexX = vertex_x->at(0);
      ttreeAnalysis.pvVertexY = vertex_y->at(0);
      ttreeAnalysis.pvVertexZ = vertex_z->at(0);
      ttreeAnalysis.pvVertexR = sqrt(vertex_x->at(0) * vertex_x->at(0) + vertex_y->at(0) * vertex_y->at(0));
      ttreeAnalysis.pvVertexPhi = atan2(vertex_y->at(0), vertex_x->at(0));

      for (auto evtVz : *vertex_z){
	ttreeAnalysis.allVertexZ.push_back(evtVz);
      }

      if(strcmp(mode, "MC_Muon")==0 || strcmp(mode, "MC_Electron")==0){
	ttreeAnalysis.nGenLeptons = genleptons_pt->size();
	ttreeAnalysis.nGenParticles = nGenParticles;
	ttreeAnalysis.nGenElectrons = nGenElectrons;
	ttreeAnalysis.nGenMuons = nGenMuons;
	if(genleptons_pt->size()>0){
	  ttreeAnalysis.genleadingLeptonPDGId = genleptons_pdgid->at(0);
	  ttreeAnalysis.genleadingLeptonEnergy = genleptons_energy->at(0);
	  ttreeAnalysis.genleadingLeptonPx = genleptons_px->at(0);
	  ttreeAnalysis.genleadingLeptonPy = genleptons_py->at(0);
	  ttreeAnalysis.genleadingLeptonPz = genleptons_pz->at(0);
	  ttreeAnalysis.genleadingLeptonPt = genleptons_pt->at(0);
	  ttreeAnalysis.genleadingLeptonEta = genleptons_eta->at(0);
	  ttreeAnalysis.genleadingLeptonPhi = genleptons_phi->at(0);
	  ttreeAnalysis.genleadingLeptonCharge = genleptons_charge->at(0);
	  ttreeAnalysis.genleadingLeptonVx = genleptons_vx->at(0);
	  ttreeAnalysis.genleadingLeptonVy = genleptons_vy->at(0);
	  ttreeAnalysis.genleadingLeptonVz = genleptons_vz->at(0);
	  ttreeAnalysis.genleadingLeptonVr = sqrt(genleptons_vx->at(0) * genleptons_vx->at(0) + genleptons_vy->at(0) * genleptons_vy->at(0));
	  ttreeAnalysis.genleadingLeptonVphi = atan2(genleptons_vy->at(0), genleptons_vx->at(0));
	}

	if(genleptons_pt->size()>1){
	  ttreeAnalysis.gensecondLeptonPDGId = genleptons_pdgid->at(1);
	  ttreeAnalysis.gensecondLeptonEnergy = genleptons_energy->at(1);
	  ttreeAnalysis.gensecondLeptonPx = genleptons_px->at(1);
	  ttreeAnalysis.gensecondLeptonPy = genleptons_py->at(1);
	  ttreeAnalysis.gensecondLeptonPz = genleptons_pz->at(1);
	  ttreeAnalysis.gensecondLeptonPt = genleptons_pt->at(1);
	  ttreeAnalysis.gensecondLeptonEta = genleptons_eta->at(1);
	  ttreeAnalysis.gensecondLeptonPhi = genleptons_phi->at(1);
	  ttreeAnalysis.gensecondLeptonCharge = genleptons_charge->at(1);
	  ttreeAnalysis.gensecondLeptonVx = genleptons_vx->at(1);
	  ttreeAnalysis.gensecondLeptonVy = genleptons_vy->at(1);
	  ttreeAnalysis.gensecondLeptonVz = genleptons_vz->at(1);
	  ttreeAnalysis.gensecondLeptonVr = sqrt(genleptons_vx->at(1) * genleptons_vx->at(1) + genleptons_vy->at(1) * genleptons_vy->at(1));
	  ttreeAnalysis.gensecondLeptonVphi = atan2(genleptons_vy->at(1), genleptons_vx->at(1));
	}

	if(genleptons_pt->size()>1){
	  ttreeAnalysis.gendileptonCharge = genleptons_charge->at(0) + genleptons_charge->at(1);
	  ttreeAnalysis.gendileptonMass = gendilepton.M();
	  ttreeAnalysis.gendileptonEta = gendilepton.Eta();
	  ttreeAnalysis.gendileptonPhi = gendilepton.Phi();
	  ttreeAnalysis.gendileptonPt = gendilepton.Pt();
	  ttreeAnalysis.gendileptonRapidity = gendilepton.Rapidity();
	}

	ttreeAnalysis.genmissEt = genmissEt;
	ttreeAnalysis.genmissEt_phi = genmissEt_phi;

	ttreeAnalysis.nGenJets = jetsak4_pt->size();
	if(genjetsak4_pt->size()>0 && jetsak4_pt->size()>0){
	  if(MCMatch(genjetsak4_phi->at(0), genjetsak4_eta->at(0), genjetsak4_pt->at(0), jetsak4_phi->at(0), jetsak4_eta->at(0), jetsak4_pt->at(0))){
	    ttreeAnalysis.genleadingJetEnergy = genjetsak4_energy->at(0);
	    ttreeAnalysis.genleadingJetPx = genjetsak4_px->at(0);
	    ttreeAnalysis.genleadingJetPy = genjetsak4_py->at(0);
	    ttreeAnalysis.genleadingJetPz = genjetsak4_pz->at(0);
	    ttreeAnalysis.genleadingJetPt = genjetsak4_pt->at(0);
	    ttreeAnalysis.genleadingJetEta = genjetsak4_eta->at(0);
	    ttreeAnalysis.genleadingJetPhi = genjetsak4_phi->at(0);
	    ttreeAnalysis.genleadingJetVz = genjetsak4_vz->at(0);
	  }
	}
	if(genjetsak4_pt->size()>1 && jetsak4_pt->size()>1){
	  if(MCMatch(genjetsak4_phi->at(1), genjetsak4_eta->at(1), genjetsak4_pt->at(1), jetsak4_phi->at(1), jetsak4_eta->at(1), jetsak4_pt->at(1))){
	    ttreeAnalysis.gensecondJetEnergy = genjetsak4_energy->at(1);
	    ttreeAnalysis.gensecondJetPx = genjetsak4_px->at(1);
	    ttreeAnalysis.gensecondJetPy = genjetsak4_py->at(1);
	    ttreeAnalysis.gensecondJetPz = genjetsak4_pz->at(1);
	    ttreeAnalysis.gensecondJetPt = genjetsak4_pt->at(1);
	    ttreeAnalysis.gensecondJetEta = genjetsak4_eta->at(1);
	    ttreeAnalysis.gensecondJetPhi = genjetsak4_phi->at(1);
	    ttreeAnalysis.gensecondJetVz = genjetsak4_vz->at(1);
	  }
	}

	ttreeAnalysis.gendijetMass = gendijet.M();
	ttreeAnalysis.gendijetEta = gendijet.Eta();
	ttreeAnalysis.gendijetPhi = gendijet.Phi();
	ttreeAnalysis.gendijetPt = gendijet.Pt();
	ttreeAnalysis.gendijetRapidity = gendijet.Rapidity();
      }

      ttreeAnalysis.nLeptons = leptons_pt->size();

      if(leptons_pt->size()>0){
	ttreeAnalysis.leadingLeptonPDGId = leptons_pdgid->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonEnergy = leptons_energy->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonPx = leptons_px->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonPy = leptons_py->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonPz = leptons_pz->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonPt = leptons_pt->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonEta = leptons_eta->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonPhi = leptons_phi->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonCharge = leptons_charge->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonVx = leptons_vx->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonVy = leptons_vy->at(index_leptons.at(0));	
	ttreeAnalysis.leadingLeptonVz = leptons_vz->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonVr = sqrt(leptons_vx->at(index_leptons.at(0)) * leptons_vx->at(index_leptons.at(0)) + leptons_vy->at(index_leptons.at(0)) * leptons_vy->at(index_leptons.at(0)));
	ttreeAnalysis.leadingLeptonVphi = atan2(leptons_vy->at(index_leptons.at(0)), leptons_vx->at(index_leptons.at(0)));
	ttreeAnalysis.leadingLeptonPFIso = leptons_pfIso_->at(index_leptons.at(0)); 
	ttreeAnalysis.leadingLeptonTkIso = leptons_tkIso_->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonLooseId = leptons_looseId->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonMediumId = leptons_mediumId->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonTightId = leptons_tightId->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonPfIsoMedium = leptons_pfIsoMedium_->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonMiniIsoTight = leptons_miniIsoTight_->at(index_leptons.at(0));
	ttreeAnalysis.leadingLeptonPfIsoVeryTight = leptons_pfIsoVeryTight_->at(index_leptons.at(0));
      }

      if(leptons_pt->size()>1){
	ttreeAnalysis.secondLeptonPDGId = leptons_pdgid->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonEnergy = leptons_energy->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonPx = leptons_px->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonPy = leptons_py->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonPz = leptons_pz->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonPt = leptons_pt->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonEta = leptons_eta->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonPhi = leptons_phi->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonCharge = leptons_charge->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonVx = leptons_vx->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonVy = leptons_vy->at(index_leptons.at(1)); 
	ttreeAnalysis.secondLeptonVz = leptons_vz->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonVr = sqrt(leptons_vx->at(index_leptons.at(1)) * leptons_vx->at(index_leptons.at(1)) + leptons_vy->at(index_leptons.at(1)) * leptons_vy->at(index_leptons.at(1)));
	ttreeAnalysis.secondLeptonVphi = atan2(leptons_vy->at(index_leptons.at(1)), leptons_vx->at(index_leptons.at(1)));
	ttreeAnalysis.secondLeptonPFIso = leptons_pfIso_->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonTkIso = leptons_tkIso_->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonLooseId = leptons_looseId->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonMediumId = leptons_mediumId->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonTightId = leptons_tightId->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonPfIsoMedium = leptons_pfIsoMedium_->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonMiniIsoTight = leptons_miniIsoTight_->at(index_leptons.at(1));
	ttreeAnalysis.secondLeptonPfIsoVeryTight = leptons_pfIsoVeryTight_->at(index_leptons.at(1));
      }     

      ttreeAnalysis.nJets = jetsak4_pt->size();
      if(jetsak4_pt->size()>0){
	ttreeAnalysis.leadingJetEnergy = jetsak4_energy->at(0);
	ttreeAnalysis.leadingJetPx = jetsak4_px->at(0);
	ttreeAnalysis.leadingJetPy = jetsak4_py->at(0);
	ttreeAnalysis.leadingJetPz = jetsak4_pz->at(0);
	ttreeAnalysis.leadingJetPt = jetsak4_pt->at(0);
	ttreeAnalysis.leadingJetEta = jetsak4_eta->at(0);
	ttreeAnalysis.leadingJetPhi = jetsak4_phi->at(0);
	ttreeAnalysis.leadingJetVx = jetsak4_vx->at(0);
	ttreeAnalysis.leadingJetVy = jetsak4_vy->at(0);
	ttreeAnalysis.leadingJetVz = jetsak4_vz->at(0);
	ttreeAnalysis.leadingJetVr = sqrt(jetsak4_vx->at(0) * jetsak4_vx->at(0) + jetsak4_vy->at(0) * jetsak4_vy->at(0));
	ttreeAnalysis.leadingJetVphi = atan2(jetsak4_vy->at(0), jetsak4_vx->at(0));
	ttreeAnalysis.leadingJetBtag = jetsak4_bdis->at(0);
	ttreeAnalysis.leadingJetQGdis = jetsak4_qgdis->at(0);
	ttreeAnalysis.leadingJetNeutralEmFrac = jetsak4_neutralemfrac->at(0);
	ttreeAnalysis.leadingJetNeutralHadFrac = jetsak4_neutralhadfrac->at(0);
	ttreeAnalysis.leadingJetChargedEmFrac = jetsak4_chargedemfrac->at(0);
	ttreeAnalysis.leadingJetChargedHadFrac = jetsak4_chargedhadfrac->at(0);
	ttreeAnalysis.leadingJetMuonFrac = jetsak4_muonfrac->at(0);
	ttreeAnalysis.leadingJetNeutralMulti = jetsak4_neutralmulti->at(0);
	ttreeAnalysis.leadingJetChargedMulti = jetsak4_chargedmulti->at(0);
	ttreeAnalysis.leadingJetpuIdfdisc = jetsak4_puIdfdisc->at(0);
	ttreeAnalysis.leadingJetpuIdcbased = jetsak4_puIdcbased->at(0);
	ttreeAnalysis.leadingJetpuIdfid = jetsak4_puIdfid->at(0);
	ttreeAnalysis.leadingJetLooseId = jetsak4_looseId->at(0);
	ttreeAnalysis.leadingJetTightId = jetsak4_tightId->at(0);
	ttreeAnalysis.leadingJetLepVeto = jetsak4_lepVeto->at(0);
      }
      if(jetsak4_pt->size()>1){
	ttreeAnalysis.secondJetEnergy = jetsak4_energy->at(1);
	ttreeAnalysis.secondJetPx = jetsak4_px->at(1);
	ttreeAnalysis.secondJetPy = jetsak4_py->at(1);
	ttreeAnalysis.secondJetPz = jetsak4_pz->at(1);
	ttreeAnalysis.secondJetPt = jetsak4_pt->at(1);
	ttreeAnalysis.secondJetEta = jetsak4_eta->at(1);
	ttreeAnalysis.secondJetPhi = jetsak4_phi->at(1);
	ttreeAnalysis.secondJetVx = jetsak4_vx->at(1);
	ttreeAnalysis.secondJetVy = jetsak4_vy->at(1);
	ttreeAnalysis.secondJetVz = jetsak4_vz->at(1);
	ttreeAnalysis.secondJetVr = sqrt(jetsak4_vx->at(1) * jetsak4_vx->at(1) + jetsak4_vy->at(1) * jetsak4_vy->at(1));
	ttreeAnalysis.secondJetVphi = atan2(jetsak4_vy->at(1), jetsak4_vx->at(1));
	ttreeAnalysis.secondJetBtag = jetsak4_bdis->at(1);
	ttreeAnalysis.secondJetQGdis = jetsak4_qgdis->at(1);
	ttreeAnalysis.secondJetNeutralEmFrac = jetsak4_neutralemfrac->at(1);
	ttreeAnalysis.secondJetNeutralHadFrac = jetsak4_neutralhadfrac->at(1);
	ttreeAnalysis.secondJetChargedEmFrac = jetsak4_chargedemfrac->at(1);
	ttreeAnalysis.secondJetChargedHadFrac = jetsak4_chargedhadfrac->at(1);
	ttreeAnalysis.secondJetMuonFrac = jetsak4_muonfrac->at(1);
	ttreeAnalysis.secondJetNeutralMulti = jetsak4_neutralmulti->at(1);
	ttreeAnalysis.secondJetChargedMulti = jetsak4_chargedmulti->at(1);
	ttreeAnalysis.secondJetpuIdfdisc = jetsak4_puIdfdisc->at(1);
	ttreeAnalysis.secondJetpuIdcbased = jetsak4_puIdcbased->at(1);
	ttreeAnalysis.secondJetpuIdfid = jetsak4_puIdfid->at(1);
	ttreeAnalysis.secondJetLooseId = jetsak4_looseId->at(1);
	ttreeAnalysis.secondJetTightId = jetsak4_tightId->at(1);
	ttreeAnalysis.secondJetLepVeto = jetsak4_lepVeto->at(1);
      }

      ttreeAnalysis.dijetMass = dijet.M();
      ttreeAnalysis.dijetEta = dijet.Eta();
      ttreeAnalysis.dijetPhi = dijet.Phi();
      ttreeAnalysis.dijetPt = dijet.Pt();
      ttreeAnalysis.dijetRapidity = dijet.Rapidity();

      ttreeAnalysis.jetsystemMass = jetsystem.M();
      ttreeAnalysis.jetsystemEta = jetsystem.Eta();
      ttreeAnalysis.jetsystemPhi = jetsystem.Phi();
      ttreeAnalysis.jetsystemPt = jetsystem.Pt();
      ttreeAnalysis.jetsystemRapidity = jetsystem.Rapidity();

      ttreeAnalysis.jetcandidatesystemMass = jetcandidatesystem.M();
      ttreeAnalysis.jetcandidatesystemEta = jetcandidatesystem.Eta();
      ttreeAnalysis.jetcandidatesystemPhi = jetcandidatesystem.Phi();
      ttreeAnalysis.jetcandidatesystemPt = jetcandidatesystem.Pt();
      ttreeAnalysis.jetcandidatesystemRapidity = jetcandidatesystem.Rapidity();

      if(leptons_pt->size()>0){
	ttreeAnalysis.leptonmetCharge = leptons_charge->at(index_leptons.at(0));
	ttreeAnalysis.leptonmetMassT = leptonmet.M();
	ttreeAnalysis.leptonmetEta = leptonmet.Eta();
	ttreeAnalysis.leptonmetPhi = leptonmet.Phi();
	ttreeAnalysis.leptonmetPt = leptonmet.Pt();
	ttreeAnalysis.leptonmetRapidity = leptonmet.Rapidity();
      }

      if(leptons_pt->size()>1){
	ttreeAnalysis.dileptonCharge = leptons_charge->at(index_leptons.at(0)) + leptons_charge->at(index_leptons.at(1));
	ttreeAnalysis.dileptonMass = dilepton.M();
	ttreeAnalysis.dileptonEta = dilepton.Eta();
	ttreeAnalysis.dileptonPhi = dilepton.Phi();
	ttreeAnalysis.dileptonPt = dilepton.Pt();
	ttreeAnalysis.dileptonRapidity = dilepton.Rapidity();
      }

      ttreeAnalysis.leptonsystemMass = leptonsystem.M();
      ttreeAnalysis.leptonsystemEta = leptonsystem.Eta();
      ttreeAnalysis.leptonsystemPhi = leptonsystem.Phi();
      ttreeAnalysis.leptonsystemPt = leptonsystem.Pt();
      ttreeAnalysis.leptonsystemRapidity = leptonsystem.Rapidity();

      ttreeAnalysis.missingMassDijetRP210 = X_dijetRP210.M();
      ttreeAnalysis.missingEtaDijetRP210 = X_dijetRP210.Eta();
      ttreeAnalysis.missingPhiDijetRP210 = X_dijetRP210.Phi();
      ttreeAnalysis.missingPtDijetRP210 = X_dijetRP210.Pt();
      ttreeAnalysis.missingRapidityDijetRP210 = X_dijetRP210.Rapidity();

      ttreeAnalysis.missingMassDijetRP220 = X_dijetRP220.M();
      ttreeAnalysis.missingEtaDijetRP220 = X_dijetRP220.Eta();
      ttreeAnalysis.missingPhiDijetRP220 = X_dijetRP220.Phi();
      ttreeAnalysis.missingPtDijetRP220 = X_dijetRP220.Pt();
      ttreeAnalysis.missingRapidityDijetRP220 = X_dijetRP220.Rapidity();

      ttreeAnalysis.missingMassDijetMulti = X_dijetMulti.M();
      ttreeAnalysis.missingEtaDijetMulti = X_dijetMulti.Eta();
      ttreeAnalysis.missingPhiDijetMulti = X_dijetMulti.Phi();
      ttreeAnalysis.missingPtDijetMulti = X_dijetMulti.Pt();
      ttreeAnalysis.missingRapidityDijetMulti = X_dijetMulti.Rapidity();

      ttreeAnalysis.missingMassDileptonRP210 = X_dileptonRP210.M();
      ttreeAnalysis.missingEtaDileptonRP210 = X_dileptonRP210.Eta();
      ttreeAnalysis.missingPhiDileptonRP210 = X_dileptonRP210.Phi();
      ttreeAnalysis.missingPtDileptonRP210 = X_dileptonRP210.Pt();
      ttreeAnalysis.missingRapidityDileptonRP210 = X_dileptonRP210.Rapidity();

      ttreeAnalysis.missingMassDileptonRP220 = X_dileptonRP220.M();
      ttreeAnalysis.missingEtaDileptonRP220 = X_dileptonRP220.Eta();
      ttreeAnalysis.missingPhiDileptonRP220 = X_dileptonRP220.Phi();
      ttreeAnalysis.missingPtDileptonRP220 = X_dileptonRP220.Pt();
      ttreeAnalysis.missingRapidityDileptonRP220 = X_dileptonRP220.Rapidity();

      ttreeAnalysis.missingMassDileptonMulti = X_dileptonMulti.M();
      ttreeAnalysis.missingEtaDileptonMulti = X_dileptonMulti.Eta();
      ttreeAnalysis.missingPhiDileptonMulti = X_dileptonMulti.Phi();
      ttreeAnalysis.missingPtDileptonMulti = X_dileptonMulti.Pt();
      ttreeAnalysis.missingRapidityDileptonMulti = X_dileptonMulti.Rapidity();

      ttreeAnalysis.missingMassDileptonJet1RP220 = X_dileptonjet1RP220.M();
      ttreeAnalysis.missingEtaDileptonJet1RP220 = X_dileptonjet1RP220.Eta();
      ttreeAnalysis.missingPhiDileptonJet1RP220 = X_dileptonjet1RP220.Phi();
      ttreeAnalysis.missingPtDileptonJet1RP220 = X_dileptonjet1RP220.Pt();
      ttreeAnalysis.missingRapidityDileptonJet1RP220 = X_dileptonjet1RP220.Rapidity();

      ttreeAnalysis.missingMassDileptonJet1Multi = X_dileptonjet1Multi.M();
      ttreeAnalysis.missingEtaDileptonJet1Multi = X_dileptonjet1Multi.Eta();
      ttreeAnalysis.missingPhiDileptonJet1Multi = X_dileptonjet1Multi.Phi();
      ttreeAnalysis.missingPtDileptonJet1Multi = X_dileptonjet1Multi.Pt();
      ttreeAnalysis.missingRapidityDileptonJet1Multi = X_dileptonjet1Multi.Rapidity();

      ttreeAnalysis.missingMassDileptonDijetRP210 = X_dileptondijetRP210.M();
      ttreeAnalysis.missingEtaDileptonDijetRP210 = X_dileptondijetRP210.Eta();
      ttreeAnalysis.missingPhiDileptonDijetRP210 = X_dileptondijetRP210.Phi();
      ttreeAnalysis.missingPtDileptonDijetRP210 = X_dileptondijetRP210.Pt();
      ttreeAnalysis.missingRapidityDileptonDijetRP210 = X_dileptondijetRP210.Rapidity();

      ttreeAnalysis.missingMassDileptonDijetRP220 = X_dileptondijetRP220.M();
      ttreeAnalysis.missingEtaDileptonDijetRP220 = X_dileptondijetRP220.Eta();
      ttreeAnalysis.missingPhiDileptonDijetRP220 = X_dileptondijetRP220.Phi();
      ttreeAnalysis.missingPtDileptonDijetRP220 = X_dileptondijetRP220.Pt();
      ttreeAnalysis.missingRapidityDileptonDijetRP220 = X_dileptondijetRP220.Rapidity();

      ttreeAnalysis.missingMassDileptonDijetMulti = X_dileptondijetMulti.M();
      ttreeAnalysis.missingEtaDileptonDijetMulti = X_dileptondijetMulti.Eta();
      ttreeAnalysis.missingPhiDileptonDijetMulti = X_dileptondijetMulti.Phi();
      ttreeAnalysis.missingPtDileptonDijetMulti = X_dileptondijetMulti.Pt();
      ttreeAnalysis.missingRapidityDileptonDijetMulti = X_dileptondijetMulti.Rapidity();

      ttreeAnalysis.missEt = missEt;
      ttreeAnalysis.missEt_phi = missEt_phi;

      ttreeAnalysis.nElectrons = nElectrons;
      ttreeAnalysis.nMuons = nMuons;
      ttreeAnalysis.nElectrons_looseId = nElectrons_looseId;
      ttreeAnalysis.nElectrons_mediumId = nElectrons_mediumId;
      ttreeAnalysis.nElectrons_tightId = nElectrons_tightId;
      ttreeAnalysis.nMuons_looseId = nMuons_looseId;
      ttreeAnalysis.nMuons_mediumId = nMuons_mediumId;
      ttreeAnalysis.nMuons_tightId = nMuons_tightId; 
      ttreeAnalysis.nChargedPFMultiPV_Loose = nChargedPFMultiPV_Loose;
      ttreeAnalysis.nChargedPFMultiPV_Tight = nChargedPFMultiPV_Tight;
      ttreeAnalysis.nChargedPFMultiPV_UsedInFit = nChargedPFMultiPV_UsedInFit;
      ttreeAnalysis.nChargedPFMultiPV_Tight_Fit = nChargedPFMultiPV_Tight_Fit;
      ttreeAnalysis.SumChargedPFMultiPV_pt_Loose = SumChargedPFMultiPV_pt_Loose; 
      ttreeAnalysis.SumChargedPFMultiPV_pt_Tight = SumChargedPFMultiPV_pt_Tight;
      ttreeAnalysis.SumChargedPFMultiPV_pt_UsedInFit = SumChargedPFMultiPV_pt_UsedInFit;
      ttreeAnalysis.SumChargedPFMultiPV_pt_Tight_Fit = SumChargedPFMultiPV_pt_Tight_Fit;

      double dphi = fabs(dilepton.Phi()-jetsystem.Phi());
      dphi = (dphi<pi) ? dphi : 2.*pi-dphi;
      ttreeAnalysis.dphi = dphi;

      if(ttreeAnalysis.isprotonMulti &&  ttreeAnalysis.timeunc_multiArm45>0 && ttreeAnalysis.timeunc_multiArm56>0){
	ttreeAnalysis.t45andt56 = c_light*ns_to_s_*m_to_cm_*(ttreeAnalysis.time_multiArm56 + ttreeAnalysis.time_multiArm45);
	ttreeAnalysis.vzpps = (c_light/2)*ns_to_s_*m_to_cm_*(ttreeAnalysis.time_multiArm56-ttreeAnalysis.time_multiArm45);
	for (auto evtVz : *vertex_z){
	  ttreeAnalysis.allVertexZ.push_back(evtVz);
	  ttreeAnalysis.MinimumDistance.push_back(std::make_pair(fabs(ttreeAnalysis.vzpps-evtVz), evtVz));
	}
      }else{
	ttreeAnalysis.t45andt56 = -999;
	ttreeAnalysis.vzpps = -999;
      }

      // Finding the minimum distance between CMS vtx and PPS vtx
      std::sort(ttreeAnalysis.MinimumDistance.begin(), ttreeAnalysis.MinimumDistance.end());

      if(ttreeAnalysis.MinimumDistance.size()>0){
	ttreeAnalysis.distanceVertexZPPSandCMS = ttreeAnalysis.MinimumDistance[0].first;
	ttreeAnalysis.selectedVertexZ = ttreeAnalysis.MinimumDistance[0].second;
      }else{
	ttreeAnalysis.distanceVertexZPPSandCMS = -999.;
	ttreeAnalysis.selectedVertexZ = -999.;
      }

      ttreeAnalysis.Fill();

      /*
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
      */

    } // proton bracket

  } // end of the event loop

  std::cout << "\n\n<END> good luck!\n" << std::endl;

  ttreeAnalysis.Storing();
  //fileout->cd();
  //test->Write();
  //tout->Write();
  //fileout->Close();

  /*
     if(createProtonFile) {
     file_protons_reco.close();
     }
     if(createEventFile){
     file_eventlist.close();
     file_eventlist_kinematics.close();
     }
     */

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

