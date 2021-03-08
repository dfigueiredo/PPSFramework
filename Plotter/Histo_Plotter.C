#include "TFile.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TTree.h"

void Histo_Plotter(){
  
  //DY and Signal Ntuples
  TFile   *f_DY = new TFile("/eos/user/l/lpagliai/Background_Studies/DY_ntuples/Muons/Analyzer/DY_muon.root");
  TTree   *DY_ntuple = (TTree*)f_DY->Get("Events");
  TFile   *f_SD = new TFile("/eos/user/l/lpagliai/Background_Studies/SD_ntuples/Muons/Analyzer/SD_muon_new.root");
  TTree   *SD_ntuple = (TTree*)f_SD->Get("Events");
  TFile   *f_Sign = new TFile("/eos/user/l/lpagliai/Background_Studies/Toy_ntuples/Muons/Analyzer/TOY_muon.root");
  TTree   *Sign_ntuple = (TTree*)f_Sign->Get("Events");
 
  //Histograms
  TH1D *h_delta_z1_DY, *h_delta_z1_Sign, *h_delta_z1_SD;
  //TH1D *h_delta_z2_DY, *h_delta_z2_Sign, *h_delta_z2_SD;
  TH1D *h_pt_dil_DY, *h_pt_dil_Sign, *h_pt_dil_SD; 
  TH1D *h_eta_dil_DY, *h_eta_dil_Sign, *h_eta_dil_SD;
  /*TH1D *h_pt_lep1_DY, *h_pt_lep1_Sign, *h_pt_lep1_SD;
  TH1D *h_eta_lep1_DY, *h_eta_lep1_Sign, *h_eta_lep1_SD; 
  TH1D *h_pt_lep2_DY, *h_pt_lep2_Sign, *h_pt_lep2_SD;
  TH1D *h_eta_lep2_DY, *h_eta_lep2_Sign, *h_eta_lep2_SD;*/
  TH1D *h_pt_jet1_DY, *h_pt_jet1_Sign, *h_pt_jet1_SD;
  TH1D *h_eta_jet1_DY, *h_eta_jet1_Sign, *h_eta_jet1_SD; 
  /*TH1D *h_pt_jet2_DY, *h_pt_jet2_Sign, *h_pt_jet2_SD;
  TH1D *h_eta_jet2_DY, *h_eta_jet2_Sign, *h_eta_jet2_SD;*/       
  TH1D *h_met_DY, *h_met_Sign, *h_met_SD;
  //TH1D *h_rmpf_DY, *h_rmpf_Sign, *h_rmpf_SD;
  TH1D *h_miss_ma_DY, *h_miss_ma_Sign, *h_miss_ma_SD;
  /*TH1D *h_miss_ma_120_pre_DY, *h_miss_ma_120_pre_Sign, *h_miss_ma_120_pre_SD;
  TH1D *h_miss_ma_120_post_DY, *h_miss_ma_120_post_Sign, *h_miss_ma_120_post_SD;*/
  TH1D *h_miss_ma_120_DY, *h_miss_ma_120_Sign, *h_miss_ma_120_SD;
  TH1D *h_miss_ma_130_DY, *h_miss_ma_130_Sign, *h_miss_ma_130_SD;
  TH1D *h_miss_ma_140_DY, *h_miss_ma_140_Sign, *h_miss_ma_140_SD;
  TH1D *h_miss_ma_150_DY, *h_miss_ma_150_Sign, *h_miss_ma_150_SD;
  TH1D *h_jzb1_DY, *h_jzb1_Sign, *h_jzb1_SD;
  TH1D *h_jzb2_DY, *h_jzb2_Sign, *h_jzb2_SD;
  TH1D *h_phijetZ_DY, *h_phijetZ_Sign, *h_phijetZ_SD;
  TH1D *h_phimetma_DY, *h_phimetma_Sign, *h_phimetma_SD;
  TH1D *h_nchPF_DY, *h_nchPF_Sign, *h_nchPF_SD;
  TH1D *h_sumchPF_DY, *h_sumchPF_Sign, *h_sumchPF_SD;
  TH1D *h_qgdis_DY, *h_qgdis_Sign, *h_qgdis_SD;
  TH1D *h_njetslos_DY, *h_njetslos_Sign, *h_njetslos_SD;
  //TH1D *h_njetsmed_DY, *h_njetsmed_Sign, *h_njetsmed_SD;
  //TH1D *h_njetstig_DY, *h_njetstig_Sign, *h_njetstig_SD;
  TH1D *h_xi45_multi_DY, *h_xi45_multi_Sign, *h_xi45_multi_SD;
  TH1D *h_xi56_multi_DY, *h_xi56_multi_Sign, *h_xi56_multi_SD;
  TH1D *h_miss_ma_fake_SD;
  /*TH1D *h_mima_outZ_DY, *h_mima_outZ_SD, *h_mima_outZ_Sign;
  TH1D *h_mima_inZ_DY, *h_mima_inZ_SD, *h_mima_inZ_Sign;
  TH1D *h_mima_outjetpv_DY, *h_mima_outjetpv_SD, *h_mima_outjetpv_Sign;
  TH1D *h_mima_injetpv_DY, *h_mima_injetpv_SD, *h_mima_injetpv_Sign;*/
  TH1D *h_diffma_DY, *h_diffma_Sign, *h_diffma_SD;
  TH1D *h_ptvertZ_DY, *h_ptvertZ_Sign, *h_ptvertZ_SD;
  TH1D *h_ratio_metpt_DY, *h_ratio_metpt_Sign, *h_ratio_metpt_SD;

  //Cut to be applied
  //isprotonMulti is FUNDAMENTAL
  /*TCut cut_Nicola = "MCweight*PUweight*(isprotonMulti==1&&missingMassDileptonMulti>0&&dileptonCharge==0&&dileptonMass>85&&dileptonMass<97&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>20&&fabs(leadingJetVz-pvVertexZ)>0.1&&fabs(missEt_phi-missingPhiDileptonMulti)<2.5&&dileptonEta>-1.9&&dileptonEta<1.6&&leadingJetPt>0&&secondJetPt>0)";*/
  //TCut cut_Z_DY = "PUweight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))";
  TCut cut_Z_excl = "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))"; 
  //TCut cut_Z_Sign = "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))";

  TCut unblind_cut = "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex<20&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40)";

  /*TCut new_Z_DY = "MCweight*PUweight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>85&&dileptonMass<97&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex<50&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&leadingJetPt<30&&missEt>30)";
  TCut new_Z_SD = "MCweight*PUweight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>85&&dileptonMass<97&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex<50&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&leadingJetPt<30&&missEt>30&&is2PUproton==0)";
  TCut new_Z_Sign = "MCweight*PUweight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>85&&dileptonMass<97&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex<50&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&leadingJetPt<30&&missEt>30&&is45PUproton==0&&is56PUproton==0&&is2PUproton==0)";
*/
  //Drawing Distance Vertex Leading Jet Primary Vertex
  //Histograms must be defined in a range AS WIDER AS POSSIBLE (Normalization fails)
  /*Sign_ntuple->Draw("fabs(leadingJetVz-pvVertexZ)>>h_delta_z1_Sign(2000,0.,20.)",unblind_cut,"hist");
  h_delta_z1_Sign = (TH1D*)gPad->GetPrimitive("h_delta_z1_Sign");
  DY_ntuple->Draw("fabs(leadingJetVz-pvVertexZ)>>h_delta_z1_DY(2000,0.,20.)",unblind_cut,"hist");
  h_delta_z1_DY = (TH1D*)gPad->GetPrimitive("h_delta_z1_DY");
  SD_ntuple->Draw("fabs(leadingJetVz-pvVertexZ)>>h_delta_z1_SD(2000,0.,20.)",unblind_cut,"hist");
  h_delta_z1_SD = (TH1D*)gPad->GetPrimitive("h_delta_z1_SD"); */
  
  //Drawing Distance Vertex Second Jet Primary Vertex
  /*Sign_ntuple->Draw("fabs(secondJetVz-pvVertexZ)>>h_delta_z2_Sign(2000,0.,20.)", cut_Z_Sign, "hist");
  h_delta_z2_Sign = (TH1D*)gPad->GetPrimitive("h_delta_z2_Sign");
  DY_ntuple->Draw("fabs(secondJetVz-pvVertexZ)>>h_delta_z2_DY(2000,0.,20.)", cut_Z_DY, "hist");
  h_delta_z2_DY = (TH1D*)gPad->GetPrimitive("h_delta_z2_DY");
  SD_ntuple->Draw("fabs(secondJetVz-pvVertexZ)>>h_delta_z2_SD(2000,0.,20.)", cut_Z_SD, "hist");
  h_delta_z2_SD = (TH1D*)gPad->GetPrimitive("h_delta_z2_SD");
*/
  //Drawing Dilepton Pt
  Sign_ntuple->Draw("dileptonPt>>h_pt_dil_Sign(300,0.,300)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_pt_dil_Sign = (TH1D*)gPad->GetPrimitive("h_pt_dil_Sign");
  DY_ntuple->Draw("dileptonPt>>h_pt_dil_DY(300,0.,300)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_pt_dil_DY = (TH1D*)gPad->GetPrimitive("h_pt_dil_DY");
  SD_ntuple->Draw("dileptonPt>>h_pt_dil_SD(300,0.,300)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_pt_dil_SD = (TH1D*)gPad->GetPrimitive("h_pt_dil_SD");
  
  //Drawing Dilepton Pseudorapidity 
  Sign_ntuple->Draw("dileptonEta>>h_eta_dil_Sign(200,-10.,10.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_eta_dil_Sign = (TH1D*)gPad->GetPrimitive("h_eta_dil_Sign");
  DY_ntuple->Draw("dileptonEta>>h_eta_dil_DY(200,-10.,10.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_eta_dil_DY = (TH1D*)gPad->GetPrimitive("h_eta_dil_DY");
  SD_ntuple->Draw("dileptonEta>>h_eta_dil_SD(200,-10.,10.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_eta_dil_SD = (TH1D*)gPad->GetPrimitive("h_eta_dil_SD");

  //Drawing Leading Lepton Pt
  /*Sign_ntuple->Draw("leadingLeptonPt>>h_pt_lep1_Sign(300,0.,300)", cut_Z_Sign, "hist");
  h_pt_lep1_Sign = (TH1D*)gPad->GetPrimitive("h_pt_lep1_Sign");
  DY_ntuple->Draw("leadingLeptonPt>>h_pt_lep1_DY(300,0.,300)", cut_Z_DY, "hist");
  h_pt_lep1_DY = (TH1D*)gPad->GetPrimitive("h_pt_lep1_DY");
  SD_ntuple->Draw("leadingLeptonPt>>h_pt_lep1_SD(300,0.,300)", cut_Z_SD, "hist");
  h_pt_lep1_SD = (TH1D*)gPad->GetPrimitive("h_pt_lep1_SD");

  //Drawing Leading Lepton Pseudorapidity 
  Sign_ntuple->Draw("leadingLeptonEta>>h_eta_lep1_Sign(200,-10.,10.)", cut_Z_Sign, "hist");
  h_eta_lep1_Sign = (TH1D*)gPad->GetPrimitive("h_eta_lep1_Sign");
  DY_ntuple->Draw("leadingLeptonEta>>h_eta_lep1_DY(200,-10.,10.)", cut_Z_DY, "hist");
  h_eta_lep1_DY = (TH1D*)gPad->GetPrimitive("h_eta_lep1_DY");
  SD_ntuple->Draw("leadingLeptonEta>>h_eta_lep1_SD(200,-10.,10.)", cut_Z_SD, "hist");    
  h_eta_lep1_SD = (TH1D*)gPad->GetPrimitive("h_eta_lep1_SD");
 
  //Drawing Second Lepton Pt
  Sign_ntuple->Draw("secondLeptonPt>>h_pt_lep2_Sign(300,0.,300)", cut_Z_Sign, "hist");
  h_pt_lep2_Sign = (TH1D*)gPad->GetPrimitive("h_pt_lep2_Sign");
  DY_ntuple->Draw("secondLeptonPt>>h_pt_lep2_DY(300,0.,300)", cut_Z_DY, "hist");
  h_pt_lep2_DY = (TH1D*)gPad->GetPrimitive("h_pt_lep2_DY");
  SD_ntuple->Draw("secondLeptonPt>>h_pt_lep2_SD(300,0.,300)", cut_Z_SD, "hist");
  h_pt_lep2_SD = (TH1D*)gPad->GetPrimitive("h_pt_lep2_SD");

  //Drawing Second Lepton Pseudorapidity 
  Sign_ntuple->Draw("secondLeptonEta>>h_eta_lep2_Sign(200,-10.,10.)", cut_Z_Sign, "hist");
  h_eta_lep2_Sign = (TH1D*)gPad->GetPrimitive("h_eta_lep2_Sign");
  DY_ntuple->Draw("secondLeptonEta>>h_eta_lep2_DY(200,-10.,10.)", cut_Z_DY, "hist");
  h_eta_lep2_DY = (TH1D*)gPad->GetPrimitive("h_eta_lep2_DY");
  SD_ntuple->Draw("secondLeptonEta>>h_eta_lep2_SD(200,-10.,10.)", cut_Z_SD, "hist");
  h_eta_lep2_SD = (TH1D*)gPad->GetPrimitive("h_eta_lep2_SD");*/

  //Drawing Leading Jet Pt
  Sign_ntuple->Draw("leadingJetPt>>h_pt_jet1_Sign(300,0.,300)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_pt_jet1_Sign = (TH1D*)gPad->GetPrimitive("h_pt_jet1_Sign");
  DY_ntuple->Draw("leadingJetPt>>h_pt_jet1_DY(300,0.,300)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_pt_jet1_DY = (TH1D*)gPad->GetPrimitive("h_pt_jet1_DY");
  SD_ntuple->Draw("leadingJetPt>>h_pt_jet1_SD(300,0.,300)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_pt_jet1_SD = (TH1D*)gPad->GetPrimitive("h_pt_jet1_SD");

  //Drawing Leading Jet Pseudorapidity 
  Sign_ntuple->Draw("leadingJetEta>>h_eta_jet1_Sign(200,-10.,10.)", cut_Z_excl, "hist");
  h_eta_jet1_Sign = (TH1D*)gPad->GetPrimitive("h_eta_jet1_Sign");
  DY_ntuple->Draw("leadingJetEta>>h_eta_jet1_DY(200,-10.,10.)", cut_Z_excl, "hist");
  h_eta_jet1_DY = (TH1D*)gPad->GetPrimitive("h_eta_jet1_DY");
  SD_ntuple->Draw("leadingJetEta>>h_eta_jet1_SD(200,-10.,10.)", cut_Z_excl, "hist");
  h_eta_jet1_SD = (TH1D*)gPad->GetPrimitive("h_eta_jet1_SD");

  //Drawing Second Jet Pt
/*  Sign_ntuple->Draw("secondJetPt>>h_pt_jet2_Sign(300,0.,300)", cut_Z_Sign, "hist");
  h_pt_jet2_Sign = (TH1D*)gPad->GetPrimitive("h_pt_jet2_Sign");
  DY_ntuple->Draw("secondJetPt>>h_pt_jet2_DY(300,0.,300)", cut_Z_DY, "hist");
  h_pt_jet2_DY = (TH1D*)gPad->GetPrimitive("h_pt_jet2_DY");
  SD_ntuple->Draw("secondJetPt>>h_pt_jet2_SD(300,0.,300)", cut_Z_SD, "hist");
  h_pt_jet2_SD = (TH1D*)gPad->GetPrimitive("h_pt_jet2_SD");

  //Drawing Second Jet Pseudorapidity 
  Sign_ntuple->Draw("secondJetEta>>h_eta_jet2_Sign(200,-10.,10.)", cut_Z_Sign, "hist");
  h_eta_jet2_Sign = (TH1D*)gPad->GetPrimitive("h_eta_jet2_Sign");
  DY_ntuple->Draw("secondJetEta>>h_eta_jet2_DY(200,-10.,10.)", cut_Z_DY, "hist");
  h_eta_jet2_DY = (TH1D*)gPad->GetPrimitive("h_eta_jet2_DY");
  SD_ntuple->Draw("secondJetEta>>h_eta_jet2_SD(200,-10.,10.)", cut_Z_SD, "hist");
  h_eta_jet2_SD = (TH1D*)gPad->GetPrimitive("h_eta_jet2_SD");
*/
  //Drawing MET
  Sign_ntuple->Draw("missEt>>h_met_Sign(300, 0.,300.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_met_Sign = (TH1D*)gPad->GetPrimitive("h_met_Sign");
  DY_ntuple->Draw("missEt>>h_met_DY(300, 0.,300.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_met_DY = (TH1D*)gPad->GetPrimitive("h_met_DY");
  SD_ntuple->Draw("missEt>>h_met_SD(300, 0.,300.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_met_SD = (TH1D*)gPad->GetPrimitive("h_met_SD");

  //Drawing RMPF
/*  Sign_ntuple->Draw("Rmpf>>h_rmpf_Sign(200,-10.,10.)", cut_Z_Sign, "hist");
  h_rmpf_Sign = (TH1D*)gPad->GetPrimitive("h_rmpf_Sign");
  DY_ntuple->Draw("Rmpf>>h_rmpf_DY(200,-10.,10.)", cut_Z_DY, "hist");
  h_rmpf_DY = (TH1D*)gPad->GetPrimitive("h_rmpf_DY");
  SD_ntuple->Draw("Rmpf>>h_rmpf_SD(200,-10.,10.)", cut_Z_SD, "hist");
  h_rmpf_SD = (TH1D*)gPad->GetPrimitive("h_rmpf_SD");
*/
  //Drawing Missing Mass Multi
  Sign_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_Sign(2000,0.,2000.)", cut_Z_excl, "hist");
  h_miss_ma_Sign = (TH1D*)gPad->GetPrimitive("h_miss_ma_Sign");
  /*Sign_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_120_pre_Sign(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex<20&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&dileptonPt>40&&xangle==120&&era<4)", "hist");
  h_miss_ma_120_pre_Sign = (TH1D*)gPad->GetPrimitive("h_miss_ma_120_pre_Sign");
  Sign_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_120_post_Sign(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex<20&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&dileptonPt>40&&xangle==120&&era>3)", "hist");
  h_miss_ma_120_post_Sign = (TH1D*)gPad->GetPrimitive("h_miss_ma_120_post_Sign");*/
  DY_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_DY(2000,0.,2000.)", cut_Z_excl, "hist");
  h_miss_ma_DY = (TH1D*)gPad->GetPrimitive("h_miss_ma_DY");
  SD_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_SD(2000,0.,2000.)", cut_Z_excl, "hist");
  h_miss_ma_SD = (TH1D*)gPad->GetPrimitive("h_miss_ma_SD");   

  //Drawig Missing Mass Multi x-angle 120
  Sign_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_120_Sign(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==120)", "hist");
  h_miss_ma_120_Sign = (TH1D*)gPad->GetPrimitive("h_miss_ma_120_Sign");
  DY_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_120_DY(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==120)", "hist");
  h_miss_ma_120_DY = (TH1D*)gPad->GetPrimitive("h_miss_ma_120_DY");
  SD_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_120_SD(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==120)", "hist");
  h_miss_ma_120_SD = (TH1D*)gPad->GetPrimitive("h_miss_ma_120_SD");

  //Drawig Missing Mass Multi x-angle 130
  Sign_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_130_Sign(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==130)", "hist");
  h_miss_ma_130_Sign = (TH1D*)gPad->GetPrimitive("h_miss_ma_130_Sign");
  DY_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_130_DY(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==130)", "hist");
  h_miss_ma_130_DY = (TH1D*)gPad->GetPrimitive("h_miss_ma_130_DY");
  SD_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_130_SD(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==130)", "hist");
  h_miss_ma_130_SD = (TH1D*)gPad->GetPrimitive("h_miss_ma_130_SD");

  //Drawing Missing Mass Multi x-angle 140
  Sign_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_140_Sign(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==140)", "hist");
  h_miss_ma_140_Sign = (TH1D*)gPad->GetPrimitive("h_miss_ma_140_Sign");
  DY_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_140_DY(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==140)", "hist");
  h_miss_ma_140_DY = (TH1D*)gPad->GetPrimitive("h_miss_ma_140_DY");
  SD_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_140_SD(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==140)", "hist");
  h_miss_ma_140_SD = (TH1D*)gPad->GetPrimitive("h_miss_ma_140_SD");
  
  //Drawing Missing Mass Multi x-angle 150
  Sign_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_150_Sign(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==150)", "hist");
  h_miss_ma_150_Sign = (TH1D*)gPad->GetPrimitive("h_miss_ma_150_Sign");
  DY_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_150_DY(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==150)", "hist");
  h_miss_ma_150_DY = (TH1D*)gPad->GetPrimitive("h_miss_ma_150_DY");
  SD_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_150_SD(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&xangle==150)", "hist");
  h_miss_ma_150_SD = (TH1D*)gPad->GetPrimitive("h_miss_ma_150_SD");

  //Drawing Diffractive Mass Multi
  Sign_ntuple->Draw("diffMassMulti>>h_diffma_Sign(2000,0.,2000.)", cut_Z_excl, "hist");
  h_diffma_Sign = (TH1D*)gPad->GetPrimitive("h_diffma_Sign");
  DY_ntuple->Draw("diffMassMulti>>h_diffma_DY(2000,0.,2000.)", cut_Z_excl, "hist");
  h_diffma_DY = (TH1D*)gPad->GetPrimitive("h_diffma_DY");
  SD_ntuple->Draw("diffMassMulti>>h_diffma_SD(2000,0.,2000.)", cut_Z_excl, "hist");
  h_diffma_SD = (TH1D*)gPad->GetPrimitive("h_diffma_SD");

  //Drawing fake in SD
  SD_ntuple->Draw("missingMassDileptonMulti>>h_miss_ma_fake_SD(2000,0.,2000.)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4)&&is45PUproton==0&&is56PUproton==0&&is2PUproton==0)", "hist");
  h_miss_ma_fake_SD = (TH1D*)gPad->GetPrimitive("h_miss_ma_fake_SD");
 
  //Drawing Missing Mass dentro la Z
  /*Sign_ntuple->Draw("missingMassDileptonMulti>>h_mima_inZ_Sign(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&dileptonEta>-1.9&&dileptonEta<1.6&&is2PUproton==0&&is45PUproton==0&&is56PUproton==0", "hist");
  h_mima_inZ_Sign = (TH1D*)gPad->GetPrimitive("h_mima_inZ_Sign");
  DY_ntuple->Draw("missingMassDileptonMulti>>h_mima_inZ_DY(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&dileptonEta>-1.9&&dileptonEta<1.6)", "hist");
  h_mima_inZ_DY = (TH1D*)gPad->GetPrimitive("h_mima_inZ_DY");
  SD_ntuple->Draw("missingMassDileptonMulti>>h_mima_inZ_SD(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&dileptonEta>-1.9&&dileptonEta<1.6&&is2PUproton==0)", "hist");
  h_mima_inZ_SD = (TH1D*)gPad->GetPrimitive("h_mima_inZ_SD");

  //Drawing Missing Mass Eta Z < -1.9 || Eta Z >1.6
  Sign_ntuple->Draw("missingMassDileptonMulti>>h_mima_outZ_Sign(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&(dileptonEta<-1.9||dileptonEta>1.6)&&is2PUproton==0&&is45PUproton==0&&is56PUproton==0)", "hist");
  h_mima_outZ_Sign = (TH1D*)gPad->GetPrimitive("h_mima_outZ_Sign");
  DY_ntuple->Draw("missingMassDileptonMulti>>h_mima_outZ_DY(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&(dileptonEta<-1.9||dileptonEta>1.6))", "hist");
  h_mima_outZ_DY = (TH1D*)gPad->GetPrimitive("h_mima_outZ_DY");
  SD_ntuple->Draw("missingMassDileptonMulti>>h_mima_outZ_SD(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&(dileptonEta<-1.9||dileptonEta>1.6)&&is2PUproton==0)", "hist");
  h_mima_outZ_SD = (TH1D*)gPad->GetPrimitive("h_mima_outZ_SD");

  //Drawing Missing Mass Distance Leading Jet PV more than 1 mm
  Sign_ntuple->Draw("missingMassDileptonMulti>>h_mima_outjetpv_Sign(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&fabs(leadingJetVz-pvVertexZ)>0.1&&is2PUproton==0&&is45PUproton==0&&is56PUproton==0)", "hist");
  h_mima_outjetpv_Sign = (TH1D*)gPad->GetPrimitive("h_mima_outjetpv_Sign");
  DY_ntuple->Draw("missingMassDileptonMulti>>h_mima_outjetpv_DY(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&fabs(leadingJetVz-pvVertexZ)>0.1)", "hist");
  h_mima_outjetpv_DY = (TH1D*)gPad->GetPrimitive("h_mima_outjetpv_DY");
  SD_ntuple->Draw("missingMassDileptonMulti>>h_mima_outjetpv_SD(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&fabs(leadingJetVz-pvVertexZ)>0.1&&is2PUproton==0)", "hist");
  h_mima_outjetpv_SD = (TH1D*)gPad->GetPrimitive("h_mima_outjetpv_SD");

  //Drawing Missing Mass Distance Leading Jet PV less than 1 mm
  Sign_ntuple->Draw("missingMassDileptonMulti>>h_mima_injetpv_Sign(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&fabs(leadingJetVz-pvVertexZ)<0.1&&is2PUproton==0&&is45PUproton==0&&is56PUproton==0)", "hist");
  h_mima_injetpv_Sign = (TH1D*)gPad->GetPrimitive("h_mima_injetpv_Sign");
  DY_ntuple->Draw("missingMassDileptonMulti>>h_mima_injetpv_DY(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&fabs(leadingJetVz-pvVertexZ)<0.1)", "hist");
  h_mima_injetpv_DY = (TH1D*)gPad->GetPrimitive("h_mima_injetpv_DY");
  SD_ntuple->Draw("missingMassDileptonMulti>>h_mima_injetpv_SD(2000,0.,2000.)", "(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&dileptonPt>40&&leadingLeptonPt>30&&secondLeptonPt>20&&fabs(leadingJetVz-pvVertexZ)<0.1)", "hist");
  h_mima_injetpv_SD = (TH1D*)gPad->GetPrimitive("h_mima_injetpv_SD");
*/
  //Drawing JZB1
  Sign_ntuple->Draw("JZBalance>>h_jzb1_Sign(400, -200.,200.)", cut_Z_excl, "hist");
  h_jzb1_Sign = (TH1D*)gPad->GetPrimitive("h_jzb1_Sign");
  DY_ntuple->Draw("JZBalance>>h_jzb1_DY(400, -200.,200.)", cut_Z_excl, "hist");
  h_jzb1_DY = (TH1D*)gPad->GetPrimitive("h_jzb1_DY");
  SD_ntuple->Draw("JZBalance>>h_jzb1_SD(400, -200.,200.)", cut_Z_excl, "hist");
  h_jzb1_SD = (TH1D*)gPad->GetPrimitive("h_jzb1_SD");

  //Drawing JZB2
  Sign_ntuple->Draw("JZBalance2>>h_jzb2_Sign(400, -200.,200.)", cut_Z_excl, "hist");
  h_jzb2_Sign = (TH1D*)gPad->GetPrimitive("h_jzb2_Sign");
  DY_ntuple->Draw("JZBalance2>>h_jzb2_DY(400, -200.,200.)", cut_Z_excl, "hist");
  h_jzb2_DY = (TH1D*)gPad->GetPrimitive("h_jzb2_DY");
  SD_ntuple->Draw("JZBalance2>>h_jzb2_SD(400, -200.,200.)", cut_Z_excl, "hist"); 
  h_jzb2_SD = (TH1D*)gPad->GetPrimitive("h_jzb2_SD");

  //Drawing Phi Difference between Dilepton and First Jet
  Sign_ntuple->Draw("fabs(leadingJetPhi-dileptonPhi)>>h_phijetZ_Sign(63,0.,6.28)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6))", "hist");
  h_phijetZ_Sign = (TH1D*)gPad->GetPrimitive("h_phijetZ_Sign");
  DY_ntuple->Draw("fabs(leadingJetPhi-dileptonPhi)>>h_phijetZ_DY(63,0.,6.28)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6))", "hist");
  h_phijetZ_DY = (TH1D*)gPad->GetPrimitive("h_phijetZ_DY");
  SD_ntuple->Draw("fabs(leadingJetPhi-dileptonPhi)>>h_phijetZ_SD(63,0.,6.28)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6))", "hist");
  h_phijetZ_SD = (TH1D*)gPad->GetPrimitive("h_phijetZ_SD");
  
  //Drawing Phi Difference between MissingMass and MET
  Sign_ntuple->Draw("fabs(missingPhiDileptonMulti-missEt_phi)>>h_phimetma_Sign(63,0.,6.28)", cut_Z_excl, "hist");
  h_phimetma_Sign = (TH1D*)gPad->GetPrimitive("h_phimetma_Sign");
  DY_ntuple->Draw("fabs(missingPhiDileptonMulti-missEt_phi)>>h_phimetma_DY(63,0.,6.28)", cut_Z_excl, "hist");
  h_phimetma_DY = (TH1D*)gPad->GetPrimitive("h_phimetma_DY");
  SD_ntuple->Draw("fabs(missingPhiDileptonMulti-missEt_phi)>>h_phimetma_SD(63,0.,6.28)", cut_Z_excl, "hist");
  h_phimetma_SD = (TH1D*)gPad->GetPrimitive("h_phimetma_SD");

  //Drawing number of charged Particle associated to the PV
  Sign_ntuple->Draw("nChargedPFMultiPV_Tight_Fit>>h_nchPF_Sign(100,0.,100)", cut_Z_excl, "hist");
  h_nchPF_Sign = (TH1D*)gPad->GetPrimitive("h_nchPF_Sign");
  DY_ntuple->Draw("nChargedPFMultiPV_Tight_Fit>>h_nchPF_DY(100,0.,100)", cut_Z_excl, "hist");
  h_nchPF_DY = (TH1D*)gPad->GetPrimitive("h_nchPF_DY");
  SD_ntuple->Draw("nChargedPFMultiPV_Tight_Fit>>h_nchPF_SD(100,0.,100)", cut_Z_excl, "hist");
  h_nchPF_SD = (TH1D*)gPad->GetPrimitive("h_nchPF_SD");

  //Drawing pt sum of charged particles associated to the PV
  Sign_ntuple->Draw("SumChargedPFMultiPV_pt_Tight_Fit>>h_sumchPF_Sign(300,0.,300)", cut_Z_excl, "hist");
  h_sumchPF_Sign = (TH1D*)gPad->GetPrimitive("h_sumchPF_Sign");
  DY_ntuple->Draw("SumChargedPFMultiPV_pt_Tight_Fit>>h_sumchPF_DY(300,0.,300)", cut_Z_excl, "hist");
  h_sumchPF_DY = (TH1D*)gPad->GetPrimitive("h_sumchPF_DY");
  SD_ntuple->Draw("SumChargedPFMultiPV_pt_Tight_Fit>>h_sumchPF_SD(300,0.,300)", cut_Z_excl, "hist");
  h_sumchPF_SD = (TH1D*)gPad->GetPrimitive("h_sumchPF_SD");

  //Drawing ratio pt sum of charged particles associated to the PV & dilepton
  Sign_ntuple->Draw("((SumChargedPFMultiPV_pt_Tight_Fit / dileptonPt) - 1)>>h_ptvertZ_Sign(1000,-1,19)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(leadingJetPt<30&&missEt>30)&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_ptvertZ_Sign = (TH1D*)gPad->GetPrimitive("h_ptvertZ_Sign");
  DY_ntuple->Draw("((SumChargedPFMultiPV_pt_Tight_Fit / dileptonPt) - 1)>>h_ptvertZ_DY(1000,-1,19)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(leadingJetPt<30&&missEt>30)&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_ptvertZ_DY = (TH1D*)gPad->GetPrimitive("h_ptvertZ_DY");
  SD_ntuple->Draw("((SumChargedPFMultiPV_pt_Tight_Fit / dileptonPt) - 1)>>h_ptvertZ_SD(1000,-1,19)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(leadingJetPt<30&&missEt>30)&&missEt/dileptonPt>0.7&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_ptvertZ_SD = (TH1D*)gPad->GetPrimitive("h_ptvertZ_SD");

  //Drawing ratio met/dileptonpt
  Sign_ntuple->Draw("missEt/dileptonPt>>h_ratio_metpt_Sign(1000,0,20)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_ratio_metpt_Sign = (TH1D*)gPad->GetPrimitive("h_ratio_metpt_Sign");
  DY_ntuple->Draw("missEt/dileptonPt>>h_ratio_metpt_DY(1000,0,20)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_ratio_metpt_DY = (TH1D*)gPad->GetPrimitive("h_ratio_metpt_DY");
  SD_ntuple->Draw("missEt/dileptonPt>>h_ratio_metpt_SD(1000,0,20)", "weight*(isprotonMulti==1&&dileptonCharge==0&&missingMassDileptonMulti>0&&dileptonMass>81&&dileptonMass<101&&leadingLeptonPt>30&&secondLeptonPt>20&&nVertex>0&&xi_multiArm45>0.035&&xi_multiArm56>0.035&&xi_multiArm45<0.18&&xi_multiArm56<0.18&&dileptonPt>40&&(((SumChargedPFMultiPV_pt_Tight_Fit/dileptonPt)-1)>-0.2||(leadingJetPt<30&&missEt>30))&&(dileptonEta>-1.9&&dileptonEta<1.6)&&(fabs(leadingJetPhi-dileptonPhi)<2.2||fabs(leadingJetPhi-dileptonPhi)>4))", "hist");
  h_ratio_metpt_SD = (TH1D*)gPad->GetPrimitive("h_ratio_metpt_SD"); 
/*
  //Drawing Quark Gluon Discriminator for Leading Jet
  Sign_ntuple->Draw("leadingJetQGdis>>h_qgdis_Sign(200, -1., 1)", cut_Z_Sign, "hist");
  h_qgdis_Sign = (TH1D*)gPad->GetPrimitive("h_qgdis_Sign");
  DY_ntuple->Draw("leadingJetQGdis>>h_qgdis_DY(200, -1., 1)", cut_Z_DY, "hist");
  h_qgdis_DY = (TH1D*)gPad->GetPrimitive("h_qgdis_DY");
  SD_ntuple->Draw("leadingJetQGdis>>h_qgdis_SD(200, -1., 1)", cut_Z_SD, "hist");
  h_qgdis_SD = (TH1D*)gPad->GetPrimitive("h_qgdis_SD");
*/
  //Drawing Number of jets candidates (Loose criterium - 1000 um)
  Sign_ntuple->Draw("nJetsCandidatesLoose>>h_njetslos_Sign(10, 0., 10)", cut_Z_excl, "hist");
  h_njetslos_Sign = (TH1D*)gPad->GetPrimitive("h_njetslos_Sign");
  DY_ntuple->Draw("nJetsCandidatesLoose>>h_njetslos_DY(10, 0., 10)", cut_Z_excl, "hist");
  h_njetslos_DY = (TH1D*)gPad->GetPrimitive("h_njetslos_DY");
  SD_ntuple->Draw("nJetsCandidatesLoose>>h_njetslos_SD(10, 0., 10)", cut_Z_excl, "hist");
  h_njetslos_SD = (TH1D*)gPad->GetPrimitive("h_njetslos_SD"); 
/*
  //Drawing Number of jets candidates (Medium criterium - 500 um)
  Sign_ntuple->Draw("nJetsCandidatesMedium>>h_njetsmed_Sign(10, 0., 10)", cut_Z_Sign, "hist");
  h_njetsmed_Sign = (TH1D*)gPad->GetPrimitive("h_njetsmed_Sign");
  DY_ntuple->Draw("nJetsCandidatesMedium>>h_njetsmed_DY(10, 0., 10)", cut_Z_DY, "hist");
  h_njetsmed_DY = (TH1D*)gPad->GetPrimitive("h_njetsmed_DY");
  SD_ntuple->Draw("nJetsCandidatesMedium>>h_njetsmed_SD(10, 0., 10)", cut_Z_SD, "hist");
  h_njetsmed_SD = (TH1D*)gPad->GetPrimitive("h_njetsmed_SD"); 

  //Drawing Number of jets candidates (Tight criterium - 200 um)
  Sign_ntuple->Draw("nJetsCandidatesTight>>h_njetstig_Sign(10, 0., 10)", cut_Z_Sign, "hist");
  h_njetstig_Sign = (TH1D*)gPad->GetPrimitive("h_njetstig_Sign");
  DY_ntuple->Draw("nJetsCandidatesTight>>h_njetstig_DY(10, 0., 10)", cut_Z_DY, "hist");
  h_njetstig_DY = (TH1D*)gPad->GetPrimitive("h_njetstig_DY");
  SD_ntuple->Draw("nJetsCandidatesTight>>h_njetstig_SD(10, 0., 10)", cut_Z_SD, "hist");
  h_njetstig_SD = (TH1D*)gPad->GetPrimitive("h_njetstig_SD"); */

  //Drawing xi45 in Multi RP
  Sign_ntuple->Draw("xi_multiArm45>>h_xi45_multi_Sign(250, 0., 0.25)", cut_Z_excl, "hist");
  h_xi45_multi_Sign = (TH1D*)gPad->GetPrimitive("h_xi45_multi_Sign");
  DY_ntuple->Draw("xi_multiArm45>>h_xi45_multi_DY(250, 0., 0.25)", cut_Z_excl, "hist");
  h_xi45_multi_DY = (TH1D*)gPad->GetPrimitive("h_xi45_multi_DY");
  SD_ntuple->Draw("xi_multiArm45>>h_xi45_multi_SD(250, 0., 0.25)", cut_Z_excl, "hist");
  h_xi45_multi_SD = (TH1D*)gPad->GetPrimitive("h_xi45_multi_SD");

  //Drawing xi56 in Multi RP
  Sign_ntuple->Draw("xi_multiArm56>>h_xi56_multi_Sign(250, 0., 0.25)", cut_Z_excl, "hist");
  h_xi56_multi_Sign = (TH1D*)gPad->GetPrimitive("h_xi56_multi_Sign");
  DY_ntuple->Draw("xi_multiArm56>>h_xi56_multi_DY(250, 0., 0.25)", cut_Z_excl, "hist");
  h_xi56_multi_DY = (TH1D*)gPad->GetPrimitive("h_xi56_multi_DY");
  SD_ntuple->Draw("xi_multiArm56>>h_xi56_multi_SD(250, 0., 0.25)", cut_Z_excl, "hist");
  h_xi56_multi_SD = (TH1D*)gPad->GetPrimitive("h_xi56_multi_SD");

  //Output file
  TFile *outputFile;
  outputFile = new TFile("muon_plots_excl.root","RECREATE");

  h_delta_z1_DY   ->Write("h_delta_z1_DY");   
  h_delta_z1_SD   ->Write("h_delta_z1_SD");  
  h_delta_z1_Sign ->Write("h_delta_z1_Sign");
  /*h_delta_z2_DY   ->Write("h_delta_z2_DY");   
  h_delta_z2_SD   ->Write("h_delta_z2_SD");   
  h_delta_z2_Sign ->Write("h_delta_z2_Sign");*/
  h_pt_dil_DY     ->Write("h_pt_dil_DY");     
  h_pt_dil_SD     ->Write("h_pt_dil_SD");     
  h_pt_dil_Sign   ->Write("h_pt_dil_Sign");  
  h_eta_dil_DY    ->Write("h_eta_dil_DY");    
  h_eta_dil_SD    ->Write("h_eta_dil_SD");    
  h_eta_dil_Sign  ->Write("h_eta_dil_Sign");  
  /*h_pt_lep1_DY    ->Write("h_pt_lep1_DY");    
  h_pt_lep1_SD    ->Write("h_pt_lep1_SD");    
  h_pt_lep1_Sign  ->Write("h_pt_lep1_Sign");  
  h_eta_lep1_DY   ->Write("h_eta_lep1_DY");   
  h_eta_lep1_SD   ->Write("h_eta_lep1_SD");   
  h_eta_lep1_Sign ->Write("h_eta_lep1_Sign"); 
  h_pt_lep2_DY    ->Write("h_pt_lep2_DY");    
  h_pt_lep2_SD    ->Write("h_pt_lep2_SD");    
  h_pt_lep2_Sign  ->Write("h_pt_lep2_Sign");  
  h_eta_lep2_DY   ->Write("h_eta_lep2_DY");  
  h_eta_lep2_SD   ->Write("h_eta_lep2_SD");  
  h_eta_lep2_Sign ->Write("h_eta_lep2_Sign");*/
  h_pt_jet1_DY    ->Write("h_pt_jet1_DY");
  h_pt_jet1_SD    ->Write("h_pt_jet1_SD");
  h_pt_jet1_Sign  ->Write("h_pt_jet1_Sign");
  h_eta_jet1_DY   ->Write("h_eta_jet1_DY");
  h_eta_jet1_SD   ->Write("h_eta_jet1_SD");
  h_eta_jet1_Sign ->Write("h_eta_jet1_Sign");
  /*h_pt_jet2_DY    ->Write("h_pt_jet2_DY");
  h_pt_jet2_SD    ->Write("h_pt_jet2_SD");
  h_pt_jet2_Sign  ->Write("h_pt_jet2_Sign");
  h_eta_jet2_DY   ->Write("h_eta_jet2_DY");
  h_eta_jet2_SD   ->Write("h_eta_jet2_SD");
  h_eta_jet2_Sign ->Write("h_eta_jet2_Sign");*/
  h_met_DY        ->Write("h_met_DY");
  h_met_SD        ->Write("h_met_SD");
  h_met_Sign      ->Write("h_met_Sign");
  /*h_rmpf_DY       ->Write("h_rmpf_DY");
  h_rmpf_SD       ->Write("h_rmpf_SD");
  h_rmpf_Sign     ->Write("h_rmpf_Sign");*/
  h_miss_ma_120_DY->Write("h_miss_ma_120_DY");
  h_miss_ma_130_DY->Write("h_miss_ma_130_DY"); 
  h_miss_ma_140_DY->Write("h_miss_ma_140_DY"); 
  h_miss_ma_150_DY->Write("h_miss_ma_150_DY"); 
  h_miss_ma_120_SD->Write("h_miss_ma_120_SD");
  h_miss_ma_130_SD->Write("h_miss_ma_130_SD");
  h_miss_ma_140_SD->Write("h_miss_ma_140_SD");
  h_miss_ma_150_SD->Write("h_miss_ma_150_SD");
  h_miss_ma_120_Sign->Write("h_miss_ma_120_Sign");
  h_miss_ma_130_Sign->Write("h_miss_ma_130_Sign");
  h_miss_ma_140_Sign->Write("h_miss_ma_140_Sign");
  h_miss_ma_150_Sign->Write("h_miss_ma_150_Sign");
  h_miss_ma_DY    ->Write("h_miss_ma_DY");
  h_miss_ma_SD    ->Write("h_miss_ma_SD");
  h_miss_ma_Sign  ->Write("h_miss_ma_Sign");
  /*h_miss_ma_120_pre_Sign  ->Write("h_miss_ma_120_pre_Sign");
  h_miss_ma_120_post_Sign ->Write("h_miss_ma_120_post_Sign");*/
  h_diffma_DY	  ->Write("h_diffma_DY");
  h_diffma_SD     ->Write("h_diffma_SD");
  h_diffma_Sign   ->Write("h_diffma_Sign");
  h_miss_ma_fake_SD->Write("h_miss_ma_fake_SD");
  h_jzb1_DY       ->Write("h_jzb1_DY");
  h_jzb1_SD       ->Write("h_jzb1_SD");
  h_jzb1_Sign     ->Write("h_jzb1_Sign");
  h_jzb2_DY       ->Write("h_jzb2_DY");
  h_jzb2_SD       ->Write("h_jzb2_SD");
  h_jzb2_Sign     ->Write("h_jzb2_Sign");
  h_phijetZ_DY    ->Write("h_phijetZ_DY");
  h_phijetZ_SD    ->Write("h_phijetZ_SD");
  h_phijetZ_Sign  ->Write("h_phijetZ_Sign");
  h_nchPF_DY      ->Write("h_nchPF_DY");
  h_nchPF_Sign    ->Write("h_nchPF_Sign");
  h_nchPF_SD      ->Write("h_nchPF_SD");
  h_ptvertZ_Sign  ->Write("h_ptvertZ_Sign");
  h_ptvertZ_SD    ->Write("h_ptvertZ_SD");  
  h_ptvertZ_DY    ->Write("h_ptvertZ_DY");
  h_ratio_metpt_Sign  ->Write("h_ratio_metpt_Sign");
  h_ratio_metpt_SD    ->Write("h_ratio_metpt_SD");
  h_ratio_metpt_DY    ->Write("h_ratio_metpt_DY");
  h_sumchPF_DY    ->Write("h_sumchPF_DY");
  h_sumchPF_Sign  ->Write("h_sumchPF_Sign");
  h_sumchPF_SD    ->Write("h_sumchPF_SD");
  h_phimetma_DY   ->Write("h_phimetma_DY");
  h_phimetma_Sign ->Write("h_phimetma_Sign");
  h_phimetma_SD   ->Write("h_phimetma_SD");
  /*h_qgdis_DY	  ->Write("h_qgdis_DY"); 
  h_qgdis_Sign    ->Write("h_qgdis_Sign");
  h_qgdis_SD      ->Write("h_qgdis_SD");*/
  h_njetslos_DY   ->Write("h_njetslos_DY");
  h_njetslos_Sign ->Write("h_njetslos_Sign");
  h_njetslos_SD   ->Write("h_njetslos_SD");
  /*h_njetsmed_DY   ->Write("h_njetsmed_DY");
  h_njetsmed_Sign ->Write("h_njetsmed_Sign");
  h_njetsmed_SD   ->Write("h_njetsmed_SD");
  h_njetstig_DY   ->Write("h_njetstig_DY");
  h_njetstig_Sign ->Write("h_njetstig_Sign");
  h_njetstig_SD   ->Write("h_njetstig_SD"); */
  h_xi45_multi_DY  ->Write("h_xi45_multi_DY");
  h_xi45_multi_Sign->Write("h_xi45_multi_Sign");
  h_xi45_multi_SD  ->Write("h_xi45_multi_SD");
  h_xi56_multi_DY  ->Write("h_xi56_multi_DY");
  h_xi56_multi_Sign->Write("h_xi56_multi_Sign");
  h_xi56_multi_SD  ->Write("h_xi56_multi_SD");
  /*h_mima_outZ_DY      ->Write("h_mima_outZ_DY");
  h_mima_outZ_SD      ->Write("h_mima_outZ_SD");
  h_mima_outZ_Sign    ->Write("h_mima_outZ_Sign");
  h_mima_inZ_DY       ->Write("h_mima_inZ_DY");
  h_mima_inZ_SD       ->Write("h_mima_inZ_SD");
  h_mima_inZ_Sign     ->Write("h_mima_inZ_Sign");
  h_mima_outjetpv_DY  ->Write("h_mima_outjetpv_DY");
  h_mima_outjetpv_SD  ->Write("h_mima_outjetpv_SD");
  h_mima_outjetpv_Sign->Write("h_mima_outjetpv_Sign");
  h_mima_injetpv_DY   ->Write("h_mima_injetpv_DY");
  h_mima_injetpv_SD   ->Write("h_mima_injetpv_SD");
  h_mima_injetpv_Sign ->Write("h_mima_injetpv_Sign");*/

  outputFile->Close();

}
