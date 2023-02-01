#include <iostream>
#include <TH1.h>
#include <TFile.h>
#include <TObjArray.h>

using namespace std;

class ProtonsRandom {
  public:
    void createHistogram(TString);
    void loadHistogram(TString);
    void FillRp210Arm45(double);
    void FillRp210Arm56(double);
    void FillRp220Arm45(double);
    void FillRp220Arm56(double);
    void FillMultiArm45(double);
    void FillMultiArm56(double);
    double GetXiRp210Arm45();
    double GetXiRp210Arm56();
    double GetXiRp220Arm45();
    double GetXiRp220Arm56();
    double GetXiMultiArm45();
    double GetXiMultiArm56();
    void Storing();

  private:

    TH1F* hXiRp210Arm45;
    TH1F* hXiRp210Arm56;
    TH1F* hXiRp220Arm45;
    TH1F* hXiRp220Arm56;
    TH1F* hXiMultiArm45;
    TH1F* hXiMultiArm56;

    TFile* fileout;
    TFile* filein;

};

void ProtonsRandom::createHistogram(TString filename){
  fileout = new TFile(filename, "RECREATE");
  hXiRp210Arm45 = new TH1F("hXiRp210Arm45","hXiRp210Arm45;#xi;entries",500,0,0.3);
  hXiRp210Arm56 = new TH1F("hXiRp210Arm56","hXiRp210Arm56;#xi;entries",500,0,0.3);
  hXiRp220Arm45 = new TH1F("hXiRp220Arm45","hXiRp220Arm45;#xi;entries",500,0,0.3);
  hXiRp220Arm56 = new TH1F("hXiRp220Arm56","hXiRp220Arm56;#xi;entries",500,0,0.3);
  hXiMultiArm45 = new TH1F("hXiMultiArm45","hXiMultiArm45;#xi;entries",500,0,0.3);
  hXiMultiArm56 = new TH1F("hXiMultiArm56","hXiMultiArm56;#xi;entries",500,0,0.3);
}

void ProtonsRandom::loadHistogram(TString filename){
  filein = new TFile(filename);
  try{
    hXiRp210Arm45 = (TH1F*)filein->Get("hXiRp210Arm45");
    hXiRp210Arm56 = (TH1F*)filein->Get("hXiRp210Arm56");
    hXiRp220Arm45 = (TH1F*)filein->Get("hXiRp220Arm45");
    hXiRp220Arm56 = (TH1F*)filein->Get("hXiRp220Arm56");
    hXiMultiArm45 = (TH1F*)filein->Get("hXiMultiArm45");
    hXiMultiArm56 = (TH1F*)filein->Get("hXiMultiArm56");
  }
  catch(...){
    std::cout << "\n\t ---> [ProtonsRandom] Corrupted file! Please try another file!\n" << std::endl; 
    exit(EXIT_FAILURE);
  }
}

void ProtonsRandom::FillRp210Arm45(double var){
  hXiRp210Arm45->Fill(var);
}

void ProtonsRandom::FillRp210Arm56(double var){
  hXiRp210Arm56->Fill(var);
}

void ProtonsRandom::FillRp220Arm45(double var){
  hXiRp220Arm45->Fill(var);
}

void ProtonsRandom::FillRp220Arm56(double var){
  hXiRp220Arm56->Fill(var);
}

void ProtonsRandom::FillMultiArm45(double var){
  hXiMultiArm45->Fill(var);
}

void ProtonsRandom::FillMultiArm56(double var){
  hXiMultiArm56->Fill(var);
}

double ProtonsRandom::GetXiRp210Arm45(){
  return hXiRp210Arm45->GetRandom();
}

double ProtonsRandom::GetXiRp210Arm56(){
  return hXiRp210Arm56->GetRandom();
}

double ProtonsRandom::GetXiRp220Arm45(){
  return hXiRp220Arm45->GetRandom();
}

double ProtonsRandom::GetXiRp220Arm56(){
  return hXiRp220Arm56->GetRandom();
}

double ProtonsRandom::GetXiMultiArm45(){
  return hXiMultiArm45->GetRandom();
}

double ProtonsRandom::GetXiMultiArm56(){
  return hXiMultiArm56->GetRandom();
}

void ProtonsRandom::Storing(){

  fileout->cd(); 
  hXiRp210Arm45->Write();
  hXiRp210Arm56->Write();
  hXiRp220Arm45->Write();
  hXiRp220Arm56->Write();
  hXiMultiArm45->Write();
  hXiMultiArm56->Write();
  fileout->Close(); 

}
