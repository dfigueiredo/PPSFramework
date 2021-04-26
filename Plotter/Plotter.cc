#include <iostream> 
#include <iomanip>
#include <ctime>
#include <algorithm>

//ROOT
#include "TFile.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TBox.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TColor.h"

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

bool sigmaXi(double nsigma, double x, double y){

  double sqx=(0.055*x)*(0.055*x);
  double sqy=(0.039*y)*(0.039*y);
  double sum=sqx+sqy;
  double combinedError=sqrt(sum);

  long double delta = fabs(x-y);

  if( (delta/combinedError) <= nsigma ){ return 1;}
  return 0;

}

void set_plot_style(){
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void DrawHistoTH1(TTree *T1, TTree *T2, const char *variable, TString title, int bin, int xmin, int xmax, TCut cutsig, TCut cutbkg, bool hidden, const char *name){

  TString varsig;
  varsig.Form("%s>>%s%s_sig", variable, variable, name);

  TString varbkg;
  varbkg.Form("%s>>%s%s_bkg", variable, variable, name);

  TString hsig;
  hsig.Form("%s%s_sig", variable, name);

  TString hbkg;
  hbkg.Form("%s%s_bkg", variable, name);

  TH1F * htemp_sig = new TH1F(hsig, title, bin, xmin, xmax);
  TH1F * htemp_bkg = new TH1F(hbkg, title, bin, xmin, xmax);

  T1->Draw(varsig,cutsig,"goff");
  T2->Draw(varbkg,cutbkg,"goff");

  htemp_bkg->SetLineColor(kRed);
  htemp_sig->SetLineColor(kBlack);
  htemp_sig->SetMarkerStyle(20);
  htemp_bkg->SetMarkerSize(0.8);
  htemp_sig->Sumw2();
  htemp_bkg->Sumw2();
  htemp_sig->Scale(1/htemp_sig->GetSumOfWeights());
  htemp_bkg->Scale(1/htemp_bkg->GetSumOfWeights());
  htemp_bkg->GetYaxis()->SetRangeUser(0.,htemp_sig->GetMaximum()+0.1*htemp_sig->GetMaximum());
  htemp_sig->GetYaxis()->SetRangeUser(0.,htemp_bkg->GetMaximum()+0.1*htemp_bkg->GetMaximum());
  htemp_sig->Draw("ep");
  htemp_bkg->Draw("histosame");
  htemp_sig->GetYaxis()->SetTitleOffset(1.62);
  htemp_bkg->GetYaxis()->SetTitleOffset(1.62);

  // Blinding
  if(hidden){
    TBox *box = new TBox(1000,0,1250,htemp_bkg->GetMaximum()+0.1*htemp_bkg->GetMaximum());
    box->SetFillColor(kYellow);
    box->Draw();
    gPad->RedrawAxis();
  }

}

void DrawHistoTH1(TTree *T1, const char *variable, TString title, int bin, int xmin, int xmax, TCut cutsig, TCut cutbkg, bool hidden, const char *name){

  TString varsig;
  varsig.Form("%s>>%s%s_sig", variable, variable, name);

  TString varbkg;
  varbkg.Form("%s>>%s%s_bkg", variable, variable, name);

  TString hsig;
  hsig.Form("%s%s_sig", variable, name);

  TString hbkg;
  hbkg.Form("%s%s_bkg", variable, name);

  TH1F * htemp_sig = new TH1F(hsig, title, bin, xmin, xmax);
  TH1F * htemp_bkg = new TH1F(hbkg, title, bin, xmin, xmax);

  T1->Draw(varsig,cutsig,"goff");
  T1->Draw(varbkg,cutbkg,"goff");

  htemp_bkg->SetLineColor(kRed);
  htemp_sig->SetLineColor(kBlack);
  htemp_sig->SetMarkerStyle(20);
  htemp_bkg->SetMarkerSize(0.8);
  htemp_sig->Sumw2();
  htemp_bkg->Sumw2();
  htemp_sig->Scale(1/htemp_sig->GetSumOfWeights());
  htemp_bkg->Scale(1/htemp_bkg->GetSumOfWeights());
  htemp_bkg->GetYaxis()->SetRangeUser(0.,htemp_sig->GetMaximum()+0.1*htemp_sig->GetMaximum());
  htemp_sig->GetYaxis()->SetRangeUser(0.,htemp_bkg->GetMaximum()+0.1*htemp_bkg->GetMaximum());
  htemp_sig->Draw("ep");
  htemp_bkg->Draw("histosame");
  htemp_sig->GetYaxis()->SetTitleOffset(1.62);
  htemp_bkg->GetYaxis()->SetTitleOffset(1.62);

  // Blinding
  if(hidden){
    TBox *box = new TBox(1000,0,1250,htemp_bkg->GetMaximum()+0.1*htemp_bkg->GetMaximum());
    box->SetFillColor(kYellow);
    box->Draw();
    gPad->RedrawAxis();
  }

}

void DrawHistoTH2(TTree *T, const char *variable1, const char *variable2, TString title, int bin, int xmin, int xmax, TCut cut, const char *name){

  TString varsig;
  varsig.Form("%s:%s>>%s%s%s_h", variable1, variable2, variable1, variable2, name);

  TString h_name;
  h_name.Form("%s%s%s_h", variable1, variable2, name);

  TH2F * htemp = new TH2F(h_name, title, bin, xmin, xmax, bin, xmin, xmax);
  T->Draw(varsig,cut,"goff");

  htemp->Draw("COLZ");
  htemp->GetYaxis()->SetTitleOffset(1.62);
  htemp->GetYaxis()->SetTitleOffset(1.62);

}

void DrawHistoTH2(TTree *T, const char *variable1, const char *variable2, TString title, int bin1, int xmin1, int xmax1, int bin2, int xmin2, int xmax2, TCut cut, const char *name){

  TString varsig;
  varsig.Form("%s:%s>>%s%s%s_h", variable1, variable2, variable1, variable2, name);

  TString h_name;
  h_name.Form("%s%s%s_h", variable1, variable2, name);

  TH2F * htemp = new TH2F(h_name, title, bin1, xmin1, xmax1, bin2, xmin2, xmax2);
  T->Draw(varsig,cut,"goff");

  htemp->Draw("COLZ");
  htemp->GetYaxis()->SetTitleOffset(1.62);
  htemp->GetYaxis()->SetTitleOffset(1.62);

}

int main(int argc, char **argv){

  char * variable = getCmdOption(argv, argv + argc, "--variable");
  char * filesignal = getCmdOption(argv, argv + argc, "--fsignal");
  char * filebkg = getCmdOption(argv, argv + argc, "--fbkg");
  char * signalcut = getCmdOption(argv, argv + argc, "--signalcuts");
  char * bkgcut = getCmdOption(argv, argv + argc, "--bkgcuts");
  char * tagname = getCmdOption(argv, argv + argc, "--tagname");
  char * firstbin = getCmdOption(argv, argv + argc, "--first_bin");
  char * lastbin = getCmdOption(argv, argv + argc, "--last_bin");
  char * binsize = getCmdOption(argv, argv + argc, "--binsize");
  char * title = getCmdOption(argv, argv + argc, "--title");

  if(cmdOptionExists(argv, argv+argc, "--h")||cmdOptionExists(argv, argv+argc, "--help"))
  {
    std::cout << "\n== Help ==\n" << std::endl;
    std::cout << "\t --y year (--y 2018)\n" << std::endl;
    return 0;
  }

  clock_t begin = clock();

  gStyle->SetOptStat(1001100);
  //gStyle->SetOptStat(0);
  gStyle->SetPalette(kDarkBodyRadiator);
  set_plot_style();

  TCut signal = signalcut;
  TCut bkg = bkgcut;

  TFile *f1 = new TFile(filesignal);
  TTree *T1 = (TTree*)f1->Get("Events");

  // Missing Mass
  Double_t wc1 = 1500;
  Double_t hc1 = 800;
  TCanvas * c1 = new TCanvas(tagname, tagname, wc1, hc1);
  TString filename = TString(tagname) + ".pdf";

  c1->SetWindowSize(wc1 + (wc1 - c1->GetWw()), hc1 + (hc1 - c1->GetWh()));
  c1->cd();

  if((std::string)filebkg=="none"){
    TFile *f1 = new TFile(filesignal);
    TTree *T1 = (TTree*)f1->Get("Events");
    DrawHistoTH1(T1, variable, title, std::stoi(binsize), std::stoi(firstbin), std::stoi(lastbin), signal, bkg, false, tagname);
  }else{
    TFile *f1 = new TFile(filesignal);
    TFile *f2 = new TFile(filebkg);
    TTree *T1 = (TTree*)f1->Get("Events");
    TTree *T2 = (TTree*)f2->Get("Events");
    DrawHistoTH1(T1, T2, variable, title, std::stoi(binsize), std::stoi(firstbin), std::stoi(lastbin), signal, bkg, false, tagname);
  }
  c1->SaveAs(filename);

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << "Time: " << elapsed_secs << std::endl;

  return 0;

}

