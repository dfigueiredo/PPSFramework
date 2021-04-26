void PlotterExample(){

  gStyle->SetOptStat(1001100);
  Double_t w = 1200;
  Double_t h = 350;
  TCanvas * c1 = new TCanvas("plotter", "plotter", w, h);
  c1->SetWindowSize(w + (w - c1->GetWw()), h + (h - c1->GetWh()));
  c1->Divide(2);

  TCut * ncut = new TCut("nLeptons>10");

  TFile *f = new TFile("/eos/cms/store/group/phys_pps/MissingMassSearch/NtuplesAOD/Ntuples_aod-v5/data_muon_aod_v5.root");
  TTree *T = (TTree*)f->Get("Events");

  // Case 1
  c1->cd(1);
  TH1F * hnLeptons = new TH1F("hnLeptons", "nLeptons", 100, 0, 100);
  T->Project("hnLeptons","nLeptons",*ncut);
  hnLeptons->SetLineColor(kRed);
  hnLeptons->GetXaxis()->SetTitle("nLeptons");
  hnLeptons->GetYaxis()->SetTitle("Multiplicity");
  hnLeptons->Draw(); 

  // Case 2
  c1->cd(2);
  T->Draw("nLeptons>>h2(60,0,60)",*ncut);
  TH1 *h2 = (TH1*)gPad->GetPrimitive("h2");
  h2->SetLineColor(kBlue);
  h2->GetXaxis()->SetTitle("nLeptons");
  h2->GetYaxis()->SetTitle("Multiplicity");

}
