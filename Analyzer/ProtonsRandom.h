#include <iostream>
#include <TH1.h>
#include <TFile.h>
#include <TObjArray.h>

using namespace std;

class ProtonsRandom {
  public:
    void createHistogram();
};

void ProtonsRandom::createHistogram(){

  char name[10], title[20];

  // Create an array of histograms.
  TObjArray Hlist(0);

  // Create a pointer to a histogram.
  TH1F* h;

  // Make and fill 15 histograms and add them to the object array.
  for (Int_t i = 0; i < 15; i++) {
    sprintf(name,"h%d",i);
    sprintf(title,"histo nr:%d",i);
    h = new TH1F(name,title,100,-4,4);
    Hlist.Add(h);
    h->FillRandom("gaus",1000);
  }

  // Open a ROOT file and write the array to the ROOT file.
  TFile f("demo.root","RECREATE");
  Hlist.Write();

  // Closing the ROOT file.
  f.Close(); 

}
