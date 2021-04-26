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
#include "TTreePlayer.h"

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


int main(int argc, char **argv){

  char * file = getCmdOption(argv, argv + argc, "--file");

  if(cmdOptionExists(argv, argv+argc, "--h")||cmdOptionExists(argv, argv+argc, "--help"))
  {
    std::cout << "\n== Help ==\n" << std::endl;
    std::cout << "\t --y year (--y 2018)\n" << std::endl;
    return 0;
  }

  clock_t begin = clock();

  TFile *f1 = new TFile(file);
  TTree *T1 = (TTree*)f1->Get("Events");

  ((TTreePlayer*)(T1->GetPlayer()))->SetScanRedirect(true);
  ((TTreePlayer*)(T1->GetPlayer()))->SetScanFileName("event_list.txt"); 

  T1->Scan("run.AsInt():event.AsInt():lumiblock.AsInt()", "leadingLeptonPt>30&&secondLeptonPt>20&&fabs(dileptonMass-91)<10&&nVertex<20&&dileptonPt>40&&isprotonRP220==1&&xi_rp220_Arm45>0.035&&xi_rp220_Arm56>0.035","");

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << "Time: " << elapsed_secs << std::endl;

  return 0;

}

