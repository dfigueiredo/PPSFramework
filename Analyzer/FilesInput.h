#include <iostream>
#include <list>
#include <TH1.h>
#include <TFile.h>
#include <TString.h>
#include <TObjArray.h>

class FilesInput {
  public:
    TString GetFileName(std::list<TString>&);
};

TString FilesInput::GetFileName(std::list<TString>& ids){

  TString filename;

  for (TString x : ids) {
    filename += x;
  }

  return filename;

}
