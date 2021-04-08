/* Routine for the conversion of starlight data format
 * into something readable by CMSSW_1_4_5
 *
 * Modification by X. Rouby on a routine kindly offered by J. Hollar
 * Sept. 28, 2007
 */

#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;

void makeEventsFile(int num=0)
{
  string filename = "jhforlhesuperchic7tevall.txt"; //input
  ifstream infile(filename.c_str());
  char outfilename[100];
  sprintf(outfilename,"chibtomumugamma.lhe",num);
  ofstream output(outfilename);
  if (! infile.is_open()) { cout << "\t ERROR: I can not open \"" << filename << "\"" << endl; return; }

  string temp_string, temp;
  istringstream curstring;
  const unsigned int M = 10000; // N_events
  const double MU = 0.105658369; // muon mass [GeV]
  double charge = 0.0;
  int particlecounter = 0;
  int nn = 0;

  //  Float_t px[5], py[5], pz[5], e[5], m[5];

  getline(infile,temp_string); // The very first line is useless
  output << "<LesHouchesEvents version=\"1.0\">"  << endl;
  output << "<header>" << endl;
  output << "This file was created from the output of the FPMC generator" << endl;
  output << "</header>" << endl;

  output << "<init>" << endl;
  output << "2212  2212  0.35000000000E+04  0.35000000000E+04 0 0 10042 10042 2  1" << endl;
  output << "0.10508723460E+01  0.96530000000E-02  0.26731120000E-03   0" << endl;
  output << "</init>" << endl;

  while (getline(infile,temp_string)) {

    curstring.clear(); // needed when using several tims istringstream::str(string)
    curstring.str(temp_string);

    if(particlecounter == 5) {
      double px, py, pz, energy, mass;
      curstring >> px >> py >> pz >> energy >> mass;
      output << "-13 1 1 2 0 0 " << px << " " << py << " " << pz << " " << energy << " " << mass << " 0. " << -1 << endl;
      particlecounter = 0;
    }
    if(particlecounter == 4) {
      double px, py, pz, energy, mass;
      curstring >> px >> py >> pz >> energy >> mass;
      output << "13 1 1 2 0 0 " << px << " " << py << " " << pz << " " << energy << " " << mass << " 0. " << 1 << endl;
      particlecounter = 5;
    }
    if(particlecounter == 3) {
      double px, py, pz, energy, mass;
      curstring >> px >> py >> pz >> energy >> mass;
      output << "22 1 1 2 0 0 " << px << " " << py << " " << pz << " " << energy << " " << mass << " 0. " << 1 << endl;
      particlecounter = 4;
    }
    if(particlecounter == 2) {
      double px, py, pz, energy, mass;
      curstring >> px >> py >> pz >> energy >> mass;
      output << "2212 1 1 2 0 0 " << px << " " << py << " " << pz << " " << energy << " " << mass << " 0. " << 1 << endl;
      particlecounter = 3;
    }
    if(particlecounter == 1) {
      double px, py, pz, energy, mass;
      curstring >> px >> py >> pz >> energy >> mass;
      output << "2212 1 1 2 0 0 " << px << " " << py << " " << pz << " " << energy << " " << mass << " 0. " << 1 << endl;
      particlecounter = 2;
    }

    if(strstr(temp_string.c_str(),"<event>")){
      /*
      TLorentzVector m1;
      TLorentzVector m2;
      TLorentzVector gam;
      m1.SetXYZM(px[2],py[2],pz[2],m[2]);
      m2.SetXYZM(px[3],py[3],pz[3],m[3]);
      gam.SetXYZM(px[4],py[4],pz[4],m[4]);
      TLorentzVector chi = m1 + m2 + gam;
      TLorentzVector ups = m1 + m2;
      Float_t mass = chi.M();
      */

      output << "<event>" << endl;
      //      output << "7   0  0.2983460E-04  0.9118800E+02  0.7546772E-02  0.1300000E+00" << endl;
      output << "7   0  0.2983460E-04  0.9118800E+02  0.7546772E-02  0.1300000E+00" << endl;      
      output << "22   -1    0    0    0    0  0.00000000000E+00  0.00000000000E+00 0.00000000000E+02  0.10000000000E+02  0.00000000000E+00 0.  1." << endl;
      output << "22   -1    0    0    0    0  0.00000000000E+00  0.00000000000E+00 0.00000000000E+00  0.10000000000E+02  0.00000000000E+00 0. -1." << endl;
      particlecounter = 1;
    }
    if(strstr(temp_string.c_str(),"</event>")){
      output << "</event>" << endl;
      particlecounter = 0;
      nn++;
      if (nn%10 == 0)
	cout << "Event " << nn << endl;
    }

  } // reading loop of the input file
  output << "</LesHouchesEvents>" << endl;
  infile.close();
  output.close();
  cout << nn << " events written in " << outfilename << endl;
  return;
}

void makefiles(int number_of_files=5) {
  for(int i=0; i<number_of_files; i++)
    makeEventsFile(i);
}
