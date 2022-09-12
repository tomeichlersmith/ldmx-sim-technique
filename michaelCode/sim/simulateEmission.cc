//To compile : g++ simulateEmission.cc -o simulateDbrem `root-config --cflags --glibs` DarkPhotons.cc -lgsl -lgslcblas -lm

#include "DarkPhotons.hh"
//#include "combsim.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>
#include <string>

//./whatever A'mass ebeam randomSeed scalefile outputfilename
int main(int argc, char* argv[])
{
   srand48(std::atof(argv[3]));
   std::string scalefile;
   double map, ebeam;
   ebeam = std::atof(argv[2]);
   map = std::atof(argv[1]); 
   scalefile = argv[4];
   std::string mname = argv[1];
   std::string ename = argv[2];
   std::string jnum = argv[3];
   std::string ofname = argv[5];
   double e_m = 0.1056;
   std::string fname = ofname + ".root";
   std::unique_ptr<TFile> f(TFile::Open(fname.c_str(),"recreate"));
   f->cd();
   TTree * tree = new TTree("ForwardEvents", "Tree containing fw only Lorentz Vectors from Geant");
   TLorentzVector * evec = new TLorentzVector();
   TLorentzVector * avec = new TLorentzVector();

   tree->Branch("IncidentParticle","TLorentzVector",evec);
   tree->Branch("APrime","TLorentzVector",avec);

   DarkPhotons* dphoton = new DarkPhotons(map, 0, 1, 28, 14, 2.32, 1, scalefile);
   dphoton->PrepareTable();
   TLorentzVector* pchange = new TLorentzVector();
   double e_mom, a_z, a_t, e_z, e_t, e_x, e_y, a_E, a_y, a_x;
/*
   for(int i=0; i<100; i++)
   {
      printf("Cross section at %d Gev is %e\n", i, dphoton->TotalCrossSectionCalc(i));
   }
*/
   for(int i=0;i<1000000;i++)
   {
      TLorentzVector* pchange = dphoton->SimulateEmission(ebeam, "forward_only");
      evec->SetPxPyPzE(pchange->X(),pchange->Y(),pchange->Z(),pchange->E());
      a_z = sqrt(ebeam*ebeam - e_m*e_m) - evec->Z();
      a_x = evec->X();
      a_y = evec->Y();
      a_E = sqrt(a_x*a_x+a_y*a_y+a_z*a_z+map*map);
      avec->SetPxPyPzE(a_x,a_y,a_z,a_E);
      tree->Fill();
   }
   f->cd();
   tree->Write("forward_only");
   TTree * cmtree = new TTree("cmEvents", "Tree containing cm scaled Lorentz Vectors from Geant");
   TLorentzVector * cmevec = new TLorentzVector();
   TLorentzVector * cmavec = new TLorentzVector();

   cmtree->Branch("IncidentParticle","TLorentzVector",cmevec);
   cmtree->Branch("APrime","TLorentzVector",cmavec);

   TLorentzVector* cmpchange = new TLorentzVector();
/*
   for(int i=0; i<100; i++)
   {
      printf("Cross section at %d Gev is %e\n", i, dphoton->TotalCrossSectionCalc(i));
   }
*/
   for(int i=0;i<10000;i++)
   {
      TLorentzVector* cmpchange = dphoton->SimulateEmission(ebeam, "cm_scaling");
      cmevec->SetPxPyPzE(cmpchange->X(),cmpchange->Y(),cmpchange->Z(),cmpchange->E());
      a_z = sqrt(ebeam*ebeam - e_m*e_m) - cmevec->Z();
      a_x = cmevec->X();
      a_y = cmevec->Y();
      a_E = sqrt(a_x*a_x+a_y*a_y+a_z*a_z+map*map);
      cmavec->SetPxPyPzE(a_x,a_y,a_z,a_E);
      cmtree->Fill();
   }
   f->cd();
   cmtree->Write("cm_scaling");
//   f->Write();
   f->Close();

   return 0;
}
