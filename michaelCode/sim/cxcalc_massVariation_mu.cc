//To compile : g++ cxcalc_mu.cc -o simulate `root-config --cflags --glibs` DarkPhotons.cc -lgsl -lgslcblas -lm

#include "DarkPhotons.hh"
//#include "combsim.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>
#include <string>
#include <sstream>
#include "TGraph.h"
#include <dirent.h>
#include <map>

std::pair<double, double> lheparse(std::string madfname)
{
   madfname = "/local/cms/user/revering/madgraph/Repository/ebeam_120_mu_fe/" + madfname; 
   std::ifstream ifile;
   ifile.open(madfname.c_str());
   std::string s;
   double ebeam, madcx;
   for(int j=0; j<432; j++)
   {
      std::getline(ifile, s);
      if(j==183)
      {
         std::istringstream iss(s);
         std::string temp, apmass;
         iss >> temp >> ebeam >> temp >> apmass;
         if(apmass!="APMASS")
         {
            printf("Incorrect ebeam line found, bad filename is %s.\n", madfname.c_str());
            printf("%s\n",apmass.c_str());
         }
      }
      if(j==430)
      {
         std::istringstream iss(s);
         std::string temp, integrated;
         iss >> temp >> integrated >> temp >> temp >> temp >> madcx;
         if(integrated!="Integrated"){printf("Cross section line did not match expectation\n");}
      }
   }
   printf("found mass %f and cross section %f.\n",ebeam,madcx);
   return std::make_pair(ebeam, madcx);
}

int main(int argc, char* argv[])
{
   srand48(std::atof(argv[3]));
   double mAp=0.2; 
   std::string fname = "ebeam120_fe_WWcx.root";
   TFile *f = new TFile(fname.c_str(),"recreate");
  
   int n=40;
   double x[n], y[n], z[n]; 
   for(int i=0;i<n;i++)
   {
      DarkPhotons* dphoton = new DarkPhotons(i/20.+0.2, 0, 1, 56, 26, 8.96, 1, "temp.txt");
      dphoton->PrepareTable();
 
      x[i] = i/20.+0.2;
      y[i] = dphoton->TotalMuCrossSectionCalc(120);
   } 
   TGraph *gr1 = new TGraph(n,x,y);
   gr1->Write("WW Approx");

   std::map<double, std::pair<double, int>> maddata;
   DIR *pDIR;
   struct dirent *entry;
   if( pDIR=opendir("/local/cms/user/revering/madgraph/Repository/ebeam_120_mu_fe/") )
   {
      while(entry = readdir(pDIR))
      {
         if(strcmp(entry->d_name,".") != 0 && strcmp(entry->d_name, "..") != 0)
         {
            std::pair<double, double> madvals = lheparse(entry->d_name);
            if(maddata.find(madvals.first)==maddata.end())
            {
                maddata[madvals.first] = std::make_pair(madvals.second,1);
            }
            else
            {
                maddata[madvals.first].first += madvals.second;
                maddata[madvals.first].second++;
            }
         }
      }
   }
   
   double madx[maddata.size()], mady[maddata.size()], ratx[maddata.size()], raty[maddata.size()], ratz[maddata.size()];
   int i=0;
   for(auto const& imap: maddata)
   {
      madx[i] = imap.first;
      mady[i] = imap.second.first/imap.second.second;
//      raty[i] = dphoton->TotalMuCrossSectionCalc(imap.first)/(imap.second.first/imap.second.second);
      i++;
   }
   TGraph * gr2 = new TGraph(maddata.size(), madx, mady);
   gr2->Write("Madgraph data");

//   TGraph * gr3 = new TGraph(maddata.size(), madx, raty);
//   gr3->Write("WWOverMG");

   f->Write();
   f->Close();

   return 0;
}

