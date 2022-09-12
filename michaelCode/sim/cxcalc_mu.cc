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

std::pair<double, double> lheparse(std::string filebase, std::string madfname)
{
   madfname = filebase + madfname; 
   std::ifstream ifile;
   ifile.open(madfname.c_str());
   std::string s;
   double ebeam, madcx;
   for(int j=0; j<432; j++)
   {
      std::getline(ifile, s);
      if(j==244)
      {
         std::istringstream iss(s);
         std::string temp, ebeam1;
         iss >> ebeam >> temp >> ebeam1;
         if(ebeam1!="ebeam1")
         {
            printf("Incorrect ebeam line found, bad filename is %s.\n", madfname.c_str());
            printf("%s\n",ebeam1.c_str());
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
   //printf("found energy %f and cross section %f.\n",ebeam,madcx);
   return std::make_pair(ebeam, madcx);
}

int main(int argc, char* argv[])
{
   srand48(std::atof(argv[3]));
   double mAp=0.2; 
   std::string fname = "map0p2_WWcx.root";
   TFile *f = new TFile(fname.c_str(),"recreate");
   
   DarkPhotons* dphoton = new DarkPhotons(mAp, 0, 1, 64, 29, 8.96, 1, "temp.txt");
   dphoton->PrepareTable();
 
   int n = 500;
   double x[n], y[n], z[n], a[n], b[n], c[n];
   for(int i=0; i<n; i++)
   {
      x[i] = 2*i;
      y[i] = dphoton->TotalMuCrossSectionCalc(2*i);
      z[i] = dphoton->DMG4CrossSectionCalc(2*i,1);
      a[i] = dphoton->DMG4CrossSectionCalc(2*i,2);
      b[i] = dphoton->DMG4CrossSectionCalc(2*i,3);
      c[i] = dphoton->DMG4CrossSectionCalc(2*i,4);
  
   }
  
   TGraph *gr1 = new TGraph(n,x,y);
   gr1->Write("WW Approx");

   TGraph *gr4 = new TGraph(n,x,z);
   gr4->Write("DMG4 WW Approx");

   TGraph *gr5 = new TGraph(n,x,a);
   gr5->Write("DMG4 WW Approx 2");
   TGraph *gr6 = new TGraph(n,x,b);
   gr6->Write("DMG4 WW Approx 3");
   TGraph *gr7 = new TGraph(n,x,c);
   gr7->Write("DMG4 WW Approx 4");

   std::map<double, std::pair<double, int>> maddata;
   DIR *pDIR;
   struct dirent *entry;
   if( pDIR=opendir("/local/cms/user/revering/madgraph/Repository/text/AN_doc/map_0p2_Cu/") )
   {
      while(entry = readdir(pDIR))
      {
         if(strcmp(entry->d_name,".") != 0 && strcmp(entry->d_name, "..") != 0)
         {
            std::pair<double, double> madvals = lheparse("/local/cms/user/revering/madgraph/Repository/text/AN_doc/map_0p2_Cu/",entry->d_name);
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
   
   std::map<double, std::pair<double, int>> noinelmaddata;
   DIR *pnoinDIR;
   if( pnoinDIR=opendir("/local/cms/user/revering/madgraph/Repository/muon/noInelTerm/") )
   {
      while(entry = readdir(pnoinDIR))
      {
         if(strcmp(entry->d_name,".") != 0 && strcmp(entry->d_name, "..") != 0)
         {
            std::pair<double, double> madvals = lheparse("/local/cms/user/revering/madgraph/Repository/muon/noInelTerm/",entry->d_name);
            if(noinelmaddata.find(madvals.first)==noinelmaddata.end())
            {
                noinelmaddata[madvals.first] = std::make_pair(madvals.second,1);
                printf("Found energy %f.\n",madvals.first);
            }
            else
            {
                noinelmaddata[madvals.first].first += madvals.second;
                noinelmaddata[madvals.first].second++;
            }
         }
      }
   }
   
   double madx[maddata.size()], mady[maddata.size()], ratx[maddata.size()], raty[maddata.size()], ratz[maddata.size()], rata[maddata.size()],ratb[maddata.size()],ratc[maddata.size()];
   double noinelmadx[noinelmaddata.size()], noinelmady[noinelmaddata.size()];
   int i=0;
   for(auto const& imap: maddata)
   {
      madx[i] = imap.first;
      mady[i] = imap.second.first/imap.second.second;
      raty[i] = dphoton->TotalMuCrossSectionCalc(imap.first)/(imap.second.first/imap.second.second);
      ratz[i] = dphoton->DMG4CrossSectionCalc(imap.first,1)/(imap.second.first/imap.second.second);
      rata[i] = dphoton->DMG4CrossSectionCalc(imap.first,2)/(imap.second.first/imap.second.second);
      ratb[i] = dphoton->DMG4CrossSectionCalc(imap.first,3)/(imap.second.first/imap.second.second);
      ratc[i] = dphoton->DMG4CrossSectionCalc(imap.first,4)/(imap.second.first/imap.second.second);
      i++;
   }
   int j=0;
   for(auto const& imap: noinelmaddata)
   {
      noinelmadx[j] = imap.first;
      noinelmady[j] = imap.second.first/imap.second.second;
      j++;
   }

   TGraph * gr2 = new TGraph(maddata.size(), madx, mady);
   gr2->Write("Madgraph data");

   TGraph * gr3 = new TGraph(maddata.size(), madx, raty);
   gr3->Write("WWOverMG");

   TGraph * gr8 = new TGraph(maddata.size(), madx, ratz);
   gr8->Write("DMG4OverMG");
   TGraph * gr9 = new TGraph(maddata.size(), madx, rata);
   gr9->Write("DMG42OverMG");
   TGraph * gr10 = new TGraph(maddata.size(), madx, ratb);
   gr10->Write("DMG43OverMG");
   TGraph * gr11 = new TGraph(maddata.size(), madx, ratc);
   gr11->Write("DMG44OverMG");
   TGraph * gr12 = new TGraph(noinelmaddata.size(), noinelmadx, noinelmady);
   gr12->Write("Madgraph, no Inelastic Term");

   f->Write();
   f->Close();

   return 0;
}

