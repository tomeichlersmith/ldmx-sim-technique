 #!/usr/bin/env python

import argparse
from time import strftime
import os
import math
import ROOT
from array import array   
import plotter

def angle_comp_multi(rescaledFiles,madgraphFile):
   names = []
   
   for f in rescaledFiles:
      infileType = "."+f.split(".")[-1]
   
      if not os.path.isfile(f):
         print("Input file "+f+" does not exist!")
         quit()
      if infileType != ".root":
         print("Input file is of incorrect type \"%s\"!"%(infileType))
         quit()
      names.append(f)  
 
   if not os.path.isfile(madgraphFile):
        print("Input file does not exist!")
        quit()
   
   infileType = "."+madgraphFile.split(".")[-1]
   
   if infileType != ".root":
        print("Input file is of incorrect type \"%s\"!"%(infileType))
        quit()
   
   mghist = ROOT.TH1F("MadgraphData","MadgraphData",50,0,3.14)
   
   hists = []
   ebeam=-1
   mass=-1
   mu=False
   for name in names:
     rfile = ROOT.TFile(name,"READ")
     tree = rfile.Get("ForwardEvents")
     evec = ROOT.TLorentzVector()
     avec = ROOT.TLorentzVector()
     scaleDistance = name.split("/")[-1].split(".")[-2]
     if(scaleDistance.split("_")[0]=="mu"):
        scaleDistance=scaleDistance.split("_")[-1]
        mu=True
     histname = "Scale to "+scaleDistance
     x = ROOT.TH1F(histname,histname,50,0,3.14)
     tree.SetBranchAddress("IncidentParticle",evec)
     tree.SetBranchAddress("APrime",avec)
     entries = tree.GetEntries()
     tree.GetEntry(0)
     ebeam = round(evec.E()+avec.E())
     mass = round(avec.M(),2)
     for i in range(entries):
       tree.GetEntry(i)
       x.Fill(evec.Theta())
     x.SetDirectory(0)
     hists.append(x)
   
   mgfile = ROOT.TFile(madgraphFile,"READ")
   mgtree = mgfile.Get("Events")
   mgevec = ROOT.TLorentzVector()
   mgtree.SetBranchAddress("IncidentParticle",mgevec)
   mgentries = mgtree.GetEntries()
   mgtree.GetEntry(0)

   for i in range(mgentries):
      mgtree.GetEntry(i)
      mghist.Fill(mgevec.Theta())
   
   if not mu: 
      d = plotter.plot(hists,mghist, xtitle="Outgoing e- Angle", ytitle="Normalized Rate (1/#sigma d#sigma/d#theta)",imageTitle="scaleDistanceAngle",log=True,APMassTag=str(mass)+" GeV",y_min=0.00011,y_max=1.,lab_x0=0.68,lab_y0=0.88,normalize=True,x_max=3.14,leg_y0=0.53,leg_x0=0.65,ratMax=1.24,rebin=2)
   else:
      d = plotter.plot(hists,mghist, xtitle="Outgoing #mu- Angle", ytitle="Normalized Rate (1/#sigma d#sigma/d#theta)",imageTitle="scaleDistanceAngleMu",log=True,APMassTag=str(mass)+" GeV",y_min=0.000011,y_max=1.,lab_x0=0.58,lab_y0=0.8,normalize=True,x_max=3.14,leg_y0=0.45,ratMax=3.99, ratMin=0.51,rebin=2)

   rfile.Close()

if __name__=="__main__":
   parser = argparse.ArgumentParser(description="Parse Root Files to create angle and energy distributions")
   parser.add_argument("madgraphFile", help="name of files in directory")
   parser.add_argument("rescaledFiles", help="name of rescaled files in directory", nargs='+')
   parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd())
   arg = parser.parse_args()
   
   outDir = arg.outDir
   nevts = 0
   #Check for trailing slash on ouput dir and delete
   if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
   
   if not os.path.isdir(outDir):
      print "Output directory does not exist!"
      quit()

   angle_comp_multi(arg.rescaledFiles,arg.madgraphFile)
