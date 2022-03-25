 #!/usr/bin/env python
#Quick script to make root files with ttree of events from lhe for plotting purposes
#Usage: python LHE_Root.py [file.lhe] --o outputDir
import argparse
from time import strftime
import os
import math
import ROOT
from array import array

parser = argparse.ArgumentParser(description="Parse LHE Files to create angle and energy distributions")
parser.add_argument("lhefile", help="name of LHE file to convert to root")
parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd())
arg = parser.parse_args()

outDir = arg.outDir
nevts = 0
#Check for trailing slash on ouput dir and delete
if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
infileType = "."+arg.lhefile.split(".")[-1]

if not os.path.isdir(outDir):
     print "Output directory does not exist!"
     quit()

if not os.path.isfile(arg.lhefile):
     print "Input file does not exist!"
     quit()

if infileType != ".lhe":
     print "Input file is of incorrect type \"%s\"!"%(infileType)
     quit()

f = open(arg.lhefile, "r")

linesList = f.readlines()
outFilename = outDir+"/"+arg.lhefile.split("/")[-1].split(".lhe")[0]

outfile = ROOT.TFile("%s.root"%(outFilename), "recreate")
tree = ROOT.TTree("Events", "Tree containing Lorentz Vectors of Dark Brem Products")
evec = ROOT.TLorentzVector()
filled = ROOT.TLorentzVector()
avec = ROOT.TLorentzVector()
nvec = ROOT.TLorentzVector()
tree.Branch( "IncidentParticle", "TLorentzVector", evec)
tree.Branch( "APrime", "TLorentzVector", avec)
tree.Branch( "Nucleus", "TLorentzVector", nvec)
i = 0
for j in range(0,len(linesList)):
    parsedLine = linesList[j].split(" ")
    parsedLine = filter(None, parsedLine)
    if parsedLine[0] == "11" and parsedLine[1] == "1":            #Find a final state "electron" (mu sim uses electrons with changed mass).
       i = i+1
       incomingline = linesList[j-2].split(" ")                   #Get the initial "electon".
       incomingline = filter(None, incomingline)
       evec.SetPxPyPzE(float(parsedLine[6]),float(parsedLine[7]),
                                      float(parsedLine[8]),float(parsedLine[9]))
       nucleusline = linesList[j+1].split(" ")
       nucleusline = filter(None, nucleusline)
       if nucleusline[0] == "-623" and nucleusline[1] == "1":
          nvec.SetPxPyPzE(float(nucleusline[6]),float(nucleusline[7]),
	                                float(nucleusline[8]),float(nucleusline[9]))
       else:
          nvec.SetPxPyPzE(0,0,0,0)
       dphotonline = linesList[j+2].split(" ")
       dphotonline = filter(None, dphotonline)
       if dphotonline[0] == "622" and dphotonline[1] == "1":
          avec.SetPxPyPzE(float(dphotonline[6]),float(dphotonline[7]),
	                                float(dphotonline[8]),float(dphotonline[9]))
       else:
          avec.SetPxPyPzE(0,0,0,0)
       tree.Fill()

outfile.Write()
outfile.Close()
f.close()

