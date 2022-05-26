 #!/usr/bin/env python

import argparse
from time import strftime
import os
import math
import ROOT

def cx_comp(rootfile):
   if not os.path.isfile(rootfile):
        print "Input file does not exist!"
        quit()
   
   if infileType != ".root":
        print "Input file is of incorrect type \"%s\"!"%(infileType)
        quit()
   
   rfile = ROOT.TFile(rootfile,"READ")
   wwgraph = rfile.Get("WW Approx")
   madgraph = rfile.Get("Madgraph data")
   multgraph = ROOT.TMultiGraph()
   
   plot = ROOT.TCanvas("Cx comparison","Cx comparison", 2)
   ROOT.gStyle.SetPadTickY(1)
   ROOT.gStyle.SetPadTickX(1)
   pad1 = ROOT.TPad("pad1","pad1",0,0,1.0,1.0)
   pad1.SetTopMargin(0.05)
   pad1.SetBottomMargin(0.1)
   pad1.SetRightMargin(0.03)
   pad1.SetLeftMargin(0.13)
   pad1.Draw()
   pad1.cd()
   wwgraph.SetTitle("")
   wwgraph.SetLineColor(2)
   wwgraph.SetMarkerColor(2)
   wwgraph.SetMarkerSize(0.2)
   wwgraph.SetMarkerStyle(21)
   multgraph.Add(wwgraph)
   
   madgraph.SetTitle("")
   #madgraph.SetLineColor(2)
   #madgraph.SetMarkerColor(2)
   madgraph.SetMarkerSize(1.0)
   madgraph.SetMarkerStyle(22)
   multgraph.Add(madgraph)
   multgraph.Draw("AP")
   multgraph.GetYaxis().SetTitle("Cross section/#epsilon^{2} (pb)")
   multgraph.GetYaxis().SetTitleSize(23)
   multgraph.GetYaxis().SetTitleFont(43)
   multgraph.GetYaxis().SetTitleOffset(1.25)
   multgraph.GetYaxis().SetLabelFont(43)
   multgraph.GetYaxis().SetLabelSize(18)
   multgraph.GetYaxis().SetMaxDigits(3)
   
   leg_x0 = 0.55
   leg_y0 = 0.25
   leg = ROOT.TLegend(leg_x0,leg_y0,leg_x0+0.35,leg_y0+0.15)
   leg.SetBorderSize(0)
   leg.AddEntry(wwgraph, "WW Approximation","l")
   leg.AddEntry(madgraph, "Madgraph Simulation","p")
   leg.Draw()
   ROOT.gStyle.SetTextFont(43)
   
   lab_x0 = 0.17
   lab_y0 = 0.85
   tag1 = ROOT.TLatex(lab_x0,lab_y0,"LDMX")
   tag1.SetNDC()
   tag1.SetTextFont(62)
   tag2 = ROOT.TLatex(lab_x0+0.15,lab_y0,"Internal")
   tag2.SetNDC()
   tag2.SetTextFont(52)
   tag1.SetTextSize(0.055)
   tag2.SetTextSize(0.045)
   tag1.Draw()
   tag2.Draw()
   tag3 = ROOT.TLatex(lab_x0,lab_y0-0.06,"m_{A'} = 100 MeV")
   tag3.SetNDC()
   tag3.SetTextFont(52)
   tag3.SetTextSize(0.05)
   tag3.Draw()
   
   xlabel = ROOT.TPaveText(0.4,0.0,0.99,0.06,"brNDC")
   xlabel.SetTextSize(23)
   xlabel.AddText("e^{-} Beam Energy (GeV)")
   xlabel.SetFillColor(0)
   xlabel.SetBorderSize(0)
   xlabel.Draw()
   img = ROOT.TImage.Create()
   img.FromPad(plot)
   img.WriteImage("cx_comp_0p1.png")
   plot.Print("cx_comp_0p1.pdf")
   plot.Print("cx_comp_0p1.C")
   
   rfile.Close()

if __name__=="__main__":
   parser = argparse.ArgumentParser(description="Parse Root Files to create pretty cross section comps")
   parser.add_argument("rootfile", help="name of files in directory")
   parser.add_argument("--o", "--outputdir"  , dest="outDir"    , help="output directory"          , default=os.getcwd())
   arg = parser.parse_args()
   
   outDir = arg.outDir
   #Check for trailing slash on ouput dir and delete
   if arg.outDir.split("/")[-1] == "": outDir = arg.outDir[:-1]
   infileType = "."+arg.rootfile.split(".")[-1]
   
   if not os.path.isdir(outDir):
        print "Output directory does not exist!"
        quit()

   cx_comp(arg.rootfile)   
