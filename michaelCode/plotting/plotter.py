import ROOT

def integrate_hist(histin,histout):
   nbins = histin.GetNbinsX()
   TotalInt = histin.Integral(1,nbins+1)
   for i in range(nbins+1):
      HistInt = histin.Integral(i+1,nbins+1)
      if(TotalInt>0):
         histout.SetBinContent(i,1-HistInt/TotalInt)
   return 

#When ratio is enabled "mainHist" is used as the common divisor.
def plot(hists, mainHist, ytitle="Set me", y_max=-1,y_min=-1,imageTitle="image",ratio=True, xtitle="Set me",leg_x0=0.55,leg_y0=0.25, log=True, normalize=False,ratMin=0.81,ratMax=1.09,ratioTitle="Rescaled/Madgraph",APMassTag="SET ME",lab_x0=0.725,lab_y0=0.8,x_max=1.,ebeam=100, rebin=-1):
   plot = ROOT.TCanvas(imageTitle+"_canvas","Plot comparison",2) 
   div_line = 0.35
   ROOT.gStyle.SetPadTickY(1)
   ROOT.gStyle.SetPadTickX(1)
   if(ratio):
      pad1 = ROOT.TPad("pad1","pad1",0,div_line,1.0,1.0)
      pad1.SetBottomMargin(0.)
      pad1.SetRightMargin(0.03)
   else:
      pad1 = ROOT.TPad("pad1","pad1",0,1,1.0,1.0)
      pad1.SetBottomMargin(0.1)
      pad1.SetRightMargin(0.03)
     
   pad1.Draw()
   pad1.cd()
   i=0
   mainHist.SetTitle("")
   mainHist.SetLineColor(1)
   mainHist.SetLineWidth(3)
   mainHist.SetLineStyle(7)
   if(rebin>0):
      mainHist.Rebin(rebin)
   if(normalize):
      mainHist.Scale(1/mainHist.GetEntries())
   mainHist.SetStats(0)
   mainHist.Draw("sameHist")
   mainHist.GetYaxis().SetTitle(ytitle)
   mainHist.GetYaxis().SetTitleSize(16)
   mainHist.GetYaxis().SetTitleFont(43)
   mainHist.GetYaxis().SetTitleOffset(1.2)
   mainHist.GetYaxis().SetLabelFont(43)
   mainHist.GetYaxis().SetLabelSize(15)
   if(y_min>0):
      mainHist.SetMinimum(y_min)
   if(y_max>0):
      mainHist.SetMaximum(y_max)
   if(log):
      pad1.SetLogy()
   leg = ROOT.TLegend(leg_x0,leg_y0,leg_x0+0.3,leg_y0+0.25)
   leg.SetBorderSize(0)
   leg.AddEntry(mainHist, "Madgraph","l")
   goodColors = [1,2,3,4,6]
   for hist in hists:
      hist.SetTitle("")
      hist.SetStats(0)
      if(rebin>0):
         hist.Rebin(rebin)
      #Use known-good colors, just increment from the last one once done
      if(i<len(goodColors)):
         hist.SetLineColor(goodColors[i])
      else:
         hist.SetLineColor(goodColors[-1]+i-len(goodColors)+1)
      i=i+1
      hist.SetLineWidth(2)
      if(normalize):
         hist.Scale(1/hist.GetEntries())
      hist.Draw("sameHist")
      leg.AddEntry(hist,hist.GetName(),"l")
    
   leg.Draw()
   ROOT.gStyle.SetTextFont(43)
#   label = ROOT.TPaveText(0.4,0.9,0.99,0.96,"brNDC")
#   label.SetTextSize(20)
#   label.AddText("A' mass = "+str(mA)+" GeV, Beam energy = "+str(ebeam)+" GeV")
#   label.SetFillColor(0)
#   label.SetBorderSize(0)
#   label.Draw() 
   #tag1 = ROOT.TLatex(lab_x0,lab_y0,"LDMX")
   #tag1.SetNDC()
   #tag1.SetTextFont(62)
   #tag2 = ROOT.TLatex(lab_x0+0.1,lab_y0,"Internal")
   #tag2.SetNDC()
   #tag2.SetTextFont(52)
   #tag1.SetTextSize(0.055)
   #tag2.SetTextSize(0.045)
   tag2 = ROOT.TLatex(lab_x0-0.08,lab_y0-0.42,"Madgraph sim at "+str(ebeam)+" GeV")
   tag2.SetNDC()
   tag2.SetTextFont(52)
   tag2.SetTextSize(0.05)
  
   #tag1.Draw()
   tag2.Draw()
   tag3 = ROOT.TLatex(lab_x0-0.02,lab_y0-0.06,"m_{A'} = "+APMassTag)
   tag3.SetNDC()
   tag3.SetTextFont(52)
   tag3.SetTextSize(0.05)
   tag3.Draw()
   plot.cd()
   if(ratio):
      pad2 = ROOT.TPad("pad2","pad2",0.0, 0.06, 1, div_line)
      pad2.SetTopMargin(0.0)
      pad2.SetBottomMargin(0.1)
      pad2.SetRightMargin(0.03)
      pad2.Draw()
      pad2.cd()
      temprat = mainHist.Clone()
      temprat.SetStats(0)
      temprat.Divide(mainHist*2)
      temprat.SetLineWidth(0)
      temprat.Draw("Hist")
      divs = []
      for hist in hists:
         ratio = hist.Clone("h3")
         ratio.SetDirectory(0)
         ratio.SetStats(0)
         ratio.Divide(mainHist)
         divs.append(ratio)
      j=0
      for div in divs:
         if(j<len(goodColors)):
            div.SetLineColor(goodColors[j])
         else:
            div.SetLineColor(goodColors[-1]+i+1-len(goodColors))
         j=j+1
         div.SetLineWidth(2)
         div.Draw("sameHist")
      temprat.SetMinimum(ratMin)
      temprat.SetMaximum(ratMax)
      temprat.SetTitle("")
      temprat.GetYaxis().SetTitle(ratioTitle)
      temprat.GetYaxis().CenterTitle()
      temprat.GetYaxis().SetTitleSize(13)
      temprat.GetYaxis().SetTitleFont(43)
      temprat.GetYaxis().SetTitleOffset(1.62)
      temprat.GetYaxis().SetLabelFont(43)
      temprat.GetYaxis().SetLabelSize(15)

   #temprat.GetXaxis().SetTitle(xtitle)
   #temprat.GetXaxis().SetTitleSize(20)
   #temprat.GetXaxis().SetTitleFont(43)
   #temprat.GetXaxis().SetTitleOffset(3.)
   temprat.GetXaxis().SetLabelFont(43)
   temprat.GetXaxis().SetLabelSize(15)
   temprat.GetXaxis().SetTickSize(0.07)
   
   line = ROOT.TLine(0,1,x_max,1)
   line.SetLineStyle(7)
   line.Draw()
#   ratio.SetLineWidth(2)
  
   plot.cd()
   xlabel = ROOT.TPaveText(0.4,0.0,0.99,0.06,"brNDC")
   xlabel.SetTextSize(16)
   xlabel.AddText(xtitle)
   xlabel.SetFillColor(0)
   xlabel.SetBorderSize(0)
   xlabel.Draw()
   img = ROOT.TImage.Create()
   img.FromPad(plot)
   img.WriteImage(imageTitle+".png")
   plot.Print(imageTitle+".pdf")
   plot.Print(imageTitle+".C")
   return plot
