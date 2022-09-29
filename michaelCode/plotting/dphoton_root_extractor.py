import uproot
import pandas
import numpy

def extractData(rootfile, treeName):
    inputFile=uproot.open(rootfile)
    events = inputFile[treeName]
    dataArray = events.arrays("IncidentParticle",library="pd")
    x=dataArray['IncidentParticle']['fP']['fX']
    y=dataArray['IncidentParticle']['fP']['fY']
    z=dataArray['IncidentParticle']['fP']['fZ']
    E=dataArray['IncidentParticle']['fE']
    df=pandas.DataFrame()
    df['Energy']=E
    pt=numpy.sqrt(x*x+y*y)
    angle=numpy.arctan(pt/z)
    df['Angle']=angle
    return df

def extractRescaledData(rootfile):
    return extractData(rootfile, "forward_only")

def extractMadgraphData(rootfile):
    return extractData(rootfile, "Events")