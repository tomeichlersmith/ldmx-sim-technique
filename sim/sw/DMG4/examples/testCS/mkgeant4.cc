#include <math.h>
#include "globals.hh"
#include "DarkMatter.hh"
#include "DarkPhotons.hh"
//#include "DM_TYPE_NAME.hh"

#define nPointsE0 2
#define nPointsMass 3
#define RELATIVE_DELTA 0.01

int main() {

/*GeV*/double EThresh = 0.01; // for total CS testing applicable to mass range mA<1MeV
/*GeV*/double valuesE0[nPointsE0]={20.,100.};
/*MeV*/double testMassValues[nPointsMass]={16.7, 100., 1000.};
/*GeV*/double E0, massTested;
       double csCalcResult, csRefResult; //, epsilBench = 0.0001;
       int countNotEqual = 0;

       // include file with ETL reference values
       #include "DarkPhotonsSigmaTotETL.inc"
       //#include "DM_TYPE_NAMESigmaTotETL.inc"        
  
       for(int ii=0; ii<nPointsE0; ii++){
         E0 = valuesE0[ii];
         double* pETLii = pdataETL[ii];
         //std::out<<" "<<E0<<std::endl;
         for(int jj=0; jj<nPointsMass; jj++) {
           massTested = testMassValues[jj]/1000.; //convert to GeV    
           //std::cout<<std::setw(7)<<massTested;
           DarkMatter* myDarkMatter = new DarkPhotons(massTested, EThresh); // Initialize DM by default for Pb with eps=0.0001
           //DarkMatter* myDarkMatter = new DM_TYPE_NAME(massTested, EThresh); // Initialize DM by default for Pb with eps=0.0001
           //get calculated ETL cs from DMG4       
           csCalcResult = myDarkMatter->TotalCrossSectionCalc(E0);
           //get reference ETL cs from include file
           //csRefResult = GeVtoPb*epsilBench*epsilBench*(*pETLii); pETLii++;
           csRefResult = (*pETLii); pETLii++;
          
           //check if calculated tot cs is equal to reference cs 
           if (fabs(csCalcResult-csRefResult)/csRefResult<RELATIVE_DELTA) {std::cout<<" OK ";} 
           else {std::cout<<" != ";countNotEqual++;}
           std::cout<<std::setw(7)<<massTested<<" d = "<<fabs(csCalcResult-csRefResult)/csRefResult*100.<<"  REF: "<<csRefResult<<" vs "<<csCalcResult<<" :CALC"<<std::endl;
         }
       }
       if(countNotEqual) {std::cout<<" Error: Found "<<countNotEqual<<" deviation(s) from reference, exiting..."<<std::endl; exit(1);}
       return 0;
}
