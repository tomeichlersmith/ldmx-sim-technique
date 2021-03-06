#include "DarkPhotons.hh"

#include <sstream>

#define EPSPARINV 1.e-8
#define  nMA 24
#define  nE0 8

template <int NY>
double BilinearInterpolation(double X, double Y, double ArgX[], double ArgY[], double Func[][NY], int NX, int iprint)
{
  // searching for cell which bounds the point (MA, E0). Cell is numerated by its upper left index (i,j)
  int ii=-1, jj=-1;
  for (int i2 = 1; i2 < NX; i2++) {
    for (int j2 = 1; j2 < NY; j2++) {
      if( X <= ArgX[i2] && X > ArgX[i2-1] && Y <= ArgY[j2] && Y > ArgY[j2-1] ) {
         ii=i2;
         jj=j2;
         if(iprint) printf("i=%i j=%i , X=%.3f, Y=%.2f \n", ii, jj, X, Y );
      }
    }
  }
  if(ii < 0 || jj < 0) {std::cout << "Error in bilinear interpolation" << std::endl; exit(1);}
  // calculate interpolated value of K-factor
  double tangentX1=(X-ArgX[ii-1])/(ArgX[ii]-ArgX[ii-1]);
  double tangentY1=(Y-ArgY[jj-1])/(ArgY[jj]-ArgY[jj-1]);
  double tangentX2=(X-ArgX[ii])/(ArgX[ii-1]-ArgX[ii]);
  double tangentY2=(Y-ArgY[jj])/(ArgY[jj-1]-ArgY[jj]);
  double Result = Func[ii-1][jj-1]*tangentX2*tangentY2 +
      Func[ii-1][jj]*tangentX2*tangentY1   +
      Func[ii][jj-1]*tangentX1*tangentY2   +
      Func[ii][jj]*tangentX1*tangentY1;
  if(iprint) printf("ResultBilinear=%.4f \n", Result);
  return Result;
}
double KfactorApproximate(double MAtest, double E0test)
{
  double KKffactor[nMA][nE0];
  
  // mass of A' in GeV 
  double  MMAA[nMA] = {0.001, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.3, 0.4, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}; 
  // Beam energy of electrons in GeV
  double  EE00[nE0] = {20, 40, 60, 80, 100, 120, 140, 160}; 
  double KfMA1MeV[nE0] = {0.606286, 0.605283, 0.604929, 0.60474, 0.604626, 0.604549, 0.604485, 0.604527};
  double KfMA2MeV[nE0] = {0.837599, 0.83442, 0.833288, 0.832704, 0.832335, 0.832088, 0.831908, 0.831281};
  double KfMA4MeV[nE0] = {0.969624, 0.960278, 0.956935, 0.955195, 0.954116, 0.95338, 0.952846, 0.952438};
  double KfMA6MeV[nE0] = {1.04025, 1.02264, 1.01628, 1.01296, 1.0109, 1.00949, 1.00846, 1.00767};
  double KfMA8MeV[nE0] = {1.09678, 1.06925, 1.0592, 1.05392, 1.05064, 1.04838, 1.04674, 1.04549};
  double KfMA10MeV[nE0] = {1.14791, 1.10876, 1.09496, 1.0874, 1.08268, 1.07945, 1.07708, 1.07527};
  double KfMA20MeV[nE0] = {1.37278, 1.27614, 1.23623, 1.21422, 1.20018, 1.19038, 1.18313, 1.17754};
  double KfMA30MeV[nE0] = {1.54771, 1.41733, 1.3544, 1.31719, 1.29257, 1.27496, 1.26174, 1.25144};
  double KfMA40MeV[nE0] = {1.67437, 1.53055, 1.45506, 1.40683, 1.37324, 1.34847, 1.32944, 1.31436};
  double KfMA50MeV[nE0] = {1.77288, 1.61665, 1.53582, 1.48201, 1.44284, 1.4129, 1.38925, 1.37014};
  double KfMA60MeV[nE0] = {1.85602, 1.68401, 1.59904, 1.54266, 1.50082, 1.46799, 1.4415, 1.4195};
  double KfMA70MeV[nE0] = {1.92991, 1.74007, 1.64994, 1.5915, 1.5483, 1.51409, 1.486, 1.46248};
  double KfMA80MeV[nE0] = {1.99782, 1.78938, 1.69307, 1.63207, 1.58764, 1.5526, 1.52374, 1.49944};
  double KfMA90MeV[nE0] = {2.0617, 1.83451, 1.73141, 1.66729, 1.6213, 1.58541, 1.55599, 1.53126};
  double KfMA100MeV[nE0] = {2.12278, 1.87693, 1.76675, 1.69911, 1.65122, 1.61423, 1.58415, 1.55903};
  double KfMA200MeV[nE0] = {2.06543, 1.99617, 1.97534, 1.96104, 1.88813, 1.83371, 1.79096, 1.75693};
  double KfMA300MeV[nE0] = {2.11679, 1.99508, 1.94777, 1.92431, 1.91109, 1.90242, 1.89753, 1.8933};
  double KfMA400MeV[nE0] = {2.20752, 2.03888, 1.96707, 1.92846, 1.90484, 1.8891, 1.87848, 1.8697};
  double KfMA500MeV[nE0] = {2.31917, 2.10361, 2.00917, 1.95564, 1.92204, 1.89876, 1.88227, 1.86935};
  double KfMA600MeV[nE0] = {2.44447, 2.18024, 2.06445, 1.99649, 1.95258, 1.92186, 1.89907, 1.88231};
  double KfMA700MeV[nE0] = {2.5793, 2.26485, 2.1282, 2.04634, 1.9922, 1.9539, 1.92531, 1.90398};
  double KfMA800MeV[nE0] = {2.72069, 2.35541, 2.19782, 2.10255, 2.03836, 1.9924, 1.95799, 1.93186};
  double KfMA900MeV[nE0] = {2.86596, 2.45069, 2.27181, 2.16346, 2.08944, 2.03584, 1.99582, 1.96446};
  double KfMA1000MeV[nE0] = {3.01242, 2.54971, 2.34923, 2.22797, 2.14436, 2.08319, 2.03674, 2.00084};
      
  // initialize K-factor matrix
  for (int j = 0; j < nE0; j++) {
     KKffactor[0][j] = KfMA1MeV[j];
     KKffactor[1][j] = KfMA2MeV[j];
     KKffactor[2][j] = KfMA4MeV[j];
     KKffactor[3][j] = KfMA6MeV[j];
     KKffactor[4][j] = KfMA8MeV[j];
     KKffactor[5][j] = KfMA10MeV[j];
     KKffactor[6][j] = KfMA20MeV[j];
     KKffactor[7][j] = KfMA30MeV[j];
     KKffactor[8][j] = KfMA40MeV[j];
     KKffactor[9][j] = KfMA50MeV[j];
     KKffactor[10][j] = KfMA60MeV[j];
     KKffactor[11][j] = KfMA70MeV[j];
     KKffactor[12][j] = KfMA80MeV[j];
     KKffactor[13][j] = KfMA90MeV[j];
     KKffactor[14][j] = KfMA100MeV[j];
     KKffactor[15][j] = KfMA200MeV[j];
     KKffactor[16][j] = KfMA300MeV[j];
     KKffactor[17][j] = KfMA400MeV[j];
     KKffactor[18][j] = KfMA500MeV[j];
     KKffactor[19][j] = KfMA600MeV[j];
     KKffactor[20][j] = KfMA700MeV[j];
     KKffactor[21][j] = KfMA800MeV[j];
     KKffactor[22][j] = KfMA900MeV[j];
     KKffactor[23][j] = KfMA1000MeV[j];        
  }
  
  // bounds for energy E0
  if (E0test <= EE00[0])
     E0test = 1.005 * EE00[0]; //  Lower limit for energy
  
  if (E0test >= EE00[nE0-1] )
     E0test  = 0.995 * EE00[nE0-1]; // Upper limit for energy
  
  // bounds for mass of Dark state MA
  if (MAtest <= MMAA[0])
     MAtest = 1.005 * MMAA[0]; //  Lower limit for mass
 
  if (MAtest >= MMAA[nMA-1] )
     MAtest  = 0.995 * MMAA[nMA-1]; // Upper limit for mass
                               
  return BilinearInterpolation<nE0> (MAtest, E0test, MMAA, EE00, KKffactor, nMA, 1);
}

double parinv(double x, double a[], double f[], int n)
{
//
//    Interpolation at the point x. Function f(a) is tabulated
//    in arrays a, f with dimension n.
//
  int k1, k2, k3;

  if(n < 3) {std::cerr << "parinv: insufficient number of points" << std::endl; exit(1);}
  if(x < a[0]) {
    double c = fabs(x - a[0]);
    if(c < EPSPARINV*fabs(a[1]-a[0])) return a[0];
    k1 = 0;
  }
  else if(x > a[n-1]) {
    double c = fabs(x - a[n-1]);
    if(c < EPSPARINV*fabs(a[n-1]-a[n-2])) return a[n-1];
    k1 = n-3;
  }
  else {
    k1 = 0;
    k2 = n-1;
    k3 = k2 - k1;
    while(k3 > 1) {
      k3 = k1 + k3/2;
      if( a[k3]-x == 0 ) return f[k3];
      if( a[k3]-x < 0 ) k1 = k3;
      if( a[k3]-x > 0 ) k2 = k3;
      k3 = k2 - k1;
    }
    if(k2 == n-1) k1 = n - 3;
  }
  if(k1 < 0 || k1 > n-3) {std::cerr << "parinv: wrong index found" << std::endl; exit(1);}
  double b1 = a[k1];
  double b2 = a[k1+1];
  double b3 = a[k1+2];
  double b4 = f[k1];
  double b5 = f[k1+1];
  double b6 = f[k1+2];
  return b4 * ((x-b2)*(x-b3))/((b1-b2)*(b1-b3)) +
         b5 * ((x-b1)*(x-b3))/((b2-b1)*(b2-b3)) +
         b6 * ((x-b1)*(x-b2))/((b3-b1)*(b3-b2));
}


DarkPhotons::DarkPhotons(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                         double epsilIn, std::string fname)
:MA(MAIn), EThresh(EThreshIn), SigmaNorm(SigmaNormIn),
ANucl(ANuclIn), ZNucl(ZNuclIn), Density(DensityIn), epsilBench(1), epsil(epsilIn),
AccumulatedProbability(0.)
{
   nptable = NPTAB;
   double epi[NPTAB]={MA+Mmu, MA+Mmu+.1,MA+Mmu+2., MA+Mmu+10., MA+Mmu+20., MA+Mmu+50., MA+Mmu+100., MA+Mmu+150., MA+Mmu+250., MA+Mmu+500., MA+Mmu+800., MA+Mmu+1500., MA+Mmu+2000., MA+Mmu+2500., MA+Mmu+4000.};
   for(int ip=0; ip < nptable; ip++) {ep[ip] = epi[ip];}
   if(fname.find("lhe")!=std::string::npos){ParseLHE(fname);}
   else if(fname.find("root")!=std::string::npos){ParseROOT(fname);}
}


DarkPhotons::~DarkPhotons()
{
}

double DsigmaDx(double x, void * pp) 
{
   ParamsForChi* params = (ParamsForChi*)pp;
   
   double beta = sqrt(1 - (params->MMA)*(params->MMA)/(params->EE0)/(params->EE0));
   double num = 1.-x+x*x/3.;
   double denom = (params->MMA)*(params->MMA)*(1.-x)/x+Mel*Mel*x;
   double DsDx = beta*num/denom;

   return DsDx;
}

void DarkPhotons::ParseLHE(std::string fname)
{
   std::ifstream ifile;
   ifile.open(fname.c_str());
   if (!ifile)
   {
      std::cout << "Unable to open LHE file\n";
      exit(1);
   }   
   std::string line;
   int n=0;
   while(std::getline(ifile, line))
   {
      std::istringstream iss(line);
      int ptype, state;
      double skip, px, py, pz, E, pt, efrac, M;
      if (iss >> ptype >> state >> skip >> skip >> skip >> skip >> px >> py >> pz >> E >> M ) 
      {
         if((ptype==11)&&(state==1))
	 {
           double ebeam = E;
           double e_px, e_py, e_pz, a_px, a_py,  a_pz, e_E, a_E, e_M, a_M;
            if(mgdata.count(ebeam) == 0) {mgdata[ebeam];}
            for(int i=0;i<2;i++) {std::getline(ifile,line);}
            std::istringstream jss(line);
            jss >> ptype >> state >> skip >> skip >> skip >> skip >> e_px >> e_py >> e_pz >> e_E >> e_M;
            if((ptype==11)&&(state==1))
            {
               for(int i=0;i<2;i++) {std::getline(ifile,line);}
               std::istringstream kss(line);
               kss >> ptype >> state >> skip >> skip >> skip >> skip >> a_px >> a_py >> a_pz >> a_E >> a_M;
               if((ptype==622)&&(state==1))
               {
                  frame evnt;
                  double cmpx = a_px+e_px;
                  double cmpy = a_py+e_py;
                  double cmpz = a_pz+e_pz;
                  double cmE = a_E+e_E;
                  evnt.fEl = new TLorentzVector(e_px,e_py,e_pz,e_E);
                  evnt.cm = new TLorentzVector(cmpx,cmpy,cmpz,cmE);
                  evnt.E = ebeam;
                  mgdata[ebeam].push_back(evnt);
                  n++;
               }
            }
         }
      }
   }
   printf("Number of events: %i.\n", n);
   ifile.close();
   MakePlaceholders();
}

void DarkPhotons::ParseROOT(std::string fname)
{
   TFile *f = TFile::Open(fname.c_str());
   TTree *tree = (TTree*)f->Get("Events");
   TLorentzVector* mvec = new TLorentzVector();
   TLorentzVector* avec = new TLorentzVector();
   TLorentzVector* nvec = new TLorentzVector();
   tree->SetBranchAddress("IncidentParticle",&mvec);
   tree->SetBranchAddress("APrime",&avec);
   tree->SetBranchAddress("Nucleus",&nvec);
   int entries = tree->GetEntries();
   for(int i=0; i<entries; i++)
   {
      if(i<entries){tree->GetEntry(i);}
      else{tree->GetEntry(i-entries);}
      frame evnt;
      evnt.fEl = (TLorentzVector*)mvec->Clone();
      evnt.cm = (TLorentzVector*)avec->Clone();
      *evnt.cm = *evnt.cm+*evnt.fEl;
      TLorentzVector* ebeam = (TLorentzVector*)nvec->Clone();
      *ebeam = *ebeam+*evnt.cm;
      evnt.E = round(ebeam->Z()*10.)/10.;
      if(mgdata.count(evnt.E)==0){mgdata[evnt.E];}
      mgdata[evnt.E].push_back(evnt);
   }
   f->Close();
   MakePlaceholders();
}

void DarkPhotons::MakePlaceholders()
{
   for ( const auto &iter : mgdata )
   {                                                
      energies.push_back(std::make_pair(iter.first,iter.second.size()));
   } 

   for(uint64_t i=0;i<energies.size();i++)
   {  
      energies[i].second=int(drand48()*mgdata[energies[i].first].size());
   }
}

double chi (double t, void * pp) 
{
  ParamsForChi* params = (ParamsForChi*)pp;

/* Reminder II:
   params->AA;
   params->ZZ;
   params->MMA;
   params->EE0;
*/

  double d = 0.164/pow((params->AA),2./3.);
  double ap = 773.0/Mel/pow((params->ZZ),2./3.);
  double a = 111.0/Mel/pow((params->ZZ),1./3.);
  double G2el = (params->ZZ)*(params->ZZ)*a*a*a*a*t*t/(1.0+a*a*t)/(1.0+a*a*t)/(1.0+t/d)/(1.0+t/d);
  double G2in = (params->ZZ)*ap*ap*ap*ap*t*t/(1.0+ap*ap*t)/(1.0+ap*ap*t)/(1.0+t/0.71)/(1.0+t/0.71)
    /(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)
    *(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr)*(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr);
  double G2 = G2el+G2in;
  double ttmin = (params->MMA)*(params->MMA)*(params->MMA)*(params->MMA)/4.0/(params->EE0)/(params->EE0);
  //double ttmin = lowerLimit(x,theta,p);
  double Under = G2*(t-ttmin)/t/t;
//  std::cout << "Under: " << Under << " d: " << d << " AA " << params->AA << " a: " << a << std::endl;
  return Under;
}

double chi_inel (double t, void * pp) 
{
  ParamsForChi* params = (ParamsForChi*)pp;

/* Reminder II:
   params->AA;
   params->ZZ;
   params->MMA;
   params->EE0;
*/

  double d = 0.164/pow((params->AA),2./3.);
  double ap = 773.0/Mel/pow((params->ZZ),2./3.);
  double a = 111.0/Mel/pow((params->ZZ),1./3.);
  double G2el = (params->ZZ)*(params->ZZ)*a*a*a*a*t*t/(1.0+a*a*t)/(1.0+a*a*t)/(1.0+t/d)/(1.0+t/d);
//  double G2in = (params->ZZ)*ap*ap*ap*ap*t*t/(1.0+ap*ap*t)/(1.0+ap*ap*t)/(1.0+t/0.71)/(1.0+t/0.71)
//    /(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)
//    *(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr)*(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr);
  double G2in=0;
  double G2 = G2el+G2in;
  double ttmin = (params->MMA)*(params->MMA)*(params->MMA)*(params->MMA)/4.0/(params->EE0)/(params->EE0);
  //double ttmin = lowerLimit(x,theta,p);
  double Under = G2*(t-ttmin)/t/t;
//  std::cout << "Under: " << Under << " d: " << d << " AA " << params->AA << " a: " << a << std::endl;
  return Under;
}


double DarkPhotons::TotalCrossSectionCalc(double E0, bool KFactor)
{
  double Xmax;
  double sigmaTot;

    if(E0 < 2.*MA) return 0.;

    //begin: chi-formfactor calculation

    gsl_integration_workspace * w
      = gsl_integration_workspace_alloc (1000);

    double result, error;
    double tmin = MA*MA*MA*MA/(4.*E0*E0);
    double tmax = MA*MA;

    gsl_function F;
    ParamsForChi alpha = {1.0, 1.0, 1.0, 1.0};
    if(!KFactor){F.function = &chi;}
    else{F.function = &chi_inel;}
    F.params = &alpha;

    alpha.AA = ANucl;
    alpha.ZZ = ZNucl;
    alpha.MMA = MA;
    alpha.EE0 = E0;

    gsl_integration_qags (&F, tmin, tmax, 0, 1e-7, 1000,
                          w, &result, &error);

    //printf ("chi/Z^2 = % .18f\n", result/(ZNucl*ZNucl));
    //printf ("result    = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);

    double ChiRes = result;
//    std::cout << "Chi: " << result << " E0: " << E0 << " MA: " << MA << std::endl;
    if(!KFactor)
    {
       gsl_integration_workspace_free (w);
   
       gsl_integration_workspace * s 
          = gsl_integration_workspace_alloc (1000);
       gsl_function G;
       G.function = &DsigmaDx;
       G.params = &alpha;
       double xmin = 0;
       double xmax = 1;
       if((Mel/E0)>(MA/E0)) xmax = 1-Mel/E0;
       else xmax = 1-MA/E0;
       double res, err;
   
       gsl_integration_qags (&G, xmin, xmax, 0, 1e-7, 1000,
                             s, &res, &err);
       double DsDx = res;
   
       gsl_integration_workspace_free(s);
       
   //    end: chi-formfactor calculation
   
       sigmaTot= GeVtoPb*4.*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*ChiRes*DsDx;
    }
    else
    {
       double beta = sqrt(1. - MA*MA/(E0*E0));
       double cutoff1 = Mel/MA;
       double cutoff2 = MA/E0;
       double cutoff = cutoff2;
       if(cutoff1 > cutoff2) cutoff = cutoff1;
       sigmaTot= GeVtoPb*(4./3.)*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*ChiRes*beta*log(1./(cutoff*cutoff))/(MA*MA);
       double KFactor = KfactorApproximate(MA,E0);
       KFactor /= 1.049;
       sigmaTot /= KFactor;
    }

    if(sigmaTot < 0.) sigmaTot=0.;

    return sigmaTot;
}

double DsigmaDxmu(double x, void * pp)
{
   ParamsForChi* params = (ParamsForChi*)pp;

   double MMu = 105.658/1000.;
   double beta = sqrt(1- (params->MMA)*(params->MMA)/(params->EE0)/(params->EE0));
   double num = 1.-x+x*x/3.;
   double denom = (params->MMA)*(params->MMA)*(1.-x)/x+MMu*MMu*x;
   double DsDx = beta*num/denom;

   return DsDx;
}

double chimu (double t, void * pp) {
  ParamsForChi* params = (ParamsForChi*)pp;

/* Reminder II:
   params->AA;
   params->ZZ;
   params->MMA;
   params->EE0;
*/
//Originally used Mmu instead of Mel, check for differences and verify that Mel was used.
  double d = 0.164/pow((params->AA),2./3.);
  double ap = 773.0/Mel/pow((params->ZZ),2./3.);
  double a = 111.0/Mel/pow((params->ZZ),1./3.);
  double G2el = (params->ZZ)*(params->ZZ)*a*a*a*a*t*t/(1.0+a*a*t)/(1.0+a*a*t)/(1.0+t/d)/(1.0+t/d);
  double G2in = (params->ZZ)*ap*ap*ap*ap*t*t/(1.0+ap*ap*t)/(1.0+ap*ap*t)/(1.0+t/0.71)/(1.0+t/0.71)
    /(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)
    *(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr)*(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr);
  double G2 = G2el+G2in;
  double ttmin = (params->MMA)*(params->MMA)*(params->MMA)*(params->MMA)/4.0/(params->EE0)/(params->EE0);
  //double ttmin = lowerLimit(x,theta,p);
  double Under = G2*(t-ttmin)/t/t;
//  std::cout << "Under: " << Under << " MMA: " << params->MMA << " EE0: " << params->EE0 << std::endl;
  return Under;
}

double DarkPhotons::TotalMuCrossSectionCalc(double E0)
{
  double Xmin;
  double Xmax;
  double sigmaTot;

    if(E0 < 2.*MA) return 0.;

    Xmin = MA/E0;
    Xmax = 1.0-Xmin;

    //begin: chi-formfactor calculation
    gsl_integration_workspace * w
      = gsl_integration_workspace_alloc (1000);

    double result, error;
    double tmin = MA*MA*MA*MA/(4.*E0*E0);
    double tmax = MA*MA;

    gsl_function F;
    ParamsForChi alpha = {1.0, 1.0, 1.0, 1.0};
    F.function = &chimu;
    F.params = &alpha;

    alpha.AA = ANucl;
    alpha.ZZ = ZNucl;
    alpha.MMA = MA;
    alpha.EE0 = E0;

    gsl_integration_qags (&F, tmin, tmax, 0, 1e-7, 1000,
                          w, &result, &error);

    //printf ("chi/Z^2 = % .18f\n", result/(ZNucl*ZNucl));
    //printf ("result    = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);

   double ChiRes = result;
   gsl_integration_workspace_free (w);
//    end: chi-formfactor calculation

   gsl_integration_workspace * dxspace = gsl_integration_workspace_alloc (1000);
   gsl_function G;
   G.function = &DsigmaDxmu;
   G.params = &alpha;
   double xmin = 0;
   double xmax = 1;
   if((Mmu/E0)>(MA/E0)) xmax = 1-Mmu/E0;
   else xmax = 1-MA/E0;
   double res, err;

   gsl_integration_qags (&G, xmin, xmax, 0, 1e-7, 1000, dxspace, &res, &err);
   
   double DsDx = res;
   gsl_integration_workspace_free(dxspace);

   sigmaTot = GeVtoPb*4.*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*ChiRes*DsDx;
   if(sigmaTot < 0.)
   {
      sigmaTot = 0.;
   }

   return sigmaTot;
}
double DarkPhotons::MaxCrossSectionCalc(double E0)
{
  double Xmax, Xmin;
  double sigmaMax;

    if(E0 < 2.*MA) return 0.;

    Xmin = MA/E0;
    Xmax = 0.998;

    double UxthetaMax = MA*MA*(1. - Xmax)/Xmax + Mel*Mel*Xmax;
    double AAMax = (1. - Xmax + Xmax*Xmax/2.0) / (UxthetaMax*UxthetaMax);
    double BBMax = (1. - Xmax)*(1. - Xmax)*MA*MA / (UxthetaMax*UxthetaMax*UxthetaMax*UxthetaMax);
    double CCMax = MA*MA - UxthetaMax*Xmax/(1. - Xmax);
    sigmaMax = Xmax * (AAMax + BBMax*CCMax);

    return sigmaMax;
}


void DarkPhotons::PrepareTable()
{
  for(int ip=0; ip < nptable; ip++) {
    sigmap[ip] = TotalMuCrossSectionCalc(ep[ip]);
    sigmax[ip] = MaxCrossSectionCalc(ep[ip]);
  }
}


double DarkPhotons::GetsigmaTot(double E0)
{
  if(E0<(MA+Mmu))
  {
     return 0;
  }
  double st = parinv(E0, ep, sigmap, nptable);
  if(st<0) {return TotalMuCrossSectionCalc(E0);}
  return parinv(E0, ep, sigmap, nptable);
}


double DarkPhotons::GetsigmaMax(double E0)
{
  return parinv(E0, ep, sigmax,	nptable);
}


bool DarkPhotons::Emission(double E0, double DensityMat, double StepLength)
{
  if(E0 < EThresh) return false;
  if(fabs(DensityMat - Density) > 0.1) return false;
  double prob = SigmaNorm*GetsigmaTot(E0)*StepLength;
  AccumulatedProbability += prob;
  if(drand48() < prob) return true;
  return false;
}

TLorentzVector* DarkPhotons::SimulateEmission(double E0, std::string type)
{
   frame data = GetMadgraphData(E0);
   double EAcc, Pt, P, PhiAcc;
   TLorentzVector* fel;
   if(type == "forward_only")
   {
      EAcc = (data.fEl->E()-Mel)/(data.E-Mel-MA)*(E0-Mel-MA);
      Pt = data.fEl->Pt();
      P = sqrt(EAcc*EAcc-Mel*Mel);
      PhiAcc = data.fEl->Phi();
      int i = 0;
      while(Pt*Pt+Mel*Mel>EAcc*EAcc) //Skip events until the Pt is less than the energy.
      {
         i++;
         data = GetMadgraphData(E0);
         EAcc = (data.fEl->E()-Mel)/(data.E-Mel-MA)*(E0-Mel-MA);
         Pt = data.fEl->Pt();
         P = sqrt(EAcc*EAcc-Mel*Mel);
         PhiAcc = data.fEl->Phi();
         if(i>10000)
         {
            printf("Did not manage to simulate. E0 = %e, EAcc = %e\n", E0, EAcc);
         }  
      }
   }
   else if(type == "cm_scaling")
   {
      TLorentzVector* el = new TLorentzVector(data.fEl->X(),data.fEl->Y(),data.fEl->Z(),data.fEl->E());
      double ediff = data.E-E0;
      TLorentzVector* newcm = new TLorentzVector(data.cm->X(),data.cm->Y(),data.cm->Z()-ediff,data.cm->E()-ediff);
      el->Boost(-1.*data.cm->BoostVector());
      el->Boost(newcm->BoostVector());
      double newE = (data.fEl->E()-Mel)/(data.E-Mel-MA)*(E0-Mel-MA);
      el->SetE(newE);
      fel = el;
      EAcc = el->E();
      Pt = el->Pt();
      P = el->P();
   }
   TLorentzVector* fParticle = new TLorentzVector;
   if (type == "forward_only")
   {
      double Theta = asin(Pt/P);
      double Eta = -log(tan(Theta/2.));
      double Phi = drand48()*2.*3.14159;
      fParticle->SetPtEtaPhiE(Pt,Eta,Phi,EAcc);
   }
   else if(type == "cm_scaling")
   {
      fParticle = fel;
   }
   return fParticle;
}

frame DarkPhotons::GetMadgraphData(double E0)
{
   double samplingE = energies[0].first;
   frame cmdata;
   int i=-1;
   bool pass = false;
   while(!pass)
   {
      i++;
      samplingE = energies[i].first;
      if(E0<=samplingE||i>=energies.size()){pass=true;}
   }
   if(i==energies.size()) {i=i-1;}
   if(energies[i].second>=double(mgdata[energies[i].first].size())) {energies[i].second = 0;}
   cmdata = mgdata[energies[i].first].at(energies[i].second);
   energies[i].second=energies[i].second+1;

   return cmdata;
}

TLorentzVector* DarkPhotons::MuSimulateEmission(double E0, std::string type)
{
   frame data = GetMadgraphData(E0);
   double EAcc, Pt, P, PhiAcc;
   TLorentzVector* fel; 
   if(type == "forward_only")
   {
      EAcc = (data.fEl->E()-Mmu)/(data.E-Mmu-MA)*(E0-Mmu-MA);
      EAcc = Mmu+EAcc;
      Pt = data.fEl->Pt();
      P = sqrt(EAcc*EAcc-Mmu*Mmu);
      PhiAcc = data.fEl->Phi();
      int i = 0;
      while(Pt*Pt+Mmu*Mmu>EAcc*EAcc) //Skip events until the Pt is less than the energy.
      {
         i++; 
         data = GetMadgraphData(E0);
         EAcc = (data.fEl->E()-Mmu)/(data.E-Mmu-MA)*(E0-Mmu-MA);
         EAcc = Mmu+EAcc;
         Pt = data.fEl->Pt();
         P = sqrt(EAcc*EAcc-Mmu*Mmu);
         PhiAcc = data.fEl->Phi();

         if(i>10000)
         {
            printf("Did not manage to simulate. E0 = %e, EAcc = %e\n", E0, EAcc);
            printf("Pt is %e, Data E is %e, Mmu is %e, MA is %e\n", data.fEl->Pt(),data.E,Mmu,MA);
         }
      }
   }
   else if(type == "cm_scaling")
   {
      TLorentzVector* el = new TLorentzVector(data.fEl->X(),data.fEl->Y(),data.fEl->Z(),data.fEl->E());
      double ediff = data.E-E0;
      TLorentzVector* newcm = new TLorentzVector(data.cm->X(),data.cm->Y(),data.cm->Z()-ediff,data.cm->E()-ediff);
      el->Boost(-1.*data.cm->BoostVector());
      el->Boost(newcm->BoostVector());
      double newE = (data.fEl->E()-Mmu)/(data.E-Mmu-MA)*(E0-Mmu-MA);
      el->SetE(newE); 
      fel=el;
      EAcc = el->E();
      Pt = el->Pt();
      P = el->P();
      PhiAcc=el->Phi();
   }
   else
   {
      std::cout << "Method not recognized. Skipping Event.\n";
      EAcc = E0;
   }
   TLorentzVector* fParticle = new TLorentzVector;
   if (type == "forward_only")
   {
      double Theta = asin(Pt/P);
      double Eta = -log(tan(Theta/2.));
      double Phi = drand48()*2.*3.14159;
      fParticle->SetPtEtaPhiE(Pt,Eta,Phi,EAcc);
   }
   else if(type == "cm_scaling")
   {
      fParticle = fel;
   }
   return fParticle;
}
