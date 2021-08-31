#include "mRICH.h"

using namespace std;


mRICH::mRICH(double trackResolution, double timePrecision, double pix, double p)
{
	
  // In this arameterization we assume that the particle enters the mRICH perpendicular to its front face. More realistic case will be implemented in the subsequent updates.


  fTrackResolution = trackResolution;
  fTimePrecision = timePrecision;
  mom     = p;
  pLow    = 3;
  pHigh   = 10.1;
  c       = 0.0299792458;  // cm/picosecond
  n       = 1.03; //Aerogel
  a       = pix;// pixel size 3.0; // mm -- one side
  f       = 152.4; //focal length mm =  6"
  N_gam   = 10;
  mPion   = 0.13957018;    //GeV/c^2
  mKaon   = 0.493677;      //GeV/c^2
  mProton = 0.93827208816; //GeV/c^2
  pi      = 3.14159;
  alpha   = 0.0072973525693; // hyperfine const
  L       = 3.0;//Aerogel block thickness in cm

  //===============
  th0 = 0.; //incidence angle in radians

  //Angle difference 
  double dth = getAng(mPion) - getAng(mKaon);
  //Detector uncertainty 
  double sigTh = sqrt(pow(getdAng(mPion),2)+pow(getdAng(mKaon),2));
  //Global uncertainty
  double sigThTrk = sqrt(pow(getdAngTrk(mPion),2)+pow(getdAngTrk(mKaon),2));
  double sigThc = sqrt(pow(sigTh/sqrt(getNgamma(L,mKaon)),2)+pow(sigThTrk,2));
  Nsigma = dth/sigThc;

#if 0
  //From simulations
  ReadMap("mRICHresol.root");
  Nsigma = fpixMap->GetBinContent(p,pix);
#endif

  //cout<<getNgamma(mPion)<<"\t"<<getNgamma(mKaon)<<endl;

  char nameit[500];
  sprintf(nameit,"mRICH Pixel Size =%dx%d (mm^{2}) / p =%.f (GeV/c)",pix,pix,p);
  myName = nameit;

  //cout<<".................\t"<<Nsigma<<endl;
}

double mRICH::maxP(double eta, double p, double numSigma, PID::type PID)
{
  if (valid(eta,p))
    {
      return Nsigma;
    }
  else
    {
      cout << "mRICH.C:  Invalid (eta) for this detector." <<endl;
      return 0.0;
    }
}

//Angle exiting the Aerogel
double mRICH::getAng(double mass)
{

  double beta = mom/sqrt(pow(mom,2)+pow(mass,2));
  double thc = acos(1./n/beta);
  double th0p = asin(sin(th0)/n);
  double dthc = thc - th0p;
  double theta = asin(n*sin(dthc));

  return theta;
}

//Uncertainty due to detector effects
double mRICH::getdAng(double mass)
{

  double beta = mom/sqrt(pow(mom,2)+pow(mass,2));
  double thc = acos(1./n/beta);
  double th0p = asin(sin(th0)/n);
  double dthc = thc - th0p;
  double theta = asin(n*sin(dthc));

  double sig_ep = 0;//Emission point error
  double sig_chro = 0; //Chromatic dispersion error
  double sig_pix = a*pow(cos(theta),2)/f/sqrt(12.);;

  double sigTh = sqrt(pow(sig_ep,2)+pow(sig_chro,2)+pow(sig_pix,2));

  return sigTh;
}

//Uncertainty due to tracking resolution
double mRICH::getdAngTrk(double mass)
{

  double beta = mom/sqrt(pow(mom,2)+pow(mass,2));
  double thc = acos(1./n/beta);
  double th0p = asin(sin(th0)/n);
  double dthc = thc - th0p;
  double theta = asin(n*sin(dthc));

  double sig_trk = (cos(dthc)/cos(theta))*(cos(th0)/cos(th0p))*fTrackResolution;;

  return sig_trk;
}

//no. of gamms
double mRICH::getNgamma(double t, double mass)
{
  int tot = 10000;
  double beta = mom/sqrt(pow(mom,2)+pow(mass,2));
  double fact = 2.*pi*alpha*t*(1.-1./pow(n*beta,2));
  double T_lensWin = 0.92*0.92;
  double xmin = 300.e-7;
  double xmax = 650.e-7;
  double dx = (xmax-xmin)/tot;
  double sum = 0;
  for(int j=0; j<tot; j++){
    double x = xmin + j*dx + dx/2.;
    sum += T_QE(x)*T_Aer(t,x)/pow(x,2);
  }
  return fact*T_lensWin*sum*dx;
}

//Quantum efficiency
double mRICH::T_QE(double lam)
{
  return 0.34*exp(-1.*pow(lam-345.e-7,2)/(2.*pow(119.e-7,2)));
}

//Transmissions of the radiator block
double mRICH::T_Aer(double t, double lam)
{
  return 0.83*exp(-1.*t*56.29e-20/pow(lam,4));
}


//Read histogram from MC parametrization
void  mRICH::ReadMap(TString name){
  TFile* file = TFile::Open(name);
  fpixMap = (TH2F *) file->Get("h2pix");
}


void mRICH::description()
{
  //  Here one should describe what this detector is and what assumptions are 
  //  made in calculating the performance specifications.

  cout << "My name is \"" << myName << "\" and I am described as follows:" <<endl;
  cout << "    Momentum coverage for k/pi separation=  [" << 3.0 << " to " << 10.0 <<" (GeV/c) ]"<<endl;
  cout << "    Since mRICH detector system consists of array of identical shoebox-like modules, in principle, " << endl;
  cout << "    the mRICH pid performance is invariant relative to the tracking polar angle but the momentum of the " << endl;
  cout << "    charge particle and the entrance point location on the front face of a given mRICH module." << endl;
  cout << "    Assumed time precision = " << fTimePrecision << " ns" <<endl;
  cout << "    Assumed track resolution = "<< fTrackResolution << " mrad" <<endl;
  //cout << "    Assumed quantum efficiency of the photosensors = "<< fSensorQMefficiency << "%" <<endl;
  cout << endl;

}
