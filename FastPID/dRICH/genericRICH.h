/// \author R+Preghenella
/// \email  preghenella@bo.infn.it
/// \date   March 2020

#ifndef __GENERICRICH_H__
#define __GENERICRICH_H__
	
#include "genericDetector.h"

class genericRICH : public genericDetector
{
 public:
  genericRICH() = default;
  virtual ~genericRICH() = default;

  /** setters **/
  void setIndex(double val) { mIndex = val; };
  void setEfficiency(double val) { mEfficiency = val; };
  void setMinPhotons(double val) { mMinPhotons = val; };
  void setThresholdMode(bool val) { mThresholdMode = val; };
  
  void setChromaticSigma(int n, double *valx, double *valy) {
    if (mChromaticSigma) delete mChromaticSigma;
    mChromaticSigma = new TGraph(n, valx, valy);
  }
  void setPositionSigma(int n, double *valx, double *valy) {
    if (mPositionSigma) delete mPositionSigma;
    mPositionSigma = new TGraph(n, valx, valy);
  }
  void setEmissionSigma(int n, double *valx, double *valy) {
    if (mEmissionSigma) delete mEmissionSigma;
    mEmissionSigma = new TGraph(n, valx, valy);
  }
  void setFieldSigma(int n, double *valx, double *valy) {
    if (mFieldSigma) delete mFieldSigma;
    mFieldSigma = new TGraph(n, valx, valy);
  }
  void setTrackingSigma(int n, double *valx, double *valy) {
    if (mTrackingSigma) delete mTrackingSigma;
    mTrackingSigma = new TGraph(n, valx, valy);
  }
  
  /** methods to override **/
  double numSigma (double eta, double p, PID::type PID) override;
  double maxP (double eta, double nsigma, PID::type PID) override;
  double minP (double eta, double nsigma, PID::type PID) override;

  double cherenkovAngle(double p, double m) const { return acos( sqrt( m * m + p * p ) / ( mIndex * p ) ); };
  double cherenkovThreshold(double m) const { return m / sqrt(mIndex * mIndex - 1.); };
  double numberOfPhotons(double angle) const { return 490. * sin(angle) * sin(angle) * mLength; };
  double numberOfDetectedPhotons(double angle) const { return numberOfPhotons(angle) * mEfficiency; };
  double cherenkovAngleSigma(double eta, double p, double m) const;
  
 protected:
  
  // RICH parameters
  double mIndex = 1.0014;    // refractive index
  double mEfficiency = 0.25; // overall photon detection efficiency 
  double mMinPhotons = 3.;   // minimum number of detected photons

  // contributions to resolution
  TGraph *mChromaticSigma = nullptr; // chromatic resolution vs. polar angle [rad]
  TGraph *mPositionSigma = nullptr; // position resolution vs. polar angle [rad]
  TGraph *mEmissionSigma = nullptr; // emission resolution vs. polar angle [rad]
  TGraph *mFieldSigma = nullptr; // field resolution vs. polar angle [rad]
  TGraph *mTrackingSigma = nullptr; // tracking resolution vs. polar angle [rad]

  // threshold mode
  bool mThresholdMode = true;
  
};
	
double genericRICH::cherenkovAngleSigma(double eta, double p, double m) const
{

  auto theta = 2. * atan(exp(-eta)) * 57.295780;
  auto chromatic = mChromaticSigma ? mChromaticSigma->Eval(theta) : 0.;
  auto position = mPositionSigma ? mPositionSigma->Eval(theta) : 0.;
  auto emission = mEmissionSigma ? mEmissionSigma->Eval(theta) : 0.;
  auto field = mFieldSigma ? mFieldSigma->Eval(theta) : 0.;
  auto tracking = mTrackingSigma ? mTrackingSigma->Eval(theta) : 0.;
  
  // contributions that scale with number of detected photons
  auto ndet = numberOfDetectedPhotons(cherenkovAngle(p, m));
  auto sigma1 = sqrt(chromatic * chromatic +
		     position * position +
		     emission * emission +
		     field * field +
		     tracking * tracking);
  // contributions that do not
  auto sigma2 = 0.;
  // 
  return sqrt(sigma1 * sigma1 / ndet + sigma2 * sigma2);
};


double genericRICH::numSigma(double eta, double p, PID::type PID)
{
  double mass1, mass2;
  switch (PID) {
  case pi_k:
    mass1 = mMassPion;
    mass2 = mMassKaon;
    break;
  case k_p:
    mass1 = mMassKaon;
    mass2 = mMassProton;
    break;
  }

  double thr1 = cherenkovThreshold(mass1);
  double thr2 = cherenkovThreshold(mass2);

  /** both particles are above threshold **/
  if (p > thr1 && p > thr2)
    return (cherenkovAngle(p, mass1) - cherenkovAngle(p, mass2)) / cherenkovAngleSigma(eta, p, mass1);
    
  /** lightest particle above threshold **/
  if (mThresholdMode && p > thr1)
    return (cherenkovAngle(thr2 + 0.001, mass1) - cherenkovAngle(thr2 + 0.001, mass2)) / cherenkovAngleSigma(eta, thr2 + 0.001, mass1);
    
  /** none above threshold **/
  return 0.;
}

double genericRICH::maxP(double eta, double nsigma, PID::type PID)
{
  double mass1, mass2;
  switch (PID) {
  case pi_k:
    mass1 = mMassPion;
    mass2 = mMassKaon;
    break;
  case k_p:
    mass1 = mMassKaon;
    mass2 = mMassProton;
    break;
  }

  /** let's do it numerically, starting from the threshold **/
  double p = minP(eta, nsigma, PID) + 0.001;
  while (numSigma(eta, p, PID) > nsigma) p += 0.001;
  return p;
}

double genericRICH::minP(double eta, double nsigma, PID::type PID)
{
  double mass;
  switch (PID) {
  case pi_k:
    mass = mThresholdMode ? mMassPion : mMassKaon;
    break;
  case k_p:
    mass = mThresholdMode ? mMassKaon : mMassProton;
    break;
  }

  /** let's do it numerically, starting from the threshold **/
  double p = cherenkovThreshold(mass);
  while (numberOfPhotons(cherenkovAngle(p, mass)) * mEfficiency < mMinPhotons) p += 0.001;
  return std::max(p, pMin(eta));
}

#endif /* __GENERICRICH_H__ */
