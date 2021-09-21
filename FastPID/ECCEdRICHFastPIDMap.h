// $Id: $

/*!
 * \file ECCEdRICHFastPIDMap.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef ECCEdRICHFastPIDMap_H_
#define ECCEdRICHFastPIDMap_H_

#include "ECCEFastPIDMap.h"

#include <string>

#include <cmath>

class TH2F;
class TF1;

#include <TGraph.h>
/*!
 * \brief ECCEdRICHFastPIDMap
 * Import from dRICH/genericRICH
 */
class ECCEdRICHFastPIDMap : public ECCEFastPIDMap {
public:
  ECCEdRICHFastPIDMap();
  virtual ~ECCEdRICHFastPIDMap();

  //! set one of these two first
  void dualRICH_aerogel();
  void dualRICH_C2F6();

  PIDCandidate_LogLikelihood_map
  getFastSmearLogLikelihood(int truth_pid, const double momentum,
                            const double theta) const override;

  enum type { pi_k, k_p };

  enum EDetector_t { kBarrel, kForward };

  /** setters **/
  void setIndex(double val) { mIndex = val; };
  void setEfficiency(double val) { mEfficiency = val; };
  void setMinPhotons(double val) { mMinPhotons = val; };
  void setThresholdMode(bool val) { mThresholdMode = val; };

  void setChromaticSigma(int n, double *valx, double *valy) {
    if (mChromaticSigma)
      delete mChromaticSigma;
    mChromaticSigma = new TGraph(n, valx, valy);
  }
  void setPositionSigma(int n, double *valx, double *valy) {
    if (mPositionSigma)
      delete mPositionSigma;
    mPositionSigma = new TGraph(n, valx, valy);
  }
  void setEmissionSigma(int n, double *valx, double *valy) {
    if (mEmissionSigma)
      delete mEmissionSigma;
    mEmissionSigma = new TGraph(n, valx, valy);
  }
  void setFieldSigma(int n, double *valx, double *valy) {
    if (mFieldSigma)
      delete mFieldSigma;
    mFieldSigma = new TGraph(n, valx, valy);
  }
  void setTrackingSigma(int n, double *valx, double *valy) {
    if (mTrackingSigma)
      delete mTrackingSigma;
    mTrackingSigma = new TGraph(n, valx, valy);
  }

  /** methods to override **/
  double numSigma(double eta, double p, type PID) const;

  double cherenkovAngle(double p, double m) const {
    return acos(sqrt(m * m + p * p) / (mIndex * p));
  };
  double cherenkovThreshold(double m) const {
    return m / sqrt(mIndex * mIndex - 1.);
  };
  double numberOfPhotons(double angle) const {
    return 490. * sin(angle) * sin(angle) * mLength;
  };
  double numberOfDetectedPhotons(double angle) const {
    return numberOfPhotons(angle) * mEfficiency;
  };
  double cherenkovAngleSigma(double eta, double p, double m) const;

  double etaMin() const;
  double etaMax() const;
  /** setters **/
  void setType(EDetector_t val) { mType = val; };
  void setName(const std::string &val) { mName = val; };
  void setLength(double val) { mLength = val; };
  void setRadius(double val) { mRadius = val; };
  void setPositionZ(double val) { mPositionZ = val; };
  void setRadiusIn(double val) { mRadiusIn = val; };
  void setRadiusOut(double val) { mRadiusOut = val; };
  void setMagneticField(double val) { mMagneticField = val; };

protected:
  // RICH parameters
  double mIndex = 1.0014;    // refractive index
  double mEfficiency = 0.25; // overall photon detection efficiency
  double mMinPhotons = 3.;   // minimum number of detected photons

  // contributions to resolution
  TGraph *mChromaticSigma =
      nullptr; // chromatic resolution vs. polar angle [rad]
  TGraph *mPositionSigma = nullptr; // position resolution vs. polar angle [rad]
  TGraph *mEmissionSigma = nullptr; // emission resolution vs. polar angle [rad]
  TGraph *mFieldSigma = nullptr;    // field resolution vs. polar angle [rad]
  TGraph *mTrackingSigma = nullptr; // tracking resolution vs. polar angle [rad]

  // threshold mode
  bool mThresholdMode = true;

  bool initialized = false;
  std::string mName = "genericDetector";
  std::string mDescription = "Detector description";
  EDetector_t mType = kBarrel;
  double mLength = 200.;      // [cm]
  double mRadius = 200.;      // [cm]
  double mPositionZ = 200.;   // [cm]
  double mRadiusIn = 20.;     // [cm]
  double mRadiusOut = 200.;   // [cm]
  double mMagneticField = 2.; // [T]

  const double mLightSpeed = 29.9792458;      // speed of light [cm/ns]
  const double mMassElectron = 0.00051099891; // electron mass [GeV]
  const double mMassPion = 0.13957018;        // pion mass [GeV]
  const double mMassKaon = 0.493677;          // kaon mass [GeV]
  const double mMassProton = 0.93827208816;   // proton mass [GeV]
};

#endif /* ECCEdRICHFastPIDMap_H_ */
