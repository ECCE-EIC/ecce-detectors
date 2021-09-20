// $Id: $

/*!
 * \file ECCEmRICHFastPIDMap.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef ECCEmRICHFastPIDMap_H_
#define ECCEmRICHFastPIDMap_H_

#include "ECCEFastPIDMap.h"

#include <string>

class TH2F;
class TF1;

/*!
 * \brief ECCEmRICHFastPIDMap
 * Import from ./mRICH/mRICH
 */
class ECCEmRICHFastPIDMap : public ECCEFastPIDMap
{
 public:

  //   Detectors.push_back( new mRICH(0.00175, 1, 3, mom) ); // 20 psec @ 100 cm
  ECCEmRICHFastPIDMap(double trackResolution = 0.00175, double timePrecision = 1.0, double pixS = 3);
  virtual ~ECCEmRICHFastPIDMap();

  PIDCandidate_LogLikelihood_map getFastSmearLogLikelihood(int truth_pid, const double momentum, const double theta) const override;

  void setThetaAcceptanceMinMax(double min, double max)
  {
    m_acceptanceThetaMin = min;
    m_acceptanceThetaMax = max;
  }

  /// Sets the verbosity of this module (0 by default=quiet).
  virtual void Verbosity(const int ival) { m_Verbosity = ival; }

  /// Gets the verbosity of this module.
  virtual int Verbosity() const { return m_Verbosity; }

 private:
  double m_acceptanceThetaMin = 2.658;  // eta = -1.4;
  double m_acceptanceThetaMax = 3.04;   // eta = -3

  double getAng(double mass, double mom) const;
  double getdAng(double mass, double mom) const;
  double getdAngTrk(double mass, double mom) const;
  double getNgamma(double t, double mass, double mom) const;
  double T_Aer(double t, double lam) const;
  double T_QE(double lam) const;

  // Physical constants (should come from elsewhere!)
  double mPion;    // GeV/c^2
  double mKaon;    // GeV/c^2
  double mProton;  // GeV/c^2
  double c;        // cm/picosecond;
  double n;
  double f;  //mm
  double a;  //mm
  double N_gam;
  double pi;
  double alpha;
  double L;
  double th0;

  double fTrackResolution;
  double fTimePrecision;

  double pLow;
  double pHigh;

  /// The verbosity level. 0 means not verbose at all.
  int m_Verbosity = 0;
};

#endif /* ECCEmRICHFastPIDMap_H_ */
