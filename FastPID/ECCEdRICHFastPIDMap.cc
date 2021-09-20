// $Id: $

/*!
 * \file ECCEdRICHFastPIDMap.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "ECCEdRICHFastPIDMap.h"

#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

#include <cassert>
#include <cmath>
#include <iostream>

ECCEdRICHFastPIDMap::ECCEdRICHFastPIDMap()
{
}

ECCEdRICHFastPIDMap::~ECCEdRICHFastPIDMap()
{
}

ECCEdRICHFastPIDMap::PIDCandidate_LogLikelihood_map
ECCEdRICHFastPIDMap::getFastSmearLogLikelihood(int truth_pid, const double momentum, const double theta_rad) const
{
  assert(initialized);

  PIDCandidate_LogLikelihood_map ll_map;

  const int abs_truth_pid = abs(truth_pid);
  const double eta = -log(tan(0.5 * theta_rad));

  if (eta < etaMin() or eta > etaMax())
  {
    // not processing out of acceptance tracks
    if (Verbosity())
      std::cout << __PRETTY_FUNCTION__ << " " << mName << ": not processing out of acceptance tracks eta = " << eta
                << "  etaMin()  = " << etaMin()
                << "  etaMax()  = " << etaMax()
                << std::endl;

    return ll_map;
  }

  const double Nsigma_piK = numSigma(eta, momentum, pi_k);
  const double Nsigma_Kp = numSigma(eta, momentum, k_p);

  if (Verbosity())
    std::cout << __PRETTY_FUNCTION__ << " " << mName << ": processing tracks momentum = " << momentum
              << " eta = " << eta
              << " Nsigma_piK = " << Nsigma_piK << " Nsigma_Kp = " << Nsigma_Kp
              << std::endl;

  const double pion_sigma_space_ring_radius = +Nsigma_piK;
  const double kaon_sigma_space_ring_radius = 0;
  const double proton_sigma_space_ring_radius = -Nsigma_Kp;

  double sigma_space_ring_radius = gRandom->Gaus(0, 1);
  if (abs_truth_pid == EICPIDDefs::PionCandiate)
    sigma_space_ring_radius += pion_sigma_space_ring_radius;
  else if (abs_truth_pid == EICPIDDefs::KaonCandiate)
    sigma_space_ring_radius += kaon_sigma_space_ring_radius;
  else if (abs_truth_pid == EICPIDDefs::ProtonCandiate)
    sigma_space_ring_radius += proton_sigma_space_ring_radius;

  ll_map[EICPIDDefs::PionCandiate] = -0.5 * pow(sigma_space_ring_radius - pion_sigma_space_ring_radius, 2);
  ll_map[EICPIDDefs::KaonCandiate] = -0.5 * pow(sigma_space_ring_radius - kaon_sigma_space_ring_radius, 2);
  ll_map[EICPIDDefs::ProtonCandiate] = -0.5 * pow(sigma_space_ring_radius - proton_sigma_space_ring_radius, 2);

  return ll_map;
}

double
ECCEdRICHFastPIDMap::etaMin() const
{
  switch (mType)
  {
  case kBarrel:
    return -log(tan(atan2(mRadius, -mLength) * 0.5));
  case kForward:
    return -log(tan(atan2(mRadiusOut, mPositionZ) * 0.5));
  }
  return 0.;
}

double
ECCEdRICHFastPIDMap::etaMax() const
{
  switch (mType)
  {
  case kBarrel:
    return -log(tan(atan2(mRadius, mLength) * 0.5));
  case kForward:
    return -log(tan(atan2(mRadiusIn, mPositionZ) * 0.5));
  }
  return 0.;
}

void ECCEdRICHFastPIDMap::dualRICH_aerogel()
{
  initialized = true;
  setName("aerogel");
  /** geometry **/
  setType(kForward);
  setRadiusIn(10.);    // [cm]
  setRadiusOut(120.);  // [cm]
  setPositionZ(250.);  // [cm]
  /** radiator **/
  setLength(4.);  // [cm]
  setIndex(1.02);
  /** overall photon-detection efficiency **/
  setEfficiency(0.08);
  /** single-photon angular resolution **/
  double angle[5] = {5., 10., 15., 20., 25.};                                              // [deg]
  double chromatic[5] = {0.00260572, 0.00223447, 0.00229996, 0.00237615, 0.00245689};      // [rad] from actual file
  double emission[5] = {0.000658453, 0.000297004, 0.00014763, 0.000196477, 0.000596087};   // [rad] from actual file
  double pixel[5] = {0.000502646, 0.000575427, 0.000551095, 0.000555055, 0.000564831};     // [rad] from actual file
  double field[5] = {8.13634e-05, 6.41901e-05, 3.92289e-05, 9.76800e-05, 2.58328e-05};     // [rad] from actual file
  double tracking[5] = {0.000350351, 0.000306691, 0.000376006, 0.000401814, 0.000389742};  // [rad] from actual file
  setChromaticSigma(5, angle, chromatic);
  setPositionSigma(5, angle, pixel);
  setEmissionSigma(5, angle, emission);
  setFieldSigma(5, angle, field);
  setTrackingSigma(5, angle, tracking);
}

void ECCEdRICHFastPIDMap::dualRICH_C2F6()
{
  initialized = true;
  setName("C2F6");
  /** geometry **/
  setType(kForward);
  setRadiusIn(10.);    // [cm]
  setRadiusOut(120.);  // [cm]
  setPositionZ(250.);  // [cm]
  /** radiator **/
  setLength(160.);  // [cm]
  setIndex(1.0008);
  /** overall photon-detection efficiency **/
  setEfficiency(0.15);
  /** single-photon angular resolution **/
  double angle[5] = {5., 10., 15., 20., 25.};                                               // [deg]
  double chromatic[5] = {0.000516327, 0.000527914, 0.000525467, 0.000515349, 0.000489377};  // [rad] from actual file
  double emission[5] = {0.001439090, 0.000718037, 0.000656786, 0.000946782, 0.001404630};   // [rad] from actual file
  double pixel[5] = {0.000480520, 0.000533282, 0.000564187, 0.000577872, 0.000605236};      // [rad] from actual file
  double field[5] = {8.60521e-05, 7.64798e-05, 0.000167358, 0.000475598, 0.000629863};      // [rad] from actual file
  double tracking[5] = {0.000389136, 0.000328530, 0.000402517, 0.000417901, 0.000393391};   // [rad] from actual file
  setChromaticSigma(5, angle, chromatic);
  setPositionSigma(5, angle, pixel);
  setEmissionSigma(5, angle, emission);
  setFieldSigma(5, angle, field);
  setTrackingSigma(5, angle, tracking);
}

double ECCEdRICHFastPIDMap::cherenkovAngleSigma(double eta, double p, double m) const
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

double ECCEdRICHFastPIDMap::numSigma(double eta, double p, ECCEdRICHFastPIDMap::type PID) const
{
  double mass1(0), mass2(0);
  switch (PID)
  {
  case pi_k:
    mass1 = mMassPion;
    mass2 = mMassKaon;
    break;
  case k_p:
    mass1 = mMassKaon;
    mass2 = mMassProton;
    break;
  default:
    assert(0);  // exit the code for logical error
    exit(1);
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
