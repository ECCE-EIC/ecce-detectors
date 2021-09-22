// $Id: $

/*!
 * \file ECCEmRICHFastPIDMap.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "ECCEmRICHFastPIDMap.h"

#include <TF1.h>
#include <TFile.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

#include <cassert>
#include <cmath>
#include <iostream>

ECCEmRICHFastPIDMap::ECCEmRICHFastPIDMap(double trackResolution,
                                         double incidentAngle, double pixS) {
  fTrackResolution = trackResolution;
  pLow = 0.6;
  pHigh = 20;
  c = 0.0299792458; // cm/picosecond
  n = 1.03;         // Aerogel
  a = pixS;         // pixel size 3.0; // mm -- one side
  f = 152.4;        // focal length mm =  6"
  N_gam = 10;
  mElectron = 0.00051099895;    //GeV/c^2
  mPion = 0.13957018;      // GeV/c^2
  mKaon = 0.493677;        // GeV/c^2
  mProton = 0.93827208816; // GeV/c^2
  pi = 3.14159;
  alpha = 0.0072973525693; // hyperfine const
  L = 3.0;                 // Aerogel block thickness in cm

  //===============
  th0 = incidentAngle; // incidence angle in radians
}

ECCEmRICHFastPIDMap::~ECCEmRICHFastPIDMap() {}

ECCEmRICHFastPIDMap::PIDCandidate_LogLikelihood_map
ECCEmRICHFastPIDMap::getFastSmearLogLikelihood(int truth_pid,
                                               const double momentum,
                                               const double theta_rad) const {
  const int abs_truth_pid = abs(truth_pid);

  PIDCandidate_LogLikelihood_map ll_map;

  if (abs_truth_pid != EICPIDDefs::ElectronCandiate and
      abs_truth_pid != EICPIDDefs::PionCandiate and
      abs_truth_pid != EICPIDDefs::KaonCandiate and
      abs_truth_pid != EICPIDDefs::ProtonCandiate) {
    // not processing non-hadronic tracks
    if (Verbosity())
      std::cout << __PRETTY_FUNCTION__
                << ":  not processing non-hadronic tracks " << truth_pid
                << std::endl;

    return ll_map;
  }
  if (theta_rad < m_acceptanceThetaMin or theta_rad > m_acceptanceThetaMax) {
    // not processing out of acceptance tracks
    if (Verbosity())
      std::cout << __PRETTY_FUNCTION__
                << ": not processing out of acceptance tracks theta_rad = "
                << theta_rad << std::endl;

    return ll_map;
  }
  if (momentum < pLow or momentum > pHigh) {
    // not processing out of acceptance tracks
    if (Verbosity())
      std::cout << __PRETTY_FUNCTION__
                << ": not processing out of acceptance tracks momentum = "
                << momentum << std::endl;

    return ll_map;
  }

  double Nsigma_epi = 0;
  double Nsigma_piK = 0;
  double Nsigma_Kp = 0;

  //electron-pion case
  if(momentum>0.59){ // pion cerenkov's threhsold at n=1.03 is ~0.59 GeV/c
    // Angle difference
    double dth = getAng(mElectron) - getAng(mPion);
    // Detector uncertainty
    double sigTh = sqrt(pow(getdAng(mPion, momentum), 2) +
                        pow(getdAng(mElectron, momentum), 2));
    // Global uncertainty
    double sigThTrk = getdAngTrk(mElectron, momentum);
    double sigThc = sqrt(pow(sigTh / sqrt(getNgamma(L, mElectron, momentum)), 2) +
                         pow(sigThTrk, 2));
    Nsigma_epi = dth / sigThc;
    if (isnan(Nsigma_epi)) Nsigma_epi = 0;
  }
  //pion-Kaon case
  if(momentum>2.0){ // kaon cerenkov's threhsold at n=1.03 is ~2 GeV/c
    // Angle difference
    double dth = getAng(mPion, momentum) - getAng(mKaon, momentum);
    // Detector uncertainty
    double sigTh = sqrt(pow(getdAng(mPion, momentum), 2) +
                        pow(getdAng(mKaon, momentum), 2));
    // Global uncertainty
    double sigThTrk = getdAngTrk(mPion, momentum);
    double sigThc = sqrt(pow(sigTh / sqrt(getNgamma(L, mPion, momentum)), 2) +
                         pow(sigThTrk, 2));
    Nsigma_piK = dth / sigThc;
    if (isnan(Nsigma_piK)) Nsigma_piK = 0;
  }
  //kaon-proton case
  if(momentum>3.8){ // Proton cerenkov's threhsold at n=1.03 is ~3.8 GeV/c
    // Angle difference
    double dth = getAng(mKaon, momentum) - getAng(mProton, momentum);
    // Detector uncertainty
    double sigTh = sqrt(pow(getdAng(mKaon, momentum), 2) +
                        pow(getdAng(mProton, momentum), 2));
    // Global uncertainty
    double sigThTrk = getdAngTrk(mKaon, momentum);
    double sigThc = sqrt(pow(sigTh / sqrt(getNgamma(L, mKaon, momentum)), 2) +
                         pow(sigThTrk, 2));
    Nsigma_Kp = dth / sigThc;
    if (isnan(Nsigma_Kp)) Nsigma_Kp = 0;
  }

  if (Verbosity())
    std::cout << __PRETTY_FUNCTION__
              << ": processing tracks momentum = " << momentum
              << " Nsigma_epi = " << Nsigma_epi
              << " Nsigma_piK = " << Nsigma_piK
              << " Nsigma_Kp = " << Nsigma_Kp
              << std::endl;

  const double electron_sigma_space_ring_radius = +Nsigma_piK + Nsigma_epi;
  const double pion_sigma_space_ring_radius = +Nsigma_piK;
  const double kaon_sigma_space_ring_radius = 0;
  const double proton_sigma_space_ring_radius = -Nsigma_Kp;

  double sigma_space_ring_radius = gRandom->Gaus(0, 1);
  if (abs_truth_pid == EICPIDDefs::ElectronCandiate)
    sigma_space_ring_radius += electron_sigma_space_ring_radius;
  else if (abs_truth_pid == EICPIDDefs::PionCandiate)
    sigma_space_ring_radius += pion_sigma_space_ring_radius;
  else if (abs_truth_pid == EICPIDDefs::KaonCandiate)
    sigma_space_ring_radius += kaon_sigma_space_ring_radius;
  else if (abs_truth_pid == EICPIDDefs::ProtonCandiate)
    sigma_space_ring_radius += proton_sigma_space_ring_radius;

  ll_map[EICPIDDefs::ElectronCandiate] =
      -0.5 * pow(sigma_space_ring_radius - electron_sigma_space_ring_radius, 2);
  ll_map[EICPIDDefs::PionCandiate] =
      -0.5 * pow(sigma_space_ring_radius - pion_sigma_space_ring_radius, 2);
  ll_map[EICPIDDefs::KaonCandiate] =
      -0.5 * pow(sigma_space_ring_radius - kaon_sigma_space_ring_radius, 2);
  ll_map[EICPIDDefs::ProtonCandiate] =
      -0.5 * pow(sigma_space_ring_radius - proton_sigma_space_ring_radius, 2);

  return ll_map;
}

// Angle exiting the Aerogel
double ECCEmRICHFastPIDMap::getAng(double mass, double mom) const {
  double beta = mom / sqrt(pow(mom, 2) + pow(mass, 2));
  double thc = acos(1. / n / beta);
  double th0p = asin(sin(th0) / n);
  double dthc = thc - th0p;
  double theta = asin(n * sin(dthc));

  return theta;
}

// Uncertainty due to detector effects
double ECCEmRICHFastPIDMap::getdAng(double mass, double mom) const {
  double beta = mom / sqrt(pow(mom, 2) + pow(mass, 2));
  double thc = acos(1. / n / beta);
  double th0p = asin(sin(th0) / n);
  double dthc = thc - th0p;
  double theta = asin(n * sin(dthc));

  double sig_ep = 0;   // Emission point error
  double sig_chro = 0; // Chromatic dispersion error
  double sig_pix = a * pow(cos(theta), 2) / f / sqrt(12.);

  double sigTh = sqrt(pow(sig_ep, 2) + pow(sig_chro, 2) + pow(sig_pix, 2));

  return sigTh;
}

// Uncertainty due to tracking resolution
double ECCEmRICHFastPIDMap::getdAngTrk(double mass, double mom) const {
  double beta = mom / sqrt(pow(mom, 2) + pow(mass, 2));
  double thc = acos(1. / n / beta);
  double th0p = asin(sin(th0) / n);
  double dthc = thc - th0p;
  double theta = asin(n * sin(dthc));

  double sig_trk =
      (cos(dthc) / cos(theta)) * (cos(th0) / cos(th0p)) * fTrackResolution;

  return sig_trk;
}

// no. of gamms
double ECCEmRICHFastPIDMap::getNgamma(double t, double mass, double mom) const {
  int tot = 10000;
  double beta = mom / sqrt(pow(mom, 2) + pow(mass, 2));
  double fact = 2. * pi * alpha * t * (1. - 1. / pow(n * beta, 2));
  double T_lensWin = 0.92 * 0.92;
  double xmin = 300.e-7;
  double xmax = 650.e-7;
  double dx = (xmax - xmin) / tot;
  double sum = 0;
  for (int j = 0; j < tot; j++) {
    double x = xmin + j * dx + dx / 2.;
    sum += T_QE(x) * T_Aer(t, x) / pow(x, 2);
  }
  return fact * T_lensWin * sum * dx;
}

// Quantum efficiency
double ECCEmRICHFastPIDMap::T_QE(double lam) const {
  return 0.34 * exp(-1. * pow(lam - 345.e-7, 2) / (2. * pow(119.e-7, 2)));
}

// Transmissions of the radiator block
double ECCEmRICHFastPIDMap::T_Aer(double t, double lam) const {
  return 0.83 * exp(-1. * t * 56.29e-20 / pow(lam, 4));
}
