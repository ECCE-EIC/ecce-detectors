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
#include <TH2F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

#include <cassert>
#include <cmath>

ECCEdRICHFastPIDMap::ECCEdRICHFastPIDMap()
{
  int barid = 0;

  fMass[0] = 0.000511;
  fMass[1] = 0.105658;
  fMass[2] = 0.139570;
  fMass[3] = 0.49368;
  fMass[4] = 0.938272;

  // multiple scattering for 17 mm thick radiator at 30 deg
  fMs_mom = new TF1("", "expo(0)+expo(2)+expo(4)");
  fMs_mom->SetParameters(4.40541e+00, -5.52436e+00, 2.35058e+00, -1.02703e+00, 9.55032e-01,
                         -1.48500e-01);
  // fMs_mom->SetParameters(9.39815e-01, -1.48243e-01, 4.69733e+00, -4.33960e+00, 2.19745e+00,
  //                       -9.68617e-01);

  fMs_thickness = new TF1("", "pol1");

  if (barid == 1)
    fMs_thickness->SetParameters(3.5, 0.0214286);  // 10 mm bar
  else
    fMs_thickness->SetParameters(4.5, 0.0357143);  // 17 mm bar

  TF1 *fMs_thickness_17 = new TF1("", "pol1");
  fMs_thickness_17->SetParameters(4.5, 0.0357143);  // 17 mm bar
  fMs_thickness_max = fMs_thickness_17->Eval(70);
}

ECCEdRICHFastPIDMap::~ECCEdRICHFastPIDMap()
{
}

ECCEdRICHFastPIDMap::PIDCandidate_LogLikelihood_map
ECCEdRICHFastPIDMap::getFastSmearLogLikelihood(int truth_pid, const double momentum, const double theta_rad) const
{
  assert(fTrrMap);

  // preprocessing
  double theta = theta_rad * TMath::RadToDeg();
  truth_pid = abs(truth_pid);
  // track_err - error assosiated with track direction [mrad]
  double track_err = 0.326;
  double p = momentum;

  PIDCandidate_LogLikelihood_map ll_map;

  // copy from DrcPidFast::
  // pdg - Particle Data Group code of the particle
  // mom - 3-momentum of the particle [GeV/c]
  // track_err - error assosiated with track direction [mrad]
  //  DrcPidInfo GetInfo(int pdg, TVector3 mom, double track_err = 0);

  const int max = 5;
  DrcPidInfo info;
  int pid = get_pid(truth_pid);

  if (pid == 0)
  {
    return ll_map;
  }

  // set default values
  for (int i = 0; i < max; i++)
  {
    info.probability[i] = 0.25;
    info.sigma[i] = 100;
  }
  info.cangle = 0;
  info.cctr = 0;

  // check range
  //  if (theta < 19.99 || theta > 160.01){
  //    std::cout<<"theta out of [20,160] deg range: "<<theta<<std::endl;
  //  }

  double ms_mom_err = fMs_mom->Eval(p);  // vector deviation after radiator

  double alpha = (theta < 90) ? 90 - theta : theta - 90;
  double ms_thick_frac = fMs_thickness->Eval(alpha) / fMs_thickness_max;

  // 0.31 for averaging direction vector over the radiator thickness
  double ms_err = 0.31 * ms_mom_err * ms_thick_frac;

  // ctr map is for theta = [25,153] and p = [0,10] GeV/c
  if (theta < 25) theta = 25;
  if (theta > 153) theta = 153;
  if (p > 10) p = 10;

  int bin = fTrrMap->FindBin(theta, p);
  double ctr = fTrrMap->GetBinContent(bin);  // Cherenkov track resolution [mrad]
  double cctr = sqrt(ctr * ctr + track_err * track_err + ms_err * ms_err) *
                0.001;  // combined Cherenkov track resolution[rad]

  // 1.46907 - fused silica
  double true_cangle = acos(sqrt(p * p + fMass[pid] * fMass[pid]) / p / 1.46907);
  true_cangle += gRandom->Gaus(0, cctr);

  // return default values if momentum below Cherenkov threshold (true_cangle is NaN)
  if (isnan(true_cangle)) return ll_map;

  double cangle, sum = 0, fsum = 0;
  double delta[max] = {0};  //, probability[max] = {0};

  for (int i = 0; i < max; i++)
  {
    cangle = acos(sqrt(p * p + fMass[i] * fMass[i]) / p / 1.46907);
    if (isnan(cangle))
    {
      ll_map[get_PIDCandidate(i)] = -100;  // set non-firing particle candidate to low probability
      continue;
    }
    delta[i] = fabs(cangle - true_cangle);
    sum += delta[i];
    info.sigma[i] = (cangle - true_cangle) / cctr;
    if (i == pid) info.cangle = cangle;

    ll_map[get_PIDCandidate(i)] = -0.5 * info.sigma[i] * info.sigma[i];
  }
  // normalization
  for (int i = 0; i < max; i++)
  {
    if (delta[i] > 0) info.probability[i] = sum / delta[i];
    fsum += info.probability[i];
  }
  for (int i = 0; i < max; i++) info.probability[i] /= fsum;
  info.cctr = cctr;

  return ll_map;
}

void ECCEdRICHFastPIDMap::ReadMap(const std::string &name)
{
  TFile *file = TFile::Open(name.c_str());
  assert(file);
  //  fTrrMap = new TH2F();
  file->GetObject("htrr", fTrrMap);
  assert(fTrrMap);
}

int ECCEdRICHFastPIDMap::get_pid(int pdg)
{
  int pid = 0;
  if (pdg == 11) pid = 0;    // e
  if (pdg == 13) pid = 1;    // mu
  if (pdg == 211) pid = 2;   // pi
  if (pdg == 321) pid = 3;   // K
  if (pdg == 2212) pid = 4;  // p
  return pid;
}

EICPIDDefs::PIDCandidate ECCEdRICHFastPIDMap::get_PIDCandidate(int pid)
{
  EICPIDDefs::PIDCandidate id = EICPIDDefs::InvalidCandiate;
  if (pid == 0) id = EICPIDDefs::ElectronCandiate;  // e
  if (pid == 1) id = EICPIDDefs::MuonCandiate;      // mu
  if (pid == 2) id = EICPIDDefs::PionCandiate;      // pi
  if (pid == 3) id = EICPIDDefs::KaonCandiate;      // K
  if (pid == 4) id = EICPIDDefs::ProtonCandiate;    // p
  assert(id != EICPIDDefs::InvalidCandiate);

  return id;
}
