// $Id: $

/*!
 * \file Test_ECCEFastPIDMap.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include <eccefastpidreco/ECCEFastPIDReco.h>
#include <eccefastpidreco/ECCEdRICHFastPIDMap.h>
#include <eccefastpidreco/ECCEhpDIRCFastPIDMap.h>
#include <eccefastpidreco/ECCEmRICHFastPIDMap.h>
#include <eicpidbase/EICPIDDefs.h>

#include <iostream>
#include <map>
#include <utility>

R__LOAD_LIBRARY(libECCEFastPIDReco.so)

void Test_ECCEFastPIDMap(const double truth_pid = 211,   //
                         const double momentum_GeV = 8,  //
                         const double theta_rad = 3,     //
                         const double LogLikelihoodCut = 0)
{
  ECCEmRICHFastPIDMap* pidmap = new ECCEmRICHFastPIDMap();

  //  ECCEhpDIRCFastPIDMap * pidmap = new ECCEhpDIRCFastPIDMap();
  //  pidmap->ReadMap( string(getenv("CALIBRATIONROOT")) + string("/hpDIRC/FastPID/ctr_map_p1_0.95.root") );

  //  ECCEdRICHFastPIDMap *pidmap = new ECCEdRICHFastPIDMap();
  //  pidmap->dualRICH_aerogel(); // pick one
  //  pidmap->dualRICH_C2F6();// pick one

  pidmap->Verbosity(1);

  const int n_trial = 1000;
  int n_electronID = 0;
  int n_pionID = 0;
  int n_kaonID = 0;
  int n_protonID = 0;

  for (int i = 0; i < n_trial; ++i)
  {
    std::map<EICPIDDefs::PIDCandidate, float> ll_map =
        pidmap->getFastSmearLogLikelihood(truth_pid, momentum_GeV,
                                          theta_rad);

    if (ll_map[EICPIDDefs::PionCandiate] > ll_map[EICPIDDefs::KaonCandiate] + LogLikelihoodCut         //
        and ll_map[EICPIDDefs::PionCandiate] > ll_map[EICPIDDefs::ProtonCandiate] + LogLikelihoodCut)  //
      ++n_pionID;
    if (ll_map[EICPIDDefs::KaonCandiate] > ll_map[EICPIDDefs::PionCandiate] + LogLikelihoodCut         //
        and ll_map[EICPIDDefs::KaonCandiate] > ll_map[EICPIDDefs::ProtonCandiate] + LogLikelihoodCut)  //
      ++n_kaonID;
    if (ll_map[EICPIDDefs::ProtonCandiate] > ll_map[EICPIDDefs::PionCandiate] + LogLikelihoodCut       //
        and ll_map[EICPIDDefs::ProtonCandiate] > ll_map[EICPIDDefs::KaonCandiate] + LogLikelihoodCut)  //
      ++n_protonID;

    // a simple electron ID
    if (ll_map[EICPIDDefs::ElectronCandiate] > ll_map[EICPIDDefs::PionCandiate] + LogLikelihoodCut)  //
      ++n_electronID;
  }

  std::cout << "Probability for truth_pid = " << truth_pid << " to be identified as following candidate with LogLikelihoodCut = " << LogLikelihoodCut << std::endl;

  std::cout << "\t"
            << "Pion"
            << " = " << (double) n_pionID / n_trial << std::endl;
  std::cout << "\t"
            << "Kaon"
            << " = " << (double) n_kaonID / n_trial << std::endl;
  std::cout << "\t"
            << "Proton"
            << " = " << (double) n_protonID / n_trial << std::endl;
  std::cout << "\t"
            << "ElectronID"
            << " = " << (double) n_electronID / n_trial << std::endl;
}
