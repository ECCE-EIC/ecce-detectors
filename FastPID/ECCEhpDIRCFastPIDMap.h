// $Id: $

/*!
 * \file ECCEhpDIRCFastPIDMap.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef ECCEHPDIRCFASTPIDMAP_H_
#define ECCEHPDIRCFASTPIDMAP_H_

#include "ECCEFastPIDMap.h"

#include <string>

class TH2F;
class TF1;

/*!
 * \brief ECCEhpDIRCFastPIDMap
 * Import from DrcPidFast
 */
class ECCEhpDIRCFastPIDMap : public ECCEFastPIDMap {
public:
  ECCEhpDIRCFastPIDMap();
  virtual ~ECCEhpDIRCFastPIDMap();

  PIDCandidate_LogLikelihood_map
  getFastSmearLogLikelihood(int truth_pid, const double momentum,
                            const double theta) const override;

  //! read Cherenkov track resolution map from a file
  void ReadMap(const std::string &name);

  TH2F *GetTrrMap() { return fTrrMap; }

private:
  //! probability - normalized to 1 probability for e,mu,pi,k,p
  //! sigma - deviation of the determined Cherenkov angle from expected in terms
  //! of Cherenkov track resolution cangle - Cherenkov angle cctr -  combined
  //! Cherenkov track resolution
  struct DrcPidInfo {
    double probability[5] = {0};
    double sigma[5] = {0};
    double cangle = 0;
    double cctr = 0;
  };

  static int get_pid(int pdg);
  static EICPIDDefs::PIDCandidate get_PIDCandidate(int id);

  TH2F *fTrrMap = nullptr;
  double fMass[5] = {0};
  TF1 *fMs_mom = nullptr;
  TF1 *fMs_thickness = nullptr;
  double fMs_thickness_max = (0);
};

#endif /* ECCEHPDIRCFASTPIDMAP_H_ */
