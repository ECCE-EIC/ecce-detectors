// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ECCEFastPIDMap_H
#define ECCEFastPIDMap_H

#include <eicpidbase/EICPIDDefs.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class PHCompositeNode;

class ECCEFastPIDMap
{
 public:
  ECCEFastPIDMap(const std::string &name = "ECCEFastPIDMap");

  virtual ~ECCEFastPIDMap();

  typedef std::map<EICPIDDefs::PIDCandidate, float> PIDCandidate_LogLikelihood_map;

  virtual PIDCandidate_LogLikelihood_map
    getFastSmearLogLikelihood(int truth_pid, const double momentum, const double theta) const = 0;

//  EICPIDDefs::PIDCandidate getFastSmearPID(int truth_pid, const double momentum);

 private:
};

#endif  // ECCEFastPIDMap_H
