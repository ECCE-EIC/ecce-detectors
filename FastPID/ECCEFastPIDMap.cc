

#include "ECCEFastPIDMap.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

//____________________________________________________________________________..
ECCEFastPIDMap::ECCEFastPIDMap(const std::string &name) {}

//____________________________________________________________________________..
ECCEFastPIDMap::~ECCEFastPIDMap() {}

EICPIDDefs::PIDCandidate getFastSmearPID(int truth_pid, const double momentum) {
  return EICPIDDefs::InvalidCandiate;
}
