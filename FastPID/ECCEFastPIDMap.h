// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ECCEFastPIDMap_H
#define ECCEFastPIDMap_H

#include <eicpidbase/EICPIDDefs.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>

class PHCompositeNode;

class ECCEFastPIDMap {
public:
  ECCEFastPIDMap(const std::string &name = "ECCEFastPIDMap");

  virtual ~ECCEFastPIDMap();

  typedef std::map<EICPIDDefs::PIDCandidate, float>
      PIDCandidate_LogLikelihood_map;

  virtual PIDCandidate_LogLikelihood_map
  getFastSmearLogLikelihood(int truth_pid, const double momentum,
                            const double theta) const = 0;

  //  EICPIDDefs::PIDCandidate getFastSmearPID(int truth_pid, const double
  //  momentum);

  /// Sets the verbosity of this module (0 by default=quiet).
  virtual void Verbosity(const int ival) { m_Verbosity = ival; }

  /// Gets the verbosity of this module.
  virtual int Verbosity() const { return m_Verbosity; }

private:
  /// The verbosity level. 0 means not verbose at all.
  int m_Verbosity = 0;
};

#endif // ECCEFastPIDMap_H
