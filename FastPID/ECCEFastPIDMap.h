// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ECCEFastPIDMap_H
#define ECCEFastPIDMap_H

#include <eicpidbase/EICPIDDefs.h>

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class ECCEFastPIDMap
{
 public:

  ECCEFastPIDMap(const std::string &name = "ECCEFastPIDMap");

  virtual ~ECCEFastPIDMap();


 private:
};

#endif // ECCEFastPIDMap_H
