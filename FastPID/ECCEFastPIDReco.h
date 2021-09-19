// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ECCEFASTPIDRECO_H
#define ECCEFASTPIDRECO_H

#include <fun4all/SubsysReco.h>

#include <eicpidbase/EICPIDDefs.h>

#include <string>

class PHCompositeNode;
class ECCEFastPIDMap;
class SvtxTrackMap;
class EICPIDParticleContainer;

class ECCEFastPIDReco : public SubsysReco
{
 public:
  ECCEFastPIDReco(const std::string &name = "ECCEFastPIDReco");

  virtual ~ECCEFastPIDReco();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void setTrackmapNodeName(const std::string &name) { m_TrackmapNodeName = name; }

  void setEICPIDParticleMapNodeName(const std::string &name) { m_EICPIDParticleMapNodeName = name; }

 private:
  bool m_matchG4Hit = false;

  ECCEFastPIDMap *m_pidmap = nullptr;

  std::string m_TrackmapNodeName = "SvtxTrackMap";
  std::string m_EICPIDParticleMapNodeName = "EICPIDParticleMap";

  SvtxTrackMap * m_SvtxTrackMap = nullptr;
  EICPIDParticleContainer * m_EICPIDParticleContainer = nullptr;
};

#endif  // ECCEFASTPIDRECO_H
