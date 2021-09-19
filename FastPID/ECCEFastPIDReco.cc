
#include "ECCEFastPIDReco.h"
#include "ECCEFastPIDMap.h"

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>

#include <eicpidbase/EICPIDParticle.h>
#include <eicpidbase/EICPIDParticleContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cassert>
#include <iostream>

using namespace std;

//____________________________________________________________________________..
ECCEFastPIDReco::ECCEFastPIDReco(const std::string &name)
  : SubsysReco(name)
{
  std::cout << "ECCEFastPIDReco::ECCEFastPIDReco(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
ECCEFastPIDReco::~ECCEFastPIDReco()
{
  if (m_pidmap) delete m_pidmap;
}

//____________________________________________________________________________..
int ECCEFastPIDReco::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCEFastPIDReco::InitRun(PHCompositeNode *topNode)
{
  m_SvtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, m_TrackmapNodeName);
  if (!m_SvtxTrackMap)
  {
    cout << __PRETTY_FUNCTION__ << " fatal error missing node " << m_TrackmapNodeName << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_EICPIDParticleContainer = findNode::getClass<EICPIDParticleContainer>(topNode, m_EICPIDParticleMapNodeName);
  if (!m_EICPIDParticleContainer)
  {
    m_EICPIDParticleContainer = new EICPIDParticleContainer;

    PHIODataNode<PHObject> *pid_node = new PHIODataNode<PHObject>(m_EICPIDParticleContainer, m_EICPIDParticleMapNodeName, "PHObject");
    topNode->addNode(pid_node);
    if (Verbosity() > 0)
    {
      cout << m_EICPIDParticleMapNodeName << " node added" << endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCEFastPIDReco::process_event(PHCompositeNode *topNode)
{
  assert(m_SvtxTrackMap);
  assert(m_EICPIDParticleContainer);

  for (const auto &track_pair : *m_SvtxTrackMap)
  {
    const SvtxTrack *track = track_pair.second;
    assert(track);

    const SvtxTrack_FastSim *fasttrack = dynamic_cast<const SvtxTrack_FastSim *>(track);

    if (!fasttrack)
    {
      if (Verbosity())
      {
        cout << __PRETTY_FUNCTION__ << " : ignore none SvtxTrack_FastSim track: ";
        track->identify();
      }
    }
    else
    {  // process track
      auto iter = m_EICPIDParticleContainer->findOrAddPIDParticle(track->get_id());

      EICPIDParticle *pidparticle = iter->second;

      assert(pidparticle);
    }
  }

  if (Verbosity())
  {
    cout << __PRETTY_FUNCTION__ << " : done processing from trackmap ";
    m_SvtxTrackMap->identify();
    cout << __PRETTY_FUNCTION__ << " : produced EICPIDParticleContainer ";
    m_EICPIDParticleContainer->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCEFastPIDReco::End(PHCompositeNode *topNode)
{
  std::cout << "ECCEFastPIDReco::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
