
#include "ECCEFastPIDReco.h"
#include "ECCEFastPIDMap.h"

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>

#include <eicpidbase/EICPIDParticle.h>
#include <eicpidbase/EICPIDParticleContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
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

#include <CLHEP/Vector/ThreeVector.h>

#include <cassert>
#include <iostream>

using namespace std;

//____________________________________________________________________________..
ECCEFastPIDReco::ECCEFastPIDReco(ECCEFastPIDMap *map,
                                 const std::string &name)
  : SubsysReco(name)
  , m_pidmap(map)
{
  if (!map)
  {
    std::cout << "ECCEFastPIDReco::ECCEFastPIDReco(): Fatal Error missing ECCEFastPIDMap" << std::endl;
    exit(1);
  }
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

  /// G4 truth particle node
  m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthInfo)
  {
    cout << __PRETTY_FUNCTION__
         << ": PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_matchG4Hit)
  {
    m_g4hits = findNode::getClass<PHG4HitContainer>(topNode, m_G4HitNodeName);
    if (!m_g4hits)
    {
      cout << __PRETTY_FUNCTION__
           << ": PHG4HitContainer " << m_G4HitNodeName << " node is missing, can't match to g4hits"
           << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  m_EICPIDParticleContainer = findNode::getClass<EICPIDParticleContainer>(topNode, m_EICPIDParticleMapNodeName);
  if (!m_EICPIDParticleContainer)
  {
    m_EICPIDParticleContainer = new EICPIDParticleContainer;

    PHIODataNode<PHObject> *pid_node = new PHIODataNode<PHObject>(m_EICPIDParticleContainer, m_EICPIDParticleMapNodeName, "PHObject");

    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    assert(dstNode);
    dstNode->addNode(pid_node);
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

  if (Verbosity() >= 1)
  {
    cout << __PRETTY_FUNCTION__ << ": m_SvtxTrackMap = ";
    m_SvtxTrackMap->identify();
  }

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

      assert(m_truthInfo);

      PHG4Particle *g4particle =
          m_truthInfo->GetParticle(fasttrack->get_truth_track_id());
      assert(g4particle);

      const int truth_pid = g4particle->get_pid();

      CLHEP::Hep3Vector momentum;
      if (m_matchG4Hit)
      {
        // find hit in detector vol.
        assert(m_g4hits);

        auto hit_range = m_g4hits->getHits();
        for (auto hititer = hit_range.first; hititer != hit_range.second; ++hititer)
        {
          const PHG4Hit *hit = hititer->second;
          assert(hit);

          // match first hit of track in PID detector:
          if (hit->get_trkid() == g4particle->get_track_id())
          {
            momentum.set(
                hit->get_px(0),
                hit->get_py(0),
                hit->get_pz(0));

            if (Verbosity() >= 2)
            {
              cout << __PRETTY_FUNCTION__ << " Named " << Name() << " with hits " << m_G4HitNodeName << ": matching track ";
              fasttrack->identify();
              cout << "with particle: ";
              g4particle->identify();
              cout << "with hit: ";
              hit->identify();
              cout << "Result in momentum of " << momentum.mag() << "GeV/c at eta = " << momentum.eta() << endl;
            }

            break;
          }
        }  // for

        if (momentum.mag() == 0)
        {
          if (Verbosity() >= 2)
          {
            cout << __PRETTY_FUNCTION__ << " Named " << Name() << " did NOT match hits in "
                 << m_G4HitNodeName << " with track ";
            fasttrack->identify();
            cout << "with particle: ";
            g4particle->identify();
          }
        }
      }
      else
      {
        // use vertex kinematics
        momentum.set(
            g4particle->get_px(),
            g4particle->get_py(),
            g4particle->get_pz());
      }

      auto pid_iter = m_EICPIDParticleContainer->findOrAddPIDParticle(fasttrack->get_id());
      EICPIDParticle *pidparticle = pid_iter->second;
      assert(pidparticle);

      pidparticle->set_property(EICPIDParticle::Truth_PID, truth_pid);

      assert(m_pidmap);

      if (momentum.mag() > 0)
      {
        ECCEFastPIDMap::PIDCandidate_LogLikelihood_map ll_map =
            m_pidmap->getFastSmearLogLikelihood(truth_pid, momentum.mag(), momentum.theta());

        for (const auto &pair : ll_map)
          pidparticle->set_LogLikelyhood(pair.first, EICPIDDefs::DIRC, pair.second);
      }

      if (Verbosity() >= 2)
      {
        pidparticle->identify();
      }
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
