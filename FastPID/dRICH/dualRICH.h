/// \author R+Preghenella
/// \email  preghenella@bo.infn.it
/// \date   April 2020

#ifndef __DUALRICH_H__
#define __DUALRICH_H__

#include "genericRICH.h"

/** dualRICH aerogel **/

class dualRICH_aerogel : public genericRICH
{
 public:
  dualRICH_aerogel();
  virtual ~dualRICH_aerogel() = default;
};

dualRICH_aerogel::dualRICH_aerogel()
{  
  setName("aerogel");
  /** geometry **/
  setType(genericDetector::kForward);
  setRadiusIn(10.); // [cm]
  setRadiusOut(120.); // [cm]
  setPositionZ(250.); // [cm]
  /** radiator **/
  setLength(4.); // [cm]
  setIndex(1.02);
  /** overall photon-detection efficiency **/
  setEfficiency(0.08);
  /** single-photon angular resolution **/
  double angle[5] = {5., 10., 15., 20., 25.}; // [deg]
  double chromatic[5] = {0.00260572, 0.00223447, 0.00229996, 0.00237615, 0.00245689}; // [rad] from actual file
  double emission[5] = {0.000658453, 0.000297004, 0.00014763, 0.000196477, 0.000596087}; // [rad] from actual file
  double pixel[5] = {0.000502646, 0.000575427, 0.000551095, 0.000555055, 0.000564831}; // [rad] from actual file
  double field[5] = {8.13634e-05, 6.41901e-05, 3.92289e-05, 9.76800e-05, 2.58328e-05}; // [rad] from actual file
  double tracking[5] = {0.000350351, 0.000306691, 0.000376006, 0.000401814, 0.000389742}; // [rad] from actual file
  setChromaticSigma(5, angle, chromatic);
  setPositionSigma(5, angle, pixel);
  setEmissionSigma(5, angle, emission);
  setFieldSigma(5, angle, field);
  setTrackingSigma(5, angle, tracking);
}

/** dualRICH C2F6 **/

class dualRICH_C2F6 : public genericRICH
{
 public:
  dualRICH_C2F6();
  virtual ~dualRICH_C2F6() = default;
};

dualRICH_C2F6::dualRICH_C2F6()
{
  setName("C2F6");
  /** geometry **/
  setType(genericDetector::kForward);
  setRadiusIn(10.); // [cm]
  setRadiusOut(120.); // [cm]
  setPositionZ(250.); // [cm]
  /** radiator **/
  setLength(160.); // [cm]
  setIndex(1.0008);
  /** overall photon-detection efficiency **/
  setEfficiency(0.15);
  /** single-photon angular resolution **/
  double angle[5] = {5., 10., 15., 20., 25.}; // [deg]
  double chromatic[5] = {0.000516327, 0.000527914, 0.000525467, 0.000515349, 0.000489377}; // [rad] from actual file
  double emission[5] = {0.001439090, 0.000718037, 0.000656786, 0.000946782, 0.001404630}; // [rad] from actual file
  double pixel[5] = {0.000480520, 0.000533282, 0.000564187, 0.000577872, 0.000605236}; // [rad] from actual file
  double field[5] = {8.60521e-05, 7.64798e-05, 0.000167358, 0.000475598, 0.000629863}; // [rad] from actual file
  double tracking[5] = {0.000389136, 0.000328530, 0.000402517, 0.000417901, 0.000393391}; // [rad] from actual file
  setChromaticSigma(5, angle, chromatic);
  setPositionSigma(5, angle, pixel);
  setEmissionSigma(5, angle, emission);
  setFieldSigma(5, angle, field);
  setTrackingSigma(5, angle, tracking);
}

#endif /* __DUALRICH_H__ */
