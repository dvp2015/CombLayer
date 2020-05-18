/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   LinacInc/L2SPFsegment5.h
 *
 * Copyright (c) 2004-2020 by Stuart Ansell
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 *
 ****************************************************************************/
#ifndef tdcSystem_L2SPFsegment5_h
#define tdcSystem_L2SPFsegment5_h

namespace constructSystem
{
  class VacuumPipe;
  class Bellows;
  class portItem;
  class BlankTube;
  class PipeTube;
}

/*!
  \namespace xraySystem
  \brief General xray optics system
  \version 1.0
  \date January 2018
  \author S. Ansell
*/

namespace tdcSystem
{
  class CorrectorMag;
  class LQuad;
  class LSexupole;
  class BPM;

  
  /*!
    \class L2SPFsegment5
    \version 1.0
    \author S. Ansell
    \date April 2020
    \brief Second segment in the TDC from the linac
  */

class L2SPFsegment5 :
  public TDCsegment
{
 private:

  /// first pipe 
  std::shared_ptr<constructSystem::VacuumPipe> pipeA;   

  /// Beam position monitor
  std::shared_ptr<tdcSystem::BPM> bpmA;

  // first pipe in dipole
  std::shared_ptr<constructSystem::VacuumPipe> pipeB;   
  
  /// Quad begining QSQ 
  std::shared_ptr<tdcSystem::LQuad> QuadA;

  /// sexupole
  std::shared_ptr<tdcSystem::LSexupole> SexuA;

  /// Quad endng QSQ 
  std::shared_ptr<tdcSystem::LQuad> QuadB;
  
  /// yag station
  std::shared_ptr<tdcSystem::YagUnit> yagUnit;

  /// exit pipe for corrector mags
  std::shared_ptr<constructSystem::Bellows> bellowA;   

  /// exit pipe for corrector mags
  std::shared_ptr<constructSystem::VacuumPipe> pipeC;   
  
  std::shared_ptr<tdcSystem::CorrectorMag> cMagHorC;    ///< corrector mag
  
  std::shared_ptr<tdcSystem::CorrectorMag> cMagVertC;   ///< corrector mag

  void buildObjects(Simulation&);
  void createLinks();
  
 public:
  
  L2SPFsegment5(const std::string&);
  L2SPFsegment5(const L2SPFsegment5&);
  L2SPFsegment5& operator=(const L2SPFsegment5&);
  ~L2SPFsegment5();

  
  void createAll(Simulation&,const attachSystem::FixedComp&,
		 const long int);

};

}

#endif