/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   LinacInc/Segment8.h
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
#ifndef tdcSystem_Segment8_h
#define tdcSystem_Segment8_h

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
  class LQuadF;
  class CorrectorMag;

  /*!
    \class Segment8
    \version 1.0
    \author S. Ansell
    \date May 2020
    \brief Seventh segment
  */

class Segment8 :
  public TDCsegment
{
 private:

  /// first bellow
  std::shared_ptr<constructSystem::Bellows> bellowA;   
  /// Beam stop
  std::shared_ptr<tdcSystem::EBeamStop> eBeamStop;   
  /// second bellow
  std::shared_ptr<constructSystem::Bellows> bellowB;   
  /// first pipe
  std::shared_ptr<constructSystem::VacuumPipe> pipeA;   
  
  void buildObjects(Simulation&);
  void createLinks();
  
 public:
  
  Segment8(const std::string&);
  Segment8(const Segment8&);
  Segment8& operator=(const Segment8&);
  ~Segment8();


  using FixedComp::createAll;
  virtual void createAll(Simulation&,const attachSystem::FixedComp&,
		 const long int);

};

}

#endif
