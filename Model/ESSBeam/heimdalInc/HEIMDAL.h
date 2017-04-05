/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   essBuildInc/HEIMDAL.h
 *
 * Copyright (c) 2004-2017 by Stuart Ansell
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
#ifndef essSystem_HEIMDAL_h
#define essSystem_HEIMDAL_h

namespace attachSystem
{
  class FixedComp;
  class TwinComp;
  class CellMap;
}

namespace instrumentSystem
{
  class CylSample;
}

namespace constructSystem
{  
  class ChopperPit;
  class Cryostat;   
  class DiskChopper;
  class Jaws;
  class JawSet;
  class LineShield;
  class RotaryCollimator;
  class VacuumBox;
  class VacuumPipe;
  class VacuumWindow;
  class ChopperUnit;
  class TwinChopper;
  class HoleShape;
  class CrystalMount;
  class TubeDetBox;  
}

namespace essSystem
{  
  class GuideItem;
  class HeimdalHut;
  class DetectorTank;

  /*!
    \class HEIMDAL
    \version 1.0
    \author S. Ansell
    \date September 2015
    \brief HEIMDAL beamline constructor for the ESS
  */
  
class HEIMDAL : public attachSystem::CopiedComp
{
 private:

  /// Start at [0:Complete / 1:Cave]
  int startPoint;  
  /// Stop at [0:Complete / 1:Mono Wall / 2:Inner Bunker / 3:Outer Bunker ]
  int stopPoint;  

  /// Main Beam Axis [for construction]
  std::shared_ptr<attachSystem::FixedOffset> heimdalAxis;

  /// Elliptic focus in bulkshield [m5]
  std::shared_ptr<beamlineSystem::GuideLine> FocusTA;
  /// Elliptic focus in bulkshield [m5]
  std::shared_ptr<beamlineSystem::GuideLine> FocusCA;

  /// Vac pipe in gamma shutter
  std::shared_ptr<constructSystem::VacuumPipe> VPipeB;
  /// Tapered guide from 5.5 to 6metre
  std::shared_ptr<beamlineSystem::GuideLine> FocusTB;
  /// Tapered guide from 5.5 to 6metre
  std::shared_ptr<beamlineSystem::GuideLine> FocusCB;

  /// 6.0 - 6.5m Vac piper
  std::shared_ptr<constructSystem::VacuumPipe> VPipeC;
  /// Cold guide from 6.0 to 6.5metre
  std::shared_ptr<beamlineSystem::GuideLine> FocusTC;
  /// Thermal guide from 6.0 to 6metre
  std::shared_ptr<beamlineSystem::GuideLine> FocusCC;

  /// First single chopper pair
  std::shared_ptr<constructSystem::ChopperUnit> ChopA;
  /// Top twin disk
  std::shared_ptr<constructSystem::DiskChopper> ADiskOne;
  /// Lower twin disk
  std::shared_ptr<constructSystem::DiskChopper> ADiskTwo;

  
  void setBeamAxis(const FuncDataBase&,const GuideItem&,const bool);

  void buildBunkerUnits(Simulation&,const attachSystem::FixedComp&,
			const long int,const attachSystem::FixedComp&,
			const long int,const int);
  
  void buildOutGuide(Simulation&,const attachSystem::FixedComp&,
		     const long int,const int);
  void buildHut(Simulation&,const attachSystem::FixedComp&,
		const long int,const int);
  void buildDetectorArray(Simulation&,const attachSystem::FixedComp&,
			  const long int,const int);
  
 public:
  
  HEIMDAL(const std::string&);
  HEIMDAL(const HEIMDAL&);
  HEIMDAL& operator=(const HEIMDAL&);
  ~HEIMDAL();
  
  void buildIsolated(Simulation&,const int);
  void build(Simulation&,const GuideItem&,
	     const Bunker&,const int);

};

}

#endif
