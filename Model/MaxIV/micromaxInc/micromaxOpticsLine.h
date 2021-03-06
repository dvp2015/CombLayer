/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   micromaxInc/micromaxOpticsLine.h
 *
 * Copyright (c) 2004-2019 by Stuart Ansell
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>. 
 *
 ****************************************************************************/
#ifndef xraySystem_micromaxOpticsLineX_h
#define xraySystem_micromaxOpticsLineX_h

namespace insertSystem
{
  class insertPlate;
}

namespace constructSystem
{
  class SupplyPipe;
  class CrossPipe;
  class VacuumPipe;
  class Bellows;
  class VacuumBox;
  class portItem;
  class PortTube;
  class PipeTube;
  class GateValveCube;
  class JawValveCube;
  class JawFlange;

}

namespace xraySystem
{
  class OpticsHutch;
  class BremColl;
  class FlangeMount;
  class Mirror;
  class MonoCrystals;
  class MonoBox;
  class MonoShutter;
  class PipeShield;
  class ShutterUnit;
    
  /*!
    \class micromaxOpticsLine
    \version 1.0
    \author S. Ansell
    \date January 2019
    \brief Constructor for the micromax optics components
  */

class micromaxOpticsLine :
  public attachSystem::CopiedComp,
  public attachSystem::ContainedComp,
  public attachSystem::FixedOffset,
  public attachSystem::ExternalCut,
  public attachSystem::CellMap
{
 private:

  /// construction space for main object
  attachSystem::InnerZone buildZone;
  
  /// Shared point to use for last component:
  std::shared_ptr<attachSystem::FixedComp> lastComp;

  /// Inital bellow
  std::shared_ptr<constructSystem::Bellows> pipeInit;
  /// vauucm trigger system
  std::shared_ptr<constructSystem::CrossPipe> triggerPipe;
  /// first ion pump
  std::shared_ptr<constructSystem::CrossPipe> gaugeA;
  /// bellows after ion pump to filter
  std::shared_ptr<constructSystem::Bellows> bellowA;
  /// Vacuum pipe for collimator
  std::shared_ptr<xraySystem::BremColl> bremCollA;
  /// Filter tube
  std::shared_ptr<constructSystem::PortTube> filterBoxA;
  /// Filter stick [only one blade type -- fix]
  std::shared_ptr<xraySystem::FlangeMount> filterStick;
  /// First gate valve
  std::shared_ptr<constructSystem::GateValveCube> gateA;
  /// bellows after gateA ->view
  std::shared_ptr<constructSystem::Bellows> bellowB;
  
  /// diamon screen(?)
  std::shared_ptr<constructSystem::PipeTube> screenPipeA;
  /// View/something(?)
  std::shared_ptr<constructSystem::PipeTube> screenPipeB;
  /// Primary jaw (Box)
  std::shared_ptr<constructSystem::VacuumBox> primeJawBox;
  /// Bellow to gate on mono
  std::shared_ptr<constructSystem::Bellows> bellowC;
  /// First gate valve
  std::shared_ptr<constructSystem::GateValveCube> gateB;
  /// Mono box
  std::shared_ptr<xraySystem::MonoBox> monoBox;
  /// Mono Xstal 
  std::shared_ptr<xraySystem::MonoCrystals> monoXtal;
  // Gate to isolate mono
  std::shared_ptr<constructSystem::GateValveCube> gateC;
  /// Bellow to diagnositics
  std::shared_ptr<constructSystem::Bellows> bellowD;
  /// Diagnostic unit 1:
  std::shared_ptr<constructSystem::PortTube> diagBoxA;
  /// Bellow from diagnositics
  std::shared_ptr<constructSystem::Bellows> bellowE;
  // Gate fro first mirror
  std::shared_ptr<constructSystem::GateValveCube> gateD;

  /// Mirror box 
  std::shared_ptr<constructSystem::VacuumBox> mirrorA;
  // Gate fro first mirror
  std::shared_ptr<constructSystem::GateValveCube> gateE;
  /// Bellow to diagnositics
  std::shared_ptr<constructSystem::Bellows> bellowF;
  /// Diagnostic unit 2:
  std::shared_ptr<constructSystem::PortTube> diagBoxB;
  /// Diag Box B :: Jaw units
  std::array<std::shared_ptr<constructSystem::JawFlange>,2> jawCompB;

  /// Bellow to mirror B
  std::shared_ptr<constructSystem::Bellows> bellowG;
  // Gate valve
  std::shared_ptr<constructSystem::GateValveCube> gateF;
  /// Mirror box B
  std::shared_ptr<constructSystem::VacuumBox> mirrorB;
  // Gate valve
  std::shared_ptr<constructSystem::GateValveCube> gateG;
  /// Bellow to mirror B
  std::shared_ptr<constructSystem::Bellows> bellowH;

  /// Diagnostic unit 3:
  std::shared_ptr<constructSystem::PortTube> diagBoxC;
  /// Diag Box C :: Jaw units
  std::array<std::shared_ptr<constructSystem::JawFlange>,2> jawCompC;
  /// Monochromatic shutter
  std::shared_ptr<xraySystem::MonoShutter> monoShutter;
  /// Screen at end of hut
  std::shared_ptr<xraySystem::PipeShield> screenA;

  /// Bellow to end station
  std::shared_ptr<constructSystem::Bellows> bellowI;

  double outerLeft;    ///< Left Width for cut rectangle
  double outerRight;   ///< Right width for cut rectangle
  double outerTop;     ///< Top lift for cut rectangle
  
  int constructDiag
    (Simulation&,
     MonteCarlo::Object**,
     constructSystem::PortTube&,
     std::array<std::shared_ptr<constructSystem::JawFlange>,2>&,
     const attachSystem::FixedComp&,
     const long int);
  
  void populate(const FuncDataBase&);
  void createSurfaces();
  void buildObjects(Simulation&);
  void createLinks();
  
 public:
  
  micromaxOpticsLine(const std::string&);
  micromaxOpticsLine(const micromaxOpticsLine&);
  micromaxOpticsLine& operator=(const micromaxOpticsLine&);
  ~micromaxOpticsLine();
  
  void createAll(Simulation&,const attachSystem::FixedComp&,
		 const long int);

};

}

#endif
