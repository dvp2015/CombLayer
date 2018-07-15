/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File: balder/BALDER.cxx
 *
 * Copyright (c) 2004-2018 by Stuart Ansell
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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <iterator>
#include <memory>

#include "Exception.h"
#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "GTKreport.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "inputParam.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "Rules.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "Qhull.h"
#include "Simulation.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "ContainedComp.h"
#include "SpaceCut.h"
#include "ContainedSpace.h"
#include "ContainedGroup.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "SurfMap.h"
#include "FrontBackCut.h"
#include "CopiedComp.h"
#include "ExternalCut.h"
#include "World.h"
#include "AttachSupport.h"

#include "VacuumPipe.h"
#include "SplitFlangePipe.h"
#include "Bellows.h"
#include "LeadPipe.h"
#include "VacuumBox.h"
#include "portItem.h"
#include "PortTube.h"

#include "OpticsHutch.h"
#include "ExperimentalHutch.h"
#include "CrossPipe.h"
#include "MonoVessel.h"
#include "MonoCrystals.h"
#include "GateValve.h"
#include "JawUnit.h"
#include "JawValve.h"
#include "FlangeMount.h"
#include "FrontEndCave.h"
#include "FrontEnd.h"
#include "OpticsBeamline.h"
#include "ConnectZone.h"
#include "PipeShield.h"
#include "BALDER.h"

namespace xraySystem
{

BALDER::BALDER(const std::string& KN) :
  attachSystem::CopiedComp("Balder",KN),
  ringCaveA(new FrontEndCave(newName+"RingCaveA")),
  ringCaveB(new FrontEndCave(newName+"RingCaveB")),
  frontBeam(new FrontEnd(newName+"FrontBeam")),
  joinPipe(new constructSystem::VacuumPipe(newName+"JoinPipe")),
  opticsHut(new OpticsHutch(newName+"OpticsHut")),
  opticsBeam(new OpticsBeamline(newName+"OpticsLine")),
  joinPipeB(new constructSystem::LeadPipe(newName+"JoinPipeB")),
  pShield(new xraySystem::PipeShield(newName+"PShield")),
  nShield(new xraySystem::PipeShield(newName+"NShield")),
  outerShield(new xraySystem::PipeShield(newName+"OuterShield")),
  connectZone(new ConnectZone(newName+"Connect")),
  joinPipeC(new constructSystem::LeadPipe(newName+"JoinPipeC")),
  exptHut(new ExperimentalHutch(newName+"Expt"))
  /*!
    Constructor
    \param KN :: Keyname
  */
{
  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();
  
  OR.addObject(ringCaveA);
  OR.addObject(ringCaveB);
  OR.addObject(frontBeam);
  OR.addObject(joinPipe);
  
  OR.addObject(opticsHut);
  OR.addObject(opticsBeam);
  OR.addObject(joinPipeB);
  OR.addObject(pShield);
  OR.addObject(nShield);
  OR.addObject(outerShield);
  OR.addObject(joinPipeC);
  OR.addObject(exptHut);
}

BALDER::~BALDER()
  /*!
    Destructor
   */
{}

void 
BALDER::build(Simulation& System,
	      const attachSystem::FixedComp& FCOrigin,
	      const long int sideIndex)
  /*!
    Carry out the full build
    \param System :: Simulation system
    \param FCOrigin :: Start origin
    \param sideIndex :: link point for origin
   */
{
  // For output stream
  ELog::RegMethod RControl("BALDER","build");

  int voidCell(74123);
 
  ringCaveA->addInsertCell(voidCell);
  ringCaveA->createAll(System,FCOrigin,sideIndex);

  ringCaveB->addInsertCell(voidCell);
  ringCaveB->setCutSurf("front",*ringCaveA,"connectPt");
  ringCaveB->createAll(System,*ringCaveA,
		       ringCaveA->getSideIndex("connectPt"));


  const HeadRule caveVoid=ringCaveA->getCellHR(System,"Void");
  frontBeam->addInsertCell(ringCaveA->getCell("Void"));
  frontBeam->createAll(System,*ringCaveA,-1);

  if (stopPoint=="frontEnd") return;

  opticsHut->addInsertCell(voidCell);
  opticsHut->setCutSurf("ringWall",*ringCaveB,"outerWall");
  opticsHut->createAll(System,*ringCaveA,2);

  // UgLY HACK to get the two objects to merge
  const std::string OH=opticsHut->SurfMap::getSurfString("ringFlat");
  ringCaveB->insertComponent
    (System,"OuterWall",*opticsHut,opticsHut->getSideIndex("frontCut"));
  ringCaveB->insertComponent
    (System,"FloorA",*opticsHut,opticsHut->getSideIndex("floorCut"));
  ringCaveB->insertComponent
    (System,"RoofA",*opticsHut,opticsHut->getSideIndex("roofCut"));


  joinPipe->addInsertCell(ringCaveA->getCell("Void"));
  joinPipe->addInsertCell(ringCaveA->getCell("FrontWallHole"));
  joinPipe->addInsertCell(opticsHut->getCell("Inlet"));
  joinPipe->addInsertCell(opticsHut->getCell("Void"));

  joinPipe->setPrimaryCell(opticsHut->getCell("Void"));
  joinPipe->setFront(*frontBeam,2);
  joinPipe->setSpaceLinkCopy(0,*opticsHut,-12);
  joinPipe->registerSpaceCut(0,2);
  joinPipe->createAll(System,*frontBeam,2);


  joinPipe->clear();
  joinPipe->setPrimaryCell(caveVoid);
  joinPipe->registerSpaceCut(1,0);
  joinPipe->insertObjects(System);

  System.removeCell(ringCaveA->getCell("Void"));

  opticsBeam->addInsertCell(opticsHut->getCell("Void"));
  opticsBeam->createAll(System,*joinPipe,2);

  joinPipeB->addInsertCell(opticsHut->getCell("ExitHole"));
  joinPipeB->setPrimaryCell(opticsHut->getCell("Void"));
  joinPipeB->setFront(*opticsBeam,2);
  joinPipeB->setSpaceLinkCopy(1,*opticsHut,
			      opticsHut->getSideIndex("-innerBack"));
  joinPipeB->registerSpaceCut(1,0);
  joinPipeB->createAll(System,*opticsBeam,2);

  // pipe shield goes around joinPipeB:
  pShield->addInsertCell(joinPipeB->getCell("OuterSpace"));
  pShield->setCutSurf("inner",*joinPipeB,"outerPipe");
  pShield->setCutSurf("front",*opticsHut,"innerBack");
  pShield->createAll(System,*opticsHut,opticsHut->getSideIndex("exitHole"));

  // pipe shield goes around joinPipeB:
  nShield->addInsertCell(joinPipeB->getCell("OuterSpace"));
  nShield->setCutSurf("inner",*joinPipeB,"outerPipe");
  nShield->createAll(System,*opticsHut,opticsHut->getSideIndex("exitHole"));

  System.removeCell(opticsHut->getCell("Void"));

  exptHut->addInsertCell(voidCell);
  exptHut->createAll(System,*ringCaveA,2);

  connectZone->addInsertCell(voidCell);
  connectZone->setFront(*opticsHut,2);
  connectZone->setBack(*exptHut,1);
  connectZone->createAll(System,*joinPipeB,2);

  // horrid way to create a SECOND register space [MAKE INTERNAL]
  joinPipeB->clear();  // reinitialize
  joinPipeB->setSpaceLinkCopy(0,*opticsHut,-2);
  //  joinPipeB->registerSpaceIsolation(0,2);
  joinPipeB->registerSpaceCut(0,2);
  joinPipeB->setPrimaryCell(connectZone->getCell("OuterVoid"));
  joinPipeB->insertObjects(System);

  connectZone->insertFirstRegion
    (System,joinPipeB->getCell("OuterSpace",1));
  

  
  // pipe shield goes around joinPipeB:

  outerShield->addInsertCell(joinPipeB->getCell("OuterSpace",1));
  outerShield->setCutSurf("inner",*joinPipeB,"outerPipe");
  outerShield->setCutSurf("front",*opticsHut,"back");
  outerShield->createAll(System,*opticsHut,
			 opticsHut->getSideIndex("-exitHole"));

  joinPipeC->addInsertCell(connectZone->getCell("OuterVoid"));
  joinPipeC->addInsertCell(exptHut->getCell("Void"));
  joinPipeC->addInsertCell(exptHut->getCell("EnteranceHole"));
  joinPipeC->setFront(*connectZone,2);
  joinPipeC->registerSpaceCut(1,0);
  joinPipeC->setSpaceLinkCopy(1,*exptHut,-1);
  joinPipeC->setPrimaryCell(connectZone->getCell("OuterVoid"));
  joinPipeC->createAll(System,*connectZone,2);

  System.removeCell(connectZone->getCell("OuterVoid"));

  return;
}


}   // NAMESPACE xraySystem

