/*********************************************************************
  CombLayer : MCNP(X) Input builder

 * File: Linac/Segment38.cxx
 *
 * Copyright (c) 2004-2020 by Konstantin Batkov
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

#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "Vec3D.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "groupRange.h"
#include "objectGroups.h"
#include "Simulation.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "FixedRotate.h"
#include "ContainedComp.h"
#include "ContainedGroup.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "SurfMap.h"
#include "ExternalCut.h"
#include "FrontBackCut.h"
#include "InnerZone.h"
#include "generalConstruct.h"

#include "SplitFlangePipe.h"
#include "Bellows.h"
#include "VacuumPipe.h"
#include "portItem.h"
#include "VirtualTube.h"
#include "BlankTube.h"

#include "TDCsegment.h"
#include "Segment38.h"

namespace tdcSystem
{

// Note currently uncopied:

Segment38::Segment38(const std::string& Key) :
  TDCsegment(Key,2),
  pipeA(new constructSystem::VacuumPipe(keyName+"PipeA")),
  ionPump(new constructSystem::BlankTube(keyName+"IonPump")),
  pipeB(new constructSystem::VacuumPipe(keyName+"PipeB")),
  bellowA(new constructSystem::Bellows(keyName+"BellowA"))
  /*!
    Constructor
    \param Key :: Name of construction key
  */
{
  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();

  OR.addObject(pipeA);
  OR.addObject(ionPump);
  OR.addObject(pipeB);
  OR.addObject(bellowA);

  setFirstItems(pipeA);
}

Segment38::~Segment38()
  /*!
    Destructor
   */
{}

void
Segment38::buildObjects(Simulation& System)
  /*!
    Build all the objects relative to the main FC
    point.
    \param System :: Simulation to use
  */
{
  ELog::RegMethod RegA("Segment38","buildObjects");

  int outerCell;
  MonteCarlo::Object* masterCell=buildZone->getMaster();

  pipeA->createAll(System,*this,0);
  if (!masterCell)
    masterCell=buildZone->constructMasterCell(System,*pipeA,-1);
  outerCell=buildZone->createOuterVoidUnit(System,masterCell,*pipeA,2);
  pipeA->insertInCell(System,outerCell);

  const constructSystem::portItem& ionPumpBackPort =
    buildIonPump2Port(System,*buildZone,masterCell,*pipeA,"back",*ionPump);

  constructSystem::constructUnit
    (System,*buildZone,masterCell,ionPumpBackPort,"OuterPlate",*pipeB);

  constructSystem::constructUnit
    (System,*buildZone,masterCell,*pipeB,"back",*bellowA);

  buildZone->removeLastMaster(System);

  return;
}

void
Segment38::createLinks()
  /*!
    Create a front/back link
   */
{
  ELog::RegMethod RegA("Segment38","createLinks");

  setLinkSignedCopy(0,*pipeA,1);
  setLinkSignedCopy(1,*bellowA,2);

  joinItems.push_back(FixedComp::getFullRule(2));
  return;
}

void
Segment38::createAll(Simulation& System,
		       const attachSystem::FixedComp& FC,
		       const long int sideIndex)
  /*!
    Carry out the full build
    \param System :: Simulation system
    \param FC :: Fixed component
    \param sideIndex :: link point
   */
{
  // For output stream
  ELog::RegMethod RControl("Segment38","build");

  FixedRotate::populate(System.getDataBase());
  createUnitVector(FC,sideIndex);

  buildObjects(System);
  createLinks();
  return;
}


}   // NAMESPACE tdcSystem
