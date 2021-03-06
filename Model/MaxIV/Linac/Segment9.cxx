
/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File: Linac/Segment9.cxx
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
#include "groupRange.h"
#include "objectGroups.h"
#include "Simulation.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "FixedGroup.h"
#include "FixedRotate.h"
#include "ContainedComp.h"
#include "ContainedGroup.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "SurfMap.h"
#include "ExternalCut.h"
#include "FrontBackCut.h"
#include "InnerZone.h"
#include "AttachSupport.h"
#include "generateSurf.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generalConstruct.h"

#include "VacuumPipe.h"
#include "SplitFlangePipe.h"
#include "portItem.h"
#include "VirtualTube.h"
#include "BlankTube.h"
#include "Bellows.h"
#include "LQuadF.h"
#include "BPM.h"
#include "CorrectorMag.h"
#include "CeramicSep.h"

#include "LObjectSupport.h"
#include "TDCsegment.h"
#include "Segment9.h"

namespace tdcSystem
{

// Note currently uncopied:

  
Segment9::Segment9(const std::string& Key) :
  TDCsegment(Key,2),

  ceramicBellowA(new tdcSystem::CeramicSep(keyName+"CeramicBellowA")),
  pumpA(new constructSystem::BlankTube(keyName+"PumpA")),
  pipeA(new constructSystem::VacuumPipe(keyName+"PipeA")),
  cMagVertA(new tdcSystem::CorrectorMag(keyName+"CMagVertA")),
  cMagHorA(new tdcSystem::CorrectorMag(keyName+"CMagHorA")),

  bellowB(new constructSystem::Bellows(keyName+"BellowB")),
  bpm(new tdcSystem::BPM(keyName+"BPM")),
  pipeB(new constructSystem::VacuumPipe(keyName+"PipeB")),
  QuadA(new tdcSystem::LQuadF(keyName+"QuadA")),
  
  bellowC(new constructSystem::Bellows(keyName+"BellowC"))

  /*!
    Constructor
    \param Key :: Name of construction key
  */
{
  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();

  OR.addObject(ceramicBellowA);
  OR.addObject(pumpA);
  OR.addObject(pipeA);
  OR.addObject(cMagVertA);
  OR.addObject(cMagHorA);
  OR.addObject(bellowB);
  OR.addObject(bpm);
  OR.addObject(pipeB);
  OR.addObject(QuadA);
  OR.addObject(bellowC);

  setFirstItems(ceramicBellowA);
}
  
Segment9::~Segment9()
  /*!
    Destructor
   */
{}

void
Segment9::buildObjects(Simulation& System)
  /*!
    Build all the objects relative to the main FC
    point.
    \param System :: Simulation to use
  */
{
  ELog::RegMethod RegA("Segment9","buildObjects");

  int outerCell;

  MonteCarlo::Object* masterCell=buildZone->getMaster();
  if (!masterCell)
    masterCell=buildZone->constructMasterCell(System);

  if (isActive("front"))
    ceramicBellowA->copyCutSurf("front",*this,"front");
  ceramicBellowA->createAll(System,*this,0);
  outerCell=buildZone->createOuterVoidUnit(System,masterCell,*ceramicBellowA,2);
  ceramicBellowA->insertInCell(System,outerCell);


  // FAKE INSERT REQUIRED
  pumpA->addAllInsertCell(masterCell->getName());
  pumpA->setPortRotation(3,Geometry::Vec3D(1,0,0));
  pumpA->createAll(System,*ceramicBellowA,"back");

  const constructSystem::portItem& VPB=pumpA->getPort(1);
  outerCell=buildZone->createOuterVoidUnit
    (System,masterCell,VPB,VPB.getSideIndex("OuterPlate"));
  pumpA->insertAllInCell(System,outerCell);

  pipeA->createAll(System,VPB,"OuterPlate");
  pipeMagUnit(System,*buildZone,pipeA,"#front","outerPipe",cMagVertA);
  pipeMagUnit(System,*buildZone,pipeA,"#front","outerPipe",cMagHorA);
  pipeTerminate(System,*buildZone,pipeA);
  
  constructSystem::constructUnit
    (System,*buildZone,masterCell,*pipeA,"back",*bellowB);

  constructSystem::constructUnit
    (System,*buildZone,masterCell,*bellowB,"back",*bpm);
  
  pipeB->createAll(System,*bpm,"back");
  pipeMagUnit(System,*buildZone,pipeB,"#front","outerPipe",QuadA);
  pipeTerminate(System,*buildZone,pipeB);

  constructSystem::constructUnit
    (System,*buildZone,masterCell,*pipeB,"back",*bellowC);


  buildZone->removeLastMaster(System);  
  return;
}

void
Segment9::createLinks()
  /*!
    Create a front/back link
   */
{
  setLinkSignedCopy(0,*ceramicBellowA,1);
  setLinkSignedCopy(1,*bellowC,2);

  joinItems.push_back(FixedComp::getFullRule(2));
  return;
}

void 
Segment9::createAll(Simulation& System,
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
  ELog::RegMethod RControl("Segment9","createAll");

  FixedRotate::populate(System.getDataBase());
  createUnitVector(FC,sideIndex);
  buildObjects(System);
  createLinks();
  
  return;
}


}   // NAMESPACE tdcSystem

