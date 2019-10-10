/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File: danmax/danmaxOpticsLine.cxx
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
#include "FixedRotate.h"
#include "FixedGroup.h"
#include "FixedOffsetGroup.h"
#include "ContainedComp.h"
#include "ContainedGroup.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "SurfMap.h"
#include "ExternalCut.h"
#include "InnerZone.h"
#include "FrontBackCut.h"
#include "CopiedComp.h"
#include "World.h"
#include "AttachSupport.h"
#include "ModelSupport.h"
#include "generateSurf.h"

#include "insertObject.h"
#include "insertPlate.h"
#include "VacuumPipe.h"
#include "SplitFlangePipe.h"
#include "Bellows.h"
#include "VacuumBox.h"
#include "portItem.h"
#include "PipeTube.h"
#include "PortTube.h"

#include "CrossPipe.h"
#include "BremBlock.h"
#include "BremMonoColl.h"
#include "MonoVessel.h"
#include "MonoCrystals.h"
#include "GateValveCube.h"
#include "GateValveCylinder.h"
#include "JawUnit.h"
#include "JawFlange.h"
#include "FlangeMount.h"
#include "Mirror.h"
#include "MonoBox.h"
#include "MonoShutter.h"
#include "BeamMount.h"
#include "BeamPair.h"
#include "DCMTank.h"
#include "MonoBlockXstals.h"
#include "MLMono.h"
#include "danmaxOpticsLine.h"

namespace xraySystem
{

// Note currently uncopied:
  
danmaxOpticsLine::danmaxOpticsLine(const std::string& Key) :
  attachSystem::CopiedComp(Key,Key),
  attachSystem::ContainedComp(),
  attachSystem::FixedOffset(newName,2),
  attachSystem::ExternalCut(),
  attachSystem::CellMap(),

  buildZone(*this,cellIndex),
  
  pipeInit(new constructSystem::Bellows(newName+"InitBellow")),
  triggerPipe(new constructSystem::PipeTube(newName+"TriggerPipe")),
  gateTubeA(new constructSystem::PipeTube(newName+"GateTubeA")),
  gateTubeAItem(new xraySystem::FlangeMount(newName+"GateTubeAItem")),
  bellowA(new constructSystem::Bellows(newName+"BellowA")),
  pipeA(new constructSystem::VacuumPipe(newName+"PipeA")),    
  bellowB(new constructSystem::Bellows(newName+"BellowB")),
  collTubeA(new constructSystem::PipeTube(newName+"CollTubeA")),
  bremColl(new xraySystem::BremBlock(newName+"BremColl")),
  filterPipe(new constructSystem::VacuumPipe(newName+"FilterPipe")),
  gateA(new constructSystem::GateValveCylinder(newName+"GateA")),
  bellowC(new constructSystem::Bellows(newName+"BellowC")),
  lauePipe(new constructSystem::VacuumPipe(newName+"LauePipe")),    
  bellowD(new constructSystem::Bellows(newName+"BellowD")),
  slitTube(new constructSystem::PortTube(newName+"SlitTube")),
  jaws({
      std::make_shared<xraySystem::BeamPair>(newName+"JawX"),
      std::make_shared<xraySystem::BeamPair>(newName+"JawZ")
	}),
  gateB(new constructSystem::GateValveCylinder(newName+"GateB")),
  bellowE(new constructSystem::Bellows(newName+"BellowE")),
  monoVessel(new xraySystem::DCMTank(newName+"MonoVessel")),
  mbXstals(new xraySystem::MonoBlockXstals(newName+"MBXstals")),
  gateC(new constructSystem::GateValveCylinder(newName+"GateC")),
  viewTube(new constructSystem::PipeTube(newName+"ViewTube")),
  viewTubeScreen(new xraySystem::FlangeMount(newName+"ViewTubeScreen")),
  bellowF(new constructSystem::Bellows(newName+"BellowF")),
  MLMVessel(new constructSystem::VacuumBox(newName+"MLMVessel")),
  MLM(new xraySystem::MLMono(newName+"MLM"))

  /*!
    Constructor
    \param Key :: Name of construction key
  */
{
  ELog::RegMethod RegA("danmaxOpticsLine","danmaxOpticsLine(constructor)");

  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();
  
  OR.addObject(pipeInit);
  OR.addObject(triggerPipe);
  OR.addObject(gateTubeA);
  OR.addObject(gateTubeAItem);

  OR.addObject(bellowA);
  OR.addObject(pipeA);
  OR.addObject(bellowB);
  OR.addObject(collTubeA);
  OR.addObject(bremColl);
  OR.addObject(filterPipe);
  OR.addObject(gateA);
  OR.addObject(bellowC);
  OR.addObject(lauePipe);
  OR.addObject(bellowD);
  OR.addObject(slitTube);
  OR.addObject(jaws[0]);
  OR.addObject(jaws[1]);
  
  OR.addObject(gateB);
  OR.addObject(bellowE);
  OR.addObject(monoVessel);
  OR.addObject(mbXstals);
  OR.addObject(gateC);
  OR.addObject(viewTube);
  OR.addObject(viewTubeScreen);
  OR.addObject(bellowF);
  
}
  
danmaxOpticsLine::~danmaxOpticsLine()
  /*!
    Destructor
   */
{}

void
danmaxOpticsLine::populate(const FuncDataBase& Control)
  /*!
    Populate the intial values [movement]
   */
{
  FixedOffset::populate(Control);

  outerLeft=Control.EvalDefVar<double>(keyName+"OuterLeft",0.0);
  outerRight=Control.EvalDefVar<double>(keyName+"OuterRight",outerLeft);
  outerTop=Control.EvalDefVar<double>(keyName+"OuterTop",outerLeft);
  
  return;
}


void
danmaxOpticsLine::createSurfaces()
  /*!
    Create surfaces for outer void
    Mainly for the masterCell 
  */
{
  ELog::RegMethod RegA("danmaxOpticsLine","createSurface");

  if (outerLeft>Geometry::zeroTol &&  isActive("floor"))
    {
      std::string Out;
      ModelSupport::buildPlane
	(SMap,buildIndex+3,Origin-X*outerLeft,X);
      ModelSupport::buildPlane
	(SMap,buildIndex+4,Origin+X*outerRight,X);
      ModelSupport::buildPlane
	(SMap,buildIndex+6,Origin+Z*outerTop,Z);
      Out=ModelSupport::getComposite(SMap,buildIndex," 3 -4 -6");
      const HeadRule HR(Out+getRuleStr("floor"));
      buildZone.setSurround(HR);
    }
  return;
}

void
danmaxOpticsLine::constructViewScreen(Simulation& System,
				      MonteCarlo::Object* masterCell,
				      const attachSystem::FixedComp& initFC, 
				      const std::string& sideName)
  /*!
    Sub build of the slit package unit
    \param System :: Simulation to use
    \param masterCell :: Main master volume
    \param initFC :: Start point
    \param sideName :: start link point
  */
{
  ELog::RegMethod RegA("danmaxOpticsLine","constructViewScreen");

  // FAKE insertcell: required
  viewTube->addAllInsertCell(masterCell->getName());
  viewTube->setPortRotation(3,Geometry::Vec3D(1,0,0));
  viewTube->createAll(System,initFC,sideName);
  viewTube->intersectPorts(System,1,2);

  const constructSystem::portItem& VPA=viewTube->getPort(0);
  const constructSystem::portItem& VPB=viewTube->getPort(1);
  const constructSystem::portItem& VPC=viewTube->getPort(2); // screen)
  
  int outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,VPB,VPB.getSideIndex("OuterPlate"));
  const Geometry::Vec3D  Axis=viewTube->getY()*(VPB.getY()+VPC.getY())/2.0;
  this->splitObjectAbsolute(System,1501,outerCell,
			      viewTube->getCentre(),VPB.getY());
  this->splitObjectAbsolute(System,1502,outerCell+1,
			      viewTube->getCentre(),Axis);
  
  const std::vector<int> cellUnit=this->getCells("OuterVoid");
  viewTube->insertMainInCell(System,cellUnit);
  VPA.insertInCell(System,this->getCell("OuterVoid"));
  
  //  viewTube->insertPortInCell(System,{cellN,cellM,cellX});
  viewTube->insertPortInCell
    (System,{{outerCell},{outerCell+1},{outerCell+2}});
  cellIndex+=2;

  
  viewTubeScreen->addInsertCell("Body",viewTube->getCell("Void"));
  viewTubeScreen->addInsertCell("Body",VPC.getCell("Void"));
  viewTubeScreen->setBladeCentre(*viewTube,0);
  viewTubeScreen->createAll(System,VPC,"-InnerPlate");


  bellowF->FrontBackCut::setFront(VPB,"OuterPlate");
  bellowF->createAll(System,VPB,"OuterPlate");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*bellowF,2);
  bellowF->insertInCell(System,outerCell);
  VPC.insertInCell(System,outerCell);
  

  return;
}

void
danmaxOpticsLine::constructMono(Simulation& System,
				MonteCarlo::Object* masterCell,
				const attachSystem::FixedComp& initFC, 
				const std::string& sideName)
  /*!
    Sub build of the slit package unit
    \param System :: Simulation to use
    \param masterCell :: Main master volume
    \param initFC :: Start point
    \param sideName :: start link point
  */
{
  ELog::RegMethod RegA("danmaxOpticsLine","buildMono");

  int outerCell;

  // FAKE insertcell: required
  monoVessel->addInsertCell(masterCell->getName());
  monoVessel->createAll(System,initFC,sideName);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*monoVessel,2);
  monoVessel->insertInCell(System,outerCell);

  
  mbXstals->addInsertCell(monoVessel->getCell("Void"));
  //  mbXstals->copyCutSurf("innerCylinder",*monoVessel,"innerRadius");
  mbXstals->createAll(System,*monoVessel,0);


  return;
}
void
danmaxOpticsLine::constructMirrorMono(Simulation& System,
				      MonteCarlo::Object* masterCell,
				      const attachSystem::FixedComp& initFC, 
				      const std::string& sideName)
  /*!
    Sub build of the slit package unit
    \param System :: Simulation to use
    \param masterCell :: Main master volume
    \param initFC :: Start point
    \param sideName :: start link point
  */
{
  ELog::RegMethod RegA("danmaxOpticsLine","buildMono");

  int outerCell;

  // FAKE insertcell: required
  MLMVessel->addInsertCell(masterCell->getName());
  MLMVessel->createAll(System,initFC,sideName);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*MLMVessel,2);
  MLMVessel->insertInCell(System,outerCell);

  
  MLM->addInsertCell(MLMVessel->getCell("Void"));
  //  MLM->copyCutSurf("innerCylinder",*MLMVessel,"innerRadius");
  MLM->createAll(System,*MLMVessel,0);


  return;
}

void
danmaxOpticsLine::constructSlitTube(Simulation& System,
				    MonteCarlo::Object* masterCell,
				    const attachSystem::FixedComp& initFC, 
				    const std::string& sideName)
  /*!
    Build the DM2 split package
    \param System :: Simuation to use
   */
{
  ELog::RegMethod RegA("danmaxOpticsLine","buildSlitTube");

  int outerCell;
  
  // FAKE insertcell: required
  slitTube->addAllInsertCell(masterCell->getName());
  slitTube->createAll(System,initFC,sideName);
  slitTube->intersectPorts(System,0,1);
  slitTube->intersectPorts(System,1,2);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*slitTube,2);
  slitTube->insertAllInCell(System,outerCell);

  slitTube->splitVoidPorts(System,"SplitVoid",1001,
   			   slitTube->getCell("Void"),
   			   Geometry::Vec3D(0,1,0));


  slitTube->splitObject(System,1501,outerCell,
			Geometry::Vec3D(0,0,0),
			Geometry::Vec3D(0,0,1));
  cellIndex++;  // remember creates an extra cell in  primary


  for(size_t i=0;i<jaws.size();i++)
    {
      const constructSystem::portItem& PI=slitTube->getPort(i);
      jaws[i]->addInsertCell("SupportA",PI.getCell("Void"));
      jaws[i]->addInsertCell("SupportA",slitTube->getCell("SplitVoid",i));
      jaws[i]->addInsertCell("SupportB",PI.getCell("Void"));
      jaws[i]->addInsertCell("SupportB",slitTube->getCell("SplitVoid",i));
      jaws[i]->addInsertCell("BlockA",slitTube->getCell("SplitVoid",i));
      jaws[i]->addInsertCell("BlockB",slitTube->getCell("SplitVoid",i));
      jaws[i]->createAll(System,*slitTube,0,
			 PI,PI.getSideIndex("InnerPlate"));
    }  

  return;
}

  

void
danmaxOpticsLine::buildObjects(Simulation& System)
  /*!
    Build all the objects relative to the main FC
    point.
    \param System :: Simulation to use
  */
{
  ELog::RegMethod RegA("danmaxOpticsLine","buildObjects");

  int outerCell;
  
  buildZone.setFront(getRule("front"));
  buildZone.setBack(getRule("back"));

  MonteCarlo::Object* masterCell=
    buildZone.constructMasterCell(System,*this);


  // dummy space for first item
  // This is a mess but want to preserve insert items already
  // in the hut beam port
  pipeInit->createAll(System,*this,0);
  // dump cell for joinPipe
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*pipeInit,-1);
  // real cell for initPipe
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*pipeInit,2);
  pipeInit->insertInCell(System,outerCell);

  // FAKE insertcell: required due to rotation ::
  triggerPipe->addAllInsertCell(masterCell->getName());
  triggerPipe->setPortRotation(3,Geometry::Vec3D(1,0,0));
  triggerPipe->createAll(System,*pipeInit,2);  

  const constructSystem::portItem& TPI=triggerPipe->getPort(1);
  outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,TPI,TPI.getSideIndex("OuterPlate"));
  triggerPipe->insertAllInCell(System,outerCell);

  // FAKE insertcell: required
  gateTubeA->addAllInsertCell(masterCell->getName());
  gateTubeA->setPortRotation(3,Geometry::Vec3D(1,0,0));
  gateTubeA->createAll(System,TPI,TPI.getSideIndex("OuterPlate"));
  
  const constructSystem::portItem& GPI=gateTubeA->getPort(1);
  outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,GPI,GPI.getSideIndex("OuterPlate"));
  gateTubeA->insertAllInCell(System,outerCell);

  gateTubeAItem->addInsertCell("Body",gateTubeA->getCell("Void"));
  gateTubeAItem->setBladeCentre(*gateTubeA,0);
  gateTubeAItem->createAll(System,*gateTubeA,std::string("InnerBack"));

  // after gate value
  bellowA->setFront(GPI,GPI.getSideIndex("OuterPlate"));
  bellowA->createAll(System,GPI,"OuterPlate");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*bellowA,2);
  bellowA->insertInCell(System,outerCell);

  pipeA->setFront(*bellowA,2);
  pipeA->createAll(System,*bellowA,"back");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*pipeA,2);
  pipeA->insertInCell(System,outerCell);

  bellowB->setFront(*pipeA,2);
  bellowB->createAll(System,*pipeA,"back");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*bellowB,2);
  bellowB->insertInCell(System,outerCell);

  // brem:
  // FAKE insertcell: required
  collTubeA->addAllInsertCell(masterCell->getName());
  collTubeA->setPortRotation(3,Geometry::Vec3D(1,0,0));
  collTubeA->createAll(System,*bellowB,"back");
  
  const constructSystem::portItem& CPI=collTubeA->getPort(1);
  outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,CPI,CPI.getSideIndex("OuterPlate"));
  collTubeA->insertAllInCell(System,outerCell);

  bremColl->addInsertCell(collTubeA->getCell("Void"));
  bremColl->createAll(System,*collTubeA,"OrgOrigin");

  filterPipe->setFront(CPI,CPI.getSideIndex("OuterPlate"));
  filterPipe->createAll(System,CPI,"OuterPlate");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*filterPipe,2);
  filterPipe->insertInCell(System,outerCell);

  gateA->setFront(*filterPipe,filterPipe->getSideIndex("back"));
  gateA->createAll(System,*filterPipe,"back");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*gateA,2);
  gateA->insertInCell(System,outerCell);

    // after gate value
  bellowC->setFront(*gateA,gateA->getSideIndex("back"));
  bellowC->createAll(System,*gateA,"back");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*bellowC,2);
  bellowC->insertInCell(System,outerCell);

  lauePipe->setFront(*bellowC,2);
  lauePipe->createAll(System,*bellowC,"back");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*lauePipe,2);
  lauePipe->insertInCell(System,outerCell);

  bellowD->setFront(*lauePipe,2);
  bellowD->createAll(System,*lauePipe,"back");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*bellowD,2);
  bellowD->insertInCell(System,outerCell);


  constructSlitTube(System,masterCell,*bellowD,"back");


  gateB->setFront(*slitTube,slitTube->getSideIndex("back"));
  gateB->createAll(System,*slitTube,"back");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*gateB,2);
  gateB->insertInCell(System,outerCell);

    // after gate value
  bellowE->setFront(*gateB,gateB->getSideIndex("back"));
  bellowE->createAll(System,*gateB,"back");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*bellowE,2);
  bellowE->insertInCell(System,outerCell);

  constructMono(System,masterCell,*bellowE,"back");

  gateC->setFront(*monoVessel,"back");
  gateC->createAll(System,*monoVessel,"back");
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*gateC,2);
  gateC->insertInCell(System,outerCell);

  constructViewScreen(System,masterCell,*gateC,"back");

  constructMirrorMono(System,masterCell,*bellowF,"back");
  lastComp=bellowF;
  return;

  /*
  setCell("LastVoid",masterCell->getName());
  lastComp=gateJ;
  */
  return;
}

void
danmaxOpticsLine::createLinks()
  /*!
    Create a front/back link
   */
{
  ELog::RegMethod RControl("danmaxOpticsLine","createLinks");
  
  setLinkSignedCopy(0,*pipeInit,1);
  setLinkSignedCopy(1,*lastComp,2);
  return;
}
  
  
void 
danmaxOpticsLine::createAll(Simulation& System,
			    const attachSystem::FixedComp& FC,
			    const long int sideIndex)
  /*!
    Carry out the full build
    \param System :: Simulation system
    \param FC :: Fixed component
    \param sideIndex :: link point
  */
{
  ELog::RegMethod RControl("danmaxOpticsLine","createAll");

  populate(System.getDataBase());
  createUnitVector(FC,sideIndex);
  createSurfaces();
  
  buildObjects(System);
  createLinks();
  return;
}


}   // NAMESPACE xraySystem

