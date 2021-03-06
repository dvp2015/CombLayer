/*********************************************************************
  CombLayer : MCNP(X) Input builder

 * File: softimax/softimaxOpticsLine.cxx
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
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "Vec3D.h"
#include "surfRegister.h"
#include "objectRegister.h"
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
#include "ModelSupport.h"
#include "generateSurf.h"
#include "generalConstruct.h"
#include "insertObject.h"
#include "insertPlate.h"

#include "VacuumPipe.h"
#include "SplitFlangePipe.h"
#include "OffsetFlangePipe.h"
#include "Bellows.h"
#include "FlangeMount.h"
#include "BremOpticsColl.h"
#include "portItem.h"
#include "VirtualTube.h"
#include "PipeTube.h"
#include "PortTube.h"

#include "BlockStand.h"
#include "GateValveCube.h"
#include "JawUnit.h"
#include "JawFlange.h"
#include "JawValveBase.h"
#include "JawValveTube.h"
#include "TankMonoVessel.h"
#include "GratingUnit.h"
#include "Mirror.h"
#include "BeamPair.h"
#include "TwinPipe.h"
#include "BiPortTube.h"
#include "PipeShield.h"
#include "softimaxOpticsLine.h"

namespace xraySystem
{

// Note currently uncopied:

softimaxOpticsLine::softimaxOpticsLine(const std::string& Key) :
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
  pumpM1(new constructSystem::PipeTube(newName+"PumpM1")),
  gateA(new constructSystem::GateValveCube(newName+"GateA")),
  bellowB(new constructSystem::Bellows(newName+"BellowB")),
  M1TubeFront(new constructSystem::OffsetFlangePipe(newName+"M1TubeFront")),
  M1Tube(new constructSystem::PipeTube(newName+"M1Tube")),
  M1TubeBack(new constructSystem::OffsetFlangePipe(newName+"M1TubeBack")),
  M1Mirror(new xraySystem::Mirror(newName+"M1Mirror")),
  M1Stand(new xraySystem::BlockStand(newName+"M1Stand")),
  bellowC(new constructSystem::Bellows(newName+"BellowC")),
  pumpTubeA(new constructSystem::PipeTube(newName+"PumpTubeA")),
  bremCollA(new xraySystem::BremOpticsColl(newName+"BremCollA")),
  gateB(new constructSystem::GateValveCube(newName+"GateB")),
  bellowD(new constructSystem::Bellows(newName+"BellowD")),
  slitTube(new constructSystem::PortTube(newName+"SlitTube")),
  jaws({
	std::make_shared<xraySystem::BeamPair>(newName+"JawX"),
	std::make_shared<xraySystem::BeamPair>(newName+"JawZ")
    }),
  monoVessel(new xraySystem::TankMonoVessel(newName+"MonoVessel")),
  grating(new xraySystem::GratingUnit(newName+"Grating")),
  zeroOrderBlock(new xraySystem::FlangeMount(newName+"ZeroOrderBlock")),
  gateC(new constructSystem::GateValveCube(newName+"GateC")),
  bellowE(new constructSystem::Bellows(newName+"BellowE")),
  pumpTubeB(new constructSystem::PipeTube(newName+"PumpTubeB")),
  gateD(new constructSystem::GateValveCube(newName+"GateD")),
  joinPipeA(new constructSystem::VacuumPipe(newName+"JoinPipeA")),
  bellowF(new constructSystem::Bellows(newName+"BellowF")),
  slitsA(new constructSystem::JawValveTube(newName+"SlitsA")),
  pumpTubeM3(new constructSystem::PipeTube(newName+"PumpTubeM3")),
  pumpTubeM3Baffle(new xraySystem::FlangeMount(newName+"PumpTubeM3Baffle")),
  bellowG(new constructSystem::Bellows(newName+"BellowG")),
  M3Front(new constructSystem::VacuumPipe(newName+"M3Front")),
  M3Tube(new constructSystem::PipeTube(newName+"M3Tube")),
  M3Mirror(new xraySystem::Mirror(newName+"M3Mirror")),
  M3Stand(new xraySystem::BlockStand(newName+"M3Stand")),
  M3Back(new constructSystem::VacuumPipe(newName+"M3Back")),
  bellowH(new constructSystem::Bellows(newName+"BellowH")),
  gateE(new constructSystem::GateValveCube(newName+"GateE")),
  joinPipeB(new constructSystem::VacuumPipe(newName+"JoinPipeB")),
  pumpTubeC(new constructSystem::PipeTube(newName+"PumpTubeC")),
  bellowI(new constructSystem::Bellows(newName+"BellowI")),
  joinPipeC(new constructSystem::VacuumPipe(newName+"JoinPipeC")),
  gateF(new constructSystem::GateValveCube(newName+"GateF")),
  bellowJ(new constructSystem::Bellows(newName+"BellowJ")),
  M3STXMFront(new constructSystem::VacuumPipe(newName+"M3STXMFront")),
  M3STXMTube(new constructSystem::PipeTube(newName+"M3STXMTube")),
  splitter(new xraySystem::TwinPipe(newName+"Splitter")),
  bellowAA(new constructSystem::Bellows(newName+"BellowAA")),
  bellowBA(new constructSystem::Bellows(newName+"BellowBA")),
  M3Pump(new constructSystem::BiPortTube(newName+"M3Pump")),
  bellowAB(new constructSystem::Bellows(newName+"BellowAB")),
  joinPipeAA(new constructSystem::VacuumPipe(newName+"JoinPipeAA")),
  bremCollAA(new xraySystem::BremOpticsColl(newName+"BremCollAA")),
  joinPipeAB(new constructSystem::VacuumPipe(newName+"JoinPipeAB")),
  bellowBB(new constructSystem::Bellows(newName+"BellowBB")),
  joinPipeBA(new constructSystem::VacuumPipe(newName+"JoinPipeBA")),
  bremCollBA(new xraySystem::BremOpticsColl(newName+"BremCollBA")),
  joinPipeBB(new constructSystem::VacuumPipe(newName+"JoinPipeBB")),

  screenA(new xraySystem::PipeShield(newName+"ScreenA")),
  innerScreen(new insertSystem::insertPlate(newName+"InnerScreen"))

  /*!
    Constructor
    \param Key :: Name of construction key
  */
{
  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();

  OR.addObject(pipeInit);
  OR.addObject(triggerPipe);
  OR.addObject(gateTubeA);
  OR.addObject(gateTubeAItem);
  OR.addObject(bellowA);
  OR.addObject(pipeA);
  OR.addObject(pumpM1);
  OR.addObject(gateA);
  OR.addObject(bellowB);
  OR.addObject(M1TubeFront);
  OR.addObject(M1Tube);
  OR.addObject(M1TubeBack);
  OR.addObject(M1Mirror);
  OR.addObject(M1Stand);
  OR.addObject(bellowC);
  OR.addObject(pumpTubeA);
  OR.addObject(bremCollA);
  OR.addObject(gateB);
  OR.addObject(bellowD);
  OR.addObject(slitTube);
  OR.addObject(jaws[0]);
  OR.addObject(jaws[1]);
  OR.addObject(monoVessel);
  OR.addObject(grating);
  OR.addObject(zeroOrderBlock);
  OR.addObject(gateC);
  OR.addObject(bellowE);
  OR.addObject(pumpTubeB);
  OR.addObject(gateD);
  OR.addObject(joinPipeA);
  OR.addObject(bellowF);
  OR.addObject(slitsA);
  OR.addObject(pumpTubeM3);
  OR.addObject(pumpTubeM3Baffle);
  OR.addObject(bellowG);
  OR.addObject(M3Front);
  OR.addObject(M3Tube);
  OR.addObject(M3Mirror);
  OR.addObject(M3Stand);
  OR.addObject(M3Back);
  OR.addObject(bellowH);
  OR.addObject(gateE);
  OR.addObject(joinPipeB);
  OR.addObject(pumpTubeC);
  OR.addObject(joinPipeC);
  OR.addObject(gateF);
  OR.addObject(bellowJ);
  OR.addObject(M3STXMFront);
  OR.addObject(M3STXMTube);
  OR.addObject(splitter);
  OR.addObject(bellowAA);
  OR.addObject(bellowBA);
  OR.addObject(M3Pump);
  OR.addObject(bellowAB);
  OR.addObject(joinPipeAA);
  OR.addObject(bremCollAA);
  OR.addObject(joinPipeAB);
  OR.addObject(bellowBB);
  OR.addObject(joinPipeBA);
  OR.addObject(bremCollBA);
  OR.addObject(joinPipeBB);
  OR.addObject(screenA);
  OR.addObject(innerScreen);
}

softimaxOpticsLine::~softimaxOpticsLine()
  /*!
    Destructor
   */
{}

void
softimaxOpticsLine::populate(const FuncDataBase& Control)
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
softimaxOpticsLine::createSurfaces()
  /*!
    Create surfaces for outer void
  */
{
  ELog::RegMethod RegA("softimaxOpticsLine","createSurface");

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

int
softimaxOpticsLine::constructDiag
  (Simulation& System,
   MonteCarlo::Object** masterCellPtr,
   constructSystem::PortTube& diagBoxItem,
   std::array<std::shared_ptr<constructSystem::JawFlange>,2>& jawComp,
   const attachSystem::FixedComp& FC,const long int linkPt)
/*!
    Construct a diagnostic box
    \param System :: Simulation for building
    \param diagBoxItem :: Diagnostic box item
    \param jawComp :: Jaw componets to build in diagnostic box
    \param FC :: FixedComp for start point
    \param linkPt :: side index
    \return outerCell
   */
{
  ELog::RegMethod RegA("softimaxOpticsLine","constructDiag");

  int outerCell;

  // fake insert

  diagBoxItem.addAllInsertCell((*masterCellPtr)->getName());
  diagBoxItem.setFront(FC,linkPt);
  diagBoxItem.createAll(System,FC,linkPt);
  outerCell=buildZone.createOuterVoidUnit(System,*masterCellPtr,diagBoxItem,2);
  diagBoxItem.insertAllInCell(System,outerCell);


  for(size_t index=0;index<2;index++)
    {
      const constructSystem::portItem& DPI=diagBoxItem.getPort(index);
      jawComp[index]->setFillRadius
	(DPI,DPI.getSideIndex("InnerRadius"),DPI.getCell("Void"));

      jawComp[index]->addInsertCell(diagBoxItem.getCell("Void"));
      if (index)
	jawComp[index]->addInsertCell(jawComp[index-1]->getCell("Void"));
      jawComp[index]->createAll
	(System,DPI,DPI.getSideIndex("InnerPlate"),diagBoxItem,0);
    }

  diagBoxItem.splitVoidPorts(System,"SplitOuter",2001,
			     diagBoxItem.getCell("Void"),{0,2});
  diagBoxItem.splitObject(System,-11,outerCell);
  diagBoxItem.splitObject(System,12,outerCell);
  diagBoxItem.splitObject(System,2001,outerCell);
  cellIndex+=3;

  return outerCell;
}

void
softimaxOpticsLine::buildM1Mirror(Simulation& System,
				  MonteCarlo::Object* masterCell,
				  const attachSystem::FixedComp& initFC,
				  const std::string& side)
  /*!
    Sub build of the m1-mirror package
    \param System :: Simulation to use
    \param masterCell :: Main master volume
    \param initFC :: Start point
    \param sideIndex :: start link point
  */
{
  ELog::RegMethod RegA("softimaxOpticsBeamline","buildM1Mirror");

  int outerCell;

  constructSystem::constructUnit
    (System,buildZone,masterCell,initFC,side,*M1TubeFront);

  M1Tube->setFront(*M1TubeFront,11);
  M1Tube->createAll(System,*M1TubeFront,11);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*M1Tube,2);
  M1Tube->insertAllInCell(System,outerCell);

  M1Mirror->addInsertCell(M1Tube->getCell("Void"));
  M1Mirror->createAll(System,*M1Tube,0);

  M1Stand->setCutSurf("floor",this->getRule("floor"));
  M1Stand->setCutSurf("front",*M1Tube,-1);
  M1Stand->setCutSurf("back",*M1Tube,-2);
  M1Stand->addInsertCell(outerCell);
  M1Stand->createAll(System,*M1Tube,0);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*M1Tube,"back",*M1TubeBack);

  return;
}

void
softimaxOpticsLine::buildM3Mirror(Simulation& System,
				  MonteCarlo::Object* masterCell,
				  const attachSystem::FixedComp& initFC,
				  const std::string& side)
  /*!
    Sub build of the m1-mirror package
    \param System :: Simulation to use
    \param masterCell :: Main master volume
    \param initFC :: Start point
    \param side :: start link point
  */
{
  ELog::RegMethod RegA("softimaxOpticsBeamline","buildM3Mirror");

  int outerCell;

  constructSystem::constructUnit
    (System,buildZone,masterCell,initFC,side,*M3Front);

  M3Tube->setFront(*M3Front,2);
  M3Tube->createAll(System,*M3Front,2);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*M3Tube,2);
  M3Tube->insertAllInCell(System,outerCell);

  M3Mirror->addInsertCell(M3Tube->getCell("Void"));
  M3Mirror->createAll(System,*M3Tube,0);

  M3Stand->setCutSurf("floor",this->getRule("floor"));
  M3Stand->setCutSurf("front",*M3Tube,-1);
  M3Stand->setCutSurf("back",*M3Tube,-2);
  M3Stand->addInsertCell(outerCell);
  M3Stand->createAll(System,*M3Tube,0);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*M3Tube,"back",*M3Back);

  return;
}

void
softimaxOpticsLine::buildM3STXMMirror(Simulation& System,
				      MonteCarlo::Object* masterCell,
				      const attachSystem::FixedComp& initFC,
				      const std::string& side)
  /*!
    Sub build of the m1-mirror package
    \param System :: Simulation to use
    \param masterCell :: Main master volume
    \param initFC :: Start point
    \param sideIndex :: start link point
  */
{
  ELog::RegMethod RegA("softimaxOpticsBeamline","buildM3STXMMirror");

  int outerCell;

  constructSystem::constructUnit
    (System,buildZone,masterCell,initFC,side,*M3STXMFront);

  M3STXMTube->setFront(*M3STXMFront,2);
  M3STXMTube->createAll(System,*M3STXMFront,2);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*M3STXMTube,2);
  M3STXMTube->insertAllInCell(System,outerCell);

  return;
}

void
softimaxOpticsLine::constructSlitTube(Simulation& System,
				      MonteCarlo::Object* masterCell,
				      const attachSystem::FixedComp& initFC,
				      const std::string& sideName)
  /*!
    Build the Baffle before Monochromator
    \param System :: Simuation to use
   */
{
  ELog::RegMethod RegA("softimaxOpticsLine","buildSlitTube");

  int outerCell;

  // FAKE insertcell: required
  slitTube->addAllInsertCell(masterCell->getName());
  slitTube->createAll(System,initFC,sideName);
  slitTube->intersectPorts(System,0,2);
  slitTube->intersectPorts(System,0,1);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*slitTube,2);
  // split front/back ends off
  const int newCellB=
    slitTube->splitObject(System,"PortBCut",outerCell).back();
  this->addCell("OuterVoid",newCellB);
  const int newCellA=
    slitTube->splitObject(System,"-PortACut",outerCell).back();
  this->addCell("OuterVoid",newCellA);
  cellIndex+=2;  // remember creates an extra cell in  primary

  slitTube->insertInCell("FlangeA",System,newCellA);
  slitTube->insertInCell("FlangeB",System,newCellB);
  slitTube->insertInCell("Main",System,outerCell);
  slitTube->insertPortInCell(System,outerCell);

  slitTube->splitVoidPorts(System,"SplitVoid",1001,
			   slitTube->getCell("Void"),
			   Geometry::Vec3D(0,1,0));

  std::vector<int> cellItems=
    slitTube->splitObject(System,1501,outerCell,
			  Geometry::Vec3D(0,0,0),
			  Geometry::Vec3D(0,0,1));
  this->addCell("OuterVoid",cellItems[1]);

  // Note need cell number later
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
softimaxOpticsLine::buildMono(Simulation& System,
				 MonteCarlo::Object* masterCell,
				 const attachSystem::FixedComp& initFC,
				 const long int sideIndex)
  /*!
    Sub build of the slit package unit
    \param System :: Simulation to use
    \param masterCell :: Main master volume
    \param initFC :: Start point
    \param sideIndex :: start link point
  */
{
  ELog::RegMethod RegA("softimaxOpticsBeamline","buildMono");

  int outerCell;

  // offPipeA->createAll(System,initFC,sideIndex);
  // outerCell=buildZone.createOuterVoidUnit(System,masterCell,*offPipeA,2);
  // offPipeA->insertInCell(System,outerCell);

  // FAKE insertcell: required
  monoVessel->addInsertCell(masterCell->getName());
  monoVessel->createAll(System,initFC,sideIndex);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*monoVessel,2);
  monoVessel->insertInCell(System,outerCell);

  grating->addInsertCell(monoVessel->getCell("Void"));
  grating->copyCutSurf("innerCylinder",*monoVessel,"innerRadius");
  grating->createAll(System,*monoVessel,0);

  zeroOrderBlock->addInsertCell("Body", monoVessel->getCell("Void"));
  //zeroOrderBlock->setBladeCentre(*monoVessel,0);
  zeroOrderBlock->setFront(*monoVessel,3);
  zeroOrderBlock->createAll(System,*monoVessel,3);

  // offPipeB->createAll(System,*monoVessel,2);
  // outerCell=buildZone.createOuterVoidUnit(System,masterCell,*offPipeB,2);
  // offPipeB->insertInCell(System,outerCell);

  return;
}

void
softimaxOpticsLine::buildSplitter(Simulation& System,
				     MonteCarlo::Object* masterCellA,
				     MonteCarlo::Object* masterCellB,
				     const attachSystem::FixedComp& initFC,
				     const long int sideIndex)
  /*!
    Sub build of the spliter package
    \param System :: Simulation to use
    \param masterCellA :: Current master cell
    \param masterCellB :: Secondary master cell
    \param initFC :: Start point
    \param sideIndex :: start link point
  */

{
  ELog::RegMethod RegA("softimaxOpticsBeamLine","buildSplitter");

  int cellA(0),cellB(0);

  // /////////  1: build splitter without creating two outer void units
  // splitter->createAll(System,*M3STXMBack,2);
  // cellA=buildZone.createOuterVoidUnit(System,masterCellA,*splitter,2);
  // splitter->insertInCell("Flange",System,cellA);
  // splitter->insertInCell("PipeA",System,cellA);
  // splitter->insertInCell("Flange",System,cellA);
  // splitter->insertInCell("PipeB",System,cellA);


  //////////////////// 2: build splitter with creating two outer void units
  ////////////////////////////////////////////////////////////////////////////////////
  //  buildZone.constructMiddleSurface(SMap,buildIndex+10,*M3STXMBack,2);

  // No need for insert -- note removal of old master cell

  const int deadCell=masterCellA->getName();

  splitter->createAll(System,initFC,sideIndex);

  //  buildZone.constructMiddleSurface(SMap,buildIndex+10,*M3STXMBack,2);
  buildZone.constructMiddleSurface(SMap,buildIndex+10,*splitter,2,*splitter,3);

  attachSystem::InnerZone leftZone=buildZone.buildMiddleZone(-1);
  attachSystem::InnerZone rightZone=buildZone.buildMiddleZone(1);

  masterCellA=leftZone.constructMasterCell(System);
  masterCellB=rightZone.constructMasterCell(System);

  cellA=leftZone.createOuterVoidUnit(System,masterCellA,*splitter,2);
  cellB=rightZone.createOuterVoidUnit(System,masterCellB,*splitter,3);
  System.removeCell(deadCell);

  splitter->insertInCell("Flange",System,cellA);
  splitter->insertInCell("PipeA",System,cellA);

  splitter->insertInCell("Flange",System,cellB);
  splitter->insertInCell("PipeB",System,cellB);
  ////////////////////////////////////////////////////////////////////////////////////

  leftZone.setExtra();
  constructSystem::constructUnit
    (System,leftZone,masterCellA,*splitter,"back",*bellowAA);
  leftZone.removeExtra();

    bellowBA->setFront(*splitter,3);
  bellowBA->createAll(System,*splitter,3);
  int outerCell=rightZone.createOuterVoidUnit(System,masterCellB,*bellowBA,2);
  bellowBA->insertInCell(System,outerCell);

  // TODO: optimize
  M3Pump->addAllInsertCell(masterCellA->getName());
  M3Pump->addAllInsertCell(masterCellB->getName());
  M3Pump->setPortRotation(3,Geometry::Vec3D(1,0,0));
  M3Pump->createAll(System,*bellowAA,2);

  const constructSystem::portItem& CPI2=M3Pump->getPort(2);
  cellA=leftZone.createOuterVoidUnit(System,masterCellA,CPI2,CPI2.getSideIndex("OuterPlate"));
  const constructSystem::portItem& CPI3=M3Pump->getPort(3);
  cellB=rightZone.createOuterVoidUnit(System,masterCellB,CPI3,CPI3.getSideIndex("OuterPlate"));
  M3Pump->insertAllInCell(System,cellA);
  M3Pump->insertAllInCell(System,cellB);

  // now build left/ right
  // LEFT
  constructSystem::constructUnit
    (System,leftZone,masterCellA,CPI2,"OuterPlate",*bellowAB);
  constructSystem::constructUnit
    (System,leftZone,masterCellA,*bellowAB,"back",*joinPipeAA);

  constructSystem::constructUnit
    (System,leftZone,masterCellA,*joinPipeAA,"back",*bremCollAA);

  // constructSystem::constructUnit
  //   (System,leftZone,masterCellA,*bremCollAA,"back",*joinPipeAB);



  // RIGHT
  constructSystem::constructUnit
    (System,rightZone,masterCellB,CPI3,"OuterPlate",*bellowBB);

  constructSystem::constructUnit
    (System,rightZone,masterCellB,*bellowBB,"back",*joinPipeBA);

  constructSystem::constructUnit
    (System,rightZone,masterCellB,*joinPipeBA,"back",*bremCollBA);

  // constructSystem::constructUnit
  //   (System,rightZone,masterCellB,*bremCollBA,"back",*joinPipeBB);


  // gateAA->createAll(System,*bellowAB,2);
  // cellA=leftZone.createOuterVoidUnit(System,masterCellA,*gateAA,2);
  // gateAA->insertInCell(System,cellA);

  // // make build necessary
  // pumpTubeAA->addAllInsertCell(masterCellA->getName());
  // pumpTubeAA->createAll(System,*gateAA,2);
  // cellA=leftZone.createOuterVoidUnit(System,masterCellA,*pumpTubeAA,2);
  // pumpTubeAA->insertAllInCell(System,cellA);


  // // RIGHT
  // bellowBB->createAll(System,*splitter,3);
  // cellB=rightZone.createOuterVoidUnit(System,masterCellB,*bellowBB,2);
  // bellowBB->insertInCell(System,cellB);

  // gateBA->createAll(System,*bellowBB,2);
  // cellB=rightZone.createOuterVoidUnit(System,masterCellB,*gateBA,2);
  // gateBA->insertInCell(System,cellB);

  // pumpTubeBA->addAllInsertCell(masterCellB->getName());
  // pumpTubeBA->createAll(System,*gateBA,2);
  // cellB=rightZone.createOuterVoidUnit(System,masterCellB,*pumpTubeBA,2);
  // pumpTubeBA->insertAllInCell(System,cellB);

  // Get last two cells
  setCell("LeftVoid",masterCellA->getName());
  setCell("RightVoid",masterCellB->getName());

  return;
}



void
softimaxOpticsLine::buildObjects(Simulation& System)
  /*!
    Build all the objects relative to the main FC
    point.
    \param System :: Simulation to use
  */
{
  ELog::RegMethod RegA("softimaxOpticsLine","buildObjects");

  int outerCell;

  buildZone.setFront(getRule("front"));
  buildZone.setBack(getRule("back"));
  buildZone.setInsertCells(this->getInsertCells());
  MonteCarlo::Object* masterCell=
    buildZone.constructMasterCell(System);

  // dummy space for first item
  // This is a mess but want to preserve insert items already
  // in the hut beam port
  pipeInit->createAll(System,*this,0);
  // dump cell for joinPipe
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*pipeInit,-1);
  if (preInsert)
    preInsert->insertInCell(System,outerCell);
  // real cell for initPipe
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*pipeInit,2);
  pipeInit->insertInCell(System,outerCell);

  // FAKE insertcell: required due to rotation ::
  triggerPipe->addAllInsertCell(masterCell->getName());
  triggerPipe->setPortRotation(3,Geometry::Vec3D(1,0,0));
  triggerPipe->createAll(System,*pipeInit,2);

  const constructSystem::portItem& TPI=triggerPipe->getPort(1);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,TPI,TPI.getSideIndex("OuterPlate"));
  triggerPipe->insertAllInCell(System,outerCell);

  // FAKE insertcell: required
  gateTubeA->addAllInsertCell(masterCell->getName());
  gateTubeA->setPortRotation(3,Geometry::Vec3D(1,0,0));
  gateTubeA->createAll(System,TPI,TPI.getSideIndex("OuterPlate"));

  const constructSystem::portItem& GPI1=gateTubeA->getPort(1);
  outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,GPI1,GPI1.getSideIndex("OuterPlate"));
  gateTubeA->insertAllInCell(System,outerCell);

  gateTubeAItem->addInsertCell("Body",gateTubeA->getCell("Void"));
  gateTubeAItem->setBladeCentre(*gateTubeA,0);
  gateTubeAItem->createAll(System,*gateTubeA,std::string("InnerBack"));

  constructSystem::constructUnit
    (System,buildZone,masterCell,GPI1,"OuterPlate",*bellowA);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*bellowA,"back",*pipeA);

  // FAKE insertcell: required
  pumpM1->addAllInsertCell(masterCell->getName());
  pumpM1->setPortRotation(3,Geometry::Vec3D(1,0,0));
  pumpM1->setOuterVoid();
  pumpM1->createAll(System,*pipeA,"back");
  //pumpM1->intersectPorts(System,1,2);

  ///////////// split for FLUKA
  //  const constructSystem::portItem& VP0=pumpM1->getPort(0);
  const constructSystem::portItem& VP1=pumpM1->getPort(1);
  const constructSystem::portItem& VP2=pumpM1->getPort(2);
  //  const constructSystem::portItem& VP3=pumpM1->getPort(3);
  const constructSystem::portItem& VP4=pumpM1->getPort(4);
  //  const constructSystem::portItem& VP5=pumpM1->getPort(5);
  const constructSystem::portItem& VP6=pumpM1->getPort(6);

  outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,VP1,VP1.getSideIndex("OuterPlate"));
  const Geometry::Vec3D  Axis12=pumpM1->getY()*(VP1.getY()+VP2.getY())/2.0;
  const Geometry::Vec3D  Axis26=pumpM1->getY()*(VP2.getY()+VP6.getY())/2.0;

  this->splitObjectAbsolute(System,1501,outerCell,
			    (VP1.getCentre()+VP4.getCentre())/2.0,
  			    Z);
  this->splitObjectAbsolute(System,1502,outerCell+1,
    			    pumpM1->getCentre(),
			    VP4.getY());
  this->splitObjectAbsolute(System,1503,outerCell,
			    pumpM1->getCentre(),Axis12);
  this->splitObjectAbsolute(System,1504,outerCell+3,
   			    pumpM1->getCentre(),Axis26);

  const std::vector<int> cellUnit=this->getCells("OuterVoid");
  pumpM1->insertMainInCell(System,cellUnit);
  //  VP0.insertInCell(System,this->getCell("OuterVoid"));

  pumpM1->insertPortInCell
    (System,{{outerCell+4,outerCell},{outerCell,outerCell+1,outerCell+2},{outerCell+3},{outerCell},
	     {outerCell+2},{outerCell+1},
	     {outerCell+4}});

  cellIndex+=5;
  /////////////////////////////////////////

  //  setCell("LastVoid",masterCell->getName());  lastComp=pumpM1;  return;


  constructSystem::constructUnit
    (System,buildZone,masterCell,VP1,"OuterPlate",*gateA);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*gateA,"back",*bellowB);

  buildM1Mirror(System,masterCell,*bellowB,"back");

  constructSystem::constructUnit
    (System,buildZone,masterCell,*M1TubeBack,"back",*bellowC);

  // FAKE insertcell: required
  pumpTubeA->addAllInsertCell(masterCell->getName());
  pumpTubeA->setPortRotation(3,Geometry::Vec3D(1,0,0));
  pumpTubeA->createAll(System,*bellowC,2);

  const constructSystem::portItem& CPI1=pumpTubeA->getPort(1);
  outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,CPI1,CPI1.getSideIndex("OuterPlate"));
  pumpTubeA->insertAllInCell(System,outerCell);
  //  pumpTubeA->intersectPorts(System,1,2);


  constructSystem::constructUnit
    (System,buildZone,masterCell,CPI1,"OuterPlate",*bremCollA);


  constructSystem::constructUnit
    (System,buildZone,masterCell,*bremCollA,"back",*gateB);

  //  bremCollA->createExtension(System,gateB->getCell("Void")); // !!! UGLY - it does not actually intersect gateB

  constructSystem::constructUnit
    (System,buildZone,masterCell,*gateB,"back",*bellowD);

  constructSlitTube(System,masterCell,*bellowD,"back");

  buildMono(System,masterCell,*slitTube,2);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*monoVessel,"back",*gateC);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*gateC,"back",*bellowE);

  // FAKE insertcell: required
  pumpTubeB->addAllInsertCell(masterCell->getName());
  pumpTubeB->setPortRotation(3,Geometry::Vec3D(1,0,0));
  pumpTubeB->createAll(System,*bellowE,2);

  const constructSystem::portItem& pumpTubeBCPI=pumpTubeB->getPort(1);
  outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,pumpTubeBCPI,pumpTubeBCPI.getSideIndex("OuterPlate"));
  pumpTubeB->insertAllInCell(System,outerCell);
  //  pumpTubeB->intersectPorts(System,1,2);

  constructSystem::constructUnit
    (System,buildZone,masterCell,pumpTubeBCPI,"OuterPlate",*gateD);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*gateD,"back",*joinPipeA);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*joinPipeA,"back",*bellowF);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*bellowF,"back",*slitsA);

  /////////////////// M3 Pump and baffle
  // FAKE insertcell: required
  pumpTubeM3->addAllInsertCell(masterCell->getName());
  pumpTubeM3->setPortRotation(3,Geometry::Vec3D(1,0,0));
  pumpTubeM3->createAll(System,*slitsA,"back");

  const constructSystem::portItem& GPI=pumpTubeM3->getPort(1);
  outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,GPI,GPI.getSideIndex("OuterPlate"));
  pumpTubeM3->insertAllInCell(System,outerCell);

  pumpTubeM3Baffle->addInsertCell("Body",pumpTubeM3->getCell("Void"));
  pumpTubeM3Baffle->setBladeCentre(*pumpTubeM3,0);
  pumpTubeM3Baffle->createAll(System,*pumpTubeM3,std::string("InnerBack"));

  constructSystem::constructUnit
    (System,buildZone,masterCell,GPI,"OuterPlate",*bellowG);

  buildM3Mirror(System,masterCell,*bellowG,"back");

  constructSystem::constructUnit
    (System,buildZone,masterCell,*M3Back,"back",*bellowH);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*bellowH,"back",*gateE);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*gateE,"back",*joinPipeB);

  //// pumpTubeC
  // FAKE insertcell: required
  pumpTubeC->addAllInsertCell(masterCell->getName());
  pumpTubeC->setPortRotation(3,Geometry::Vec3D(1,0,0));
  pumpTubeC->createAll(System,*joinPipeB,2);

  const constructSystem::portItem& pumpTubeCCPI=pumpTubeC->getPort(1);
  outerCell=buildZone.createOuterVoidUnit
    (System,masterCell,pumpTubeCCPI,pumpTubeCCPI.getSideIndex("OuterPlate"));
  pumpTubeC->insertAllInCell(System,outerCell);
  //  pumpTubeC->intersectPorts(System,1,2);
  ///////////////////////////////////////////////

  constructSystem::constructUnit
    (System,buildZone,masterCell,pumpTubeCCPI,"back",*bellowI);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*bellowI,"back",*joinPipeC);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*joinPipeC,"back",*gateF);

  constructSystem::constructUnit
    (System,buildZone,masterCell,*gateF,"back",*bellowJ);

  buildM3STXMMirror(System,masterCell,*bellowJ,"back");


  MonteCarlo::Object* masterCellB(0);
  buildSplitter(System,masterCell,masterCellB,*M3STXMTube,2);

  setCell("LastVoid",masterCell->getName());
  lastComp=bellowA; //gateJ;
  return;
}

void
softimaxOpticsLine::buildOutGoingPipes(Simulation& System,
				       const int leftCell,
				       const int rightCell,
				       const std::vector<int>& hutCells)
  /*!
    Construct outgoing tracks
    \param System :: Simulation
    \param leftCell :: additional left cell for insertion
    \param rightCell :: additional right cell for insertion
    \param hutCell :: Cells for construction in hut [common to both pipes]
  */
{
  ELog::RegMethod RegA("softimaxOpticsLine","buildOutgoingPipes");

  joinPipeAB->addInsertCell(hutCells);
  joinPipeAB->addInsertCell(leftCell);
  joinPipeAB->createAll(System,*bremCollAA,2);

  joinPipeBB->addInsertCell(hutCells);
  joinPipeBB->addInsertCell(rightCell);
  joinPipeBB->createAll(System,*bremCollBA,2);

  screenA->addAllInsertCell(rightCell);
  screenA->addAllInsertCell(leftCell);
  screenA->setCutSurf("inner",*joinPipeAB,"pipeOuterTop");
  screenA->setCutSurf("innerTwo",*joinPipeBA,"pipeOuterTop");
  screenA->createAll(System,*bremCollAA,2);

  // lineScreen->setNoInsert();
  // lineScreen->addInsertCell(ContainedComp::getInsertCells());
  // lineScreen->createAll(System,*gateA,"back");

  // innerScreen->setNoInsert();
  // for(size_t index=12;index<=20;index++)
  //   innerScreen->addInsertCell(this->getCell("OuterVoid",index));

  innerScreen->addInsertCell(this->getCell("OuterVoid",23));
  innerScreen->createAll(System,*gateA,"back");

  return;
}


void
softimaxOpticsLine::createLinks()
  /*!
    Create a front/back link
   */
{
  ELog::RegMethod RControl("softimaxOpticsLine","createLinks");

  setLinkSignedCopy(0,*pipeInit,1);
  setLinkSignedCopy(1,*lastComp,2);
  return;
}


void
softimaxOpticsLine::createAll(Simulation& System,
			      const attachSystem::FixedComp& FC,
			      const long int sideIndex)
  /*!
    Carry out the full build
    \param System :: Simulation system
    \param FC :: Fixed component
    \param sideIndex :: link point
  */
{
  ELog::RegMethod RControl("softimaxOpticsLine","createAll");

  populate(System.getDataBase());
  createUnitVector(FC,sideIndex);
  createSurfaces();

  buildObjects(System);
  createLinks();
  return;
}


}   // NAMESPACE xraySystem
