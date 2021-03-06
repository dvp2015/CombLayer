/*********************************************************************
  CombLayer : MCNP(X) Input builder

 * File:   Model/MaxIV/cosaxs/cosaxsTube.cxx
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
#include <memory>

#include "FileReport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "Vec3D.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "groupRange.h"
#include "objectGroups.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "generateSurf.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "FixedRotate.h"
#include "ContainedComp.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "SurfMap.h"
#include "ExternalCut.h"
#include "InnerZone.h"
#include "FrontBackCut.h"

#include "GateValveCylinder.h"
#include "cosaxsTubeNoseCone.h"
#include "cosaxsTubeStartPlate.h"
#include "cosaxsTubeBeamDump.h"
#include "cosaxsTubeWAXSDetector.h"
#include "cosaxsTubeAirBox.h"
#include "cosaxsTubeCable.h"

#include "ContainedGroup.h"
#include "VirtualTube.h"
#include "PipeTube.h"

#include "cosaxsTube.h"

namespace xraySystem
{

cosaxsTube::cosaxsTube(const std::string& Key)  :
  attachSystem::ContainedComp(),
  attachSystem::FixedOffset(Key,6),
  attachSystem::CellMap(),
  attachSystem::SurfMap(),
  attachSystem::FrontBackCut(),
  delayPortFlag(0),
  buildZone(*this,cellIndex),
  buildZoneTube(*this,cellIndex),
  noseCone(new xraySystem::cosaxsTubeNoseCone(keyName+"NoseCone")),
  gateA(new constructSystem::GateValveCylinder(keyName+"GateA")),
  startPlate(new xraySystem::cosaxsTubeStartPlate(keyName+"StartPlate")),
  beamDump(new xraySystem::cosaxsTubeBeamDump(keyName+"BeamDump")),
  waxs(new xraySystem::cosaxsTubeWAXSDetector(keyName+"WAXS")),
  airBox(new xraySystem::cosaxsTubeAirBox(keyName+"AirBox")),
  cable(new xraySystem::cosaxsTubeCable(keyName+"Cable"))
 /*!
    Constructor BUT ALL variable are left unpopulated.
    \param Key :: Name for item in search
  */
{
  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();

  OR.addObject(noseCone);
  OR.addObject(gateA);
  OR.addObject(startPlate);

  for(size_t i=0;i<8;i++)
    {
      seg[i] = std::make_shared<constructSystem::PipeTube>
	(keyName+"Segment"+std::to_string(i+1));
      OR.addObject(seg[i]);
    }

  OR.addObject(beamDump);
  OR.addObject(waxs);
  OR.addObject(airBox);
  OR.addObject(cable);
}

cosaxsTube::cosaxsTube(const cosaxsTube& A) :
  attachSystem::ContainedComp(A),
  attachSystem::FixedOffset(A),
  attachSystem::CellMap(A),
  attachSystem::SurfMap(A),
  attachSystem::FrontBackCut(A),
  delayPortFlag(A.delayPortFlag),
  outerRadius(A.outerRadius),
  outerLength(A.outerLength),
  buildZone(A.buildZone),
  buildZoneTube(A.buildZoneTube),
  noseCone(A.noseCone),
  gateA(A.gateA),
  startPlate(A.startPlate),
  seg(A.seg),
  beamDump(A.beamDump),
  waxs(A.waxs),
  airBox(A.airBox),
  cable(A.cable)
  /*!
    Copy constructor
    \param A :: cosaxsTube to copy
  */
{}

cosaxsTube&
cosaxsTube::operator=(const cosaxsTube& A)
  /*!
    Assignment operator
    \param A :: cosaxsTube to copy
    \return *this
  */
{
  if (this!=&A)
    {
      attachSystem::ContainedComp::operator=(A);
      attachSystem::FixedOffset::operator=(A);
      attachSystem::CellMap::operator=(A);
      attachSystem::SurfMap::operator=(A);
      attachSystem::FrontBackCut::operator=(A);
      delayPortFlag=A.delayPortFlag;
      outerRadius=A.outerRadius;
      outerLength=A.outerLength;
      noseCone=A.noseCone;
      gateA=A.gateA;
      startPlate=A.startPlate;
      seg=A.seg;
      beamDump=A.beamDump;
      waxs=A.waxs;
      airBox=A.airBox;
      cable=A.cable;
    }
  return *this;
}

cosaxsTube*
cosaxsTube::clone() const
/*!
  Clone self
  \return new (this)
 */
{
    return new cosaxsTube(*this);
}

cosaxsTube::~cosaxsTube()
  /*!
    Destructor
  */
{}

void
cosaxsTube::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable data base
  */
{
  ELog::RegMethod RegA("cosaxsTube","populate");

  FixedOffset::populate(Control);

  outerRadius=Control.EvalVar<double>(keyName+"OuterRadius");
  outerLength=Control.EvalVar<double>(keyName+"OuterLength");

  return;
}

void
cosaxsTube::createUnitVector(const attachSystem::FixedComp& FC,
			      const long int sideIndex)
  /*!
    Create the unit vectors
    \param FC :: object for origin
    \param sideIndex :: link point for origin
  */
{
  ELog::RegMethod RegA("cosaxsTube","createUnitVector");

  FixedComp::createUnitVector(FC,sideIndex);
  applyOffset();

  return;
}

void
cosaxsTube::createSurfaces()
  /*!
    Create All the surfaces
  */
{
  ELog::RegMethod RegA("cosaxsTube","createSurfaces");

  if (!frontActive())
    {
      ModelSupport::buildPlane(SMap,buildIndex+1,Origin,Y);
      FrontBackCut::setFront(SMap.realSurf(buildIndex+1));
    }

  if (!backActive())
    {
      ModelSupport::buildPlane(SMap,buildIndex+2,Origin+Y*(outerLength),Y);
      FrontBackCut::setBack(-SMap.realSurf(buildIndex+2));
    }

  ModelSupport::buildCylinder(SMap,buildIndex+7,Origin,Y,outerRadius);

  // // cable
  // ModelSupport::buildPlane(SMap,buildIndex+103,Origin-X*(cableWidth/2.0),X);
  // ModelSupport::buildPlane(SMap,buildIndex+104,Origin+X*(cableWidth/2.0),X);
  // ModelSupport::buildPlane(SMap,buildIndex+105,Origin-Z*(cableHeight/2.0-cableZStep),Z);
  // ModelSupport::buildPlane(SMap,buildIndex+106,Origin+Z*(cableHeight/2.0+cableZStep),Z);
  // ModelSupport::buildShiftedPlane(SMap,buildIndex+115,
  // 				  SMap.realPtr<Geometry::Plane>(buildIndex+105),
  // 				  -cableHeight);

  // const double cableTailLength(3*M_PI/4*cableTailRadius);
  // const double cableBottomLength(detYStep/2.0);
  // const double cableUpLength(cableLength-cableTailLength-cableBottomLength-detYStep);

  // ModelSupport::buildPlane(SMap,buildIndex+101,Origin+Y*(detYStep),Y);
  // // centre of cylinders
  // const double yTail(detYStep+detYStep+cableUpLength+cableTailRadius);
  // ModelSupport::buildPlane(SMap,buildIndex+102,Origin+Y*yTail,Y);
  // ModelSupport::buildShiftedPlane(SMap,buildIndex+111,
  // 				  SMap.realPtr<Geometry::Plane>(buildIndex+102),
  // 				  -cableBottomLength);

  // ModelSupport::buildCylinder(SMap,buildIndex+107,
  // 			      Origin+Y*yTail+Z*cableZStep,X,
  // 			      cableTailRadius);
  // ModelSupport::buildCylinder(SMap,buildIndex+117,
  // 			      Origin+Y*yTail+Z*cableZStep,X,
  // 			      cableTailRadius+cableHeight);


  const std::string Out=ModelSupport::getComposite(SMap,buildIndex," -7 ");
  const HeadRule HR(Out);
  buildZone.setSurround(HR);

  return;
}

void
cosaxsTube::createObjects(Simulation& System)
  /*!
    Adds the all the components
    \param System :: Simulation to create objects in
  */
{
  ELog::RegMethod RegA("cosaxsTube","createObjects");

  std::string Out;
  const std::string frontStr(frontRule());
  const std::string backStr(backRule());

  Out=ModelSupport::getComposite(SMap,buildIndex," -7 ");
  addOuterSurf(Out+frontStr+backStr);

  int outerCell;
  buildZone.setFront(getRule("front"));//HeadRule(SMap.realSurf(buildIndex+1)));
  buildZone.setBack(getRule("back"));//HeadRule(-SMap.realSurf(buildIndex+2)));
  buildZone.setInsertCells(this->getInsertCells());
  MonteCarlo::Object* masterCell=buildZone.constructMasterCell(System);

  noseCone->createAll(System, *this, 0);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*noseCone,2);
  noseCone->insertInCell(System,outerCell);

  gateA->setFront(*noseCone,2);
  gateA->createAll(System,*noseCone,2);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*gateA,2);
  gateA->insertInCell(System,outerCell);

  startPlate->setFront(*gateA,2);
  startPlate->createAll(System,*gateA,2);
  outerCell=buildZone.createOuterVoidUnit(System,masterCell,*startPlate,2);
  startPlate->insertInCell(System,outerCell);

  // tube segments
  attachSystem::FixedComp *last = startPlate.get();

  for (size_t i=0; i<8; i++)
    {
      seg[i]->setFront(*last,2);
      if (delayPortFlag) seg[i]->delayPorts();
      seg[i]->createAll(System,*last,2);

      outerCell=buildZone.createOuterVoidUnit(System,masterCell,*seg[i],2);
      seg[i]->insertAllInCell(System,outerCell);

      // delete the individual inner void cells in order to create
      // a common InnerVoid after this loop
      seg[i]->deleteCell(System,"Void");
      if (i==5)
	{
	  // paired cellVec:
	  std::vector<int> CellVec;
	  CellVec=seg[i]->splitObject
	    (System,3001,getCell("OuterVoid",i+2),
	     Geometry::Vec3D(0,0,0),Geometry::Vec3D(-1,0,0.5));
	  this->addCell("OuterVoid",CellVec.back());

	  CellVec=seg[i]->splitObject
	    (System,3002,getCell("OuterVoid",i+2),
	     Geometry::Vec3D(0,0,0),Geometry::Vec3D(1,0,0.5));
	  this->addCell("OuterVoid",CellVec.back());

	  CellVec=seg[i]->splitObject
	    (System,3003,getCell("OuterVoid",i+4),
	     Geometry::Vec3D(0,0,0),Geometry::Vec3D(1,0,0.5));
	  this->addCell("OuterVoid",CellVec.back());

	  cellIndex+=4;
	}

      last = seg[i].get();
    }

  buildZoneTube.setSurround(last->getFullRule("InnerSide"));
  buildZoneTube.setFront(seg[0]->getFullRule("InnerFront"));
  buildZoneTube.setBack(last->getFullRule("InnerBack"));
//   masterCell=buildZoneTube.constructMasterCell(System,*this);
  buildZone.setInsertCells(this->getInsertCells());

  masterCell=buildZoneTube.constructMasterCell(System);

  beamDump->setFront(*seg[0],-1);
  beamDump->createAll(System,*seg[0],-1);
  outerCell=buildZoneTube.createOuterVoidUnit(System,masterCell,*beamDump,2);
  beamDump->insertInCell(System,outerCell);

  waxs->setFront(*beamDump,2);
  waxs->createAll(System,*beamDump,2);
  outerCell=buildZoneTube.createOuterVoidUnit(System,masterCell,*waxs,2);
  waxs->insertInCell(System,outerCell);

  airBox->setFront(*waxs,2);
  airBox->createAll(System,*waxs,2);
  outerCell=buildZoneTube.createOuterVoidUnit(System,masterCell,*airBox,2);
  airBox->insertInCell(System,outerCell);

  cable->setFront(*airBox,2);
  cable->createAll(System,*airBox,2);
  outerCell=buildZoneTube.createOuterVoidUnit(System,masterCell,*cable,2);
  cable->insertInCell(System,outerCell);

  return;
}

void
cosaxsTube::createLinks()
  /*!
    Create all the linkes
  */
{
  ELog::RegMethod RegA("cosaxsTube","createLinks");

  FrontBackCut::createLinks(*this,Origin,Y);

  return;
}

void
cosaxsTube::createPorts(Simulation& System)
  /*!
    Generic function to create the ports
    \param System :: Simulation item
  */
{
  ELog::RegMethod RegA("cosaxsTube","createPorts");

  for (size_t i=0; i<8; i++)
    seg[i]->createPorts(System);
  return;
}

void
cosaxsTube::createAll(Simulation& System,
		       const attachSystem::FixedComp& FC,
		       const long int sideIndex)
  /*!
    Generic function to create everything
    \param System :: Simulation item
    \param FC :: Central origin
    \param sideIndex :: link point for origin
  */
{
  ELog::RegMethod RegA("cosaxsTube","createAll");

  populate(System.getDataBase());
  createUnitVector(FC,sideIndex);
  createSurfaces();
  createObjects(System);
  createLinks();
  insertObjects(System);

  return;
}

}  // xraySystem
