/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   Linac/EArrivalMon.cxx
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
#include <memory>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "support.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfDIter.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "surfEqual.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Line.h"
#include "Rules.h"
#include "SurInter.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "SimProcess.h"
#include "groupRange.h"
#include "objectGroups.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "LinkUnit.h"  
#include "FixedComp.h"
#include "FixedOffset.h"
#include "ContainedComp.h"
#include "ExternalCut.h"
#include "FrontBackCut.h" 
#include "BaseMap.h"
#include "SurfMap.h"
#include "CellMap.h" 

#include "EArrivalMon.h"

namespace tdcSystem
{

EArrivalMon::EArrivalMon(const std::string& Key) :
  attachSystem::FixedOffset(Key,6),
  attachSystem::ContainedComp(),
  attachSystem::FrontBackCut(),
  attachSystem::CellMap(),
  attachSystem::SurfMap()
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param Key :: KeyName
  */
{}


EArrivalMon::~EArrivalMon() 
  /*!
    Destructor
  */
{}

void
EArrivalMon::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: DataBase for variables
  */
{
  ELog::RegMethod RegA("EArrivalMon","populate");

  FixedOffset::populate(Control);

  radius=Control.EvalVar<double>(keyName+"Radius");
  length=Control.EvalVar<double>(keyName+"Length");
  thick=Control.EvalVar<double>(keyName+"Thick");
  faceThick=Control.EvalVar<double>(keyName+"FaceThick");

  frontPipeILen=Control.EvalVar<double>(keyName+"FrontPipeILen");
  frontPipeLen=Control.EvalVar<double>(keyName+"FrontPipeLen");
  frontPipeRadius=Control.EvalVar<double>(keyName+"FrontPipeRadius");
  frontPipeThick=Control.EvalVar<double>(keyName+"FrontPipeThick");

  backPipeILen=Control.EvalVar<double>(keyName+"BackPipeILen");
  backPipeLen=Control.EvalVar<double>(keyName+"BackPipeLen");
  backPipeRadius=Control.EvalVar<double>(keyName+"BackPipeRadius");
  backPipeThick=Control.EvalVar<double>(keyName+"BackPipeThick");

  flangeRadius=Control.EvalVar<double>(keyName+"FlangeRadius");
  flangeLength=Control.EvalVar<double>(keyName+"FlangeLength");

  windowRotAngle=Control.EvalVar<double>(keyName+"WindowRotAngle");
  windowRadius=Control.EvalVar<double>(keyName+"WindowRadius");
  windowThick=Control.EvalVar<double>(keyName+"WindowThick");
  

  voidMat=ModelSupport::EvalMat<int>(Control,keyName+"VoidMat");
  flangeMat=ModelSupport::EvalMat<int>(Control,keyName+"FlangeMat");
  electrodeMat=ModelSupport::EvalMat<int>(Control,keyName+"ElectrodeMat");
  outerMat=ModelSupport::EvalMat<int>(Control,keyName+"OuterMat");

  return;
}


void
EArrivalMon::createSurfaces()
  /*!
    Create All the surfaces
  */
{
  ELog::RegMethod RegA("EArrivalMon","createSurfaces");

  if (!isActive("front"))
    {
      ModelSupport::buildPlane(SMap,buildIndex+1,Origin-Y*(length/2.0),Y);
      ExternalCut::setCutSurf("front",SMap.realSurf(buildIndex+1));
    }
  if (!isActive("back"))
    {
      ModelSupport::buildPlane(SMap,buildIndex+2,Origin+Y*(length/2.0),Y);
      ExternalCut::setCutSurf("back",-SMap.realSurf(buildIndex+2));
    }

  // main pipe and thicness
  ModelSupport::buildCylinder(SMap,buildIndex+7,Origin,Y,radius);
  ModelSupport::buildCylinder(SMap,buildIndex+17,Origin,Y,radius+outerThick);
  
  ModelSupport::buildCylinder(SMap,buildIndex+107,Origin,Y,flangeARadius);
  ModelSupport::buildCylinder(SMap,buildIndex+207,Origin,Y,flangeBRadius);

  ModelSupport::buildPlane(SMap,buildIndex+101,
			   Origin-Y*(length/2.0-flangeALength),Y);
  ModelSupport::buildPlane(SMap,buildIndex+202,
			   Origin+Y*(length/2.0-flangeBLength),Y);

  // Field shaping [300]
  ModelSupport::buildCylinder(SMap,buildIndex+307,Origin,Y,innerRadius);
  ModelSupport::buildCylinder
    (SMap,buildIndex+317,Origin,Y,innerRadius+innerThick);

  const double ang(M_PI*innerAngle/(2.0*180.0));
  double midAngle(innerAngleOffset*M_PI/180.0);
  int BI(buildIndex+350);
  for(size_t i=0;i<4;i++)
    {
      const Geometry::Vec3D lowAxis(X*cos(midAngle-ang)+Z*sin(midAngle-ang));
      const Geometry::Vec3D highAxis(X*cos(midAngle+ang)+Z*sin(midAngle+ang));
      ModelSupport::buildPlane(SMap,BI+1,Origin,lowAxis);
      ModelSupport::buildPlane(SMap,BI+2,Origin,highAxis);
      midAngle+=M_PI/2.0;
      BI+=2;
    }
      
  // end point on electrons (solid)
  ModelSupport::buildPlane(SMap,buildIndex+302,
			   Origin+Y*(length/2.0-electrodeEnd),Y);

  
  // Electrode holder
  const double angleStep(M_PI/4.0);

  BI=buildIndex+411;
  for(size_t i=0;i<8;i++)
    {
      const double angle(angleStep*static_cast<double>(i));
      const Geometry::Vec3D axis(X*cos(angle)+Z*sin(angle));
      ModelSupport::buildPlane(SMap,BI,Origin+axis*electrodeRadius,axis);
      BI++;
    }

  const Geometry::Vec3D eCent(Origin+Y*(electrodeYStep-length/2.0));
  ModelSupport::buildPlane
    (SMap,buildIndex+401,eCent-Y*(electrodeThick/2.0),Y);
  ModelSupport::buildPlane
    (SMap,buildIndex+402,eCent+Y*(electrodeThick/2.0),Y);
  
  return;
}

void
EArrivalMon::createObjects(Simulation& System)
  /*!
    Builds all the objects
    \param System :: Simulation to create objects in
  */
{
  ELog::RegMethod RegA("EArrivalMon","createObjects");

  std::string Out;
  
  const std::string frontStr=getRuleStr("front");
  const std::string backStr=getRuleStr("back");

  // inner void
  Out=ModelSupport::getComposite(SMap,buildIndex," -307 ");
  makeCell("Void",System,cellIndex++,voidMat,0.0,Out+frontStr+backStr);

  
  // front void
  Out=ModelSupport::getComposite(SMap,buildIndex," -7 307 -401 ");
  makeCell("FrontVoid",System,cellIndex++,voidMat,0.0,Out+frontStr);

  // main walll
  Out=ModelSupport::getComposite(SMap,buildIndex," 101 -202 7 -17 ");
  makeCell("Outer",System,cellIndex++,outerMat,0.0,Out);

  // front flange
  Out=ModelSupport::getComposite(SMap,buildIndex," -101 7 -107 ");
  makeCell("FlangeA",System,cellIndex++,flangeMat,0.0,Out+frontStr);

  // back flange
  Out=ModelSupport::getComposite(SMap,buildIndex," 202 7 -207 ");
  makeCell("FlangeB",System,cellIndex++,flangeMat,0.0,Out+backStr);

  Out=ModelSupport::getComposite
    (SMap,buildIndex," 17 401 -402 -411 -412 -413 -414 -415 -416 -417 -418 ");
  makeCell("ElectronPlate",System,cellIndex++,electrodeMat,0.0,Out);


  Out=ModelSupport::getComposite
    (SMap,buildIndex," 101 17 -401 -411 -412 -413 -414 -415 -416 -417 -418 ");
  makeCell("OuterVoid",System,cellIndex++,voidMat,0.0,Out);
  
  Out=ModelSupport::getComposite
    (SMap,buildIndex," -202 17 402 -411 -412 -413 -414 -415 -416 -417 -418 ");
  makeCell("OuterVoid",System,cellIndex++,voidMat,0.0,Out);
  
  Out=ModelSupport::getComposite
    (SMap,buildIndex," -101 107 -411 -412 -413 -414 -415 -416 -417 -418 ");
  makeCell("OuterVoid",System,cellIndex++,voidMat,0.0,Out+frontStr);

  Out=ModelSupport::getComposite
    (SMap,buildIndex," 202 207 -411 -412 -413 -414 -415 -416 -417 -418 ");
  makeCell("OuterVoid",System,cellIndex++,voidMat,0.0,Out+backStr);

  int BI(0);
  const std::string IOut=
    ModelSupport::getComposite(SMap,buildIndex," 401 -302  307 -317 "); 
  for(size_t i=0;i<4;i++)
    {
      Out=ModelSupport::getComposite(SMap,buildIndex+BI," 351 -352");
      makeCell("Electrode",System,cellIndex++,electrodeMat,0.0,IOut+Out);
      Out=ModelSupport::getRangeComposite
	(SMap,351,358,BI,buildIndex," 352R -353R");
      makeCell("ElectrodeGap",System,cellIndex++,voidMat,0.0,IOut+Out);
      BI+=2;
    }
  // edge electrod void
  Out=ModelSupport::getComposite(SMap,buildIndex," 401 -302 317 -7");
  makeCell("EdgeVoid",System,cellIndex++,voidMat,0.0,Out);

  // edge electrod void
  Out=ModelSupport::getComposite(SMap,buildIndex," 302  307 -7");
  makeCell("ElectrodeEnd",System,cellIndex++,electrodeMat,0.0,Out+backStr);


  Out=ModelSupport::getComposite
      (SMap,buildIndex," -411 -412 -413 -414 -415 -416 -417 -418 ");
  addOuterSurf(Out);
  
  return;
}

void 
EArrivalMon::createLinks()
  /*!
    Create the linked units
   */
{
  ELog::RegMethod RegA("EArrivalMon","createLinks");

  ExternalCut::createLink("front",*this,0,Origin,Y);  //front and back
  ExternalCut::createLink("back",*this,1,Origin,Y);  //front and back
      
  return;
}

void
EArrivalMon::createAll(Simulation& System,
	       const attachSystem::FixedComp& FC,
	       const long int sideIndex)
  /*!
    Generic function to create everything
    \param System :: Simulation item
    \param FC :: Fixed point track 
    \param sideIndex :: link point
  */
{
  ELog::RegMethod RegA("EArrivalMon","createAll");
  
  populate(System.getDataBase());
  createCentredUnitVector(FC,sideIndex,length);
  createSurfaces();
  createObjects(System);
  createLinks();
  insertObjects(System);   
  
  return;
}
  
}  // NAMESPACE tdcSystem
