/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   bibBuild/RotorWheel.cxx
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
#include <memory>

#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "Surface.h"
#include "surfIndex.h"
#include "Quadratic.h"
#include "Rules.h"
#include "HeadRule.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "Object.h"
#include "groupRange.h"
#include "objectGroups.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "support.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "ContainedComp.h"
#include "ContainedGroup.h"

#include "RotorWheel.h"

namespace bibSystem
{

RotorWheel::RotorWheel(const std::string& Key) :
  attachSystem::ContainedGroup("Target","Body"),
  attachSystem::FixedOffset(Key,15)
  /*!
    Constructor
    \param Key :: Name of construction key
  */
{}

RotorWheel::RotorWheel(const RotorWheel& A) : 
  attachSystem::ContainedGroup(A),attachSystem::FixedOffset(A),
  radius(A.radius),
  waterThick(A.waterThick),wallThick(A.wallThick),
  beThick(A.beThick),beLength(A.beLength),frontVac(A.frontVac),
  backVac(A.backVac),endVac(A.endVac),outWallThick(A.outWallThick),
  boxFront(A.boxFront),boxBack(A.boxBack),boxIn(A.boxIn),
  boxOut(A.boxOut),boxTop(A.boxTop),boxBase(A.boxBase),
  boxInterior(A.boxInterior),steelThick(A.steelThick),
  pbDepth(A.pbDepth),pbThick(A.pbThick),beMat(A.beMat),
  waterMat(A.waterMat),wallMat(A.wallMat),vesselMat(A.vesselMat),
  polyMat(A.polyMat),pbMat(A.pbMat)
  /*!
    Copy constructor
    \param A :: RotorWheel to copy
  */
{}

RotorWheel&
RotorWheel::operator=(const RotorWheel& A)
  /*!
    Assignment operator
    \param A :: RotorWheel to copy
    \return *this
  */
{
  if (this!=&A)
    {
      attachSystem::ContainedGroup::operator=(A);
      attachSystem::FixedOffset::operator=(A);
      radius=A.radius;
      waterThick=A.waterThick;
      wallThick=A.wallThick;
      beThick=A.beThick;
      beLength=A.beLength;
      frontVac=A.frontVac;
      backVac=A.backVac;
      endVac=A.endVac;
      outWallThick=A.outWallThick;
      boxFront=A.boxFront;
      boxBack=A.boxBack;
      boxIn=A.boxIn;
      boxOut=A.boxOut;
      boxTop=A.boxTop;
      boxBase=A.boxBase;
      boxInterior=A.boxInterior;
      steelThick=A.steelThick;
      pbDepth=A.pbDepth;
      pbThick=A.pbThick;
      beMat=A.beMat;
      waterMat=A.waterMat;
      wallMat=A.wallMat;
      vesselMat=A.vesselMat;
      polyMat=A.polyMat;
      pbMat=A.pbMat;
    }
  return *this;
}

RotorWheel::~RotorWheel()
  /*!
    Destructor
   */
{}

void
RotorWheel::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: Variable table to use
  */
{
  ELog::RegMethod RegA("RotorWheel","populate");

  // Master object : everything references against this

  FixedOffset::populate(Control);
  
  radius=Control.EvalVar<double>(keyName+"Radius");
  waterThick=Control.EvalVar<double>(keyName+"WaterThick");
  wallThick=Control.EvalVar<double>(keyName+"WallThick");
  beThick=Control.EvalVar<double>(keyName+"BeThick");
  beLength=Control.EvalVar<double>(keyName+"BeLength");

  frontVac=Control.EvalVar<double>(keyName+"FrontVac");
  backVac=Control.EvalVar<double>(keyName+"BackVac");
  endVac=Control.EvalVar<double>(keyName+"EndVac");

  outWallThick=Control.EvalVar<double>(keyName+"OutWallThick");

  boxFront=Control.EvalVar<double>(keyName+"BoxFront");
  boxBack=Control.EvalVar<double>(keyName+"BoxBack");
  boxIn=Control.EvalVar<double>(keyName+"BoxIn");
  boxOut=Control.EvalVar<double>(keyName+"BoxOut");
  boxTop=Control.EvalVar<double>(keyName+"BoxTop");
  boxBase=Control.EvalVar<double>(keyName+"BoxBase");
  boxInterior=Control.EvalVar<double>(keyName+"BoxInterior");
  steelThick=Control.EvalVar<double>(keyName+"SteelThick");

  pbDepth=Control.EvalVar<double>(keyName+"PbDepth");
  pbThick=Control.EvalVar<double>(keyName+"PbThick");

  beMat=ModelSupport::EvalMat<int>(Control,keyName+"BeMat");  
  wallMat=ModelSupport::EvalMat<int>(Control,keyName+"WallMat");  
  waterMat=ModelSupport::EvalMat<int>(Control,keyName+"WaterMat");  

  vesselMat=ModelSupport::EvalMat<int>(Control,keyName+"VesselMat");
  polyMat=ModelSupport::EvalMat<int>(Control,keyName+"PolyMat");
  pbMat=ModelSupport::EvalMat<int>(Control,keyName+"PbMat");

  xStep+=radius-beLength/2.0;

  return;
}

void
RotorWheel::createSurfaces()
  /*!
    Create planes for the silicon and Polyethene layers
  */
{
  ELog::RegMethod RegA("RotorWheel","createSurfaces");

  ModelSupport::buildPlane(SMap,buildIndex+1,Origin-Y*waterThick,Y);  
  ModelSupport::buildPlane(SMap,buildIndex+2,Origin+Y*waterThick,Y);  
  ModelSupport::buildCylinder(SMap,buildIndex+7,Origin,Y,radius);  


  // Be Segment
  ModelSupport::buildPlane(SMap,buildIndex+101,
			   Origin-Y*(waterThick+beThick),Y);  
  ModelSupport::buildCylinder(SMap,buildIndex+8,
			      Origin,Y,radius-beLength);  

  // inner wall
  ModelSupport::buildPlane(SMap,buildIndex+11,
			   Origin-Y*(waterThick+wallThick),Y);  
  ModelSupport::buildPlane(SMap,buildIndex+12,
			   Origin+Y*(waterThick+wallThick),Y);  
  ModelSupport::buildCylinder(SMap,buildIndex+17,Origin,Y,radius+wallThick);  
  ModelSupport::buildCylinder(SMap,buildIndex+18,Origin,
			      Y,radius-(beLength+wallThick));  
  
  // VacLayer
  ModelSupport::buildPlane(SMap,buildIndex+21,
			   Origin-Y*(waterThick+wallThick+frontVac),Y);  
  ModelSupport::buildPlane(SMap,buildIndex+22,
			   Origin+Y*(waterThick+wallThick+backVac),Y);  
  ModelSupport::buildCylinder(SMap,buildIndex+27,
			      Origin,Y,radius+wallThick+endVac);  

  // Wrapper layer
  ModelSupport::buildPlane(SMap,buildIndex+31,
			   Origin-Y*(waterThick+wallThick+
				     frontVac+outWallThick),Y);  

  ModelSupport::buildPlane(SMap,buildIndex+32,
			   Origin+Y*(waterThick+wallThick+
				     backVac+outWallThick),Y);  

  ModelSupport::buildCylinder(SMap,buildIndex+37,
			      Origin,Y,radius+wallThick+
			      endVac+outWallThick);  

  ModelSupport::buildPlane(SMap,buildIndex+41,Origin-Y*(boxFront+steelThick),Y);
  ModelSupport::buildPlane(SMap,buildIndex+42,Origin+Y*(boxBack+steelThick),Y);
  ModelSupport::buildPlane(SMap,buildIndex+43,Origin-X*boxIn,X);
  ModelSupport::buildPlane(SMap,buildIndex+44,Origin+X*(boxOut+steelThick),X);
  ModelSupport::buildPlane(SMap,buildIndex+45,Origin-Z*(boxTop+steelThick),Z);
  ModelSupport::buildPlane(SMap,buildIndex+46,Origin+Z*(boxBase+steelThick),Z);
  ModelSupport::buildPlane(SMap,buildIndex+49,Origin+X*boxInterior,X);

  ModelSupport::buildPlane(SMap,buildIndex+51,Origin-Y*boxFront,Y);
  ModelSupport::buildPlane(SMap,buildIndex+52,Origin+Y*boxBack,Y);
  //  ModelSupport::buildPlane(SMap,buildIndex+53,Origin-X*boxIn,X);
  ModelSupport::buildPlane(SMap,buildIndex+54,Origin+X*boxOut,X);
  ModelSupport::buildPlane(SMap,buildIndex+55,Origin-Z*boxTop,Z);
  ModelSupport::buildPlane(SMap,buildIndex+56,Origin+Z*boxBase,Z);

  // Lead Wall

  ModelSupport::buildPlane(SMap,buildIndex+61,Origin+X*pbDepth,X);
  ModelSupport::buildPlane(SMap,buildIndex+62,Origin+X*(pbDepth+pbThick),X);


  return; 
}

void
RotorWheel::createObjects(Simulation& System)
  /*!
    Create the vaned moderator
    \param System :: Simulation to add results
   */
{
  ELog::RegMethod RegA("RotorWheel","createObjects");
  
  std::string Out;

  // Water
  Out=ModelSupport::getComposite(SMap,buildIndex,"1 -2 -7");
  System.addCell(MonteCarlo::Object(cellIndex++,waterMat,0.0,Out));

  // Be Segment
  Out=ModelSupport::getComposite(SMap,buildIndex,"-1 101 -7 8");	
  System.addCell(MonteCarlo::Object(cellIndex++,beMat,0.0,Out));

  // Wall layer [Front/Back]
  Out=ModelSupport::getComposite(SMap,buildIndex,"2 -12 -7");	
  System.addCell(MonteCarlo::Object(cellIndex++,wallMat,0.0,Out));

  // Front cut short by Be
  Out=ModelSupport::getComposite(SMap,buildIndex,"-1 11 -18");	
  System.addCell(MonteCarlo::Object(cellIndex++,wallMat,0.0,Out));
  // Be inside wall
  Out=ModelSupport::getComposite(SMap,buildIndex,"-1 101 -8 18");	
  System.addCell(MonteCarlo::Object(cellIndex++,wallMat,0.0,Out));

  // Wall Be outsde (full length)
  Out=ModelSupport::getComposite(SMap,buildIndex,"101 -12 7 -17");	
  System.addCell(MonteCarlo::Object(cellIndex++,wallMat,0.0,Out));

  // Part without anthing in middle
  Out=ModelSupport::getComposite(SMap,buildIndex,"-18 -11 101");	
  System.addCell(MonteCarlo::Object(cellIndex++,0,0.0,Out));

  // Part without anthing in middle
  Out=ModelSupport::getComposite(SMap,buildIndex,"(17:-101:12) 21 -27 -22");	
  System.addCell(MonteCarlo::Object(cellIndex++,0,0.0,Out));

  // Part without anything in middle
  Out=ModelSupport::getComposite(SMap,buildIndex,"(27:-21:22) 31 -37 -32");	
  System.addCell(MonteCarlo::Object(cellIndex++,0,0.0,Out));

  Out=ModelSupport::getComposite(SMap,buildIndex,"41 -42 49 -44 45 -46 "
				 " (-51:52:54:-55:56)");
  System.addCell(MonteCarlo::Object(cellIndex++,vesselMat,0,Out));

  Out=ModelSupport::getComposite(SMap,buildIndex,"51 -52 55 -56 61 -62");
  System.addCell(MonteCarlo::Object(cellIndex++,pbMat,0.0,Out));

  Out=ModelSupport::getComposite(SMap,buildIndex,"51 -52 55 -56 62 -54");
  System.addCell(MonteCarlo::Object(cellIndex++,polyMat,0.0,Out));

  Out=ModelSupport::getComposite(SMap,buildIndex,"51 -52 43 -61 55 -56 (37:-31:32)");
  System.addCell(MonteCarlo::Object(cellIndex++,polyMat,0.0,Out));


  Out=ModelSupport::getComposite(SMap,buildIndex,"41 -42 43 -49 45 -46 "
                                      "(-51:52:54 :-55:56)");
  System.addCell(MonteCarlo::Object(cellIndex++,polyMat,0.0,Out));


  Out=ModelSupport::getComposite(SMap,buildIndex,"-37 31 -32");
  addOuterSurf("Target",Out);
  Out=ModelSupport::getComposite(SMap,buildIndex,"41 -42 43 -44 45 -46");
  addOuterSurf("Body",Out);


  return; 
}

void
RotorWheel::createLinks()
  /*!
    Creates a full attachment set
    Surfaces pointing outwards
  */
{
  ELog::RegMethod RegA("RotorWheel","createLinks");
  // set Links :: Inner links:
  // Wrapper layer

  const double outDist(waterThick+wallThick+frontVac+outWallThick); /** Estos son variables de ayuda creadas para simplificar la definicion de los puntos. */
  const double outRadius(radius+wallThick+endVac+outWallThick);

  FixedComp::setConnect(0,Origin-Y*outDist,-Y);
  FixedComp::setLinkSurf(0,-SMap.realSurf(buildIndex+31));

  FixedComp::setConnect(1,Origin+Y*outDist,Y);
  FixedComp::setLinkSurf(1,SMap.realSurf(buildIndex+32));

  FixedComp::setConnect(2,Origin-X*outRadius,-X);
  FixedComp::setLinkSurf(2,SMap.realSurf(buildIndex+37));

  FixedComp::setConnect(3,Origin+X*outRadius,X);
  FixedComp::setLinkSurf(3,SMap.realSurf(buildIndex+37));

  FixedComp::setConnect(4,Origin-Z*outRadius,-Z);
  FixedComp::setLinkSurf(4,SMap.realSurf(buildIndex+37));

  FixedComp::setConnect(5,Origin+Z*outRadius,Z);
  FixedComp::setLinkSurf(5,SMap.realSurf(buildIndex+37));


  FixedComp::setConnect(6,Origin+Y*outDist+X*(0.1*outRadius),X);
  FixedComp::setLinkSurf(6,-SMap.realSurf(buildIndex+31));

  FixedComp::setConnect(7,Origin+Y*outDist+X*(0.5*outRadius),X);
  FixedComp::setLinkSurf(7,-SMap.realSurf(buildIndex+31));

  FixedComp::setConnect(8,Origin-Y*(waterThick+beThick)-X*(radius-beLength/2),
			-Y);
  FixedComp::setLinkSurf(8,-SMap.realSurf(buildIndex+31));



// /*
//  * CENTRO DE LA LAMINA DE BERILIO
//  *
//  *
// */
  
   FixedComp::setConnect(9,Origin+Y*outDist+X*(0.8*outRadius),X);
   FixedComp::setLinkSurf(9,-SMap.realSurf(buildIndex+31));

   FixedComp::setConnect(10,Origin+Y*outDist+X*(0.1*outRadius),X);
   FixedComp::setLinkSurf(10,-SMap.realSurf(buildIndex+31));

   FixedComp::setConnect(11,Origin+Y*outDist,Y);
   FixedComp::setLinkSurf(11,-SMap.realSurf(buildIndex+32));

   FixedComp::setConnect(12,Origin-X*(outRadius-beLength/2)+Y*outDist,-X);
   FixedComp::setLinkSurf(12,SMap.realSurf(buildIndex+32));

   FixedComp::setConnect(13,Origin-Y*(waterThick+beThick)-X*(radius-beLength/2),Y);
   FixedComp::setLinkSurf(13,-SMap.realSurf(buildIndex+31));

   FixedComp::setConnect(14,Origin+X*outRadius-Y*outDist,X);
   FixedComp::setLinkSurf(14,SMap.realSurf(buildIndex+37));

  return;
}

void
RotorWheel::createAll(Simulation& System,
		      const attachSystem::FixedComp& FC,
		      const long int sideIndex)
  /*!
    Extrenal build everything
    \param System :: Simulation
    \param FC :: FixedComponent for origin
    \param sideIndex :: linke point
   */
{
  ELog::RegMethod RegA("RotorWheel","createAll");
  populate(System.getDataBase());

  createUnitVector(FC,sideIndex);
  createSurfaces();
  createObjects(System);

  createLinks();
  insertObjects(System);       

  return;
}

}  // NAMESPACE instrumentSystem
