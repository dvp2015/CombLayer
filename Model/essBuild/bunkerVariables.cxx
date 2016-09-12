/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   essBuild/bunkerVariables.cxx
 *
 * Copyright (c) 2004-2016 by Stuart Ansell/Konstantin Batkov
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
#include "support.h"
#include "stringCombine.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Code.h"
#include "varList.h"
#include "FuncDataBase.h"
#include "variableSetup.h"
#include "essVariables.h"

namespace setVariable
{

void
setBunkerVar(FuncDataBase&,const std::string&,const std::string&);



void
setBunkerVar(FuncDataBase& Control,const std::string& AKey,
	     const std::string& BKey)
/*!
    Create all the bunker variables
    \param Control :: DataBase
    \param AKey :: Key for A/C bunker
    \param BKey :: Key for B/D bunker
*/
{
  ELog::RegMethod RegA("bunkerVariables[F]","setBunkerVar");


  Control.addVariable(AKey+"BunkerLeftPhase",-65.0);
  Control.addVariable(AKey+"BunkerRightPhase",0.0);
  Control.addVariable(AKey+"BunkerLeftAngle",0.0);
  Control.addVariable(AKey+"BunkerRightAngle",0.0);
  Control.addVariable(AKey+"BunkerNSectors",10);
  Control.addVariable(AKey+"BunkerNVert",1);
  Control.addVariable(AKey+"BunkerNLayers",1);
  Control.addVariable(AKey+"BunkerNSegment",5);
  Control.addVariable(AKey+"BunkerRoofActive",0);
  Control.addVariable(AKey+"BunkerWallRadius",2450.0);
  Control.addVariable(BKey+"BunkerWallRadius",1150.0);
  
  for(const std::string& KItem : {AKey,BKey})
    {
      Control.addVariable(KItem+"BunkerWallActive",0);

      Control.addVariable(KItem+"BunkerWallNVert",20);
      Control.addVariable(KItem+"BunkerWallNMedial",2);
      Control.addVariable(KItem+"BunkerWallNRadial",12);

      Control.addVariable(KItem+"BunkerRoofNVert",20);
      Control.addVariable(KItem+"BunkerRoofNMedial",1);
      Control.addVariable(KItem+"BunkerRoofNRadial",12);
      Control.addVariable(KItem+"BunkerRoofNBasicVert",5);
      
      
      Control.addVariable(KItem+"BunkerNSide",5);
      Control.addVariable(KItem+"BunkerNSideVert",5);
      Control.addVariable(KItem+"BunkerNSideThick",5);
      

      Control.addVariable(KItem+"BunkerFloorDepth",120.0);
      Control.addVariable(KItem+"BunkerRoofHeight",190.0);
      
      Control.addVariable(KItem+"BunkerWallThick",320.0);
      Control.addVariable(KItem+"BunkerSideThick",80.0);
      Control.addVariable(KItem+"BunkerRoofThick",195.0);
      Control.addVariable(KItem+"BunkerFloorThick",100.0);
      
      Control.addVariable(KItem+"BunkerVoidMat","Void");
      Control.addVariable(KItem+"BunkerWallMat","Steel71");
      Control.addVariable(KItem+"BunkerRoofMat","Aluminium");


      // WALL LAYERED
      Control.addVariable(KItem+"BunkerWallNBasic",5);      

      Control.addVariable(KItem+"BunkerNBasicVert",5);
      Control.addVariable(KItem+"BunkerRoofVert1",20.0);
      Control.addVariable(KItem+"BunkerRoofVert2",40.0);
      Control.addVariable(KItem+"BunkerRoofVert3",70.0);
      Control.addVariable(KItem+"BunkerRoofVert4",30.0);
  
      
      Control.addVariable(KItem+"BunkerRoofMat0","CarstonConc");
      Control.addVariable(KItem+"BunkerRoofMat1","CastIron");
      Control.addVariable(KItem+"BunkerRoofMat2","CarstonConc");
      Control.addVariable(KItem+"BunkerRoofMat3","CastIron");
      Control.addVariable(KItem+"BunkerRoofMat4","CarstonConc");

      Control.addVariable(KItem+"BunkerWallMat0","CarstonConc");
      Control.addVariable(KItem+"BunkerWallMat1","CastIron");
      Control.addVariable(KItem+"BunkerWallMat2","CarstonConc");
      Control.addVariable(KItem+"BunkerWallMat3","CastIron");
      Control.addVariable(KItem+"BunkerWallMat4","CarstonConc");
      
      Control.addVariable(KItem+"BunkerWallLen1",0.5);
      Control.addVariable(KItem+"BunkerWallLen2",10.0);
    }
  
  Control.addVariable(AKey+"BunkerLoadFile","ABunkerDef.xml");
  Control.addVariable(AKey+"BunkerOutFile","ABunker.xml");

  // PILLARS:
  
  Control.addVariable(AKey+"BunkerPillarsRadius",6.0);
  Control.addVariable(AKey+"BunkerPillarsMat","Stainless304");
  Control.addVariable(AKey+"BunkerPillarsNRadius",5);
  Control.addVariable(AKey+"BunkerPillarsNSector",11); // default
  
  Control.addVariable(AKey+"BunkerPillarsR0",600.0);
  Control.addVariable(AKey+"BunkerPillarsR1",900.0);
  Control.addVariable(AKey+"BunkerPillarsR2",1200.0);
  Control.addVariable(AKey+"BunkerPillarsR3",1500.0);
  Control.addVariable(AKey+"BunkerPillarsR4",1800.0);
  Control.addVariable(AKey+"BunkerPillarsR5",2100.0);
  
  Control.addVariable(AKey+"BunkerPillarsRS_0",-3.0);
  Control.addVariable(AKey+"BunkerPillarsRS_1",-9.0);
  Control.addVariable(AKey+"BunkerPillarsRS_2",-15.0);
  Control.addVariable(AKey+"BunkerPillarsRS_3",-21.0);
  Control.addVariable(AKey+"BunkerPillarsRS_4",-27.0);
  Control.addVariable(AKey+"BunkerPillarsRS_5",-34.0);
  Control.addVariable(AKey+"BunkerPillarsRS_6",-39.0);
  Control.addVariable(AKey+"BunkerPillarsRS_7",-46.0);
  Control.addVariable(AKey+"BunkerPillarsRS_8",-52.0);
  Control.addVariable(AKey+"BunkerPillarsRS_9",-57.0);
  Control.addVariable(AKey+"BunkerPillarsRS_10",-63.0);

  Control.addVariable(AKey+"BunkerPillarsR_0S_1Active",0);
  Control.addVariable(AKey+"BunkerPillarsR_0S_2Active",0);
  Control.addVariable(AKey+"BunkerPillarsR_0S_9Active",0);
  Control.addVariable(AKey+"BunkerPillarsR_1S_2Active",0);


  //
  // RIGHT BUNKER : B PART
  //
  Control.addVariable(BKey+"BunkerLeftPhase",0.0);
  Control.addVariable(BKey+"BunkerRightPhase",65.0);
  Control.addVariable(BKey+"BunkerLeftAngle",0.0);
  Control.addVariable(BKey+"BunkerRightAngle",0.0);
  Control.addVariable(BKey+"BunkerNSectors",9);
  Control.addVariable(BKey+"BunkerNVert",1);
  Control.addVariable(BKey+"BunkerNLayers",1);

  
  Control.addVariable(BKey+"BunkerActiveSegment",0);
  Control.addVariable(BKey+"BunkerActiveRoof",0);

  Control.addVariable(BKey+"BunkerLoadFile","BBunkerDef.xml");
  Control.addVariable(BKey+"BunkerOutFile","BBunker.xml");
  
  return;
}

void
EssBunkerVariables(FuncDataBase& Control)
  /*!
    Create all the bunker variables
    \param Control :: DataBase
  */
{
  ELog::RegMethod RegA("essVariables[F]","EssBunkerVariables");
  setBunkerVar(Control,"A","B");
  setBunkerVar(Control,"C","D");

  // BUNKER QUAKE:
  Control.addVariable("ABunkerQuakeXGap",13.0);
  Control.addVariable("ABunkerQuakeZGap",100.0);
  
  Control.addVariable("ABunkerQuakeNPoint",1);
  Control.addVariable("ABunkerQuakePtA0",Geometry::Vec3D(750,0,0));
  Control.addVariable("ABunkerQuakePtB0",Geometry::Vec3D(980,0,0));

  
  
  Control.addVariable("CurtainLeftPhase",-65.0);
  Control.addVariable("CurtainRightPhase",65.0);
  Control.addVariable("CurtainWallThick",80.0);
  Control.addVariable("CurtainInnerStep",0.0);
  Control.addVariable("CurtainTopRaise",60.0);
  Control.addVariable("CurtainHeight",250.0);
  Control.addVariable("CurtainDepth",129.0);     // fixed by Ben.
  Control.addVariable("CurtainWallMat","Concrete");
  Control.addVariable("CurtainNTopLayers",7);
  Control.addVariable("CurtainNMidLayers",2);
  Control.addVariable("CurtainNBaseLayers",7);


  Control.addVariable("BunkerFeedOffset",Geometry::Vec3D(0,0,0));
  Control.addVariable("BunkerFeedHeight",6.0);
  Control.addVariable("BunkerFeedWidth",6.0);
  Control.addVariable("BunkerFeedTrack0",Geometry::Vec3D(0,0,0));
  Control.addVariable("BunkerFeedTrack1",Geometry::Vec3D(0,0,-10));
  Control.addVariable("BunkerFeedTrack2",Geometry::Vec3D(0,100,-10));

  return;
}
  
}  // NAMESPACE setVariable
