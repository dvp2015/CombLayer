/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   chip/ChipIRFilter.cxx
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
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "support.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "localRotate.h"
#include "masterRotate.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "surfEqual.h"
#include "surfDivide.h"
#include "surfDIter.h"
#include "Quadratic.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "Object.h"
#include "groupRange.h"
#include "objectGroups.h"
#include "Simulation.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "chipDataStore.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "LinearComp.h"
#include "ContainedComp.h"
#include "ChipIRFilter.h"

namespace hutchSystem
{

ChipIRFilter::ChipIRFilter(const std::string& Key)  :
  attachSystem::ContainedComp(),attachSystem::FixedOffset(Key,6)

  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param Key :: Name for item in search
  */
{}

ChipIRFilter::ChipIRFilter(const ChipIRFilter& A) : 
  attachSystem::ContainedComp(A),attachSystem::FixedOffset(A),
  outerLen(A.outerLen),
  outerWidth(A.outerWidth),outerHeight(A.outerHeight),nLayers(A.nLayers),
  fracFrac(A.fracFrac),fracMat(A.fracMat)
  /*!
    Copy constructor
    \param A :: ChipIRFilter to copy
  */
{}

ChipIRFilter&
ChipIRFilter::operator=(const ChipIRFilter& A)
  /*!
    Assignment operator
    \param A :: ChipIRFilter to copy
    \return *this
  */
{
  if (this!=&A)
    {
      attachSystem::ContainedComp::operator=(A);
      attachSystem::FixedOffset::operator=(A);
      outerLen=A.outerLen;
      outerWidth=A.outerWidth;
      outerHeight=A.outerHeight;
      nLayers=A.nLayers;
      fracFrac=A.fracFrac;
      fracMat=A.fracMat;
    }
  return *this;
}

ChipIRFilter::~ChipIRFilter() 
 /*!
   Destructor
 */
{}

void
ChipIRFilter::populate(const FuncDataBase& Control)
 /*!
   Populate all the variables
   \param Control :: DataBase
 */
{
  ELog::RegMethod RegA("ChipIRFilter","populate");

  FixedOffset::populate(Control);
 
  outerLen=Control.EvalVar<double>(keyName+"OuterLen");
  outerWidth=Control.EvalVar<double>(keyName+"OuterWidth");
  outerHeight=Control.EvalVar<double>(keyName+"OuterHeight");

  return;
}
  
void
ChipIRFilter::createSurfaces()
  /*!
    Create All the surfaces
  */
{
  ELog::RegMethod RegA("ChipIRFilter","createSurface");


  ModelSupport::buildPlane(SMap,buildIndex+1,Origin,Y);
  ModelSupport::buildPlane(SMap,buildIndex+2,Origin+Y*outerLen,Y);
  // Widths:
  ModelSupport::buildPlane(SMap,buildIndex+3,Origin-X*outerWidth/2.0,X);
  ModelSupport::buildPlane(SMap,buildIndex+4,Origin+X*outerWidth/2.0,X);
  
  ModelSupport::buildPlane(SMap,buildIndex+5,Origin-Z*outerHeight/2.0,Z);
  ModelSupport::buildPlane(SMap,buildIndex+6,Origin+Z*outerHeight/2.0,Z);

  return;
}

void
ChipIRFilter::createObjects(Simulation& System)
  /*!
    Adds the Chip guide components
    \param System :: Simulation to create objects in
  */
{
  ELog::RegMethod RegA("ChipIRFilter","createObjects");
  
  std::string Out;
  Out=ModelSupport::getComposite(SMap,buildIndex,"1 -2 3 -4 5 -6");
  addOuterSurf(Out);

  // Outer layer
  Out=ModelSupport::getComposite(SMap,buildIndex,"1 -2 3 -4 5 -6");
  System.addCell(MonteCarlo::Object(cellIndex++,0,0.0,Out));
  return;
}

void
ChipIRFilter::createLinks()
  /*!
    Creates a full attachment set
  */
{
  // // set Links:
  // FixedComp::setConnect(0,BVec[0]+Y*(vacPosGap+alPos+terPos+
  // 				     outPos+clearNeg),Y);
  // FixedComp::setConnect(1,BVec[1]-Y*(vacNegGap+alPos+terNeg+
  // 				     outNeg+clearNeg),-Y);
  // FixedComp::setConnect(2,BVec[2]-
  // 			X*(vacSide+alSide+terSide+outSide+clearSide),-X);
  // FixedComp::setConnect(3,BVec[3]+
  // 			X*(vacSide+alSide+terSide+outSide+clearSide),X);
  // FixedComp::setConnect(4,BVec[4]-
  // 			Z*(vacBase+alBase+terBase+outBase+clearBase),-Z);
  // FixedComp::setConnect(5,BVec[5]+
  // 			Z*(vacTop+alTop+terTop+outTop+clearTop),Z);

  // // Set Connect surfaces:
  // for(int i=2;i<6;i++)
  //   FixedComp::setLinkSurf(i,SMap.realSurf(vacIndex+41+i));

  // // For Cylindrical surface must also have a divider:
  // // -- Groove:
  // FixedComp::setLinkSurf(0,SMap.realSurf(vacIndex+41));
  // FixedComp::addLinkSurf(0,SMap.realSurf(divideSurf));

  // FixedComp::setLinkSurf(1,SMap.realSurf(vacIndex+42));
  // FixedComp::addLinkSurf(1,-SMap.realSurf(divideSurf));

  return;
}

void
ChipIRFilter::createAll(Simulation& System,
			const attachSystem::FixedComp& FC,
			const long int sideIndex)
  /*!
    Global creation of the hutch
    \param System :: Simulation to add vessel to
    \param FC :: Linear object to component
    \param sideIndex :: link point
  */
{
  ELog::RegMethod RegA("ChipIRFilter","createAll");

  populate(System.getDataBase());

  createUnitVector(FC,sideIndex);
  createSurfaces();
  createObjects(System);
  createLinks();
  insertObjects(System);       
    
  return;
}

  
}  // NAMESPACE moderatorSystem
