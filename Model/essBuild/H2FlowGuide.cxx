/*********************************************************************
  CombLayer : MCNP(X) Input builder

 * File:   essBuild/H2FlowGuide.cxx
 *
 * Copyright (c) 2004-2017 by Konstantin Batkov
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
#include <complex>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <array>
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
#include "stringCombine.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "Quaternion.h"
#include "Surface.h"
#include "surfIndex.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "surfEqual.h"
#include "surfDivide.h"
#include "surfDIter.h"
#include "Quadratic.h"
#include "General.h"
#include "Plane.h"
#include "Cylinder.h"
#include "Line.h"
#include "Rules.h"
#include "varList.h"
#include "Code.h"
#include "FuncDataBase.h"
#include "HeadRule.h"
#include "RuleSupport.h"
#include "Object.h"
#include "Qhull.h"
#include "Simulation.h"
#include "Vert2D.h"
#include "Convex2D.h"
#include "ModelSupport.h"
#include "MaterialSupport.h"
#include "generateSurf.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "ContainedComp.h"
#include "LayerComp.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "AttachSupport.h"
#include "geomSupport.h"
#include "H2Wing.h"
#include "H2FlowGuide.h"

namespace essSystem
{

H2FlowGuide::H2FlowGuide(const std::string& baseKey,
			 const std::string& extraKey,
			 const std::string& finalKey ) :
  attachSystem::FixedComp(baseKey+extraKey+finalKey,6),
  baseName(baseKey),midName(extraKey),endName(finalKey),
  flowIndex(ModelSupport::objectRegister::Instance().cell(keyName)),
  cellIndex(flowIndex+1)
  /*!
    Constructor BUT ALL variable are left unpopulated.
    \param baseKey :: Butterfly main key
    \param extraKey :: H2Wing part name
    \param finalKey :: Specialized flow name
  */
{}

H2FlowGuide::H2FlowGuide(const H2FlowGuide& A) :
  attachSystem::FixedComp(A),
  baseName(A.baseName),midName(A.midName),endName(A.endName),
  flowIndex(A.flowIndex),
  cellIndex(A.cellIndex),
  wallThick(A.wallThick),len1L(A.len1L),
  baseOffset(A.baseOffset),
  len1R(A.len1R),
  angle1(A.angle1),
  dist2(A.dist2),
  len2L(A.len2L),
  len2R(A.len2R),
  dist3(A.dist3),
  len3L(A.len3L),
  len3R(A.len3R),
  sqCenterF(A.sqCenterF),
  wallMat(A.wallMat),
  wallTemp(A.wallTemp)
  /*!
    Copy constructor
    \param A :: H2FlowGuide to copy
  */
{}

H2FlowGuide&
H2FlowGuide::operator=(const H2FlowGuide& A)
  /*!
    Assignment operator
    \param A :: H2FlowGuide to copy
    \return *this
  */
{
  if (this!=&A)
    {
      attachSystem::FixedComp::operator=(A);
      cellIndex=A.cellIndex;
      wallThick=A.wallThick;
      len1L=A.len1L;
      baseOffset=A.baseOffset;
      len1R=A.len1R;
      angle1=A.angle1;
      dist2=A.dist2;
      len2L=A.len2L;
      len2R=A.len2R;
      dist3=A.dist3;
      len3L=A.len3L;
      len3R=A.len3R;
      sqCenterF=A.sqCenterF;
      wallMat=A.wallMat;
      wallTemp=A.wallTemp;
    }
  return *this;
}

H2FlowGuide*
H2FlowGuide::clone() const
  /*!
    Virtual copy constructor
    \return new (this)
   */
{
  return new H2FlowGuide(*this);
}

H2FlowGuide::~H2FlowGuide()
  /*!
    Destructor
  */
{}

void
H2FlowGuide::populate(const FuncDataBase& Control)
  /*!
    Populate all the variables
    \param Control :: DataBase for variables
  */
{
  ELog::RegMethod RegA("H2FlowGuide","populate");


  wallThick=Control.EvalPair<double>(keyName,baseName+endName,"WallThick");
  len1L=Control.EvalPair<double>(keyName,baseName+endName,"Len1L");
  len1R=Control.EvalPair<double>(keyName,baseName+endName,"Len1R");
  angle1=Control.EvalPair<double>(keyName,baseName+endName,"Angle1");
  baseOffset=Control.EvalPair<double>(keyName,baseName+endName,"BaseOffset");
  dist2=Control.EvalPair<double>(keyName,baseName+endName,"Dist2");
  len2L=Control.EvalPair<double>(keyName,baseName+endName,"Len2L");
  len2R=Control.EvalPair<double>(keyName,baseName+endName,"Len2R");
  dist3=Control.EvalPair<double>(keyName,baseName+endName,"Dist3");
  len3L=Control.EvalPair<double>(keyName,baseName+endName,"Len3L");
  len3R=Control.EvalPair<double>(keyName,baseName+endName,"Len3R");
  sqCenterF=Control.EvalPair<double>(keyName,baseName+endName,"SQCenterF");

  wallMat=ModelSupport::EvalMat<int>(Control,keyName+"WallMat",
				     baseName+endName+"WallMat");
  wallTemp=Control.EvalPair<double>(keyName,baseName+endName,"WallTemp");

  return;
}

void
H2FlowGuide::createUnitVector(const attachSystem::FixedComp& FC)
  /*!
    Create the unit vectors
    \param FC :: fixed Comp [and link comp]
  */
{
  ELog::RegMethod RegA("H2FlowGuide","createUnitVector");
  FixedComp::createUnitVector(FC);
  return;
}

void
H2FlowGuide::createSurfaces()
  /*!
    Create All the surfaces
  */
{
  ELog::RegMethod RegA("H2FlowGuide","createSurface");

  double y(baseOffset);
  ModelSupport::buildPlane(SMap,flowIndex+1,Origin+Y*(y),Y);
  y += wallThick;
  ModelSupport::buildPlane(SMap,flowIndex+2,Origin+Y*(y),Y);
  y += dist2;
  ModelSupport::buildPlane(SMap,flowIndex+11,Origin+Y*(y),Y);
  y += wallThick;
  ModelSupport::buildPlane(SMap,flowIndex+12,Origin+Y*(y),Y);
  y += dist3;
  ModelSupport::buildPlane(SMap,flowIndex+21,Origin+Y*(y),Y);
  y += wallThick;
  ModelSupport::buildPlane(SMap,flowIndex+22,Origin+Y*(y),Y);

  ModelSupport::buildPlane(SMap,flowIndex+3,Origin-X*(len1L),X);
  ModelSupport::buildPlane(SMap,flowIndex+4,Origin+X*(len1R),X);

  ModelSupport::buildPlane(SMap,flowIndex+103,Origin-X*(len2L),X);
  ModelSupport::buildPlane(SMap,flowIndex+104,Origin+X*(len2R),X);

  ModelSupport::buildPlane(SMap,flowIndex+203,Origin-X*(len3L),X);
  ModelSupport::buildPlane(SMap,flowIndex+204,Origin+X*(len3R),X);

  return;
}

void
H2FlowGuide::createObjects(Simulation& System,
			   const attachSystem::FixedComp& HW)
  /*!
    Adds the main components
    \param System :: Simulation to create objects in
    \param HW :: H2Wing object
  */
{
  ELog::RegMethod RegA("H2FlowGuide","createObjects");

  const attachSystem::CellMap* CM=
    dynamic_cast<const attachSystem::CellMap*>(&HW);
  MonteCarlo::Object* InnerObj(0);
  int innerCell(0);
  if (CM)
    {
      innerCell=CM->getCell("Inner");
      InnerObj=System.findQhull(innerCell);
    }
  if (!InnerObj)
    throw ColErr::InContainerError<int>
      (innerCell,"H2Wing inner Cell not found");

  std::string Out;
  const std::string tb(HW.getLinkString(13)+HW.getLinkString(14));
  HeadRule wallExclude;

  Out=ModelSupport::getComposite(SMap,flowIndex," 1 -2 3 -4 ");
  wallExclude.procString(Out);
  System.addCell(MonteCarlo::Qhull(cellIndex++,wallMat,wallTemp,Out+tb));

  Out=ModelSupport::getComposite(SMap,flowIndex," 11 -12 103 -104 ");
  wallExclude.addUnion(Out);
  System.addCell(MonteCarlo::Qhull(cellIndex++,wallMat,wallTemp,Out+tb));

  Out=ModelSupport::getComposite(SMap,flowIndex," 21 -22 203 -204 ");
  wallExclude.addUnion(Out);
  System.addCell(MonteCarlo::Qhull(cellIndex++,wallMat,wallTemp,Out+tb));

  wallExclude.makeComplement();
  InnerObj->addSurfString(wallExclude.display());

  return;
}


void
H2FlowGuide::createAll(Simulation& System,
		       const attachSystem::FixedComp& FC)
  /*!
    Generic function to create everything
    \param System :: Simulation item
    \param FC :: The H2Wing for the flow
  */
{
  ELog::RegMethod RegA("H2FlowGuide","createAll");

  populate(System.getDataBase());
  createUnitVector(FC);
  createSurfaces();
  createObjects(System,FC);

  return;
}

}  // NAMESPACE essSystem
