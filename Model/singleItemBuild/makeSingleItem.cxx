/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   singleItemBuild/makeSingleItem.cxx
 *
 * Copyright (c) 2004-2018 by Stuart Ansell
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
#include <utility>
#include <cmath>
#include <complex>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <array>
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
#include "stringCombine.h"
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
#include "Simulation.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "FixedGroup.h"
#include "FixedOffsetGroup.h"
#include "ContainedComp.h"
#include "ContainedSpace.h"
#include "ContainedGroup.h"
#include "FrontBackCut.h"
#include "LayerComp.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "SurfMap.h"
#include "World.h"
#include "AttachSupport.h"

#include "Cryostat.h"
#include "TwinBase.h"
#include "TwinChopper.h"
#include "TwinChopperFlat.h"
#include "SingleChopper.h"
#include "DiskChopper.h"
#include "VacuumPipe.h"

#include "CryoMagnetBase.h"

#include "makeSingleItem.h"

namespace singleItemSystem
{

makeSingleItem::makeSingleItem() 
 /*!
    Constructor
 */
{}


makeSingleItem::~makeSingleItem()
  /*!
    Destructor
  */
{}
  
void 
makeSingleItem::build(Simulation& System,
	       const mainSystem::inputParam& IParam)
  /*!
    Carry out the full build
    \param System :: Simulation system
    \param IParam :: Input parameters
   */
{
  // For output stream
  ELog::RegMethod RegA("makeSingleItem","build");

  int voidCell(74123);

  constructSystem::SingleChopper AS("singleChopper");
  AS.addInsertCell(voidCell);
  AS.createAll(System,World::masterOrigin(),0);

//  constructSystem::CryoMagnetBase A("CryoB");
//  A.addInsertCell(voidCell);
//  A.createAll(System,World::masterOrigin(),0);
}


}   // NAMESPACE singleItemSystem

