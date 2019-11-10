/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   scatMatInc/DBXrayMaterial.h
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
#ifndef scatterSystem_DBXrayMaterial_h
#define scatterSystem_DBXrayMaterial_h

namespace scatterSystem
{

/*!
  \class DBXrayMaterial
  \version 1.0
  \author S. Ansell
  \date November 2019
  \brief Storage for true xray - scattering materials
*/

class DBXrayMaterial 
{  
 private:

  /// Storage type for Materials
  typedef std::map<int,scatterSystem::xrayMaterial*> MTYPE;
  /// Store of materials
  MTYPE  MStore;

  DBXrayMaterial();

  ////\cond SINGLETON
  DBXrayMaterial(const DBXrayMaterial&);
  DBXrayMaterial& operator=(const DBXrayMaterial&);
  ////\endcond SINGLETON

  void initMaterial();
  
 public:
  
  static DBXrayMaterial& Instance();
  
  ~DBXrayMaterial() {}  ///< Destructor
  
  /// get data store
  const MTYPE& getStore() const { return MStore; }

  const scatterSystem::xrayMaterial* getMat(const int) const;

  void write(std::ostream& OX) const;
};

}

#endif
