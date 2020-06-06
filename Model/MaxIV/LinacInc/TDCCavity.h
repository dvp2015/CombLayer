/*********************************************************************
  CombLayer : MCNP(X) Input builder

 * File:   Model/MaxIV/LinacInc/TDCCavity.h
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
#ifndef tdcSystem_TDCCavity_h
#define tdcSystem_TDCCavity_h

class Simulation;

namespace tdcSystem
{

/*!
  \class TDCCavity
  \version 1.0
  \author Konstantin Batkov
  \date June 2020
  \brief TDC cavity section
*/

class TDCCavity :
    public attachSystem::ContainedComp,
    public attachSystem::FixedRotate,
    public attachSystem::CellMap,
    public attachSystem::SurfMap,
    public attachSystem::FrontBackCut
{
 private:

  int    nCells;                ///< Number of regular cells
  double cellLength;            ///< Normal cell total length (void+iris)
  double cellRadius;            ///< Normal cell inner radius
  double irisLength;            ///< Iris length
  double irisRadius;            ///< Iris inner radius
  double couplerLength;         ///< Coupler cell length
  double couplerWidth;          ///< Coupler cell width
  double wallThick;             ///< Wall thickness
  int    wallMat;               ///< Wall material
  double flangeLength;          ///< Flange length
  double flangeRadius;          ///< Flange outer radius
  int flangeMat;                ///< Flange material

  void populate(const FuncDataBase&);
  void createSurfaces();
  void createObjects(Simulation&);
  void createLinks();

 public:

  TDCCavity(const std::string&);
  TDCCavity(const TDCCavity&);
  TDCCavity& operator=(const TDCCavity&);
  virtual TDCCavity* clone() const;
  virtual ~TDCCavity();

  using FixedComp::createAll;
  void createAll(Simulation&,const attachSystem::FixedComp&,const long int);

};

}

#endif