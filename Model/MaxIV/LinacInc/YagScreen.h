/*********************************************************************
  CombLayer : MCNP(X) Input builder

 * File:   Model/MaxIV/LinacInc/YagScreen.h
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
#ifndef tdcSystem_YagScreen_h
#define tdcSystem_YagScreen_h

class Simulation;

namespace tdcSystem
{

/*!
  \class YagScreen
  \version 1.0
  \author Konstantin Batkov
  \date May 2020
  \brief Yag screen
*/

class YagScreen : public attachSystem::ContainedGroup,
		  public attachSystem::FixedRotate,
		  public attachSystem::CellMap
{
 private:

  double jbLength;              ///< electronics junction box length
  double jbWidth;               ///< electronics junction box width
  double jbHeight;              ///< electronics junction box height
  double jbWallThick;           ///< electronics junction box wall thickness
  int    jbWallMat;             ///< electronics junction box wall material
  int    jbMat;                 ///< electronics junction box cable/inside material

  double ffLength;              ///< linear pneumatics feedthrough length
  double ffInnerRadius;         ///< linear pneumatics feedthrough inner radius
  double ffWallThick;           ///< linear pneumatics feedthrough wall thickness
  double ffFlangeLen;           ///< linear pneumatics feedthrough flange length
  double ffFlangeRadius;        ///< linear pneumatics feedthrough flange radius
  int    ffWallMat;             ///< linear pneumatics feedthrough wall material

  double holderLift;            ///< screen holder lift
  double holderRad;             ///< screen holder inner radius
  int    holderMat;             ///< screen holder material
  double mirrorRadius;          ///< quartz mirror radius
  double mirrorAngle;           ///< quartz mirror inclination angle [deg]
  double mirrorThick;           ///< quartz mirror thickness
  int    mirrorMat;             ///< mirror  material

  int    voidMat;               ///< void material
  bool   closed;                ///< screen and mirror are in the beam

  void populate(const FuncDataBase&);
  void createUnitVector(const attachSystem::FixedComp&,
			const long int);

  void createSurfaces();
  void createObjects(Simulation&);
  void createLinks();

 public:

  YagScreen(const std::string&);
  YagScreen(const YagScreen&);
  YagScreen& operator=(const YagScreen&);
  virtual YagScreen* clone() const;
  virtual ~YagScreen();

  void createAll(Simulation&,const attachSystem::FixedComp&,const long int);

};

}

#endif