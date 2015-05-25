/********************************************************************* 
  CombLayer : MNCPX Input builder
 
 * File:   essBuildInc/DiskPreMod.h
 *
 * Copyright (c) 2004-2015 by Stuart Ansell
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
#ifndef essSystem_DiskPreMod_h
#define essSystem_DiskPreMod_h

class Simulation;

namespace essSystem
{
/*!
  \class DiskPreMod
  \author S. Ansell
  \version 1.0
  \date May 2015
  \brief Specialized for a cylinder pre-mod under moderator
*/

class DiskPreMod : public attachSystem::ContainedComp,
    public attachSystem::LayerComp,
    public attachSystem::FixedComp
{
 private:
  
  const int modIndex;             ///< Index of surface offset
  int cellIndex;                  ///< Cell index

  double zStep;                   ///< Step away from target
  double outerRadius;             ///< Outer radius of Be Zone
  
  std::vector<double> radius;         ///< cylinder radii [additive]
  std::vector<double> height;         ///< Full heights [additive]
  std::vector<double> depth;          ///< full depths [additive]
  std::vector<int> mat;               ///< Materials 
  std::vector<double> temp;           ///< Temperatures

  
  void populate(const FuncDataBase&,const double,const double);
  void createUnitVector(const attachSystem::FixedComp&,
			const bool);

  void createSurfaces();
  void createObjects(Simulation&);
  void createLinks();

 public:

  DiskPreMod(const std::string&);
  DiskPreMod(const DiskPreMod&);
  DiskPreMod& operator=(const DiskPreMod&);
  virtual DiskPreMod* clone() const;
  virtual ~DiskPreMod();


  virtual Geometry::Vec3D getSurfacePoint(const size_t,const size_t) const;
  virtual int getLayerSurf(const size_t,const size_t) const;
  virtual std::string getLayerString(const size_t,const size_t) const;

  /// total height of object
  double getHeight() const
    { return (depth.empty()) ? 0.0 : depth.back()+height.back(); }


  void createAll(Simulation&,const attachSystem::FixedComp&,
		 const bool,const double,const double);
  
};

}

#endif
 
