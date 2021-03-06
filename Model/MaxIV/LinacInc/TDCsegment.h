/*********************************************************************
  CombLayer : MCNP(X) Input builder

 * File:   LinacInc/TDCsegment.h
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
#ifndef tdcSystem_TDCsegment_h
#define tdcSystem_TDCsegment_h

namespace constructSystem
{
  class portItem;
  class BlankTube;
}

namespace tdcSystem
{
  /*!
    \class TDCsegment
    \version 1.0
    \author S. Ansell
    \date May 2020
    \brief Virtual top class for TDC segments
  */

class TDCsegment :
  public attachSystem::FixedRotate,
  public attachSystem::ContainedComp,
  public attachSystem::ExternalCut,
  public attachSystem::CellMap
{
 protected:

  /// System for building a divided inner
  attachSystem::InnerZone* buildZone;

  /// System for next building a divided inner
  attachSystem::InnerZone* nextZone;
  
  std::vector<HeadRule> joinItems;   ///< Stack of join items [multiface]

  /// unmanaged resource
  std::vector<attachSystem::ExternalCut*> firstItemVec;

 public:

  TDCsegment(const std::string&,const size_t);
  TDCsegment(const TDCsegment&);
  TDCsegment& operator=(const TDCsegment&);
  virtual ~TDCsegment();

  bool totalPathCheck(const FuncDataBase&,const double =0.1) const;

  /// set the current inner zone [allows joining of segments]
  void setInnerZone(attachSystem::InnerZone* IZPtr) { buildZone=IZPtr; }

  /// set the NEXT inner zone [allows joining of segments]
  void setNextZone(attachSystem::InnerZone* IZPtr)
    {  nextZone=IZPtr; }

  /// accessor to join items
  const std::vector<HeadRule>& getJoinItems() const
    { return joinItems; }
  
  const constructSystem::portItem&
  buildIonPump2Port(Simulation&,
		    attachSystem::InnerZone&,
		    MonteCarlo::Object*,
		    const attachSystem::FixedComp&,
		    const std::string&,
		    constructSystem::BlankTube&) const;

  void setFrontSurfs(const std::vector<HeadRule>&);
  void setFirstItems(const std::shared_ptr<attachSystem::FixedComp>&);
  void setFirstItems(attachSystem::FixedComp*);

  virtual void insertPrevSegment(Simulation&,const TDCsegment*) const {}
  
  virtual void createAll(Simulation&,const attachSystem::FixedComp&,
			 const long int) =0;

};

}

#endif
