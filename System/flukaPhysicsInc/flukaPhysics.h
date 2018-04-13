/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   physicsInc/flukaPhysics.h
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
#ifndef flukaSystem_flukaPhysics_h
#define flukaSystem_flukaPhysics_h

namespace flukaSystem
{
  
/*!
  \class flukaPhysics
  \version 1.0
  \date March 2018
  \author S.Ansell
  \brief Processes the physics cards in the FLUKA output
*/

class flukaPhysics 
{
 private:

  /// storage for cell system
  typedef std::tuple<int,std::string> unitTYPE;
  
  std::vector<int> cellVec;                     ///< Cell numbers in order
  std::vector<int> matVec;                      ///< Material numbers in order

  // ONE VALUE string
  std::map<std::string,strValueSet<1>> cutValue;   ///< cut values

  // ZERO VALUE
  std::map<std::string,cellValueSet<0>> flagValue;   ///< flag values
  // ONE VALUE
  std::map<std::string,cellValueSet<1>> impValue;   ///< Importance values
  // TWO VALUE
  std::map<std::string,cellValueSet<2>> emfFlag;    ///< EMF flag
  // THREE VALUE
  std::map<std::string,cellValueSet<3>> threeFlag;    ///< EMF flag

  std::map<std::string,unitTYPE> formatMap;     ///< Layout and cell/mat

 public:
   
  flukaPhysics();
  flukaPhysics(const flukaPhysics&);
  flukaPhysics& operator=(const flukaPhysics&);
  virtual ~flukaPhysics();

  void clearAll();
  
  // ALL systems setup
  void setCellNumbers(const std::vector<int>&);
  // ALL systems setup
  void setMatNumbers(const std::set<int>&);

  void setFlag(const std::string&,const int);
  void setImp(const std::string&,const int,const std::string&);
  void setEMF(const std::string&,const int,const std::string&,
	      const std::string&);
  void setTHR(const std::string&,const int,const std::string&,
	      const std::string&,const std::string&);
  void writeFLUKA(std::ostream&) const;
};

}

#endif
