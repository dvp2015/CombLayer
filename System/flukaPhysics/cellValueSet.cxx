/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   physics/cellValueSet.cxx
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
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <functional>
#include <memory>
#include <array>
#include <tuple>


#include "Exception.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "support.h"
#include "writeSupport.h"
#include "MapRange.h"
#include "Triple.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"

#include "particleConv.h"
#include "cellValueSet.h"

namespace flukaSystem
{
  
template<size_t N>
cellValueSet<N>::cellValueSet(const std::string& KN,const std::string& ON) :
  keyName(KN),outName(ON),whatValue(0.0)
  /*!
    Constructor
    \param KN :: Indentifier
    \param ON :: Output id for FLUKA
  */
{}

template<size_t N>
cellValueSet<N>::cellValueSet(const std::string& KN,
			      const std::string& ON,
			      const double WValue) :
  keyName(KN),outName(ON),whatValue(WValue)
  /*!
    Constructor
    \param KN :: Indentifier
    \param ON :: Output id for FLUKA
    \param WValue :: What value
  */
{}


template<size_t N>
cellValueSet<N>::cellValueSet(const cellValueSet& A) : 
  keyName(A.keyName),outName(A.outName),whatValue(A.whatValue)
  /*!
    Copy constructor
    \param A :: cellValueSet to copy
  */
{}

template<size_t N>  
cellValueSet<N>&
cellValueSet<N>::operator=(const cellValueSet<N>& A)
  /*!
    Assignment operator
    \param A :: cellValueSet to copy
    \return *this
  */
{
  if (this!=&A)
    {
      whatValue=A.whatValue;
      dataMap=A.dataMap;
    }
  return *this;
}

template<size_t N>
cellValueSet<N>::~cellValueSet()
  /*!
    Destructor
  */
{}

template<size_t N>  
void
cellValueSet<N>::clearAll()
  /*!
    The big reset
  */
{
  dataMap.erase(dataMap.begin(),dataMap.end());  
  return;
}

template<size_t N>  
bool
cellValueSet<N>::cellSplit(const std::vector<int>& cellN,
			   std::vector<std::tuple<int,int>>& initCell,
			   std::vector<std::array<double,N>>& outData) const
  /*!
    Process ranges to find for value
    \param cellN :: Cell vlaues
    \param initCell :: initialization range
    \param dataValue :: initialization range
    \return true if output required
   */
{
  typedef std::tuple<int,int> TITEM;
  initCell.clear();
  outData.clear();
  
  if (dataMap.empty() || cellN.empty()) return 0;
  
  size_t prev(0);
  valTYPE V;
  for(size_t i=0;i<cellN.size();i++)
    {
      const int CN=cellN[i];
      typename dataTYPE::const_iterator mc=dataMap.find(CN);
      if (mc==dataMap.end())
	{
	  if (prev)
	    {
	      initCell.push_back(TITEM(cellN[prev-1],cellN[i-1]));
	      V=mc->second;
	      prev=0;
	    }
	}
      else
	{
	  if (prev)
	    {
	      for(size_t index=0;index<N;index++)
		{
		  if (std::abs(V[index]-mc->second[index])>Geometry::zeroTol)
		    {
		      initCell.push_back(TITEM(cellN[prev-1],cellN[i-1]));
		      outData.push_back(V);
		      prev=i+1;
		      V=mc->second;
		    }
		  break;
		}
	    }
	  else if (!prev)
	    {
	      prev=i+1;
	      V=mc->second;
	    }
	}
    }
  
  if (prev)
    {
      initCell.push_back(TITEM(cellN[prev-1],cellN.back()));
      outData.push_back(V);
    }
  
  return (initCell.empty()) ? 0 : 1;
}



template<size_t N>
void
cellValueSet<N>::setValues(const int cN,const double V)
  /*!
    Set a value in the map
    \param cN :: Cell number   
    \param V :: value for cell
  */
{
  valTYPE A;
  A[0]=V;
  dataMap.emplace(cN,A);
  return;
}

template<size_t N>
void
cellValueSet<N>::setValues(const int cN,const double V,
			   const double V2)
  /*!
    Set a value in the map
    \param cN :: Cell number   
    \param V :: value for cell
    \param V2 :: value for cell
  */
{
  valTYPE A;
  A[0]=V;
  A[1]=V2;
  dataMap.emplace(cN,A);
  return;
}

template<size_t N>
void
cellValueSet<N>::writeFLUKA(std::ostream& OX,
			 const std::vector<int>& cellN,
			 const std::string& ControlStr) const 
/*!
    Process is to write keyName ControlStr units 
    \param OX :: Output stream
    \param cellN :: vector of cells
    \param ControlStr units [%0/%1/%2] for cell range/Value
  */
{
  ELog::RegMethod RegA("cellValueSet","writeFLUKA");
  
  typedef std::tuple<int,int> TITEM;

  std::ostringstream cx;
  std::vector<TITEM> Bgroup;
  std::vector<valTYPE> Bdata;

  if (cellSplit(cellN,Bgroup,Bdata))
    {

      const std::vector<std::string> Units=StrFunc::StrParts(ControlStr);
      std::vector<std::string> SArray(4);
      SArray[0]=std::to_string(whatValue);
      for(size_t index=0;index<Bgroup.size();index++)
	{
	  const TITEM& tc(Bgroup[index]);
	  const valTYPE& dArray(Bdata[index]);

	  SArray[1]="R"+std::to_string(std::get<0>(tc));
	  SArray[2]="R"+std::to_string(std::get<1>(tc));
	  for(size_t i=0;i<N;i++)
	    SArray[3+i]=std::to_string(dArray[i]);
	  cx.str("");
	  cx<<outName<<" ";
	  for(const std::string& UC : Units)
	    {
	      if (UC[0]=='%' && UC.size()==2)
		{
		  const size_t SA=(static_cast<size_t>(UC[1]-'0') % N+3); 
		  cx<<SArray[SA]<<" ";
		}
	      else
		cx<<UC<<" ";
	    }
	  StrFunc::writeFLUKA(cx.str(),OX);
	}
    }
  return;
}

///\cond TEMPLATE
template class cellValueSet<1>;
template class cellValueSet<2>;
///\endcond TEMPLATE
  
} // NAMESPACE flukaSystem
      
   
