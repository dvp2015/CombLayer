/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   scatMat/xrayMaterial.cxx
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
#include <list>
#include <vector>
#include <map>
#include <stack>
#include <string>
#include <algorithm>
#include <boost/multi_array.hpp>

#include "MersenneTwister.h"
#include "Exception.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "FileReport.h"
#include "GTKreport.h"
#include "NameStack.h"
#include "RegMethod.h"
#include "OutputLog.h"
#include "support.h"
#include "RefCon.h"
#include "MatrixBase.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "particle.h"
#include "neutron.h"
#include "Zaid.h"
#include "MXcards.h"
#include "Material.h"
#include "xrayMaterial.h"

extern MTRand RNG;

namespace scatterSystem
{

xrayMaterial::xrayMaterial() : MonteCarlo::Material()
  /*!
    Constructor
  */
{}

xrayMaterial::xrayMaterial(const std::string& N,const double M,
			   const double D) :
  Material(N,M,D)
  /*!
    Constructor for values
    \param N :: xrayMaterial name
    \param M :: Mean atomic mass
    \param D :: atomDensity [atom/angtrom^3]
  */
{}


xrayMaterial*
xrayMaterial::clone() const
  /*!
    Clone method
    \return new (this)
   */
{
  return new xrayMaterial(*this);
}

xrayMaterial::~xrayMaterial() 
  /*!
    Destructor
   */
{}

void
xrayMaterial::setScat(const double S,const double I,const double A)
  /*!
    Set the scattering factors
    \param S :: Coherrent-Scattering cross section
    \param I :: Inc-Scattering cross section
    \param A :: Absorption cross section
  */
{
  scoh=S;
  sinc=I;
  sabs=A;
  bTotal=sqrt(S+I)/(4*M_PI);
  return;
}
  
double
xrayMaterial::totalXSection(const MonteCarlo::particle& n0) const
  /*!
    Given Wavelength get the attenuation coefficient.
    \param n0 :: neutron for Wavelength [Angstrom]
    \return Attenuation (including atomDensity)
  */
{
  return atomDensity*(scoh+sinc+n0.wavelength*sabs/1.798);
}

double
xrayMaterial::scatXSection(const MonteCarlo::particle&) const
  /*!
    Given wavelength get the scattering cross section
    \param  :: particle 
    \return Attenuation (including atomDensity)
  */
{
  return atomDensity*(scoh+sinc);
}

double
xrayMaterial::calcAtten(const MonteCarlo::particle& N,
			const double Length) const
  /*!
    Calculate the attenuation factor coefficient.
    \param Wave :: Wavelength [Angstrom]
    \param Length :: Absorption length
    \return Attenuation (including atomDensity)
  */
{
  return exp(-Length*atomDensity*(scoh+sinc+N.wavelength*sabs/1.798));
}

void
xrayMaterial::scatterNeutron(MonteCarlo::particle& N) const
  /*!
    Calculate the new angle and energy of the neutron
    that is scattered. Full isotropic scattering
    Doesn't change energy:
    \param N :: neutron to scatter
  */
{
  ELog::RegMethod RegA("xrayMaterial","scatterNeutron");

  const double theta=2*M_PI*RNG.rand();
  const double phi=M_PI*RNG.rand();
  N.uVec[0]=cos(theta)*sin(phi);
  N.uVec[1]=sin(theta)*sin(phi);
  N.uVec[2]=cos(phi);

  N.addCollision();  // This also updates the ID !!!
  return;
}

double
xrayMaterial::scatTotalRatio(const MonteCarlo::particle& N) const
  /*!
    Return the ratio of scattering / total crosssection
    [for scaling by absorption coefficient]
    \param Wave :: Wavelength [Angstrom]
    \return sigma_scatter/sigma_total
  */
{
  return (atomDensity>0) ? 
    (scoh+sinc)/(scoh+sinc+N.wavelength*sabs/1.798) : 1.0;
}

double
xrayMaterial::ElasticTotalRatio(const double) const
  /*!
    Calculate the ratio of elastic/sigma_T
    against total.
    \param :: Wavelength [Angstrom]
    \return sigma_elastic/sigma_scatter
  */
{
  return (atomDensity>0) ? scoh/(scoh+sinc) : 1.0;
}


double
xrayMaterial::calcRefIndex(const double Wave) const
  /*!
    Calc refractive index 
    \param Wave :: Wavelength
    \return n - the real refractive index
  */
{
//  ELog::RegMethod RegA("xrayMaterial","calcRefIndex");
//  ELog::EM<<"Ref ["<<Name<<"]=="
//	  <<1.0-sqrt(1.0-Wave*Wave*(1e-5*atomDensity*bcoh/M_PI))<<ELog::endErr;
  return sqrt(1.0-Wave*Wave*(1e-5*atomDensity*bcoh/M_PI));
}

double
xrayMaterial::dSdOdE(const MonteCarlo::neutron& NIn,
		     const MonteCarlo::neutron& NOut) const
  /*!
    Determines the weight, relative to the total cross section
    of a neutron exiting with the channel NOut.
    Uses Ideal gas equation:
    \param NIn :: Input neutorn
    \param NOut :: Output neutorn
    \return Cross-section from A-B [Barns]
  */
{
  ELog::RegMethod RegA("xrayMaterial","dSdOdE");
  
  // First calculate dSDoDE for the given energy difference.
  // InElastic scatter:
  const double Ei=NIn.energy;
  const double Ef=NOut.energy;

  const double mu=NOut.uVec.dotProd(NIn.uVec);
  
  // Maybe a factor 2 missing?????

  const double alpha=(Ei+Ef-2.0*mu*sqrt(Ei*Ef))/
    (4*Amass*RefCon::k_bev*realTemp); 
  const double beta=(Ei-Ef)/(RefCon::k_bev*realTemp);

  // Non-Elastic scatter [neutron loss]: 
  
  const double SQW=exp(-(beta-alpha)*(beta-alpha)/(4.0*alpha)) / 
    (2*RefCon::k_bev*realTemp*sqrt(M_PI*alpha));
  return SQW*bTotal*bTotal*NIn.wavelength/NOut.wavelength;
}

void
xrayMaterial::write(std::ostream& OX) const 
  /*!
    Write out the material information for MCNPX
    \param OX :: Output stream
  */
{
  OX<<"c NeutMat:"<<std::endl;
    
  Material::write(OX);
  return;
}

} // NAMESPACE MonteCarlo

