/********************************************************************* 
  CombLayer : MCNP(X) Input builder
 
 * File:   scatMatInc/xrayMaterial.h
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
#ifndef scatterSystem_xrayMaterial_h
#define scatterSystem_xrayMaterial_h

namespace MonteCarlo
{
  class Material;
  class MXcards;
  class Zaid;
  class neutron;
}

namespace scatterSystem
{
  /*!
    \class xrayMaterial
    \brief Xray information on the material
    \author S. Ansell
    \version 1.0
    \date November 2019
    
    This processes X-ray data for both f0 and f' + f''
  */
  
class xrayMaterial : public MonteCarlo::Material
{
 protected:

  // This section is for f0: of type Sum A_i exp(-b_i K^2) + Sum c_j
  // It refers to multiple elements [all stacked toegther]
  // k = sin(theta/2)/ lambda   (reduced Q)
  std::vector<double> AFactor;      ///< A_i [5 per element]
  std::vector<double> BFactor;      ///< B_i [5 per element]
  std::vector<double> CFactor;      ///< C_j [1 per element]

  
 public:
  
  xrayMaterial();
  xrayMaterial(const std::string&,const double,const double,
	   const double,const double,const double,const double);
  xrayMaterial(const double,const double,const double,
	   const double,const double,const double);
  xrayMaterial(const xrayMaterial&);
  virtual xrayMaterial* clone() const;
  xrayMaterial& operator=(const xrayMaterial&);
  virtual ~xrayMaterial();
  
  void setScat(const double,const double,const double);

  // get Scattering prob
  virtual double scatXSection(const MonteCarlo::particle&) const;
  virtual double totalXSection(const MonteCarlo::particle&) const;

  virtual double scatTotalRatio(const MonteCarlo::particle&) const;
  virtual double ElasticTotalRatio(const double) const;


  virtual double calcRefIndex(const double) const;
  virtual double calcAtten(const MonteCarlo::particle&,const double) const;
  
  virtual void scatterNeutron(MonteCarlo::particle&) const;
  
  virtual double dSdOdE(const MonteCarlo::neutron&,
			const MonteCarlo::neutron&) const;

  virtual void write(std::ostream&) const;

};


} // NAMESPACE MonteCarlo

#endif
