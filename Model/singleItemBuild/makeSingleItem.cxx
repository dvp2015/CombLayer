/*********************************************************************
  CombLayer : MCNP(X) Input builder

 * File:   singleItemBuild/makeSingleItem.cxx
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
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
#include "OutputLog.h"
#include "BaseVisit.h"
#include "BaseModVisit.h"
#include "Vec3D.h"
#include "Line.h"
#include "inputParam.h"
#include "surfRegister.h"
#include "objectRegister.h"
#include "HeadRule.h"
#include "LinkUnit.h"
#include "FixedComp.h"
#include "FixedOffset.h"
#include "FixedRotate.h"
#include "ContainedComp.h"
#include "ContainedGroup.h"
#include "ExternalCut.h"
#include "FrontBackCut.h"
#include "BaseMap.h"
#include "CellMap.h"
#include "SurfMap.h"
#include "InnerZone.h"
#include "World.h"
#include "insertObject.h"
#include "insertSphere.h"
#include "insertCylinder.h"
#include "insertShell.h"
#include "VacuumPipe.h"
#include "Quadrupole.h"
#include "Sexupole.h"
#include "Octupole.h"
#include "LQuadF.h"
#include "LQuadH.h"
#include "LSexupole.h"
#include "CorrectorMag.h"
#include "EPSeparator.h"
#include "QuadUnit.h"
#include "DipoleChamber.h"
#include "R3ChokeChamber.h"
#include "PreDipole.h"
#include "MagnetM1.h"
#include "MagnetBlock.h"
#include "CylGateValve.h"
#include "BPM.h"
#include "BeamDivider.h"
#include "CeramicSep.h"
#include "DipoleDIBMag.h"
#include "EArrivalMon.h"
#include "EBeamStop.h"
#include "SixPortTube.h"
#include "Scrapper.h"
#include "FlatPipe.h"
#include "TriPipe.h"
#include "MultiPipe.h"
#include "YagScreen.h"
#include "YagUnit.h"
#include "YagUnitBig.h"
#include "TWCavity.h"
#include "SplitFlangePipe.h"
#include "Bellows.h"
#include "VirtualTube.h"
#include "PipeTube.h"
#include "BlankTube.h"
#include "ButtonBPM.h"


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

  std::set<std::string> validItems
    ({
      "default","CylGateValve","CorrectorMag","LQuadF","LQuadH","LSexupole",
      "MagnetBlock","Sexupole","MagnetM1","Octupole","CeramicSep",
      "EBeamStop","EPSeparator","R3ChokeChamber","QuadUnit",
      "DipoleChamber","EPSeparator","Quadrupole","TargetShield",
      "FlatPipe","TriPipe","SixPort",
      "DipoleDIBMag","EArrivalMon","YagScreen","YAG",
      "YagUnit","YagUnitBig","BPM","BeamDivider","Scrapper","TWCavity",
      "Bellow", "VacuumPipe","MultiPipe","PipeTube","BlankTube","ButtonBPM",
      "Help","help"
    });

  ModelSupport::objectRegister& OR=
    ModelSupport::objectRegister::Instance();
  const int voidCell(74123);

  const std::string item=
    IParam.getDefValue<std::string>("default","singleItem");

  if (validItems.find(item)==validItems.end())
    throw ColErr::InContainerError<std::string>
      (item,"Item no a single component");

  ELog::EM<<"Conponent == "<<item<<ELog::endDiag;
  if (item=="default" || item == "CylGateValve" )
    {
      std::shared_ptr<xraySystem::CylGateValve>
	GV(new xraySystem::CylGateValve("GV"));

      OR.addObject(GV);

      GV->addInsertCell(voidCell);
      GV->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item == "YAG" || item=="YagScreen")
    {
      std::shared_ptr<tdcSystem::YagScreen>
	YAG(new tdcSystem::YagScreen("YAG"));
      OR.addObject(YAG);

      YAG->addAllInsertCell(voidCell);
      YAG->setBeamAxis(Geometry::Vec3D(0,-10,0),
		       Geometry::Vec3D(1,0,0));
      YAG->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item == "YagUnit")
    {
      std::shared_ptr<tdcSystem::YagUnit>
	yagUnit(new tdcSystem::YagUnit("YU"));

      std::shared_ptr<tdcSystem::YagScreen>
	yagScreen(new tdcSystem::YagScreen("YAG"));

      OR.addObject(yagUnit);
      OR.addObject(yagScreen);

      yagUnit->addInsertCell(voidCell);
      yagUnit->createAll(System,World::masterOrigin(),0);

      yagScreen->setBeamAxis(*yagUnit,1);
      yagScreen->createAll(System,*yagUnit,4);
      yagScreen->insertInCell("Outer",System,voidCell);
      yagScreen->insertInCell("Connect",System,yagUnit->getCell("PlateB"));
      yagScreen->insertInCell("Connect",System,yagUnit->getCell("Void"));
      yagScreen->insertInCell("Payload",System,yagUnit->getCell("Void"));

      return;
    }
  if (item == "YagUnitBig")
    {
      std::shared_ptr<tdcSystem::YagUnitBig>
	yagUnit(new tdcSystem::YagUnitBig("YUBig"));

      std::shared_ptr<tdcSystem::YagScreen>
	yagScreen(new tdcSystem::YagScreen("YAG"));

      OR.addObject(yagUnit);
      OR.addObject(yagScreen);

      yagUnit->addInsertCell(voidCell);
      yagUnit->createAll(System,World::masterOrigin(),0);

      yagScreen->setBeamAxis(*yagUnit,1);
      yagScreen->createAll(System,*yagUnit,4);
      yagScreen->insertInCell("Outer",System,voidCell);
      yagScreen->insertInCell("Connect",System,yagUnit->getCell("Plate"));
      yagScreen->insertInCell("Connect",System,yagUnit->getCell("Void"));
      yagScreen->insertInCell("Payload",System,yagUnit->getCell("Void"));


      return;
    }
  if (item == "BPM")
    {
      std::shared_ptr<tdcSystem::BPM>
	bpm(new tdcSystem::BPM("BPM"));
      OR.addObject(bpm);

      bpm->addInsertCell(voidCell);
      bpm->createAll(System,World::masterOrigin(),0);

      return;
    }
  if (item == "EBeamStop")
    {
      std::shared_ptr<tdcSystem::EBeamStop>
	eBeam(new tdcSystem::EBeamStop("EBeam"));
      OR.addObject(eBeam);

      eBeam->addAllInsertCell(voidCell);
      eBeam->createAll(System,World::masterOrigin(),0);

      return;
    }
  if (item == "CeramicSep")
    {
      std::shared_ptr<tdcSystem::CeramicSep>
	cSep(new tdcSystem::CeramicSep("CerSep"));
      OR.addObject(cSep);

      cSep->addInsertCell(voidCell);
      cSep->createAll(System,World::masterOrigin(),0);

      return;
    }
  if (item == "BeamDivider")
    {
      std::shared_ptr<tdcSystem::BeamDivider>
	bd(new tdcSystem::BeamDivider("BeamDiv"));
      OR.addObject(bd);

      bd->addAllInsertCell(voidCell);
      bd->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item == "CorrectorMag" )
    {
      std::shared_ptr<constructSystem::VacuumPipe>
	VC(new constructSystem::VacuumPipe("CorrectorMagPipe"));
      std::shared_ptr<tdcSystem::CorrectorMag>
	CM(new tdcSystem::CorrectorMag("CM"));

      OR.addObject(VC);
      OR.addObject(CM);

      VC->addInsertCell(voidCell);
      VC->createAll(System,World::masterOrigin(),0);

      CM->setCutSurf("Inner",*VC,"outerPipe");
      CM->addInsertCell(voidCell);
      CM->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item == "EArrivalMon" )
    {
      std::shared_ptr<tdcSystem::EArrivalMon>
	EA(new tdcSystem::EArrivalMon("BeamMon"));

      OR.addObject(EA);

      EA->addInsertCell(voidCell);
      EA->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="LQuadF")
    {
      std::shared_ptr<constructSystem::VacuumPipe>
	VC(new constructSystem::VacuumPipe("QHVC"));
      std::shared_ptr<tdcSystem::LQuadF>
	QF(new tdcSystem::LQuadF("QF","QF"));

      OR.addObject(VC);
      OR.addObject(QF);

      VC->addInsertCell(voidCell);
      VC->createAll(System,World::masterOrigin(),0);

      QF->setCutSurf("Inner",*VC,"outerPipe");
      QF->addInsertCell(voidCell);
      QF->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="LQuadH")
    {
      std::shared_ptr<constructSystem::VacuumPipe>
	VC(new constructSystem::VacuumPipe("QHVC"));
      std::shared_ptr<tdcSystem::LQuadH>
	QH(new tdcSystem::LQuadH("QH","QH"));

      OR.addObject(VC);
      OR.addObject(QH);

      VC->addInsertCell(voidCell);
      VC->createAll(System,World::masterOrigin(),0);

      QH->setCutSurf("Inner",*VC,"outerPipe");
      QH->addInsertCell(voidCell);
      QH->createAll(System,World::masterOrigin(),0);


      return;
    }

  if (item=="LSexupole")
    {
      std::shared_ptr<constructSystem::VacuumPipe>
	VC(new constructSystem::VacuumPipe("QHVC"));
      std::shared_ptr<tdcSystem::LSexupole>
	LS(new tdcSystem::LSexupole("LS","LS"));

      OR.addObject(VC);
      OR.addObject(LS);

      VC->addInsertCell(voidCell);
      VC->createAll(System,World::masterOrigin(),0);

      LS->setCutSurf("Inner",*VC,"outerPipe");
      LS->addInsertCell(voidCell);
      LS->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="MagnetBlock")
    {

      std::shared_ptr<xraySystem::MagnetBlock>
	MB(new xraySystem::MagnetBlock("M1"));

      OR.addObject(MB);
      MB->addInsertCell(voidCell);
      MB->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="Sexupole")
    {
      std::shared_ptr<xraySystem::Sexupole>
	SXX(new xraySystem::Sexupole("SXX","SXX"));

      OR.addObject(SXX);

      SXX->addInsertCell(voidCell);
      SXX->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="MagnetM1")
    {
      std::shared_ptr<xraySystem::MagnetM1>
	MagBlock(new xraySystem::MagnetM1("M1Block"));

      OR.addObject(MagBlock);

      MagBlock->addAllInsertCell(voidCell);
      MagBlock->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="Octupole")
    {

      std::shared_ptr<xraySystem::Octupole>
	OXX(new xraySystem::Octupole("M1BlockOXX","M1BlockOXX"));
      OR.addObject(OXX);
      OXX->addInsertCell(voidCell);
      OXX->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="EPSeparator")
    {
      std::shared_ptr<xraySystem::EPSeparator>
	EPsep(new xraySystem::EPSeparator("EPSeparator"));
      OR.addObject(EPsep);

      EPsep->addInsertCell(voidCell);
      EPsep->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="R3ChokeChamber")
    {
      std::shared_ptr<xraySystem::R3ChokeChamber>
	CChamber(new xraySystem::R3ChokeChamber("R3Chamber"));
      OR.addObject(CChamber);
      CChamber->addAllInsertCell(voidCell);
      CChamber->createAll(System,World::masterOrigin(),0);

      return;
    }
  if (item=="")
    {
      std::shared_ptr<xraySystem::R3ChokeChamber>
	CChamber(new xraySystem::R3ChokeChamber("R3Chamber"));
      OR.addObject(CChamber);
      CChamber->addAllInsertCell(voidCell);
      CChamber->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="QuadUnit" || item=="DipoleChamber" || item=="EPSeparator")
    {
      std::shared_ptr<xraySystem::QuadUnit>
	PDipole(new xraySystem::QuadUnit("PreDipole"));
      OR.addObject(PDipole);
      PDipole->addInsertCell(voidCell);
      PDipole->createAll(System,World::masterOrigin(),0);


      std::shared_ptr<xraySystem::DipoleChamber>
	DCSep(new xraySystem::DipoleChamber("DipoleChamber"));
      OR.addObject(DCSep);
      DCSep->addAllInsertCell(voidCell);
      DCSep->createAll(System,*PDipole,2);

      std::shared_ptr<xraySystem::EPSeparator>
	EPSep(new xraySystem::EPSeparator("EPSep"));
      OR.addObject(EPSep);
      EPSep->addInsertCell(voidCell);
      EPSep->createAll(System,*PDipole,2);
      return;
    }

  if (item=="Quadrupole")
    {

      std::shared_ptr<xraySystem::Quadrupole>
	Quad(new xraySystem::Quadrupole("Quad","Quad"));
      OR.addObject(Quad);
      Quad->addInsertCell(voidCell);
      Quad->createAll(System,World::masterOrigin(),0);
      return;
    }
  if (item == "Scrapper")
    {
      std::shared_ptr<tdcSystem::Scrapper>
	sc(new tdcSystem::Scrapper("Scrapper"));
      OR.addObject(sc);

      sc->addInsertCell(voidCell);
      sc->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="TargetShield")
    {
      std::shared_ptr<insertSystem::insertSphere>
	Target(new insertSystem::insertSphere("Target"));
      std::shared_ptr<insertSystem::insertShell>
	Surround(new insertSystem::insertShell("Shield"));
      std::shared_ptr<insertSystem::insertCylinder>
	TubeA(new insertSystem::insertCylinder("TubeA"));
      std::shared_ptr<insertSystem::insertCylinder>
	TubeB(new insertSystem::insertCylinder("TubeB"));

      OR.addObject(Target);
      OR.addObject(TubeA);
      OR.addObject(TubeB);
      OR.addObject(Surround);

      TubeA->addInsertCell(voidCell);
      TubeA->createAll(System,World::masterOrigin(),0);
      TubeB->addInsertCell(voidCell);
      TubeB->createAll(System,*TubeA,2);
      return;
    }

  if (item=="SixPort")
    {
      std::shared_ptr<tdcSystem::SixPortTube>
	SP(new tdcSystem::SixPortTube("SixPort"));

      OR.addObject(SP);

      SP->addInsertCell(voidCell);
      SP->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item == "TriPipe")
    {
      std::shared_ptr<tdcSystem::TriPipe>
	tp(new tdcSystem::TriPipe("TriPipe"));
      OR.addObject(tp);

      tp->addAllInsertCell(voidCell);
      tp->createAll(System,World::masterOrigin(),0);

      return;
    }
  if (item == "MultiPipe")
    {
      std::shared_ptr<tdcSystem::MultiPipe>
	tp(new tdcSystem::MultiPipe("MultiPipe"));
      OR.addObject(tp);

      tp->addAllInsertCell(voidCell);
      tp->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item == "FlatPipe")
    {
      std::shared_ptr<tdcSystem::FlatPipe>
	tp(new tdcSystem::FlatPipe("FlatPipe"));
      OR.addObject(tp);

      tp->addAllInsertCell(voidCell);
      tp->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="DipoleDIBMag")
    {
      std::shared_ptr<constructSystem::VacuumPipe>
	VC(new constructSystem::VacuumPipe("VC"));
      std::shared_ptr<tdcSystem::DipoleDIBMag>
	DIB(new tdcSystem::DipoleDIBMag("DIB"));

      OR.addObject(VC);
      OR.addObject(DIB);

      VC->addInsertCell(voidCell);
      VC->createAll(System,World::masterOrigin(),0);

      DIB->setCutSurf("Inner",*VC,"outerPipe");
      DIB->addInsertCell(voidCell);
      DIB->createAll(System,World::masterOrigin(),0);

      return;
    }

  if (item=="TWCavity") // traveling wave cavity
    {
      std::shared_ptr<constructSystem::VacuumPipe>
	pipeA(new constructSystem::VacuumPipe("PipeA"));
      std::shared_ptr<tdcSystem::TWCavity> cavity(new tdcSystem::TWCavity("TWCavity"));
      std::shared_ptr<constructSystem::VacuumPipe>
	pipeB(new constructSystem::VacuumPipe("PipeB"));

      OR.addObject(pipeA);
      OR.addObject(cavity);
      OR.addObject(pipeB);

      pipeA->addInsertCell(voidCell);
      pipeA->createAll(System,World::masterOrigin(),0);

      cavity->addInsertCell(voidCell);
      cavity->createAll(System,*pipeA,"back");

      pipeB->addInsertCell(voidCell);
      pipeB->createAll(System,*cavity,"back");

      return;
    }

  if (item=="Bellow")
    {
      std::shared_ptr<constructSystem::Bellows>
	bellow(new constructSystem::Bellows("Bellow"));
      OR.addObject(bellow);

      bellow->addInsertCell(voidCell);
      bellow->createAll(System,World::masterOrigin(),0);

      return;
    }

    if (item == "VacuumPipe" )
    {
      std::shared_ptr<constructSystem::VacuumPipe>
	VC(new constructSystem::VacuumPipe("VC"));

      OR.addObject(VC);

      VC->addInsertCell(voidCell);
      VC->createAll(System,World::masterOrigin(),0);

      return;
    }

    if (item == "PipeTube" )
    {
      std::shared_ptr<constructSystem::PipeTube>
	pipeTube(new constructSystem::PipeTube("PipeTube"));

      OR.addObject(pipeTube);

      pipeTube->addAllInsertCell(voidCell);
      pipeTube->createAll(System,World::masterOrigin(),0);

      return;
    }
    if (item == "BlankTube" )
    {
      std::shared_ptr<constructSystem::BlankTube>
	blankTube(new constructSystem::BlankTube("BlankTube"));

      OR.addObject(blankTube);

      blankTube->addAllInsertCell(voidCell);
      blankTube->createAll(System,World::masterOrigin(),0);

      return;
    }
    if (item == "ButtonBPM" )
    {
      std::shared_ptr<tdcSystem::ButtonBPM>
	buttonBPM(new tdcSystem::ButtonBPM("ButtonBPM"));

      OR.addObject(buttonBPM);

      buttonBPM->addInsertCell(voidCell);
      buttonBPM->createAll(System,World::masterOrigin(),0);

      return;
    }


  if (item=="Help" || item=="help")

  if (item=="Help" || item=="help")
    {

      ELog::EM<<"Valid items for single selection:\n"<<ELog::endDiag;

      for(const std::string& Name : validItems)
	ELog::EM<<"Item : "<<Name<<"\n";

      ELog::EM<<"-----------"<<ELog::endDiag;
    }

  return;
}

}   // NAMESPACE singleItemSystem
