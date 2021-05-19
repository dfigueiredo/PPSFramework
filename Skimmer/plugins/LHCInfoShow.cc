/****************************************************************************
 * Authors:
 *   D. Figueiredo
 ****************************************************************************/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

#include <map>
#include <string>

//----------------------------------------------------------------------------------------------------

class LHCInfoShow : public edm::one::EDAnalyzer<>
{
  public:
    explicit LHCInfoShow(const edm::ParameterSet&);

  private:
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;

};

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

LHCInfoShow::LHCInfoShow(const edm::ParameterSet& iConfig){}

//----------------------------------------------------------------------------------------------------

void LHCInfoShow::analyze(const edm::Event& iEvent, const edm::EventSetup &iSetup)
{


  edm::ESHandle<LHCInfo> hLHCInfo;
  iSetup.get<LHCInfoRcd>().get("", hLHCInfo);

  std::cout << "Colliding bunches: " << hLHCInfo->collidingBunches() << std::endl;
  std::cout << "# bunches Beam1: " << hLHCInfo->bunchesInBeam1() << std::endl;
  std::cout << "# bunches Beam2: " << hLHCInfo->bunchesInBeam2() << std::endl;
  std::cout << "Intensity beam1: " << hLHCInfo->intensityForBeam1() << std::endl;
  std::cout << "Intensity beam2: " << hLHCInfo->intensityForBeam2() << std::endl;
  std::cout << "Inst. Luminosity: " << hLHCInfo->instLumi() << std::endl;
  std::cout << "XAngle (LHC) [murad]: " << hLHCInfo->crossingAngle() << std::endl;
  std::cout << "Beta Star: " << hLHCInfo->betaStar() << std::endl;


}

//----------------------------------------------------------------------------------------------------

void LHCInfoShow::endJob()
{
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(LHCInfoShow);
