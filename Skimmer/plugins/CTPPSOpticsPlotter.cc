/****************************************************************************
 * Authors:
 *   Jan Kašpar
 ****************************************************************************/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"

#include "CondFormats/DataRecord/interface/CTPPSInterpolatedOpticsRcd.h"
#include "CondFormats/CTPPSReadoutObjects/interface/LHCInterpolatedOpticalFunctionsSetCollection.h"

#include "TFile.h"
#include "TGraph.h"

#include <map>
#include <string>

//----------------------------------------------------------------------------------------------------

class CTPPSOpticsPlotter : public edm::one::EDAnalyzer<>
{
  public:
    explicit CTPPSOpticsPlotter(const edm::ParameterSet&);

  private:
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;

    std::string outputFile_;

    struct RPPlots
    {
      std::unique_ptr<TGraph> h_y_vs_x_disp;
      RPPlots() :
        h_y_vs_x_disp(new TGraph)
      {}

      void write() const {
        h_y_vs_x_disp->SetTitle(";x   (mm);y   (mm)");
        h_y_vs_x_disp->Write("h_y_vs_x_disp");
      }
    };

    std::map<unsigned int, RPPlots> rp_plots_;
};

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

CTPPSOpticsPlotter::CTPPSOpticsPlotter(const edm::ParameterSet& iConfig) :
  outputFile_(iConfig.getParameter<string>("outputFile"))
{}

//----------------------------------------------------------------------------------------------------

void CTPPSOpticsPlotter::analyze(const edm::Event& iEvent, const edm::EventSetup &iSetup)
{
  // stop if plots already made
  if (!rp_plots_.empty())
    return;

  // get conditions
  edm::ESHandle<LHCInterpolatedOpticalFunctionsSetCollection> hOpticalFunctions;
  iSetup.get<CTPPSInterpolatedOpticsRcd>().get(hOpticalFunctions);

  // stop if conditions invalid
  if (hOpticalFunctions->empty())
    return;

  // make plots
  for (const auto &it : *hOpticalFunctions)
  {
    CTPPSDetId rpId(it.first);
    unsigned int rpDecId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();

    auto &pl = rp_plots_[rpDecId];

    LHCInterpolatedOpticalFunctionsSet::Kinematics k_in_beam = { 0., 0., 0., 0., 0. };
    LHCInterpolatedOpticalFunctionsSet::Kinematics k_out_beam;
    it.second.transport(k_in_beam, k_out_beam);

    for (double xi = 0.; xi < 0.30001; xi += 0.001)
    {
      LHCInterpolatedOpticalFunctionsSet::Kinematics k_in = { 0., 0., 0., 0., xi };  // conversions: CMS --> LHC convention
      LHCInterpolatedOpticalFunctionsSet::Kinematics k_out;
      it.second.transport(k_in, k_out);

      const double x = (k_out.x - k_out_beam.x) * 10.;  // conversions: cm --> mm
      const double y = (k_out.y - k_out_beam.y) * 10.;

      int idx = pl.h_y_vs_x_disp->GetN();
      pl.h_y_vs_x_disp->SetPoint(idx, x, y);
    }
  }
}

//----------------------------------------------------------------------------------------------------

void CTPPSOpticsPlotter::endJob()
{
  auto f_out = std::make_unique<TFile>(outputFile_.c_str(), "recreate");

  for (const auto& p : rp_plots_) {
    gDirectory = f_out->mkdir(Form("%u", p.first));
    p.second.write();
  }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(CTPPSOpticsPlotter);
