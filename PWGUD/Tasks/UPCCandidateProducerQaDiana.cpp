// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "CommonConstants/LHCConstants.h"
#include "PWGUD/Core/UPCCutparHolder.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandProducerQaDiana {
  HistogramRegistry histRegistry{"HistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const AxisSpec axisPt{500, 0., 5., ""};
    const AxisSpec axisEta{100, -4., -2.5, ""}; // this is for muon tracks
    const AxisSpec axisPhi{628, 0., 6.28, ""};

    histRegistry.add("TracksQA/Muon/Pt/MCH_MID", "", kTH1F, {axisPt});
    histRegistry.add("TracksQA/Muon/Pt_Gen/MCH_MID", "", kTH1F, {axisPt});

    histRegistry.add("TracksQA/Muon/Eta/MCH_MID", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Muon/Phi/MCH_MID", "", kTH1F, {axisPhi});

    histRegistry.add("TracksQA/MC/PtMu", "", kTH1F, {axisPt});
  }

  using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra>;
  using FwdTracks = o2::aod::FwdTracks;

  void updateFwdTrackQA(const FwdTracks::iterator& fwdtrack)
  {
    if (fwdtrack.eta() < -4. || fwdtrack.eta() > -2.5)
      return;
    histRegistry.fill(HIST("TracksQA/Muon/Pt/MCH_MID"), fwdtrack.pt());
    histRegistry.fill(HIST("TracksQA/Muon/Eta/MCH_MID"), fwdtrack.eta());
    histRegistry.fill(HIST("TracksQA/Muon/Phi/MCH_MID"), fwdtrack.phi());
  }

  void process(FwdTracks const& fwdtracks)
  {

    for (const auto& fwdtrack : fwdtracks) {
      if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack)
        continue;
      updateFwdTrackQA(fwdtrack);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcCandProducerQaDiana>(cfgc)};
}
