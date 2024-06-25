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
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"

#include "TLorentzVector.h"
#include "TSystem.h"


namespace mpt {
  DECLARE_SOA_COLUMN(M, m, float);
  DECLARE_SOA_COLUMN(Pt, pt, float);
  //DECLARE_SOA_COLUMN(Y, y, float);
}

namespace o2::aod {
  DECLARE_SOA_TABLE(MassPt, "AOD", "MASSPT", o2::soa::Index<>,
		    mpt::M, mpt::Pt); // , mpt::Y
}

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct MuMuForward {

  Produces<o2::aod::MassPt> mptSel;
  
  HistogramRegistry hr{"hr"};

  double mMu = 0.10566;

  using CandidatesFwd = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSelsFwd>;
  using ForwardTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;

  void init(InitContext&)
  {
    const AxisSpec axisM(1000, 0., 10.);
    const AxisSpec axisPt(1000, 0., 10.);
    //const AxisSpec axisY(1000, -5., 5.);
    hr.add<TH1>("hM", "", kTH1D, {axisM});
    hr.add<TH1>("hPt", "", kTH1D, {axisPt});
    //hr.add<TH1>("hY", "", kTH1D, {axisY});
  }

  bool isMuonSelected(const ForwardTracks::iterator& fwdTrack)
  {
    float rAbsMin = 17.6;
    float rAbsMid = 26.5;
    float rAbsMax = 89.5;
    float pDca1 = 200.;
    float pDca2 = 200.;
    float etaMin = -4.0;
    float etaMax = -2.5;
    float ptMin = 0.;
    float rAbs = fwdTrack.rAtAbsorberEnd();
    float pDca = fwdTrack.pDca();
    float chi2Max = 35.;
    if (fwdTrack.chi2MatchMCHMFT() > chi2Max)
      return false;
    TLorentzVector p;
    p.SetXYZM(fwdTrack.px(), fwdTrack.py(), fwdTrack.pz(), mMu);
    float eta = p.Eta();
    float pt = p.Pt();
    float pDcaMax = rAbs < rAbsMid ? pDca1 : pDca2;
    if (eta < etaMin || eta > etaMax)
      return false;
    if (pt < ptMin)
      return false;
    if (rAbs < rAbsMin || rAbs > rAbsMax)
      return false;
    if (pDca > pDcaMax)
      return false;
    return true;
  }

/*
  bool isMCHTrackWithinAcceptance(const ForwardTracks::iterator& fwdTrack)
  {
    float etaMinMFT = -3.6;
    float etaMaxMFT = -2.45;
    //float etaMinMCH = -4.0;
    //float etaMaxMCH = -2.5;

    int mchTrackId = fwdTrack.matchMCHTrackId(); // Get the ID of the matching MCH track
    if (mchTrackId < 0) // Check if there is a valid MCH track
     { 
      return false;
     }

    TLorentzVector p;
    p.SetXYZM(fwdTrack.px(), fwdTrack.py(), fwdTrack.pz(), mMu);
    float etaMCH = p.Eta();
    
    if (etaMCH < etaMinMFT || etaMCH > etaMaxMFT)
    {
      return false;
    }
  return true;
  }

*/


  void processCandidateFwd(CandidatesFwd::iterator const& cand,
                           const ForwardTracks::iterator& tr1, const ForwardTracks::iterator& tr2)
  {
    const auto& ampsV0A = cand.amplitudesV0A();
    const auto& ampsRelBCsV0A = cand.ampRelBCsV0A();
    auto nAmpsV0A = ampsV0A.size();
    for (auto i = 0; i < nAmpsV0A; ++i) {
      if (std::abs(ampsRelBCsV0A[i]) <= 1 && ampsV0A[i] > 100.)
        return;
    }
    float m1 = mMu;
    float m2 = mMu;
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m1);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m2);
    TLorentzVector p = p1 + p2;
    
    if (!isMuonSelected(tr1))
      return;
    if (!isMuonSelected(tr2))
      return;
    int nMIDs = 0;
    if (tr1.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (tr2.chi2MatchMCHMID() > 0)
      nMIDs++;
    if (nMIDs != 2)
      return;
    bool isUS = (tr1.sign() * tr2.sign()) < 0;
    if (!isUS)
      return;
    float m = p.M();
    float pt = p.Pt();
    //float y = rapidity(p.M(), p.Px(), p.Py(), p.Pz());
    hr.fill(HIST("hM"), m);
    hr.fill(HIST("hPt"), pt);
    //hr.fill(HIST("hY"), y);
    // restric the kine
    //if (p.M() > 2 && p.M() < 6 && p.Pt() < 5) mptSel(m,pt,y);
    if (p.M() > 2 && p.M() < 6 && p.Pt() < 5) mptSel(m,pt);
  }

  template <typename TTracks>
  void collectCandIDs(std::unordered_map<int32_t, std::vector<int32_t>>& tracksPerCand, TTracks& tracks)
  {
    for (const auto& tr : tracks) {
      int32_t candId = tr.udCollisionId();
      if (candId < 0) {
        continue;
      }
      tracksPerCand[candId].push_back(tr.globalIndex());
    }
  }

  void process(CandidatesFwd const& eventCandidates,
               ForwardTracks const& fwdTracks)
  {
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      auto cand = eventCandidates.iteratorAt(candID);
      auto tr1 = fwdTracks.iteratorAt(trId1);
      auto tr2 = fwdTracks.iteratorAt(trId2);
      processCandidateFwd(cand, tr1, tr2);
    }
  }

  PROCESS_SWITCH(MuMuForward, process, "", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MuMuForward>(cfgc)};
}
