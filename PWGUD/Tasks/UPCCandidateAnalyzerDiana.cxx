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
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"

#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandAnalyzerDiana {

  float ft0DummyTime = 32.767f;
  float fT0CBBlower = -1.0; // ns
  float fT0CBBupper = 1.0;  // ns

  std::unordered_map<int32_t, float> pdgsMass;

  enum pdgs {
    kPdgElectron = 11,
    kPdgMuon = 13,
    kPdgTau = 15,
    kPdgPion = 211
  };

  // selection flags for processCandidate()
  enum selections {
    kSelUnlikeSign,
    kSelNoFT0,
    kSelFT0C,
    kSelNoFV0A,
    kSelNoFDD,
    kSelPID,
    kSelPt,
    kSelRabs,
    kNSelectors
  };

  Configurable<int32_t> fPrimaryPdg{"primaryPdg", 13, "Set 'primary' PDG code: e.g. 15 for ditau production"};
  Configurable<int32_t> fTargetPdg{"targetPdg", 13, "Target particle PDG for 'histSwitch' p_T distributions: e.g. electrons (11) from tau decays"};
  Configurable<int32_t> fHistSwitch{"histSwitch", 2, "What information to collect: 0 -- pair mass, 1 -- p_T of target particle, 2 -- both"};

  float fMinPt = 0.;
  float fMaxPt = 1.;

  static constexpr int32_t nBinsMass = 500;
  static constexpr float minMass = 0;
  static constexpr float maxMass = 5;

  static constexpr int32_t nBinsPt = 500;
  static constexpr float minPt = 0;
  static constexpr float maxPt = 5;

  HistogramRegistry registry{
    "registry",
    {
     // separate selectors stored in "SelCounter" (see init())
     {"Selection/PairMass/All", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/UnlikeSign", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/NoFT0", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/NoFV0A", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/NoFDD", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/Rabs", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"Selection/PairMass/PDca", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     //basic kinematics for all tracks
     {"ControlPlots/allTracks/Pt", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"ControlPlots/allTracks/Eta", "#eta of pair of #mu;#eta;", {HistType::kTH1D, {{100, -6., 6.}}}},
     {"ControlPlots/allTracks/Eta1", "#eta of #mu;#eta;", {HistType::kTH1D, {{100, -6., 6.}}}},
     {"ControlPlots/allTracks/Eta2", "#eta of #mu;#eta;", {HistType::kTH1D, {{100, -6., 6.}}}},
     {"ControlPlots/allTracks/Phi", "#phi of pair of #mu;#phi;", {HistType::kTH1D, {{100, -3.15, 3.15}}}},
     {"ControlPlots/allTracks/Phi1", "#phi of #mu;#phi;", {HistType::kTH1D, {{100, 0, 6.28}}}},
      {"ControlPlots/allTracks/Phi2", "#phi of #mu;#phi;", {HistType::kTH1D, {{100, 0, 6.28}}}},
     {"ControlPlots/allTracks/Chi2", ";#chi^{2};", {HistType::kTH1D, {{100, 0., 100.}}}},
     {"ControlPlots/allTracks/Charge", ";q;", {HistType::kTH1D, {{2, -1., 1.}}}},
     {"ControlPlots/allTracks/Mass", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     {"ControlPlots/allTracks/FT0Ctime", "; time, ns;", {HistType::kTH1D, {{30, -1., 1.}}}},
     {"ControlPlots/allTracks/FT0Atime", "; time, ns;", {HistType::kTH1D, {{30, -1., 1.}}}},
     //basic kinematics for selected tracks
     {"ControlPlots/selectedTracks/Pt", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"ControlPlots/selectedTracks/Eta", "#eta of pair of #mu;#eta;", {HistType::kTH1D, {{100, -6., 6.}}}},
     {"ControlPlots/selectedTracks/Phi", "#phi of pair of #mu;#phi;", {HistType::kTH1D, {{100, -3.15, 3.15}}}},
     {"ControlPlots/selectedTracks/Chi2", ";#chi^{2};", {HistType::kTH1D, {{100, 0., 100.}}}},
     {"ControlPlots/selectedTracks/Charge", ";q;", {HistType::kTH1D, {{2, -1., 1.}}}},
     {"ControlPlots/selectedTracks/Mass", ";#it{m}, GeV;", {HistType::kTH1D, {{nBinsMass, minMass, maxMass}}}},
     //
     {"Selection/TargetPt/All", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/UnlikeSign", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/NoFT0", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/NoFV0A", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/NoFDD", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/Rabs", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}},
     {"Selection/TargetPt/PDca", ";#it{p}_{T}, GeV;", {HistType::kTH1D, {{nBinsPt, minPt, maxPt}}}}}};

  using Candidates = soa::Join<o2::aod::UDCollisions, o2::aod::UDCollisionsSels>;
  using FwdTracks = soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTracksExtra>;

  void init(InitContext&)
  {    
    const AxisSpec axisSel{kNSelectors, 0., double(kNSelectors), ""};
    registry.add("Selection/SelCounter", "", kTH1F, {axisSel});
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelUnlikeSign + 1, "kSelUnlikeSign");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelNoFT0 + 1, "kSelNoFT0");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelFT0C + 1, "kSelFT0C");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelNoFV0A + 1, "kSelNoFV0A");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelNoFDD + 1, "kSelNoFDD");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelPID + 1, "kSelPID");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelPt + 1, "kSelPt");
    registry.get<TH1>(HIST("Selection/SelCounter"))->GetXaxis()->SetBinLabel(kSelRabs + 1, "kSelRabs");
    // populate "pdg->particle mass" map
    pdgsMass[kPdgElectron] = 0.000511;
    pdgsMass[kPdgMuon] = 0.10566;
    pdgsMass[kPdgTau] = 1.777;
    pdgsMass[kPdgPion] = 0.13957;
  }


  void fillMassDistr(float m, std::unordered_map<int32_t, bool>& selFlags)
  {
    registry.fill(HIST("Selection/PairMass/All"), m);
    // unlike-sign
    bool selector = selFlags[kSelPt] && selFlags[kSelUnlikeSign];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/UnlikeSign"), m);
    }
    // unlike sign + no FT0
    selector = selector && selFlags[kSelNoFT0];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/NoFT0"), m);
    }
    // unlike sign + FT0C
    selector = selector && selFlags[kSelFT0C];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/FT0C"), m);
    }
    selector = selector && selFlags[kSelNoFV0A];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/NoFV0A"), m);
    }
    selector = selector && selFlags[kSelNoFDD];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/NoFDD"), m);
    }
    selector = selector && selFlags[kSelRabs];
    if (selector) {
      registry.fill(HIST("Selection/PairMass/Rabs"), m);
    }
  }

  void fillPtDistr(float pt, std::unordered_map<int32_t, bool>& selFlags)
  {
    registry.fill(HIST("Selection/TargetPt/All"), pt);
    // unlike-sign
    bool selector = selFlags[kSelPt] && selFlags[kSelUnlikeSign];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/UnlikeSign"), pt);
    }
    // unlike sign + no FT0
    selector = selector && selFlags[kSelNoFT0];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/NoFT0"), pt);
    }
    // unlike sign + FT0C
    selector = selector && selFlags[kSelFT0C];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/FT0C"), pt);
    }
    selector = selector && selFlags[kSelNoFV0A];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/NoFV0A"), pt);
    }
    selector = selector && selFlags[kSelNoFDD];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/NoFDD"), pt);
    }
    selector = selector && selFlags[kSelRabs];
    if (selector) {
      registry.fill(HIST("Selection/TargetPt/Rabs"), pt);
    }
  }

  template <int32_t processSwitch, typename TTrack1, typename TTrack2>
  void processCandidateBasic(Candidates::iterator const& cand, TTrack1& tr1, TTrack2& tr2)
  {

    float m1 = pdgsMass[fTargetPdg];
    float m2 = pdgsMass[fTargetPdg];

    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m1);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m2);
    float eta1 = eta(tr1.px(), tr1.py(), tr1.pz());
    float eta2 = eta(tr2.px(), tr2.py(), tr2.pz());
    float phi1 = phi(tr1.px(), tr1.py());
    float phi2 = phi(tr2.px(), tr2.py());
    TLorentzVector p = p1 + p2;

    // fill control histos for all tracks (no cuts)
    registry.fill(HIST("ControlPlots/allTracks/Pt"), p.Pt());
    registry.fill(HIST("ControlPlots/allTracks/Eta"), p.Eta());
    registry.fill(HIST("ControlPlots/allTracks/Eta1"), eta1);
    registry.fill(HIST("ControlPlots/allTracks/Eta2"), eta2);
    registry.fill(HIST("ControlPlots/allTracks/Phi"), p.Phi());
    registry.fill(HIST("ControlPlots/allTracks/Phi1"), phi1);
    registry.fill(HIST("ControlPlots/allTracks/Phi2"), phi2);
    registry.fill(HIST("ControlPlots/allTracks/Charge"), tr1.sign());
    registry.fill(HIST("ControlPlots/allTracks/Charge"), tr2.sign());
    registry.fill(HIST("ControlPlots/allTracks/Chi2"), tr1.chi2());
    registry.fill(HIST("ControlPlots/allTracks/Chi2"), tr2.chi2());
    registry.fill(HIST("ControlPlots/allTracks/FT0Ctime"), cand.timeFT0C());
    registry.fill(HIST("ControlPlots/allTracks/FT0Atime"), cand.timeFT0A());
    registry.fill(HIST("ControlPlots/allTracks/Mass"), p.M());

    //  fill control histos for selected tracks
    if (tr1.sign() * tr2.sign() < 0 &&
        17.5 < tr1.rAtAbsorberEnd() < 89.5 && 17.5 < tr2.rAtAbsorberEnd() < 89.5 &&
        cand.timeFT0C() > fT0CBBlower && cand.timeFT0C() < fT0CBBupper && 
        eta1 < -2.5 && eta1 > -4.0 && eta2 < -2.5 && eta2 > -4.0 &&
        p.Pt() > fMinPt && p.Pt() < fMaxPt
        ) {
      registry.fill(HIST("ControlPlots/selectedTracks/Pt"), p.Pt());
      registry.fill(HIST("ControlPlots/selectedTracks/Eta"), p.Eta());
      registry.fill(HIST("ControlPlots/selectedTracks/Phi"), p.Phi());
      registry.fill(HIST("ControlPlots/selectedTracks/Charge"), tr1.sign());
      registry.fill(HIST("ControlPlots/selectedTracks/Charge"), tr2.sign());
      registry.fill(HIST("ControlPlots/selectedTracks/Chi2"), tr1.chi2());
      registry.fill(HIST("ControlPlots/selectedTracks/Chi2"), tr2.chi2());
      registry.fill(HIST("ControlPlots/selectedTracks/Mass"), p.M());
    }
  }
 

  template <int32_t processSwitch, typename TTrack1, typename TTrack2>
  void processCandidate(Candidates::iterator const& cand, TTrack1& tr1, TTrack2& tr2)
  {
    std::unordered_map<int32_t, bool> selFlags; // holder of selection flags
    // unlike-sign tracks requirement
    selFlags[kSelUnlikeSign] = (tr1.sign() * tr2.sign()) < 0;
    // Rabs selection
    selFlags[kSelRabs] = 17.5 < tr1.rAtAbsorberEnd() < 89.5 && 17.5 < tr2.rAtAbsorberEnd() < 89.5;
    // check FT0 signal
    bool hasNoFT0 = true;
    bool isBB = cand.bbFT0A() || cand.bbFT0C();
    bool isBG = cand.bgFT0A() || cand.bgFT0C();
    hasNoFT0 = !isBB && !isBG;
    // if there is a signal, candidate passes if timeA is dummy
    // and timeC is between +/- 1 ns
    bool checkA = std::abs(cand.timeFT0A() - ft0DummyTime) < 1e-3;
    bool checkC = cand.timeFT0C() > fT0CBBlower && cand.timeFT0C() < fT0CBBupper;
    hasNoFT0 = checkA && checkC;
    selFlags[kSelNoFT0] = hasNoFT0;
    // check FT0C signal separately
    bool hasFT0C = true;
    hasFT0C = checkC;
    selFlags[kSelFT0C] = hasFT0C;
    // check FDD signal
    bool hasNoFDD = true;
    bool isBB2 = cand.bbFDDA() || cand.bbFDDC();
    bool isBG2 = cand.bgFDDA() || cand.bgFDDC();
    hasNoFDD = !isBB2 && !isBG2;
    selFlags[kSelNoFDD] = hasNoFDD;
    float m1 = pdgsMass[fTargetPdg];
    float m2 = pdgsMass[fTargetPdg];
    TLorentzVector p1, p2;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), m1);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), m2);
    TLorentzVector p = p1 + p2;
    selFlags[kSelPt] = p.Pt() > fMinPt && p.Pt() < fMaxPt;
    // selection counters
    if (selFlags[kSelUnlikeSign]) {
      registry.fill(HIST("Selection/SelCounter"), kSelUnlikeSign, 1);
    }
    if (selFlags[kSelNoFT0]) {
      registry.fill(HIST("Selection/SelCounter"), kSelNoFT0, 1);
    }
    if (selFlags[kSelFT0C]) {
      registry.fill(HIST("Selection/SelCounter"), kSelFT0C, 1);
    }
    if (selFlags[kSelNoFV0A]) {
      registry.fill(HIST("Selection/SelCounter"), kSelNoFV0A, 1);
    }
    if (selFlags[kSelNoFDD]) {
      registry.fill(HIST("Selection/SelCounter"), kSelNoFDD, 1);
    }
    if (selFlags[kSelPID]) {
      registry.fill(HIST("Selection/SelCounter"), kSelPID, 1);
    }
    if (selFlags[kSelPt]) {
      registry.fill(HIST("Selection/SelCounter"), kSelPt, 1);
    }
    if (selFlags[kSelRabs]) {
      registry.fill(HIST("Selection/SelCounter"), kSelRabs, 1);
    }
    // collect mass distributions if needed
    if (fHistSwitch == 0 || fHistSwitch == 2) {
      float m = p.M();
      fillMassDistr(m, selFlags);
    }
    // collect pt distributions if needed
    if (fHistSwitch == 1 || fHistSwitch == 2) {
      float pt1 = p1.Pt();
      float pt2 = p2.Pt();
      fillPtDistr(pt1, selFlags);
      fillPtDistr(pt2, selFlags);
    }
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

  // process candidates with 2 muon tracks
  void processFwd(Candidates const& eventCandidates,
                  FwdTracks const& fwdTracks)
  {
    
    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    // assuming that candidates have exatly 2 muon tracks
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      const auto& cand = eventCandidates.iteratorAt(candID);
      const auto& tr1 = fwdTracks.iteratorAt(trId1);
      const auto& tr2 = fwdTracks.iteratorAt(trId2);
      processCandidate<0>(cand, tr1, tr2);
    }
  }

  // basic process with 2 muon tracks
  void processFwdBasic(Candidates const& eventCandidates,
                        FwdTracks const& fwdTracks)
  {

    std::unordered_map<int32_t, std::vector<int32_t>> tracksPerCand;
    collectCandIDs(tracksPerCand, fwdTracks);

    // assuming that candidates have exatly 2 muon tracks
    for (const auto& item : tracksPerCand) {
      int32_t trId1 = item.second[0];
      int32_t trId2 = item.second[1];
      int32_t candID = item.first;
      const auto& cand = eventCandidates.iteratorAt(candID);
      const auto& tr1 = fwdTracks.iteratorAt(trId1);
      const auto& tr2 = fwdTracks.iteratorAt(trId2);
      processCandidateBasic<0>(cand, tr1, tr2);
    }
  }

  PROCESS_SWITCH(UpcCandAnalyzerDiana, processFwdBasic, "Analyse forward candidates using basic cuts", true)
  PROCESS_SWITCH(UpcCandAnalyzerDiana, processFwd, "Analyse forward candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcCandAnalyzerDiana>(cfgc)};
}
