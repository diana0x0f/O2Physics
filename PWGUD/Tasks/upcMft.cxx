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

// c++ headers
#include <iostream>

// root headers
#include <TH1F.h>
#include <TH2F.h>
#include "TMath.h"
#include "TLorentzVector.h"

// framework headers
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

// define namespaces
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct BcTracks { 
  std::vector<aod::AmbiguousMFTTrack> tracks;
  float timeFT0C = -999; // Default is -999, because FT0 time can be negative, so the usual -1 is not ideal.
  float timeFT0A = -999;
};

//-----------------------------------------------------------------------------------------
struct UpcMftRec {

  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  // initialisation
  // -- define histograms
  void init(o2::framework::InitContext&)
  {
    // histograms for the reconstructed rho
    registry.add("hMassRhoRec", "Mass of reconstructed #rho^{0}; mass (GeV/c^{2}); Entries", kTH1F, {{150, 0., 1.5}});
    registry.add("hMassRhoRec4", "Mass of reconstructed #rho (4 tracks); mass (GeV/c^{2}); Entries", kTH1F, {{200, 0., 2.0}});
    registry.add("hPtRhoRec", "p_{T} of reconstructed #rho^{0};#it{p}_{T} (GeV/c); Entries", kTH1F, {{200, 0., 1.0}});
    registry.add("hPtRhoRec4", "p_{T} of reconstructed #rho (4 tracks);#it{p}_{T} (GeV/c); Entries", kTH1F, {{400, 0., 2.0}});
    registry.add("hMassPtRhoRec","Mass vs. p_{T};Mass of reconstructed #rho^{0} (GeV/c^{2});p_{T} of reconstructed #rho^{0} (GeV/c)", {HistType::kTH2F, {{100, 0., 2.0}, {200, 0., 2.0}}});
    registry.add("hMassPtRhoRec4","Mass vs. p_{T} (4 tracks);Mass of reconstructed #rho^{0} (GeV/c^{2});p_{T} of reconstructed #rho^{0} (GeV/c)", {HistType::kTH2F, {{100, 0., 2.0}, {200, 0., 2.0}}});
    registry.add("hMassPtRhoRecElectrons","#gamma#gamma#rightarrow#it{e}^{+}#it{e}^{-} (2 tracks);mass (GeV/c^{2}) ;#it{p}_{T} (GeV/c)", {HistType::kTH2F, {{150, 0., 1.5}, {200, 0., 1.0}}});
    registry.add("hTimeFT0C", "FT0C time (ns); Entries", kTH1F, {{200, -1.0, 1.0}});
    registry.add("hTimeFT0C4", "FT0C time (ns) (4 tracks); Entries", kTH1F, {{200, -1.0, 1.0}});
    registry.add("hTracksPerBC", "Tracks per BC; Tracks per BC; Entries", kTH1F, {{10000, 0, 10000}});
    registry.add("hChargeDistribution4Tracks", "Charge distribution of 4-track events; q_{1}+q_{2}+q_{3}+q_{4}; Entries", kTH1F, {{9, -4.5, 4.5}});
    // control plots for 2-track events
    registry.add("ControlHistos/hTrackChi2", "Track #chi^{2}; #chi^{2}; Entries", kTH1F, {{210, -0.5, 20.5}});
    registry.add("ControlHistos/hTrackNClusters", "Number of clusters per track; Number of clusters; Entries", kTH1F, {{10, 0.5, 10.5}});
    registry.add("ControlHistos/hTrackPt", "Track p_{T}; p_{T} (GeV/c); Entries", kTH1F, {{200, 0, 1.0}});
    registry.add("ControlHistos/hTrackEta", "Track #eta; #eta; Entries", kTH1F, {{50, -4, -2}});
    registry.add("ControlHistos/hTrackEtaPos", "Positive track #eta; #eta; Entries", kTH1F, {{50, -4, -2}});
    registry.add("ControlHistos/hTrackEtaNeg", "Negative track #eta; #eta; Entries", kTH1F, {{50, -4, -2}});
    registry.add("ControlHistos/hTrackPhi", "Track #phi; #phi; Entries", kTH1F, {{100, -3.2, 3.2}});
    registry.add("ControlHistos/hTrackPhiPos", "Positive track #phi; #phi; Entries", kTH1F, {{100, -3.2, 3.2}});
    registry.add("ControlHistos/hTrackPhiNeg", "Negative track #phi; #phi; Entries", kTH1F, {{100, -3.2, 3.2}});
    }

  void process(aod::BCs const& bcs,
	       aod::Collisions const& collisions,
	       aod::MFTTracks const& tracksMFT,
	       aod::AmbiguousMFTTracks const& ambiguousTracksMFT,
	       aod::FT0s const& ft0s) {
    
    // std::cout << " JGC: there are " << ambiguousTracksMFT.size() << " ambiguous MFT tracks, "
	  //     << " collisions: " << collisions.size()
	  //     << " and MFT tracks " << tracksMFT.size()
	  //     << " BCs: " << bcs.size()
	  //     << " FT0s: " << ft0s.size()
	  //     << std::endl;
    
    // define bc maps
    const int nBCs = bcs.size();
    
    std::map<int, BcTracks> bcGroups; // Map of tracks that belong to a given bunch crossing. One track can be associated with multiple BCs.

    for (auto& ambiguousTrackMFT : ambiguousTracksMFT) { 
        auto mftBCs = ambiguousTrackMFT.bc(); 
        for (auto& mftBC : mftBCs) {
          int bcId = mftBC.globalBC(); // Get global BC that will be used for grouping.
          auto element = bcGroups.find(bcId); // Find the iterator poiniting to value associated with this BCid.
          if(element == bcGroups.end()){ // If this BCid was encountered for the first time, create a vector containing the current track
            BcTracks bct;
            bct.tracks.push_back(ambiguousTrackMFT);
            bcGroups.emplace(std::make_pair(bcId, bct));
          } else { // If this BcId was already encountered, add the track to the existing vector
            element->second.tracks.push_back(ambiguousTrackMFT);
          }
        }
    }
    
    for (auto& ft0 : ft0s) {
      int collisionBC = (bcs.iteratorAt(ft0.bcId())).globalBC(); // Get global BC that will be used for grouping.
      auto element = bcGroups.find(collisionBC); // Find the iterator poiniting to value associated with this BCid.
      if(element == bcGroups.end()) { 
        //std::cout << " There is no global BC id for FT0. " << std::endl;
      } else { 
        element->second.timeFT0C = ft0.timeC();
        element->second.timeFT0A = ft0.timeA();
      }     
    }

    BcTracks defaultBcTrack; // Struct with default values used later to filter out default FT0 time
    for (auto& bcGroup : bcGroups) {
      registry.fill(HIST("hTracksPerBC"), bcGroup.second.tracks.size());
      // std::cout << " ----------DMK: Global index " << bcGroup.first << ", track count:  " << bcGroup.second.tracks.size() << ", FT0 time is: " << bcGroup.second.timeFT0C << std::endl;
      if (bcGroup.second.timeFT0C != defaultBcTrack.timeFT0C && std::abs(bcGroup.second.timeFT0C)<1 && std::abs(bcGroup.second.timeFT0A)>5 &&bcGroup.second.tracks.size() == 2) {
        // Fill control histos
        for (int i=0; i<2; i++) {
          auto mftTrack = bcGroup.second.tracks[i].mfttrack();
          registry.fill(HIST("ControlHistos/hTrackChi2"), mftTrack.chi2());
          registry.fill(HIST("ControlHistos/hTrackNClusters"), mftTrack.nClusters());
          registry.fill(HIST("ControlHistos/hTrackPt"), mftTrack.pt());
          registry.fill(HIST("ControlHistos/hTrackEta"), mftTrack.eta());
          registry.fill(HIST("ControlHistos/hTrackPhi"), mftTrack.phi());
          if (mftTrack.sign() > 0) {
            registry.fill(HIST("ControlHistos/hTrackEtaPos"), mftTrack.eta());
            registry.fill(HIST("ControlHistos/hTrackPhiPos"), mftTrack.phi());
          }
          if (mftTrack.sign() < 0) {
            registry.fill(HIST("ControlHistos/hTrackEtaNeg"), mftTrack.eta());
            registry.fill(HIST("ControlHistos/hTrackPhiNeg"), mftTrack.phi());
          }
        }
        // Get tracks for the two-track case
	      auto mftTrack1 = bcGroup.second.tracks[0].mfttrack();
        auto mftTrack2 = bcGroup.second.tracks[0].mfttrack();
	      // Make a rho candidate
	      const float mPion = 0.13957; // GeV/c2
        if (mftTrack1.sign()+mftTrack2.sign() != 0) continue; // Skip if total charge is not 0 
	      TLorentzVector pi1LV;
	      TLorentzVector pi2LV;
	      pi1LV.SetPtEtaPhiM(mftTrack1.pt(),mftTrack1.eta(),mftTrack1.phi(),mPion);
	      pi2LV.SetPtEtaPhiM(mftTrack2.pt(),mftTrack2.eta(),mftTrack2.phi(),mPion);
	      TLorentzVector rhoLV = pi1LV+pi2LV;        
	      // Fill histos
        registry.fill(HIST("hMassRhoRec"), rhoLV.M());
	      registry.fill(HIST("hPtRhoRec"), rhoLV.Pt());
        registry.fill(HIST("hMassPtRhoRec"), rhoLV.M(), rhoLV.Pt());
	      registry.fill(HIST("hTimeFT0C"), bcGroup.second.timeFT0C);
        
      }
      if (bcGroup.second.timeFT0C != defaultBcTrack.timeFT0C  && std::abs(bcGroup.second.timeFT0C)<1 && std::abs(bcGroup.second.timeFT0A)>5 && bcGroup.second.tracks.size() == 4) { //4 tracks events
        // Get tracks
	      auto mftTrack1 = bcGroup.second.tracks[0].mfttrack();
	      auto mftTrack2 = bcGroup.second.tracks[1].mfttrack();
        auto mftTrack3 = bcGroup.second.tracks[2].mfttrack();
        auto mftTrack4 = bcGroup.second.tracks[3].mfttrack();
	      // Make a rho candidate
	      const float mPion = 0.13957; // GeV/c2
        TLorentzVector pi1LV;
	      TLorentzVector pi2LV;
        TLorentzVector pi3LV;
        TLorentzVector pi4LV;
	      pi1LV.SetPtEtaPhiM(mftTrack1.pt(),mftTrack1.eta(),mftTrack1.phi(),mPion);
	      pi2LV.SetPtEtaPhiM(mftTrack2.pt(),mftTrack2.eta(),mftTrack2.phi(),mPion);
        pi3LV.SetPtEtaPhiM(mftTrack3.pt(),mftTrack3.eta(),mftTrack3.phi(),mPion);
        pi2LV.SetPtEtaPhiM(mftTrack4.pt(),mftTrack4.eta(),mftTrack4.phi(),mPion);
	      TLorentzVector rhoLV = pi1LV+pi2LV+pi3LV+pi4LV;
	      // Fill histos
        registry.fill(HIST("hChargeDistribution4Tracks"), mftTrack1.sign()+mftTrack2.sign()+mftTrack3.sign()+mftTrack4.sign());
	      if (mftTrack1.sign()+mftTrack2.sign()+mftTrack3.sign()+mftTrack4.sign() != 0) continue; // Skip if total charge is not 0 
	      registry.fill(HIST("hMassRhoRec4"), rhoLV.M());
	      registry.fill(HIST("hPtRhoRec4"), rhoLV.Pt());
        registry.fill(HIST("hMassPtRhoRec4"), rhoLV.M(), rhoLV.Pt());
        registry.fill(HIST("hTimeFT0C4"), bcGroup.second.timeFT0C);
        }
      if (bcGroup.second.timeFT0C != defaultBcTrack.timeFT0C  && std::abs(bcGroup.second.timeFT0C)<1 && std::abs(bcGroup.second.timeFT0A)>5 && bcGroup.second.tracks.size() == 2) { //electron mass hypothesis
        // Get tracks
	      auto mftTrack1 = bcGroup.second.tracks[0].mfttrack();
	      auto mftTrack2 = bcGroup.second.tracks[1].mfttrack();
	      // Make a rho candidate
	      const float mElectron = 0.000511; // GeV/c2
        if (mftTrack1.sign()+mftTrack2.sign() != 0) continue; // Skip if total charge is not 0 
	      TLorentzVector el1LV;
	      TLorentzVector el2LV;
	      el1LV.SetPtEtaPhiM(mftTrack1.pt(),mftTrack1.eta(),mftTrack1.phi(),mElectron);
	      el2LV.SetPtEtaPhiM(mftTrack2.pt(),mftTrack2.eta(),mftTrack2.phi(),mElectron);
	      TLorentzVector rhoLV = el1LV+el2LV;
	      // Fill histos
	      registry.fill(HIST("hMassPtRhoRecElectrons"), rhoLV.M(), rhoLV.Pt());
      }
    }
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
   // o2::framework::adaptAnalysisTask<UpcMftGen>(cfgc)};
   o2::framework::adaptAnalysisTask<UpcMftRec>(cfgc)};
}
